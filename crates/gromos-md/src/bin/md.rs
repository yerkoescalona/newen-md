//! md - Molecular Dynamics simulation engine (GROMOS-compatible CLI)
//!
//! Usage: md @topo <topology> @conf <coordinates> @input <parameters> [@fin <final_conf>]
//!           [@trc <trajectory>] [@tre <energies>] [@trf <forces>] [@trv <velocities>]
//!           [@verb <level>]
//!
//! The main GROMOS-RS molecular dynamics simulation program.
//! Command-line interface matches GROMOS md++ conventions.
//! All simulation parameters are read from the @input (.imd/.in) file.
//!
//! The run itself is assembled by the `gromos-run` crate (`read_imd` → `RunRecipe` → plan →
//! algorithms → `start`), the same path `py-gromos`'s `Simulation` takes; this binary owns the
//! CLI, the output files and the out-of-band GAMD/EDS blocks. `@dump` prints the recipe and the
//! plan as JSON; the effective recipe is written next to the `.tre` as `<tre>.recipe.toml`.

use gromos::{
    algorithm::SimulationState,
    configuration::BoxType,
    io::{
        coordinate::read_coordinates,
        energy::{EnergyWriter, VolumePressure},
        force::ForceWriter,
        free_energy::FreeEnergyWriter,
        topology::{build_topology, read_topology_file},
        trajectory::TrajectoryWriter,
        EdsBlock, EdsStatsWriter, EdsVrWriter, GamdBlock, GamdBoostWriter, GamdStatsWriter,
    },
    math::{truncoct_triclinic_rotmat, Vec3},
    validation::validate_energy,
};
use gromos_run::{
    build_sequence_from_imd, plan_to_json, prepare_system, start, Built, PassthroughPolicy,
    PrepareNote, Prepared, RunError, RunInputs, RunOptions,
};
use std::env;
use std::path::PathBuf;
use std::process;
use std::time::Instant;

/// Simple drop-timer that logs elapsed time at debug level.
struct Timer {
    name: &'static str,
    start: Instant,
}
impl Timer {
    fn new(name: &'static str) -> Self {
        Self {
            name,
            start: Instant::now(),
        }
    }
}
impl Drop for Timer {
    fn drop(&mut self) {
        log::debug!(
            "{} took {:.3} ms",
            self.name,
            self.start.elapsed().as_secs_f64() * 1000.0
        );
    }
}

fn print_usage() {
    eprintln!("md - Molecular Dynamics simulation (GROMOS-compatible)");
    eprintln!();
    eprintln!("Usage: md @topo <topology> @conf <coordinates> @input <parameters>");
    eprintln!("          [@fin <final_conf>] [@trc <trajectory>] [@tre <energies>]");
    eprintln!("          [@trf <forces>] [@trv <velocities>] [@verb <level>]");
    eprintln!();
    eprintln!("Input files:");
    eprintln!("  @topo       Molecular topology file (.topo/.top)");
    eprintln!("  @pttopo     Perturbation topology file (.ptp)");
    eprintln!("  @conf       Initial coordinates and restart data (.cnf/.g96)");
    eprintln!("  @input      Input parameter file (.imd/.in) — all simulation settings");
    eprintln!("  @posresspec Position restraints specification");
    eprintln!("  @refpos     Position restraints reference positions");
    eprintln!("  @distrest   Distance restraints specification");
    eprintln!("  @angrest    Angle restraints specification");
    eprintln!("  @dihrest    Dihedral restraints specification");
    eprintln!("  @colvarres  Collective variable restraints");
    eprintln!("  @jval       J-value restraints specification");
    eprintln!("  @xray       X-ray restraints specification");
    eprintln!("  @sym        Symmetry restraints specification");
    eprintln!("  @rdc        RDC restraints specification");
    eprintln!("  @order      Order-parameter restraints specification");
    eprintln!("  @lud        Local elevation umbrella database");
    eprintln!("  @led        Local elevation coordinate definition");
    eprintln!("  @bsleus     Ball & Stick Local Elevation topology");
    eprintln!("  @friction   Atomic friction coefficients");
    eprintln!("  @qmmm       QM/MM specification file");
    eprintln!("  @gamd       GaMD restraints specification");
    eprintln!();
    eprintln!("Output files:");
    eprintln!("  @fin        Final configuration output (.cnf)");
    eprintln!("  @trc        Coordinate trajectory (.trc)");
    eprintln!("  @trv        Velocity trajectory (.trv)");
    eprintln!("  @trf        Force trajectory (.trf)");
    eprintln!("  @trs        Special trajectory (.trs)");
    eprintln!("  @tre        Energy trajectory (.tre)");
    eprintln!("  @trg        Free energy trajectory (.trg)");
    eprintln!("  @bae        Block averaged energy trajectory");
    eprintln!("  @bag        Block averaged free-energy trajectory");
    eprintln!();
    eprintln!("Control:");
    eprintln!("  @verb       Verbosity level (0=quiet, 1=info, 2=debug)");
    eprintln!("  @print      Print additional information (pairlist, force)");
    eprintln!("  @version    Print version information");
    eprintln!("  @develop    Run untested development code");
    eprintln!("  @dump       Print the run recipe and algorithm plan (JSON) and exit");
    eprintln!();
    eprintln!("All simulation parameters (timestep, cutoffs, thermostat, etc.)");
    eprintln!("are specified in the @input file, following GROMOS conventions.");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  md @topo system.topo @conf initial.cnf @input run.imd");
    eprintln!("  md @topo system.topo @conf initial.cnf @input run.imd @trc out.trc @tre out.tre");
    eprintln!("  md @topo system.topo @conf initial.cnf @input run.imd @verb 1");
}

/// GROMOS-compatible command-line arguments.
/// These are file paths and control flags — all simulation parameters
/// come from the @input file.
#[derive(Debug, Default)]
struct MDArgs {
    // Required input files
    topo_file: String,
    conf_file: String,
    input_file: String,
    // Optional input files
    pttopo_file: Option<String>,
    posresspec_file: Option<String>,
    refpos_file: Option<String>,
    distrest_file: Option<String>,
    angrest_file: Option<String>,
    dihrest_file: Option<String>,
    gamd_file: Option<String>,
    // Output files
    fin_file: Option<String>,
    trc_file: Option<String>,
    trv_file: Option<String>,
    trf_file: Option<String>,
    trs_file: Option<String>,
    tre_file: Option<String>,
    trg_file: Option<String>,
    bae_file: Option<String>,
    bag_file: Option<String>,
    // Control
    verbose: usize,
    print_flags: Vec<String>,
    develop: bool,
    /// `@dump`: print the run recipe and the algorithm plan (JSON) and exit without running.
    dump: bool,
}

fn parse_args(args: Vec<String>) -> Result<MDArgs, String> {
    let mut md_args = MDArgs::default();

    let mut i = 1;
    while i < args.len() {
        let arg = args[i].as_str();

        // Helper: get the next argument value
        macro_rules! next_val {
            ($name:expr) => {{
                i += 1;
                if i >= args.len() {
                    return Err(format!("Missing value for {}", $name));
                }
                args[i].clone()
            }};
        }

        match arg {
            // Required input files
            "@topo" => md_args.topo_file = next_val!("@topo"),
            "@conf" => md_args.conf_file = next_val!("@conf"),
            "@input" | "@imd" => md_args.input_file = next_val!("@input"),
            // Optional input files
            "@pttopo" => md_args.pttopo_file = Some(next_val!("@pttopo")),
            "@posresspec" => md_args.posresspec_file = Some(next_val!("@posresspec")),
            "@refpos" => md_args.refpos_file = Some(next_val!("@refpos")),
            "@distrest" => md_args.distrest_file = Some(next_val!("@distrest")),
            "@angrest" => md_args.angrest_file = Some(next_val!("@angrest")),
            "@dihrest" => md_args.dihrest_file = Some(next_val!("@dihrest")),
            "@gamd" => md_args.gamd_file = Some(next_val!("@gamd")),
            // Output files
            "@fin" => md_args.fin_file = Some(next_val!("@fin")),
            "@trc" | "@traj" => md_args.trc_file = Some(next_val!(arg)),
            "@trv" => md_args.trv_file = Some(next_val!("@trv")),
            "@trf" => md_args.trf_file = Some(next_val!("@trf")),
            "@trs" => md_args.trs_file = Some(next_val!("@trs")),
            "@tre" | "@ene" => md_args.tre_file = Some(next_val!(arg)),
            "@trg" => md_args.trg_file = Some(next_val!("@trg")),
            "@bae" => md_args.bae_file = Some(next_val!("@bae")),
            "@bag" => md_args.bag_file = Some(next_val!("@bag")),
            // Control
            "@verb" | "@verbose" | "@v" => {
                md_args.verbose = next_val!("@verb")
                    .parse()
                    .map_err(|_| format!("Invalid value for @verb: {}", args[i]))?;
            },
            "@print" => {
                md_args.print_flags.push(next_val!("@print"));
            },
            "@version" => {
                println!("GROMOS-RS md {}", env!("CARGO_PKG_VERSION"));
                process::exit(0);
            },
            "@dump" => {
                md_args.dump = true;
            },
            "@develop" => {
                md_args.develop = true;
            },
            // Catch unrecognized
            _ if arg.starts_with('@') => {
                return Err(format!("Unknown argument: {}", arg));
            },
            _ => {
                return Err(format!("Unexpected argument: {} (did you forget @?)", arg));
            },
        }
        i += 1;
    }

    if md_args.topo_file.is_empty() {
        return Err("Missing required argument @topo".to_string());
    }
    if md_args.conf_file.is_empty() {
        return Err("Missing required argument @conf".to_string());
    }
    if md_args.input_file.is_empty() {
        return Err("Missing required argument @input".to_string());
    }

    Ok(md_args)
}

/// The `VOLUMEPRESSURE03` block of a `.tre` frame from the current state: mass, temperature,
/// volume, box, pressure, virial and kinetic tensors.
fn volume_pressure(
    state: &gromos_core::State,
    topo: &gromos_core::Topology,
    temp: f64,
) -> VolumePressure {
    let mut vp = VolumePressure {
        mass: topo.mass.iter().sum(),
        temperature: temp,
        volume: state.box_config.volume(),
        pressure: state.pressure(),
        ..Default::default()
    };
    let m = &state.box_config.vectors;
    let vir = &state.virial_tensor;
    let kin = &state.kinetic_energy_tensor;
    for a in 0..3 {
        for b in 0..3 {
            // column i of a gromos-core `Mat3` is GROMOS's row i (see `math.rs`)
            vp.box_vectors[a][b] = m.col(a)[b];
            vp.virial_tensor[a][b] = vir.col(b)[a];
            vp.kinetic_tensor[a][b] = kin.col(b)[a];
        }
    }
    vp
}

fn main() {
    let t_total = Instant::now();
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    // Parse command-line arguments (file paths only, GROMOS style)
    #[allow(unused_mut)] // mutated by the MPI setup (feature `use-mpi`)
    let mut md_args = match parse_args(args) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        },
    };

    // Set up logging
    let filter = match md_args.verbose {
        0 => "info",
        _ => "debug",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(filter)).init();
    // MPI (feature `use-mpi`, run under `mpirun`): every rank runs this whole driver on the
    // same inputs; only the nonbonded pair terms are split across ranks and summed each step
    // (`Forcefield::pair_partition` + `set_nonbonded_reducer`). Rank 0 writes the outputs.
    #[cfg(feature = "use-mpi")]
    let mpi_run = mpi_support::setup(&mut md_args);
    #[cfg(feature = "use-mpi")]
    let mpi_partition: Option<(usize, usize)> = mpi_run.partition();
    #[cfg(not(feature = "use-mpi"))]
    let mpi_partition: Option<(usize, usize)> = None;
    // Only rank 0 writes files; the other ranks write none, not even the default-named ones.
    #[cfg(feature = "use-mpi")]
    let writes_files = mpi_run.rank == 0;
    #[cfg(not(feature = "use-mpi"))]
    let writes_files = true;

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                   GROMOS-RS MD Engine                        ║");
    println!("║           Rust Implementation of GROMOS MD                   ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    log::info!("GROMOS-RS MD simulation starting");
    log::debug!("Verbose level: {}", md_args.verbose);

    // === Read simulation parameters from @input file ===
    println!("Loading input parameters: {}", md_args.input_file);
    let imd = match gromos_run::read_imd(&md_args.input_file) {
        Ok(p) => p,
        Err(e) => {
            log::error!("Failed to read input file: {}", e);
            eprintln!("Error reading input file: {}", e);
            process::exit(1);
        },
    };
    log::debug!(
        "Input parameters: steps={}, dt={}, ntb={}",
        imd.nstlim,
        imd.dt,
        imd.ntb
    );

    let n_steps = imd.nstlim;
    let dt = imd.dt;
    let nstxout = imd.ntwx;
    let nstener = imd.ntwe;
    let nstlog = imd.ntpr;
    // GROMOS convention: a write/print frequency of 0 disables that output entirely.
    let due = |step: usize, freq: usize| freq > 0 && step.is_multiple_of(freq);

    // Derive output file paths (from @args or defaults)
    let trc_file = md_args
        .trc_file
        .clone()
        .unwrap_or_else(|| "md.trc".to_string());
    let tre_file = md_args
        .tre_file
        .clone()
        .unwrap_or_else(|| "md.tre".to_string());

    // === GROMOS initialization order: parameters → topology → coordinates ===

    // === 1. Load topology (read_topology) ===
    println!("Loading topology: {}", md_args.topo_file);
    log::debug!("Reading topology file: {}", md_args.topo_file);
    let _timer = Timer::new("Topology loading");

    let topo_data = match read_topology_file(&md_args.topo_file) {
        Ok(data) => data,
        Err(e) => {
            log::error!("Failed to read topology: {}", e);
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        },
    };

    // Physical constants from the topology PHYSICALCONSTANTS block override the
    // Forcefield defaults — same as gromosXX. Carried through `prepare_system`.
    let physical_constants = topo_data.physical_constants;
    log::info!(
        "Physical constants from topo: four_pi_eps_i={}, kB={}",
        physical_constants.four_pi_eps_i,
        physical_constants.kB
    );
    let topo = build_topology(topo_data);

    // === 2. Load coordinates (read_configuration) ===
    println!("Loading coordinates: {}", md_args.conf_file);
    log::debug!("Reading coordinate file: {}", md_args.conf_file);
    let _timer = Timer::new("Coordinate loading");

    let coord_data = match read_coordinates(&md_args.conf_file) {
        Ok(data) => data,
        Err(e) => {
            log::error!("Failed to read coordinates: {}", e);
            eprintln!("Error reading coordinates: {}", e);
            process::exit(1);
        },
    };
    let n_positions = coord_data.positions.len();
    let n_velocities = coord_data.velocities.len();
    let box_type_label = match coord_data.box_type {
        0 => "vacuum",
        1 => "rectangular",
        2 => "triclinic",
        3 => "truncated octahedron",
        _ => "unknown",
    };

    // === 3. Prepare the system: perturbation topology, truncated-octahedron transform,
    //        NSM from the coordinate file, validation, initial velocities (gromos-run) ===
    let inputs = RunInputs {
        pttopo: md_args.pttopo_file.as_ref().map(PathBuf::from),
        posresspec: md_args.posresspec_file.as_ref().map(PathBuf::from),
        refpos: md_args.refpos_file.as_ref().map(PathBuf::from),
        distrest: md_args.distrest_file.as_ref().map(PathBuf::from),
    };
    let prepared = match prepare_system(&imd, topo, physical_constants, coord_data.into(), &inputs)
    {
        Ok(p) => p,
        Err(e) => {
            if let RunError::Validation { stage, report } = &e {
                report.print();
                report.print_summary();
                log::error!("Fatal errors in {} - cannot continue", stage);
            }
            eprintln!("Error: {}", e);
            process::exit(1);
        },
    };
    for note in &prepared.notes {
        match note {
            PrepareNote::PerturbationLoaded {
                atoms,
                bonds,
                angles,
                impropers,
                dihedrals,
                soft_bonds,
                soft_angles,
                soft_impropers,
            } => println!(
                "  Perturbation topology: {} perturbed atoms, {} bonds, {} angles, \
                 {} impropers, {} dihedrals, {} soft bonds, {} soft angles, {} soft impropers",
                atoms, bonds, angles, impropers, dihedrals, soft_bonds, soft_angles, soft_impropers
            ),
            PrepareNote::PerturbationIgnored { .. } => {
                log::warn!("@pttopo provided but NTG=0 — perturbation topology ignored");
            },
            PrepareNote::NsmAdjusted { imd, coordinates } => {
                println!(
                    "  Adjusting NSM: {} (imd) -> {} (from coordinates)",
                    imd, coordinates
                );
            },
            PrepareNote::EnergyGroupAdded { last_atom } => {
                println!(
                    "  FORCE block: added an additional energy group ending at atom {}",
                    last_atom
                );
            },
            PrepareNote::AtomPairLambdasDiffer { atoms } => {
                log::warn!(
                    "perturbed atom pair {}-{}: different individual lambdas for lj and crf \
                     (gromosXX: not implemented, critical message; its force uses the lj lambda)",
                    atoms.0,
                    atoms.1
                );
            },
            PrepareNote::PerturbedShakePrecalc { constraints } => {
                log::warn!(
                    "{constraints} perturbed distance constraints with PRECALCLAM: perturbed \
                     SHAKE in combination with precalculated lambdas only contributes to \
                     derivatives, not energies (BAR or reweighting in extended TI will be \
                     inappropriate)"
                );
            },
            PrepareNote::Validation { stage, report } => {
                report.print();
                if report.has_errors() {
                    report.print_summary();
                    log::warn!("{} has errors, but continuing", stage);
                } else {
                    log::debug!("{} warnings in {}", report.warnings.len(), stage);
                }
            },
        }
    }

    {
        let topo = &prepared.topology;
        println!("  Solute atoms: {}", topo.num_solute_atoms());
        if topo.num_solvent_molecules() > 0 {
            let nsm = topo.num_solvent_molecules();
            let aps = topo.atoms_per_solvent();
            println!(
                "  Solvent: {} molecules × {} atoms = {} atoms",
                nsm,
                aps,
                nsm * aps
            );
        }
        println!("  Total atoms: {}", topo.num_atoms());
        println!("  Bonds: {}", topo.moltypes[0].bonds.len());
        println!("  Angles: {}", topo.moltypes[0].angles.len());
        println!("  Dihedrals: {}", topo.moltypes[0].proper_dihedrals.len());
        println!("  Impropers: {}", topo.moltypes[0].improper_dihedrals.len());
        println!("  Chargegroups: {}", topo.chargegroups.len());
        println!();
        println!("  Positions loaded: {}", n_positions);
        if n_velocities > 0 {
            println!("  Velocities loaded: {}", n_velocities);
        }
        let box_dims = prepared.box_dims;
        println!(
            "  Box: ({:.4}, {:.4}, {:.4}) nm",
            box_dims.x, box_dims.y, box_dims.z
        );
        println!("  Box type: {}", box_type_label);
        println!();
    }

    // Parse input file for GAMD/EDS blocks
    let gamd_block = {
        log::debug!("Parsing input file for GAMD block: {}", md_args.input_file);
        match GamdBlock::parse_file(&md_args.input_file) {
            Ok(block) => block,
            Err(e) => {
                log::warn!("Failed to parse GAMD block: {}", e);
                None
            },
        }
    };

    let eds_block = {
        log::debug!("Parsing input file for EDS block: {}", md_args.input_file);
        match EdsBlock::parse_file(&md_args.input_file) {
            Ok(block) => block,
            Err(e) => {
                log::warn!("Failed to parse EDS block: {}", e);
                None
            },
        }
    };

    if let Some(ref block) = gamd_block {
        println!();
        println!("GAMD Parameters detected:");
        println!("  Search mode:  {:?}", block.search_mode);
        println!("  Boost form:   {:?}", block.boost_form);
        println!("  Threshold:    {:?}", block.threshold_type);
        println!("  Sigma0 dih:   {:.2}", block.sigma0_dih);
        println!("  Sigma0 tot:   {:.2}", block.sigma0_tot);
        if let (Some(k), Some(e)) = (block.k_tot, block.e_tot) {
            println!("  K_tot:        {:.6}", k);
            println!("  E_tot:        {:.2}", e);
        }
        println!();
    }

    if let Some(ref block) = eds_block {
        println!();
        println!("EDS Parameters detected:");
        println!("  Num states:   {}", block.num_states);
        println!("  Form:         {:?}", block.form);
        println!("  S values:     {:?}", block.s_values);
        println!("  E offsets:    {:?}", block.e_offsets);
        println!("  Temperature:  {:.1} K", block.temperature);
        if block.search_enabled {
            println!(
                "  AEDS enabled: E_max={:.2}, E_min={:.2}",
                block.e_max, block.e_min
            );
        }
        println!();
    }

    // Check for conflicting modes
    if gamd_block.is_some() && eds_block.is_some() {
        log::error!("Cannot enable both GAMD and EDS simultaneously");
        eprintln!("Error: Both GAMD and EDS blocks found in input file");
        eprintln!("       These methods cannot be used together in the same simulation");
        process::exit(1);
    }

    // Create GAMD parameters if enabled
    let mut gamd_params = gamd_block.as_ref().map(|block| {
        log::info!("Creating GAMD parameters from input block");
        block.to_parameters()
    });

    // Create EDS parameters if enabled (needs num_atoms)
    let mut eds_params = if let Some(ref block) = eds_block {
        log::info!("Creating EDS parameters from input block");
        let num_atoms = prepared.topology.num_atoms();
        if block.search_enabled {
            match block.to_aeds_parameters(num_atoms) {
                Ok(aeds) => Some(aeds),
                Err(e) => {
                    log::error!("Failed to create AEDS parameters: {}", e);
                    eprintln!("Error: Failed to create AEDS parameters: {}", e);
                    process::exit(1);
                },
            }
        } else {
            match block.to_parameters(num_atoms) {
                Ok(eds) => {
                    // Wrap in AEDS for uniform handling
                    Some(gromos::eds::AEDSParameters::new(eds, 0.0, 0.0, false))
                },
                Err(e) => {
                    log::error!("Failed to create EDS parameters: {}", e);
                    eprintln!("Error: Failed to create EDS parameters: {}", e);
                    process::exit(1);
                },
            }
        }
    } else {
        None
    };

    // === 4. Build the algorithm sequence (gromos-run: the one builder) ===
    println!("Setting up nonbonded interactions:");
    println!("  Cutoff:      {:.3} nm", imd.rcutl);
    println!("  Epsilon:     {:.2}", 1.0);
    println!("  RF epsilon:  {:.2}", imd.epsrf);
    println!("  RF kappa:    {:.3} nm^-1", imd.appak);
    println!();

    // The binary applies GAMD/EDS itself (out-of-band, after each step), so those blocks may
    // pass through the recipe; anything else unmodelled is an error (3.9 A17).
    let options = RunOptions {
        // INNERLOOP selects which inner-loop kernel gromosXX uses (special solvent loops, CUDA
        // device); it carries no physics, so it rides along instead of stopping the run.
        passthrough: PassthroughPolicy::allow(["GAMD", "EDS", "INNERLOOP"]),
        pair_partition: mpi_partition,
        ..RunOptions::default()
    };
    let built = match build_sequence_from_imd(&imd, &prepared, &inputs, &options) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        },
    };
    let Built {
        sequence: mut md_sequence,
        recipe,
        plan,
        diagnostics,
        summary,
    } = built;
    #[cfg(feature = "use-mpi")]
    mpi_run.install_reducer(&mut md_sequence);
    for note in &diagnostics.notes {
        println!("  NOTE: {}", note);
    }

    // @dump: the recipe and the plan as data, nothing else — the Rust-side A/B (3.9).
    if md_args.dump {
        match (recipe.to_json(), plan_to_json(&plan)) {
            (Ok(r), Ok(p)) => {
                println!("RECIPE");
                println!("{}", r);
                println!("END");
                println!("PLAN");
                println!("{}", p);
                println!("END");
                process::exit(0);
            },
            (Err(e), _) | (_, Err(e)) => {
                eprintln!("Error: {}", e);
                process::exit(1);
            },
        }
    }

    // Echo the effective run next to the energy output (GROMACS's `mdout.mdp` lesson): what
    // was defaulted, what passed through, and every value the engine actually used.
    if writes_files {
        let echo_path = std::path::Path::new(&tre_file).with_extension("recipe.toml");
        match recipe.to_toml() {
            Ok(toml) => {
                let mut text = String::from("# Effective run recipe written by gromos-rs md\n");
                for note in &diagnostics.notes {
                    text.push_str(&format!("# NOTE: {}\n", note));
                }
                text.push_str(&toml);
                if let Err(e) = std::fs::write(&echo_path, text) {
                    log::warn!("could not write {}: {}", echo_path.display(), e);
                }
            },
            Err(e) => log::warn!("could not render the recipe: {}", e),
        }
    }
    let Prepared {
        topology: topo,
        configuration: mut conf,
        ..
    } = prepared;

    println!("  Initial pairlist: {} pairs", summary.initial_pairs);
    println!();
    if summary.position_restraints > 0 {
        println!(
            "  Position restraints: {} atoms, CPOR={:.1} kJ/(mol·nm²)",
            summary.position_restraints, imd.cpor
        );
    }
    if imd.ntdir != 0 {
        println!(
            "  Distance restraints: {} unperturbed, {} perturbed, NTDIR={}, CDIR={:.1}, DIR0={:.3}",
            summary.distance_restraints.0,
            summary.distance_restraints.1,
            imd.ntdir,
            imd.cdir,
            imd.dir0
        );
    }

    let is_minimization = summary.is_minimization;
    if is_minimization {
        println!("Setting up algorithm sequence: Steepest Descent Energy Minimization");
        println!("  NTEM:  {} (steepest descent)", imd.ntem);
        println!("  DELE:  {:.6} kJ/mol", imd.dele);
        println!("  DX0:   {:.6} nm", imd.dx0);
        println!("  DXM:   {:.6} nm", imd.dxm);
        println!("  NMIN:  {}", imd.nmin);
        println!("  FLIM:  {:.6}", imd.flim);
        println!();
    } else {
        println!("Setting up algorithm sequence: Leap-Frog");
    }
    println!("  Sequence: {}", md_sequence.algorithm_names().join(" → "));
    println!();

    // `init()` and step 0 run through `gromos_run::start` at the top of the main loop.
    let mut sim_state = SimulationState::new(dt, n_steps);

    if let Some(t) = &summary.thermostat {
        println!("Setting up thermostat: {}", t.label);
        println!("  Target temp:   {:.1} K", t.temperature);
        println!("  Coupling time: {:.3} ps", t.tau);
        println!("  Thermostat DOF: {:.0}", summary.total_dof);
        println!();
    }
    if let Some(b) = &summary.barostat {
        println!("Setting up barostat: Berendsen");
        println!("  Target pres:   {:.1} bar", b.pressure0);
        println!("  Coupling time: {:.3} ps", b.tau);
        println!();
    }
    if let Some(s) = &summary.shake {
        println!("Setting up constraints: SHAKE");
        println!("  NTC mode:      {:?}", s.ntc);
        println!("  Tolerance:     {:.6}", s.tolerance);
        println!("  Max iter:      {}", s.max_iterations);
        println!();
    }
    if summary.constraints.settle_enabled {
        println!("Setting up constraints: SETTLE (analytical rigid water)");
        println!();
    }
    if summary.constraints.lincs_enabled {
        println!(
            "Setting up constraints: LINCS (solute={}, solvent={}, order={}/{})",
            summary.constraints.solute_lincs,
            summary.constraints.solvent_lincs,
            summary.lincs_orders.0,
            summary.lincs_orders.1
        );
        println!();
    }
    let constraints_enabled = summary.constraints.any();
    let temperature = summary
        .thermostat
        .as_ref()
        .map(|t| t.temperature)
        .unwrap_or(300.0);

    // Setup trajectory writer
    let mut traj_writer = if writes_files {
        match TrajectoryWriter::new(
            &trc_file,
            "GROMOS-RS MD trajectory",
            false, // velocities
            false, // forces
        ) {
            Ok(w) => Some(w),
            Err(e) => {
                eprintln!("Error creating trajectory file: {}", e);
                process::exit(1);
            },
        }
    } else {
        None
    };

    // Velocity trajectory (@trv, WRITETRAJ NTWV): its own file, as gromosXX writes it.
    let mut vel_writer = match (&md_args.trv_file, imd.ntwv) {
        (Some(path), n) if n > 0 => {
            match TrajectoryWriter::new(path, "GROMOS-RS MD velocity trajectory", false, false) {
                Ok(w) => Some(w),
                Err(e) => {
                    eprintln!("Error creating velocity trajectory file: {}", e);
                    process::exit(1);
                },
            }
        },
        _ => None,
    };

    // What the energy trajectories record about this run: the writer, the topology's digest
    // and the energy-group partition (`# energy-provenance` lines, read back by ene_ana).
    let provenance = gromos_io::energy_traj::Provenance {
        writer: format!("gromos-rs md {}", env!("CARGO_PKG_VERSION")),
        topology: gromos_io::energy_traj::sha256_file(&md_args.topo_file)
            .ok()
            .map(|digest| (digest, md_args.topo_file.clone())),
        energy_groups: topo
            .energy_groups
            .iter()
            .filter_map(|g| Some((g.first()? + 1, g.last()? + 1)))
            .collect(),
    };

    // Setup energy writer
    let mut ene_writer = if writes_files {
        match EnergyWriter::new(&tre_file, &imd.title)
            .map(|w| w.with_layout(imd.num_temp_baths, topo.energy_groups.len()))
            .and_then(|w| w.with_provenance(&provenance))
        {
            Ok(w) => Some(w),
            Err(e) => {
                eprintln!("Error creating energy file: {}", e);
                process::exit(1);
            },
        }
    } else {
        None
    };

    // Setup force trajectory writer (if @trf given)
    let mut force_writer = if let Some(ref trf_file) = md_args.trf_file {
        match ForceWriter::new(trf_file, "GROMOS-RS MD forces", true) {
            Ok(w) => {
                println!("  Force output:  {}", trf_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Error creating force file: {}", e);
                process::exit(1);
            },
        }
    } else {
        None
    };

    // Setup GaMD writers if enabled
    let mut gamd_stats_writer = if writes_files && gamd_params.is_some() {
        let stats_file = "gamd_stats.dat";
        match GamdStatsWriter::new(stats_file, "GROMOS-RS GaMD Statistics") {
            Ok(mut w) => {
                w.set_write_interval(10); // Write every 10 steps
                println!("  GaMD stats:   {}", stats_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Warning: Failed to create GaMD stats file: {}", e);
                None
            },
        }
    } else {
        None
    };

    let mut gamd_boost_writer = if writes_files && gamd_params.is_some() {
        let boost_file = "gamd_boost.dat";
        match GamdBoostWriter::new(boost_file, "GROMOS-RS GaMD Boost Potential") {
            Ok(w) => {
                println!("  GaMD boost:   {}", boost_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Warning: Failed to create GaMD boost file: {}", e);
                None
            },
        }
    } else {
        None
    };

    // Setup EDS writers if enabled
    let mut eds_stats_writer = if writes_files && eds_params.is_some() {
        let stats_file = "eds_stats.dat";
        match EdsStatsWriter::new(stats_file, "GROMOS-RS EDS Statistics") {
            Ok(mut w) => {
                w.set_write_interval(10); // Write every 10 steps
                println!("  EDS stats:    {}", stats_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Warning: Failed to create EDS stats file: {}", e);
                None
            },
        }
    } else {
        None
    };

    let mut eds_vr_writer = if writes_files && eds_params.is_some() {
        let vr_file = "eds_vr.dat";
        match EdsVrWriter::new(vr_file, "GROMOS-RS EDS Reference Energy") {
            Ok(w) => {
                println!("  EDS V_R:      {}", vr_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Warning: Failed to create EDS V_R file: {}", e);
                None
            },
        }
    } else {
        None
    };

    // Setup free-energy trajectory writer (@trg) — only when FEP is active (NTG != 0)
    let mut free_energy_writer: Option<FreeEnergyWriter> = if writes_files && imd.ntg != 0 {
        let trg_path = md_args
            .trg_file
            .clone()
            .unwrap_or_else(|| "md.trg".to_string());
        match FreeEnergyWriter::new(&trg_path, &imd.title)
            .map(|w| w.with_layout(imd.num_temp_baths, topo.energy_groups.len()))
            .and_then(|w| w.with_provenance(&provenance))
        {
            Ok(w) => {
                println!("  Free-energy output: {}", trg_path);
                Some(w)
            },
            Err(e) => {
                eprintln!("Error creating free-energy trajectory file: {}", e);
                process::exit(1);
            },
        }
    } else {
        None
    };

    // Parameters summary
    let energy_layout =
        gromos_io::energy_traj::builtin_schema(gromos_io::energy::ENE_VERSION, "ENERTRJ")
            .map(|s| {
                format!(
                    "ENEVERSION {} {}",
                    gromos_io::energy::ENE_VERSION,
                    s.fingerprint()
                )
            })
            .unwrap_or_default();
    if is_minimization {
        println!("Energy Minimization Parameters:");
        println!("  Max steps:     {}", n_steps);
        println!("  Tolerance:     {} kJ/mol", imd.dele);
        println!("  Step size:     {} nm (max {})", imd.dx0, imd.dxm);
        println!("  Traj output:   {}", trc_file);
        println!("  Energy output: {}", tre_file);
        println!("  Energy layout: {}", energy_layout);
        println!();

        println!("╔══════════════════════════════════════════════════════════════╗");
        println!("║              Starting Energy Minimization                    ║");
        println!("╚══════════════════════════════════════════════════════════════╝");
        println!();

        log::info!(
            "Starting energy minimization: max {} steps, dele={} kJ/mol",
            n_steps,
            imd.dele
        );
    } else {
        println!("MD Parameters:");
        println!("  Steps:         {}", n_steps);
        println!("  Time step:     {} ps", dt);
        println!("  Total time:    {} ps", n_steps as f64 * dt);
        println!("  Temperature:   {} K", temperature);
        println!("  Traj output:   {}", trc_file);
        println!("  Energy output: {}", tre_file);
        println!("  Energy layout: {}", energy_layout);
        println!();

        println!("╔══════════════════════════════════════════════════════════════╗");
        println!("║                   Starting MD Simulation                     ║");
        println!("╚══════════════════════════════════════════════════════════════╝");
        println!();

        log::info!("Starting MD simulation: {} steps, dt={} ps", n_steps, dt);
    }

    let start_time = Instant::now();
    let init_elapsed = t_total.elapsed();
    log::info!(
        "Initialization wall time: {:.3} s",
        init_elapsed.as_secs_f64()
    );
    let mut energy_history: Vec<(f64, f64, f64)> = Vec::new();
    let mut minimization_converged = false;
    let mut prev_min_energy = f64::MAX;
    let mut actual_steps: usize = 0;

    // Per-algorithm breakdown, comparable to the gromosXX TIMING block.
    // Costs two clock reads per algorithm per step, so it is opt-in via @verb.
    let timing_enabled = md_args.verbose >= 1;
    md_sequence.enable_timing(timing_enabled);

    // Main simulation loop.
    //
    // gromosXX runs exactly NSTLIM steps. Step 0 here is a *full* step — `start` is
    // `init` + `run_step` — so the range is exclusive: `0..=n_steps` ran NSTLIM+1 steps, which
    // left `@fin` (and a continuation from it) one step ahead of gromosXX. The written frames
    // were unaffected, which is why the reference tests never saw it.
    if n_steps == 0 {
        // NSTLIM=0: no dynamics, but the energies of the initial configuration are still wanted.
        if let Err(e) = start(&mut md_sequence, &topo, &mut conf, &sim_state) {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    }
    for step in 0..n_steps {
        let time = step as f64 * dt;

        log::debug!("Step {}: time = {:.6} ps", step, time);

        // Run the algorithm sequence for this step. Step 0 goes through
        // `gromos_run::start` (init + initial force evaluation) — the same entry point the
        // Python binding uses, so both start a run identically.
        let step_result = if step == 0 {
            start(&mut md_sequence, &topo, &mut conf, &sim_state)
        } else {
            md_sequence
                .run_step(&topo, &mut conf, &sim_state)
                .map_err(|e| RunError::Step { step, message: e })
        };
        if let Err(e) = step_result {
            eprintln!("Error: {}", e);
            process::exit(1);
        }

        // Debug: dump forces at step 0
        if step == 0 && md_args.verbose >= 2 {
            log::debug!("=== Forces at step 0 (old state, after exchange) ===");
            for i in 0..topo.num_atoms() {
                let f = conf.old().force[i];
                log::debug!(
                    "  Atom {:2}: ({:18.9}, {:18.9}, {:18.9})",
                    i + 1,
                    f.x,
                    f.y,
                    f.z
                );
            }
        }

        // Apply GAMD boost if enabled
        if let Some(ref mut gamd) = gamd_params {
            let dihedral_energy = conf.current().energies.dihedral_total;
            let total_potential = conf.current().energies.potential_total;

            // Update GAMD statistics
            gamd.update_statistics(dihedral_energy, total_potential);

            // Update parameters if in GaMD search mode
            if gamd.search_mode == gromos::gamd::SearchMode::GamdSearch {
                gamd.update_all_parameters();
            }

            // Apply boost potential
            // For now, use simplified version with only total potential boost
            // TODO: Separate dihedral forces for dual boost
            let dihedral_forces = vec![Vec3::ZERO; topo.num_atoms()];
            let other_forces: Vec<Vec3> = conf.current().force.clone();

            let boost = gamd.apply_boost(
                &mut conf,
                dihedral_energy,
                &dihedral_forces,
                total_potential,
                &other_forces,
            );

            // Add boost to potential energy
            conf.current_mut().energies.potential_total += boost;

            log::debug!(
                "GAMD boost applied: ΔV={:.4}, V_new={:.4}",
                boost,
                conf.current().energies.potential_total
            );

            // Write GaMD statistics
            if let Some(ref mut writer) = gamd_stats_writer {
                if let Err(e) = writer.write_frame(step, dihedral_energy, total_potential, gamd) {
                    log::warn!("Failed to write GaMD stats at step {}: {}", step, e);
                }
            }

            // Write GaMD boost potential
            if let Some(ref mut writer) = gamd_boost_writer {
                // For now, assume boost is only from total potential
                let boost_dih = 0.0;
                let boost_pot = boost;
                if let Err(e) = writer.write_frame(step, time, boost, boost_dih, boost_pot) {
                    log::warn!("Failed to write GaMD boost at step {}: {}", step, e);
                }
            }
        }

        // Apply EDS if enabled
        if let Some(ref mut aeds) = eds_params {
            let current_potential = conf.current().energies.potential_total;

            // For EDS, we need to:
            // 1. Update each state's energy (for simplicity, using current potential as state 0)
            // 2. Calculate reference energy V_R
            // 3. Apply EDS forces

            // TODO: This is a simplified implementation
            // Full EDS requires multiple state energies calculated simultaneously
            // For now, we'll just set up the framework

            let eds = &mut aeds.eds;

            // Update state 0 with current potential
            if eds.num_states > 0 {
                eds.states[0].energy = current_potential;
            }

            // Calculate reference energy based on form
            match eds.form {
                gromos::eds::EDSForm::SingleS => eds.calculate_reference_energy_single_s(),
                gromos::eds::EDSForm::MultiS => eds.calculate_reference_energy_multi_s(),
                gromos::eds::EDSForm::PairS => {
                    log::warn!("PairS EDS form not yet fully implemented");
                    eds.calculate_reference_energy_single_s();
                },
            }

            // Apply EDS forces (modifies configuration)
            eds.apply_forces(&mut conf);

            // Replace potential energy with reference energy
            conf.current_mut().energies.potential_total = eds.reference_energy;

            log::debug!(
                "EDS applied: V_R={:.4}, original V={:.4}",
                eds.reference_energy,
                current_potential
            );

            // Update AEDS parameters if search enabled
            if aeds.search_enabled {
                // TODO: Implement AEDS parameter updates
                log::debug!("AEDS search mode active");
            }

            // Write EDS statistics
            if let Some(ref mut writer) = eds_stats_writer {
                if let Err(e) = writer.write_frame(step, eds) {
                    log::warn!("Failed to write EDS stats at step {}: {}", step, e);
                }
            }

            // Write EDS reference energy
            if let Some(ref mut writer) = eds_vr_writer {
                if let Err(e) = writer.write_frame(step, time, eds.reference_energy) {
                    log::warn!("Failed to write EDS V_R at step {}: {}", step, e);
                }
            }
        }

        // Validate energy
        // GROMOS convention: after the algorithm sequence, energies are in old()
        // (Forcefield wrote to current(), exchange_state moved it to old(),
        //  Temperature_Calculation and Energy_Calculation also write to old())
        let state = conf.old();
        let n_dof = topo.inverse_mass.len() * 3;
        let temp = state.temperature(n_dof);
        let ene_validation = validate_energy(
            state.energies.kinetic_total,
            state.energies.potential_total,
            state.energies.total(),
            temp,
        );
        if ene_validation.has_errors() && md_args.verbose > 0 {
            ene_validation.print();
            log::warn!("Energy validation failed at step {}", step);
        }

        // Store energy for drift check
        energy_history.push((
            state.energies.kinetic_total,
            state.energies.potential_total,
            state.energies.total(),
        ));

        actual_steps = step;

        // Track per-step energy change for minimization display
        let min_de = if is_minimization {
            let current_pot = state.energies.potential_total + state.energies.special_total;
            let de = if prev_min_energy < f64::MAX / 2.0 {
                current_pot - prev_min_energy
            } else {
                0.0
            };
            prev_min_energy = current_pot;
            de
        } else {
            0.0
        };

        // Minimization convergence check (mirrors GROMOS: compare potential_total + special_total)
        if is_minimization && step > imd.nmin && min_de.abs() < imd.dele {
            minimization_converged = true;
            println!();
            println!("*** Energy minimization CONVERGED at step {} ***", step);
            println!(
                "    dE = {:.6e} kJ/mol (tolerance: {:.6e})",
                min_de.abs(),
                imd.dele
            );
            println!("    E_pot = {:.10e} kJ/mol", state.energies.potential_total);
            println!();
            // Write final frames before breaking. gromosXX writes the converged energies
            // once more, labelled with the previous step (its frame count is NSTEPS + 1).
            if let Some(tw) = traj_writer.as_mut() {
                if let Err(e) = tw.write_frame(step, time, &conf) {
                    eprintln!("Error writing final trajectory frame: {}", e);
                }
            }
            if let Some(ew) = ene_writer.as_mut() {
                let last = step.saturating_sub(1);
                let vp = volume_pressure(conf.old(), &topo, temp);
                let e = ew.write_energy(last, last as f64 * dt, &conf.old().energies, &vp);
                if let Err(e) = e {
                    eprintln!("Error writing final energy frame: {}", e);
                }
            }
            break;
        }

        // Log progress
        if due(step, nstlog) {
            if is_minimization {
                println!(
                    "Step {:6}  E_pot: {:18.10e}  dE: {:12.4e}",
                    step, state.energies.potential_total, min_de
                );
            } else {
                println!("Step {:6}  Time: {:8.3} ps  E_pot: {:18.10e}  E_kin: {:18.10e}  E_tot: {:18.10e}  T: {:6.1} K",
                    step, time, state.energies.potential_total, state.energies.kinetic_total,
                    state.energies.total(), temp);
            }
            log::debug!(
                "  Bond: {:.10e}  LJ: {:.10e}  CRF: {:.10e}",
                state.energies.bond_total,
                state.energies.lj_total,
                state.energies.crf_total
            );
        }

        // Write trajectory — the configuration this step's energies describe, on gromosXX's
        // schedule (frame 0 is the input configuration).
        if let (true, Some(traj_writer)) = (due(step, nstxout) && step < n_steps, &mut traj_writer)
        {
            let write = if conf.current().box_config.box_type == BoxType::TruncatedOctahedral {
                // as for the final configuration: gromosXX writes the cube frame
                let rot = truncoct_triclinic_rotmat(false);
                let pos: Vec<Vec3> = conf.old().pos.iter().map(|p| rot * *p).collect();
                let dims = conf.old().box_config.dimensions();
                traj_writer.write_trc_frame(step, time, &pos, (dims.x > 0.0).then_some(dims))
            } else {
                traj_writer.write_frame_of_step(step, time, &conf)
            };
            if let Err(e) = write {
                eprintln!("Error writing trajectory: {}", e);
            }
        }

        // Write velocities (NTWV) — the state this step's energies describe, like the positions.
        if let Some(ref mut vw) = vel_writer {
            if due(step, imd.ntwv) && step < n_steps {
                let vel = &conf.old().vel;
                let dims = conf.old().box_config.dimensions();
                if let Err(e) =
                    vw.write_velocity_frame(step, time, vel, (dims.x > 0.0).then_some(dims))
                {
                    eprintln!("Error writing velocity trajectory: {}", e);
                }
            }
        }

        // Write energies — on gromosXX's schedule: steps 0 … NSTLIM−1 (the state after the
        // last step is in the final configuration, not in the energy trajectory).
        if let (true, Some(ew)) = (due(step, nstener) && step < n_steps, &mut ene_writer) {
            let vp = volume_pressure(conf.old(), &topo, temp);
            if let Err(e) = ew.write_energy(step, time, &conf.old().energies, &vp) {
                eprintln!("Error writing energy: {}", e);
            }
        }

        // Write free-energy trajectory (dH/dλ) at the WRITETRAJ NTWG frequency (gromosXX's
        // `write.free_energy`); NTWG=0 writes no frames.
        if let Some(ref mut fw) = free_energy_writer {
            if due(step, imd.ntwg) && step < n_steps {
                let d = &conf.old().perturbed_energy_derivatives;
                if let Err(e) = fw.write_derivatives(step, time, topo.lambda, d) {
                    eprintln!("Error writing free-energy trajectory: {}", e);
                }
            }
        }

        // Write forces (same frequency as trajectory)
        if let Some(ref mut fw) = force_writer {
            if due(step, nstxout) {
                let forces = &conf.old().force;
                let constraint_forces = if constraints_enabled {
                    Some(conf.old().constraint_force.as_slice())
                } else {
                    None
                };

                // GROMOS out_configuration.cc::_print_forcered rotates
                // FREEFORCERED/CONSFORCERED back into the original (cube)
                // frame via truncoct_triclinic_rotmat(false) when
                // boundary_type==truncoct (no ROTTRANS block here, so
                // rmat(phi,theta,psi)=identity). P1.4.
                if conf.current().box_config.box_type == BoxType::TruncatedOctahedral {
                    let rot = truncoct_triclinic_rotmat(false);
                    let rotated_forces: Vec<Vec3> = forces.iter().map(|f| rot * *f).collect();
                    let rotated_constraints: Option<Vec<Vec3>> =
                        constraint_forces.map(|cf| cf.iter().map(|f| rot * *f).collect());
                    if let Err(e) =
                        fw.write_frame(step, time, &rotated_forces, rotated_constraints.as_deref())
                    {
                        eprintln!("Error writing forces: {}", e);
                    }
                } else if let Err(e) = fw.write_frame(step, time, forces, constraint_forces) {
                    eprintln!("Error writing forces: {}", e);
                }
            }
        }

        // Advance simulation state for next step
        sim_state.advance();
    }

    if is_minimization {
        if minimization_converged {
            log::info!("Energy minimization converged at step {}", actual_steps);
        } else {
            log::warn!(
                "Energy minimization did NOT converge within {} steps",
                n_steps
            );
            println!();
            println!("*** WARNING: Energy minimization did NOT converge ***");
            println!("    Max steps reached: {}", n_steps);
            println!("    Consider increasing NSTLIM or adjusting DELE");
            println!();
        }
    } else {
        log::info!("MD loop completed - {} steps", n_steps);
    }

    // Check energy drift
    log::debug!("Checking energy drift over trajectory");
    use gromos::validation::check_energy_drift;
    let drift_report = check_energy_drift(&energy_history);
    if drift_report.has_errors() || !drift_report.warnings.is_empty() {
        drift_report.print();
        drift_report.print_summary();
    } else {
        log::debug!("Energy drift check passed");
    }

    // Finalize output files
    log::debug!("Finalizing output files");
    if let Some(tw) = traj_writer.as_mut() {
        if let Err(e) = tw.flush() {
            log::error!("Failed to flush trajectory: {}", e);
            eprintln!("Error flushing trajectory: {}", e);
        } else {
            log::debug!("Trajectory file finalized: {}", trc_file);
        }
    }

    if let Some(ew) = ene_writer.as_mut() {
        if let Err(e) = ew.finalize() {
            log::error!("Failed to finalize energy file: {}", e);
            eprintln!("Error finalizing energy file: {}", e);
        } else {
            log::debug!("Energy file finalized: {}", tre_file);
        }
    }

    if let Some(ref mut fw) = force_writer {
        if let Err(e) = fw.flush() {
            log::error!("Failed to flush force file: {}", e);
            eprintln!("Error flushing force file: {}", e);
        } else if let Some(ref trf_file) = md_args.trf_file {
            log::debug!("Force file finalized: {}", trf_file);
        }
    }

    // Finalize GaMD writers if enabled
    if let Some(ref mut writer) = gamd_stats_writer {
        if let Err(e) = writer.finalize() {
            log::error!("Failed to finalize GaMD stats file: {}", e);
            eprintln!("Error finalizing GaMD stats file: {}", e);
        } else {
            log::debug!("GaMD stats file finalized");
        }
    }

    if let Some(ref mut writer) = gamd_boost_writer {
        if let Err(e) = writer.finalize() {
            log::error!("Failed to finalize GaMD boost file: {}", e);
            eprintln!("Error finalizing GaMD boost file: {}", e);
        } else {
            log::debug!("GaMD boost file finalized");
        }
    }

    // Finalize EDS writers if enabled
    if let Some(ref mut writer) = eds_stats_writer {
        if let Err(e) = writer.finalize() {
            log::error!("Failed to finalize EDS stats file: {}", e);
            eprintln!("Error finalizing EDS stats file: {}", e);
        } else {
            log::debug!("EDS stats file finalized");
        }
    }

    if let Some(ref mut writer) = eds_vr_writer {
        if let Err(e) = writer.finalize() {
            log::error!("Failed to finalize EDS V_R file: {}", e);
            eprintln!("Error finalizing EDS V_R file: {}", e);
        } else {
            log::debug!("EDS V_R file finalized");
        }
    }

    // Finalize free-energy trajectory if enabled
    if let Some(ref mut fw) = free_energy_writer {
        if let Err(e) = fw.flush() {
            log::error!("Failed to flush free-energy trajectory: {}", e);
            eprintln!("Error flushing free-energy trajectory: {}", e);
        }
    }

    // Write final configuration if @fin was specified
    if let Some(ref fin_path) = md_args.fin_file {
        // NTB = −1: the engine works in the triclinic frame `truncoct_triclinic_box` produces;
        // gromosXX rotates positions, velocities and forces back into the cube frame when it
        // writes them (`out_configuration.cc`, `truncoct_triclinic_rotmat(false)`), so a written
        // configuration is in the same frame as the input.
        let truncoct = conf.current().box_config.box_type == BoxType::TruncatedOctahedral;
        let rotated_pos: Option<Vec<Vec3>> = truncoct.then(|| {
            let rot = truncoct_triclinic_rotmat(false);
            conf.current().pos.iter().map(|p| rot * *p).collect()
        });
        let rotated_vel: Option<Vec<Vec3>> = truncoct.then(|| {
            let rot = truncoct_triclinic_rotmat(false);
            conf.current().vel.iter().map(|v| rot * *v).collect()
        });
        let positions = rotated_pos.as_ref().unwrap_or(&conf.current().pos);
        let velocities = rotated_vel.as_ref().unwrap_or(&conf.current().vel);
        let box_vec = conf.current().box_config.dimensions();
        let vels = if velocities.iter().any(|v| *v != Vec3::ZERO) {
            Some(velocities.as_slice())
        } else {
            None
        };
        let box_opt = if box_vec.x > 0.0 { Some(box_vec) } else { None };
        let title = format!("Final configuration after {} steps", actual_steps);
        // gromosXX `out_configuration.cc:540-546`: a perturbed run (NTG=1) adds PERTDATA, the λ
        // a run with NRDGL=1 continues from; time-averaged distance restraints (NTDIR < 0) add
        // their ⟨r⁻³⟩ (DISRESEXPAVE).
        let disres_averages = if imd.ntdir < 0 {
            md_sequence
                .find_mut::<gromos_integrators::algorithms::Forcefield>()
                .map(|ff| ff.distance_restraints.averages())
        } else {
            None
        };
        let write = || -> Result<(), String> {
            let mut w = gromos::io::g96::G96Writer::new(fin_path)?;
            w.write_title(&title)?;
            w.write_positions(positions, Some(&topo))?;
            if let Some(v) = vels {
                w.write_velocities(v)?;
            }
            if let Some(b) = box_opt {
                w.write_box(b)?;
            }
            if imd.ntg != 0 {
                w.write_pertdata(topo.lambda)?;
            }
            if let Some(averages) = &disres_averages {
                w.write_disresexpave(averages)?;
            }
            w.close()
        };
        if let Err(e) = write() {
            eprintln!("Error writing final configuration: {}", e);
        } else {
            log::info!("Final configuration written to: {}", fin_path);
        }
    }

    let elapsed = start_time.elapsed();
    let total_elapsed = t_total.elapsed();

    if timing_enabled {
        let total = elapsed.as_secs_f64();
        let pct = |s: f64| if total > 0.0 { s / total * 100.0 } else { 0.0 };
        let detailed = md_sequence.detailed_timings();
        let accounted: f64 = detailed.iter().map(|(_, d, _)| d.as_secs_f64()).sum();

        println!();
        println!("TIMING");
        for (name, d, subs) in &detailed {
            let s = d.as_secs_f64();
            println!("{:>38} {:>12.3} {:>9.2}%", name, s, pct(s));
            let sub_total: f64 = subs.iter().map(|(_, sd)| sd.as_secs_f64()).sum();
            for (label, sd) in subs {
                let ss = sd.as_secs_f64();
                if ss > 0.0 {
                    println!(
                        "{:>44} {:>10.3} {:>9.2}%",
                        format!("- {label}"),
                        ss,
                        pct(ss)
                    );
                }
            }
            if !subs.is_empty() && s > sub_total {
                println!(
                    "{:>44} {:>10.3} {:>9.2}%",
                    "- unaccounted",
                    s - sub_total,
                    pct(s - sub_total)
                );
            }
        }
        println!(
            "{:>38} {:>12.3} {:>9.2}%",
            "unaccounted",
            total - accounted,
            pct(total - accounted)
        );
        println!("END");
    }

    log::info!("Simulation wall time: {:.2} s", elapsed.as_secs_f64());
    log::info!(
        "Total wall time: {:.3} s (init: {:.3} s, sim: {:.3} s)",
        total_elapsed.as_secs_f64(),
        init_elapsed.as_secs_f64(),
        elapsed.as_secs_f64()
    );

    println!();
    println!("╔══════════════════════════════════════════════════════════════╗");
    if is_minimization {
        println!("║              Energy Minimization Complete                    ║");
    } else {
        println!("║                   Simulation Complete                        ║");
    }
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();
    println!("Statistics:");
    println!("  Total steps:     {}", actual_steps);
    if is_minimization {
        println!(
            "  Converged:       {}",
            if minimization_converged { "YES" } else { "NO" }
        );
    } else {
        println!("  Simulation time: {:.3} ps", n_steps as f64 * dt);
    }
    println!("  Wall time:       {:.2} s", elapsed.as_secs_f64());
    println!(
        "  Performance:     {:.1} ns/day",
        (n_steps as f64 * dt * 1e-3) / elapsed.as_secs_f64() * 86400.0
    );
    println!();
    println!("Output files:");
    if let Some(ref fin_path) = md_args.fin_file {
        println!("  Final conf: {}", fin_path);
    }
    println!("  Trajectory: {}", trc_file);
    println!("  Energies:   {}", tre_file);
    println!();
    println!("Done!");
}

/// MPI support for `md` (feature `use-mpi`): pair decomposition across ranks, rank 0 owns the
/// output. See `Forcefield::pair_partition` for the decomposition and BENCHMARKING.md §6.
#[cfg(feature = "use-mpi")]
mod mpi_support {
    use gromos_core::algorithm::AlgorithmSequence;
    use gromos_forces::nonbonded::ForceStorage;
    use gromos_integrators::algorithms::Forcefield;
    use mpi::collective::SystemOperation;
    use mpi::topology::SimpleCommunicator;
    use mpi::traits::*;

    pub struct MpiRun {
        _universe: mpi::environment::Universe,
        pub rank: usize,
        pub size: usize,
    }

    /// Initialise MPI. Ranks other than 0 write no files and send their stdout to /dev/null;
    /// their stderr (the log) stays, prefixed by mpirun.
    pub fn setup(args: &mut super::MDArgs) -> MpiRun {
        let universe = mpi::initialize().expect("MPI_Init");
        let world = universe.world();
        let (rank, size) = (world.rank() as usize, world.size() as usize);
        if rank != 0 {
            for out in [
                &mut args.fin_file,
                &mut args.trc_file,
                &mut args.trv_file,
                &mut args.trf_file,
                &mut args.trs_file,
                &mut args.tre_file,
                &mut args.trg_file,
                &mut args.bae_file,
                &mut args.bag_file,
            ] {
                *out = None;
            }
            // SAFETY: replacing fd 1 by /dev/null; the descriptor stays open and valid.
            unsafe {
                let null = std::ffi::CString::new("/dev/null").unwrap();
                let fd = libc::open(null.as_ptr(), libc::O_WRONLY);
                if fd >= 0 {
                    libc::dup2(fd, 1);
                    libc::close(fd);
                }
            }
        } else {
            println!(
                "MPI: {size} rank(s); pair terms split by first atom index, summed every step"
            );
        }
        MpiRun {
            _universe: universe,
            rank,
            size,
        }
    }

    impl MpiRun {
        pub fn partition(&self) -> Option<(usize, usize)> {
            (self.size > 1).then_some((self.rank, self.size))
        }

        /// Sum the pair-term storage over all ranks: reduce to rank 0, broadcast the sum, so
        /// every rank continues with bit-identical forces and energies.
        pub fn install_reducer(&self, sequence: &mut AlgorithmSequence) {
            if self.size <= 1 {
                return;
            }
            let ff = sequence
                .find_mut::<Forcefield>()
                .expect("the sequence carries a Forcefield");
            ff.set_nonbonded_reducer(Box::new(|s: &mut ForceStorage| {
                let world = SimpleCommunicator::world();
                let root = world.process_at_rank(0);
                let n = s.forces.len();
                // with more than one energy group the four group-pair matrices travel too
                let g = &s.groups;
                let cells = g.lj.len();
                let mut send = Vec::with_capacity(3 * n + 13 + 4 * cells);
                for f in &s.forces {
                    send.extend_from_slice(&[f.x, f.y, f.z]);
                }
                send.push(s.e_lj);
                send.push(s.e_crf);
                for row in &s.virial {
                    send.extend_from_slice(row);
                }
                send.push(s.dhdl_lj);
                send.push(s.dhdl_crf);
                for m in [&g.lj, &g.crf, &g.dhdl_lj, &g.dhdl_crf] {
                    send.extend_from_slice(m);
                }
                // and the PRECALCLAM arrays of the pairs (empty when off)
                for m in s.precalc.energy.iter().chain(s.precalc.deriv.iter()) {
                    send.extend_from_slice(m);
                }
                let mut recv = vec![0.0f64; send.len()];
                if world.rank() == 0 {
                    root.reduce_into_root(&send[..], &mut recv[..], SystemOperation::sum());
                } else {
                    root.reduce_into(&send[..], SystemOperation::sum());
                }
                root.broadcast_into(&mut recv[..]);
                for (i, f) in s.forces.iter_mut().enumerate() {
                    f.x = recv[3 * i];
                    f.y = recv[3 * i + 1];
                    f.z = recv[3 * i + 2];
                }
                let mut k = 3 * n;
                s.e_lj = recv[k];
                s.e_crf = recv[k + 1];
                k += 2;
                for row in s.virial.iter_mut() {
                    for v in row.iter_mut() {
                        *v = recv[k];
                        k += 1;
                    }
                }
                s.dhdl_lj = recv[k];
                s.dhdl_crf = recv[k + 1];
                k += 2;
                let g = &mut s.groups;
                for m in [&mut g.lj, &mut g.crf, &mut g.dhdl_lj, &mut g.dhdl_crf] {
                    m.copy_from_slice(&recv[k..k + cells]);
                    k += cells;
                }
                let pc = &mut s.precalc;
                for m in pc.energy.iter_mut().chain(pc.deriv.iter_mut()) {
                    let len = m.len();
                    m.copy_from_slice(&recv[k..k + len]);
                    k += len;
                }
            }));
        }
    }
}
