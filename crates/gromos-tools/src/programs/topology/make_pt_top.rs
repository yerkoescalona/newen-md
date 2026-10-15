//! make_pt_top - Create perturbation topology from two molecular topologies
//!
//! Usage: make_pt_top @topo <top_a> <top_b> [@softpar <alpha_lj> <alpha_crf>] [@out <output.ptp>]
//!
//! Creates a perturbation topology file (.ptp) for free energy calculations
//! by comparing two molecular topologies (states A and B).

use crate::programs::io_err;
use crate::programs::topology::pt_top::net_charge_line;
use gromos::io::{
    topology::{build_topology, read_topology_file},
    PtpWriter,
};
use gromos_core::NetCharge;
use gromos_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "make_pt_top - Write a perturbation topology from two topologies\n\nUsage: make_pt_top @topo <fileA> <fileB> [@out <file>] [@softpar <alpha_lj> <alpha_crf>]\n\nArguments:\n  @topo     State A and state B topologies\n  @out      Output perturbation topology (default perturbation.ptp; written as a file)\n  @softpar  Soft-core parameters (two values)";

pub const KNOWN: &[&str] = &["topo", "out", "softpar"];

pub const PROGRAM: Program = Program {
    name: "make_pt_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

#[allow(dead_code)]
fn print_usage() {
    eprintln!("make_pt_top - Create perturbation topology for FEP calculations");
    eprintln!();
    eprintln!("Usage: make_pt_top @topo <topology_a> <topology_b> [@softpar <alpha_lj> <alpha_crf>] [@out <output.ptp>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo       Two topology files: state A (λ=0) and state B (λ=1)");
    eprintln!("  @softpar    Softness parameters (default: 1.51 0.5)");
    eprintln!("              alpha_lj:  Lennard-Jones softness (typ. 1.51)");
    eprintln!("              alpha_crf: Coulomb reaction field softness (typ. 0.5)");
    eprintln!("  @out        Output .ptp file (default: perturbation.ptp)");
    eprintln!();
    eprintln!("Description:");
    eprintln!("  Compares two molecular topologies and generates a perturbation topology");
    eprintln!("  file for free energy perturbation (FEP) calculations. Both topologies");
    eprintln!("  must have the same number of atoms.");
    eprintln!();
    eprintln!("  The tool identifies atoms with different:");
    eprintln!("  - Integer atom codes (IAC / atom types)");
    eprintln!("  - Masses");
    eprintln!("  - Charges");
    eprintln!();
    eprintln!("  Softness parameters control the soft-core potential:");
    eprintln!("  - alpha_lj:  ~1.5 for smooth LJ transitions");
    eprintln!("  - alpha_crf: ~0.5 for smooth electrostatic transitions");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  # Create perturbation topology with default softness");
    eprintln!("  make_pt_top @topo stateA.top stateB.top");
    eprintln!();
    eprintln!("  # Specify softness parameters");
    eprintln!("  make_pt_top @topo stateA.top stateB.top @softpar 1.51 0.5");
    eprintln!();
    eprintln!("  # Custom output file");
    eprintln!("  make_pt_top @topo stateA.top stateB.top @out custom.ptp");
    eprintln!();
    eprintln!("Output:");
    eprintln!("  A .ptp file containing:");
    eprintln!("  - PERTURBEDATOM block: List of atoms that change between states");
    eprintln!("  - PERTURBATIONPARAMETERS block: Softness parameters");
    eprintln!();
    eprintln!("Reference:");
    eprintln!("  GROMOS manual section on free energy perturbation");
    eprintln!("  Liu & Berne, J. Chem. Phys. 118, 2977 (2003) - soft-core potentials");
}

#[derive(Debug)]
struct MakePtTopArgs {
    topo_a: String,
    topo_b: String,
    alpha_lj: f64,
    alpha_crf: f64,
    output: String,
}

impl Default for MakePtTopArgs {
    fn default() -> Self {
        Self {
            topo_a: String::new(),
            topo_b: String::new(),
            alpha_lj: 1.51,
            alpha_crf: 0.5,
            output: "perturbation.ptp".to_string(),
        }
    }
}

fn parse_args(args: &Arguments) -> Result<MakePtTopArgs, String> {
    let mut pt_args = MakePtTopArgs::default();

    let topo_files = args.values("topo");
    if topo_files.len() != 2 {
        return Err("Must provide exactly two topology files with @topo".to_string());
    }
    pt_args.topo_a = topo_files[0].clone();
    pt_args.topo_b = topo_files[1].clone();

    if args.has("softpar") {
        let v = args.values("softpar");
        if v.len() != 2 {
            return Err("@softpar takes alpha_lj and alpha_crf".to_string());
        }
        pt_args.alpha_lj = v[0]
            .parse()
            .map_err(|_| format!("Invalid alpha_lj: {}", v[0]))?;
        pt_args.alpha_crf = v[1]
            .parse()
            .map_err(|_| format!("Invalid alpha_crf: {}", v[1]))?;
    }
    if args.has("out") {
        pt_args.output = args.value("out")?.to_string();
    }

    Ok(pt_args)
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let pt_args = parse_args(&args).map_err(|e| format!("Error: {e}\n\n{USAGE}"))?;

    writeln!(out, "make_pt_top - Creating perturbation topology").map_err(io_err)?;
    writeln!(out, "  State A:    {}", pt_args.topo_a).map_err(io_err)?;
    writeln!(out, "  State B:    {}", pt_args.topo_b).map_err(io_err)?;
    writeln!(out, "  alpha_LJ:   {}", pt_args.alpha_lj).map_err(io_err)?;
    writeln!(out, "  alpha_CRF:  {}", pt_args.alpha_crf).map_err(io_err)?;
    writeln!(out, "  Output:     {}", pt_args.output).map_err(io_err)?;
    writeln!(out).map_err(io_err)?;

    // Read topology A
    writeln!(out, "Reading state A topology...").map_err(io_err)?;
    let blocks_a = read_topology_file(&pt_args.topo_a)
        .map_err(|e| format!("Error reading topology A: {e}"))?;

    let topo_a = build_topology(blocks_a);
    writeln!(out, "  State A: {} atoms", topo_a.num_atoms()).map_err(io_err)?;

    // Read topology B
    writeln!(out, "Reading state B topology...").map_err(io_err)?;
    let blocks_b = read_topology_file(&pt_args.topo_b)
        .map_err(|e| format!("Error reading topology B: {e}"))?;

    let topo_b = build_topology(blocks_b);
    writeln!(out, "  State B: {} atoms", topo_b.num_atoms()).map_err(io_err)?;

    // Check compatibility
    if topo_a.num_atoms() != topo_b.num_atoms() {
        return Err(format!(
            "Error: Topologies must have the same number of atoms ({} vs {})",
            topo_a.num_atoms(),
            topo_b.num_atoms()
        ));
    }

    let perturbed: Vec<bool> = (0..topo_a.num_atoms())
        .map(|i| {
            topo_a.iac[i] != topo_b.iac[i]
                || (topo_a.mass[i] - topo_b.mass[i]).abs() > 1e-6
                || (topo_a.charge[i] - topo_b.charge[i]).abs() > 1e-6
        })
        .collect();
    let n_perturbed = perturbed.iter().filter(|&&p| p).count();

    writeln!(out).map_err(io_err)?;
    writeln!(out, "Perturbation analysis:").map_err(io_err)?;
    writeln!(out, "  Total atoms:      {}", topo_a.num_atoms()).map_err(io_err)?;
    writeln!(out, "  Perturbed atoms:  {}", n_perturbed).map_err(io_err)?;
    for (m, q) in NetCharge::of_perturbed_molecules(
        &topo_a.molecules,
        &topo_a.charge,
        &topo_b.charge,
        &perturbed,
    ) {
        writeln!(out, "  {}", net_charge_line(m, &q)).map_err(io_err)?;
    }
    writeln!(out).map_err(io_err)?;

    if n_perturbed == 0 {
        eprintln!("Warning: No perturbed atoms found! States A and B appear identical.");
    }

    // Write perturbation topology
    writeln!(
        out,
        "Writing perturbation topology to {}...",
        pt_args.output
    )
    .map_err(io_err)?;
    let mut writer =
        PtpWriter::new(&pt_args.output).map_err(|e| format!("Error creating output file: {e}"))?;

    let title = format!(
        "make_pt_top generated perturbation topology\nState A: {}\nState B: {}",
        pt_args.topo_a, pt_args.topo_b
    );

    writer
        .write(
            &topo_a,
            &topo_b,
            pt_args.alpha_lj,
            pt_args.alpha_crf,
            &title,
        )
        .map_err(|e| format!("Error writing perturbation topology: {e}"))?;

    writeln!(
        out,
        "Done! Perturbation topology written to {}",
        pt_args.output
    )
    .map_err(io_err)?;
    writeln!(out).map_err(io_err)?;
    writeln!(out, "Next steps:").map_err(io_err)?;
    writeln!(
        out,
        "  1. Use this .ptp file with the md binary: md @topo ... @ptp {}",
        pt_args.output
    )
    .map_err(io_err)?;
    writeln!(
        out,
        "  2. Run FEP simulations at different λ values (0.0 to 1.0)"
    )
    .map_err(io_err)?;
    writeln!(out, "  3. Analyze with GROMOS++ bar or ext_ti_ana tools").map_err(io_err)?;
    Ok(())
}
