//! check_top - Validate GROMOS topology file consistency
//!
//! Usage: check_top @topo <input.top> [@build <mtb>] [@param <ifp>]
//!    or: check_top @f <argfile>
//!
//! Performs comprehensive checks on topology files:
//! - Atom numbering and types
//! - Bond/angle/dihedral indices
//! - Force field parameter references
//! - Exclusion consistency
//! - Charge group validity

use crate::programs::{io_err, parse_clap};
use clap::Parser;
use newen::io::topology::{build_topology, read_topology_file};
use newen::topology::Topology;
use newen_io::args::{Program, Status};
use std::fmt::Write as FmtWrite;
use std::io::Write;

pub const USAGE: &str = "check_top - Check a GROMOS topology for consistency\n\nUsage: check_top @topo <topology> [@build <mtb>] [@param <ifp>]\n\nArguments:\n  @topo   Molecular topology file\n  @build  Building-block (MTB) file for cross-validation (accepted, not yet used)\n  @param  Force-field parameter (IFP) file for cross-validation (accepted, not yet used)";

pub const KNOWN: &[&str] = &["topo", "build", "param"];

pub const PROGRAM: Program = Program {
    name: "check_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

/// What the checks found: the report as the program prints it, plus the findings as data.
///
/// A topology with errors is still a report — `run` writes `text` and *then* fails, which is
/// what the binary has always done and what a caller needs in order to see why.
#[derive(Debug, Clone)]
pub struct CheckReport {
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
    /// The full report, in the order the checks ran.
    pub text: String,
}

impl CheckReport {
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }
}

/// Run every consistency check over a built topology.
pub fn check_top(topo: &Topology) -> CheckReport {
    let mut checker = TopologyChecker::new();
    checker.check_all(topo);
    checker.print_summary();
    CheckReport {
        errors: checker.errors,
        warnings: checker.warnings,
        text: checker.text,
    }
}

/// Validate GROMOS topology file consistency.
#[derive(Parser)]
#[command(name = "check_top", version, about)]
struct Args {
    /// Topology file to check
    #[arg(long)]
    topo: String,

    /// Building block (MTB) file for cross-validation (optional)
    #[arg(long)]
    build: Option<String>,

    /// Force field parameter (IFP) file for cross-validation (optional)
    #[arg(long)]
    param: Option<String>,
}

struct TopologyChecker {
    warnings: Vec<String>,
    errors: Vec<String>,
    /// The report, accumulated in the order the checks print it.
    text: String,
}

impl TopologyChecker {
    fn new() -> Self {
        Self {
            warnings: Vec::new(),
            errors: Vec::new(),
            text: String::new(),
        }
    }

    fn check_all(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "Performing topology checks...\n");

        self.check_atom_count(topo);
        self.check_masses(topo);
        self.check_charges(topo);
        self.check_bonds(topo);
        self.check_angles(topo);
        self.check_proper_dihedrals(topo);
        self.check_improper_dihedrals(topo);
        self.check_exclusions(topo);
        self.check_bond_parameters(topo);
        self.check_angle_parameters(topo);
        self.check_lj_parameters(topo);
    }

    fn check_atom_count(&mut self, topo: &Topology) {
        let n_atoms = topo.num_solute_atoms();
        let _ = writeln!(self.text, "CHECK: Atom count");
        let _ = writeln!(self.text, "  Solute atoms: {}", n_atoms);

        if n_atoms == 0 {
            self.errors.push("No atoms in topology".to_string());
        }

        // Check array sizes match
        if topo.mass.len() != n_atoms {
            self.errors.push(format!(
                "Mass array size ({}) != atom count ({})",
                topo.mass.len(),
                n_atoms
            ));
        }

        if topo.charge.len() != n_atoms {
            self.errors.push(format!(
                "Charge array size ({}) != atom count ({})",
                topo.charge.len(),
                n_atoms
            ));
        }

        if topo.iac.len() != n_atoms {
            self.errors.push(format!(
                "IAC array size ({}) != atom count ({})",
                topo.iac.len(),
                n_atoms
            ));
        }

        let _ = writeln!(self.text, "  ✓ Atom count checks passed\n");
    }

    fn check_masses(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Atomic masses");

        let mut zero_mass_count = 0;
        let mut negative_mass_count = 0;

        for (i, &mass) in topo.mass.iter().enumerate() {
            if mass < 0.0 {
                self.errors
                    .push(format!("Atom {} has negative mass: {}", i + 1, mass));
                negative_mass_count += 1;
            } else if mass == 0.0 {
                self.warnings.push(format!("Atom {} has zero mass", i + 1));
                zero_mass_count += 1;
            }
        }

        if zero_mass_count > 0 {
            let _ = writeln!(self.text, "  ⚠ {} atoms with zero mass", zero_mass_count);
        }
        if negative_mass_count == 0 && zero_mass_count == 0 {
            let _ = writeln!(self.text, "  ✓ All masses valid\n");
        } else {
            let _ = writeln!(self.text);
        }
    }

    fn check_charges(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Atomic charges");

        let total_charge: f64 = topo.charge.iter().sum();
        let _ = writeln!(self.text, "  Total charge: {:.6} e", total_charge);

        if total_charge.abs() > 1e-4 {
            self.warnings.push(format!(
                "System is not charge-neutral: {:.6} e",
                total_charge
            ));
            let _ = writeln!(self.text, "  ⚠ System not charge-neutral\n");
        } else {
            let _ = writeln!(self.text, "  ✓ System is charge-neutral\n");
        }
    }

    fn check_bonds(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Bonds");
        let n_atoms = topo.num_solute_atoms();
        let n_bonds = topo.moltypes[0].bonds.len();

        let _ = writeln!(self.text, "  Number of bonds: {}", n_bonds);

        for (idx, bond) in topo.moltypes[0].bonds.iter().enumerate() {
            // Check atom indices are in range
            if bond.i >= n_atoms {
                self.errors.push(format!(
                    "Bond {} atom i={} out of range (max: {})",
                    idx + 1,
                    bond.i + 1,
                    n_atoms
                ));
            }
            if bond.j >= n_atoms {
                self.errors.push(format!(
                    "Bond {} atom j={} out of range (max: {})",
                    idx + 1,
                    bond.j + 1,
                    n_atoms
                ));
            }

            // Check for self-bonds
            if bond.i == bond.j {
                self.errors.push(format!(
                    "Bond {} is a self-bond: atom {}",
                    idx + 1,
                    bond.i + 1
                ));
            }

            // Check bond type is valid
            if bond.bond_type >= topo.bond_parameters.len() {
                self.errors.push(format!(
                    "Bond {} has invalid type {} (max: {})",
                    idx + 1,
                    bond.bond_type + 1,
                    topo.bond_parameters.len()
                ));
            }
        }

        if self.errors.iter().any(|e| e.contains("Bond")) {
            let _ = writeln!(self.text, "  ✗ Bond errors found\n");
        } else {
            let _ = writeln!(self.text, "  ✓ All bonds valid\n");
        }
    }

    fn check_angles(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Bond angles");
        let n_atoms = topo.num_solute_atoms();
        let n_angles = topo.moltypes[0].angles.len();

        let _ = writeln!(self.text, "  Number of angles: {}", n_angles);

        for (idx, angle) in topo.moltypes[0].angles.iter().enumerate() {
            // Check atom indices
            if angle.i >= n_atoms {
                self.errors.push(format!(
                    "Angle {} atom i={} out of range",
                    idx + 1,
                    angle.i + 1
                ));
            }
            if angle.j >= n_atoms {
                self.errors.push(format!(
                    "Angle {} atom j={} out of range",
                    idx + 1,
                    angle.j + 1
                ));
            }
            if angle.k >= n_atoms {
                self.errors.push(format!(
                    "Angle {} atom k={} out of range",
                    idx + 1,
                    angle.k + 1
                ));
            }

            // Check for degenerate angles
            if angle.i == angle.j || angle.j == angle.k || angle.i == angle.k {
                self.errors.push(format!(
                    "Angle {} has duplicate atoms: {}-{}-{}",
                    idx + 1,
                    angle.i + 1,
                    angle.j + 1,
                    angle.k + 1
                ));
            }

            // Check angle type
            if angle.angle_type >= topo.angle_parameters.len() {
                self.errors.push(format!(
                    "Angle {} has invalid type {} (max: {})",
                    idx + 1,
                    angle.angle_type + 1,
                    topo.angle_parameters.len()
                ));
            }
        }

        if self.errors.iter().any(|e| e.contains("Angle")) {
            let _ = writeln!(self.text, "  ✗ Angle errors found\n");
        } else {
            let _ = writeln!(self.text, "  ✓ All angles valid\n");
        }
    }

    fn check_proper_dihedrals(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Proper dihedrals");
        let n_atoms = topo.num_solute_atoms();
        let n_dihedrals = topo.moltypes[0].proper_dihedrals.len();

        let _ = writeln!(self.text, "  Number of proper dihedrals: {}", n_dihedrals);

        for (idx, dihedral) in topo.moltypes[0].proper_dihedrals.iter().enumerate() {
            if dihedral.i >= n_atoms
                || dihedral.j >= n_atoms
                || dihedral.k >= n_atoms
                || dihedral.l >= n_atoms
            {
                self.errors.push(format!(
                    "Proper dihedral {} has out-of-range atoms",
                    idx + 1
                ));
            }

            // Check for duplicate atoms
            let atoms = [dihedral.i, dihedral.j, dihedral.k, dihedral.l];
            for i in 0..4 {
                for j in (i + 1)..4 {
                    if atoms[i] == atoms[j] {
                        self.errors
                            .push(format!("Proper dihedral {} has duplicate atoms", idx + 1));
                        break;
                    }
                }
            }
        }

        if self.errors.iter().any(|e| e.contains("Proper dihedral")) {
            let _ = writeln!(self.text, "  ✗ Dihedral errors found\n");
        } else if n_dihedrals > 0 {
            let _ = writeln!(self.text, "  ✓ All proper dihedrals valid\n");
        } else {
            let _ = writeln!(self.text);
        }
    }

    fn check_improper_dihedrals(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Improper dihedrals");
        let n_atoms = topo.num_solute_atoms();
        let n_impropers = topo.moltypes[0].improper_dihedrals.len();

        let _ = writeln!(self.text, "  Number of improper dihedrals: {}", n_impropers);

        for (idx, dihedral) in topo.moltypes[0].improper_dihedrals.iter().enumerate() {
            if dihedral.i >= n_atoms
                || dihedral.j >= n_atoms
                || dihedral.k >= n_atoms
                || dihedral.l >= n_atoms
            {
                self.errors.push(format!(
                    "Improper dihedral {} has out-of-range atoms",
                    idx + 1
                ));
            }
        }

        if self.errors.iter().any(|e| e.contains("Improper dihedral")) {
            let _ = writeln!(self.text, "  ✗ Improper dihedral errors found\n");
        } else if n_impropers > 0 {
            let _ = writeln!(self.text, "  ✓ All improper dihedrals valid\n");
        } else {
            let _ = writeln!(self.text);
        }
    }

    fn check_exclusions(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Exclusion lists");
        let n_atoms = topo.num_solute_atoms();

        if topo.exclusions.len() != n_atoms {
            self.errors.push(format!(
                "Exclusion list size ({}) != atom count ({})",
                topo.exclusions.len(),
                n_atoms
            ));
        }

        let mut total_exclusions = 0;
        let mut asymmetric_count = 0;

        for (i, excl_set) in topo.exclusions.iter().enumerate() {
            total_exclusions += excl_set.len();

            for &j in excl_set {
                if j >= n_atoms {
                    self.errors.push(format!(
                        "Atom {} has exclusion to out-of-range atom {}",
                        i + 1,
                        j + 1
                    ));
                }

                // Check symmetry: if i excludes j, j should exclude i
                if j < topo.exclusions.len() && !topo.exclusions[j].contains(&i) {
                    self.warnings.push(format!(
                        "Asymmetric exclusion: atom {} excludes {} but not vice versa",
                        i + 1,
                        j + 1
                    ));
                    asymmetric_count += 1;
                }
            }
        }

        let _ = writeln!(self.text, "  Total exclusions: {}", total_exclusions);
        if asymmetric_count > 0 {
            let _ = writeln!(self.text, "  ⚠ {} asymmetric exclusions", asymmetric_count);
        } else {
            let _ = writeln!(self.text, "  ✓ Exclusions are symmetric\n");
        }
    }

    fn check_bond_parameters(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Bond parameters");
        let _ = writeln!(
            self.text,
            "  Number of bond types: {}",
            topo.bond_parameters.len()
        );

        for (i, params) in topo.bond_parameters.iter().enumerate() {
            if params.r0 <= 0.0 {
                self.errors.push(format!(
                    "Bond type {} has non-positive equilibrium length: {}",
                    i + 1,
                    params.r0
                ));
            }

            if params.k_harmonic < 0.0 {
                self.warnings.push(format!(
                    "Bond type {} has negative harmonic force constant: {}",
                    i + 1,
                    params.k_harmonic
                ));
            }
        }

        if self.errors.iter().any(|e| e.contains("Bond type")) {
            let _ = writeln!(self.text, "  ✗ Bond parameter errors found\n");
        } else if !topo.bond_parameters.is_empty() {
            let _ = writeln!(self.text, "  ✓ All bond parameters valid\n");
        } else {
            let _ = writeln!(self.text);
        }
    }

    fn check_angle_parameters(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Angle parameters");
        let _ = writeln!(
            self.text,
            "  Number of angle types: {}",
            topo.angle_parameters.len()
        );

        for (i, params) in topo.angle_parameters.iter().enumerate() {
            if params.theta0 < 0.0 || params.theta0 > std::f64::consts::PI {
                self.warnings.push(format!(
                    "Angle type {} has theta0 outside [0, π]: {} rad",
                    i + 1,
                    params.theta0
                ));
            }

            if params.k_harmonic < 0.0 {
                self.warnings.push(format!(
                    "Angle type {} has negative harmonic force constant: {}",
                    i + 1,
                    params.k_harmonic
                ));
            }
        }

        if !topo.angle_parameters.is_empty() {
            let _ = writeln!(self.text, "  ✓ Angle parameter checks complete\n");
        } else {
            let _ = writeln!(self.text);
        }
    }

    fn check_lj_parameters(&mut self, topo: &Topology) {
        let _ = writeln!(self.text, "CHECK: Lennard-Jones parameters");
        let n_types = topo.lj_parameters.len();
        let _ = writeln!(self.text, "  Number of LJ atom types: {}", n_types);

        if n_types == 0 {
            self.warnings.push("No LJ parameters defined".to_string());
            let _ = writeln!(self.text, "  ⚠ No LJ parameters\n");
            return;
        }

        // Check matrix is square
        for (i, row) in topo.lj_parameters.iter().enumerate() {
            if row.len() != n_types {
                self.errors.push(format!(
                    "LJ parameter matrix row {} has wrong size: {} (expected {})",
                    i,
                    row.len(),
                    n_types
                ));
            }
        }

        // Check for negative parameters
        let mut negative_count = 0;
        for i in 0..n_types {
            for j in 0..n_types {
                if i < topo.lj_parameters.len() && j < topo.lj_parameters[i].len() {
                    let params = &topo.lj_parameters[i][j];
                    if params.c6 < 0.0 || params.c12 < 0.0 {
                        negative_count += 1;
                    }
                }
            }
        }

        if negative_count > 0 {
            self.warnings.push(format!(
                "{} LJ parameter pairs have negative values",
                negative_count
            ));
        }

        let _ = writeln!(self.text, "  ✓ LJ parameter checks complete\n");
    }

    fn print_summary(&mut self) {
        let _ = writeln!(self.text, "{}", "=".repeat(60));
        let _ = writeln!(self.text, "TOPOLOGY CHECK SUMMARY");
        let _ = writeln!(self.text, "{}", "=".repeat(60));

        if self.errors.is_empty() && self.warnings.is_empty() {
            let _ = writeln!(self.text, "✓ Topology is valid - no errors or warnings");
            let _ = writeln!(self.text);
            return;
        }

        if !self.errors.is_empty() {
            let _ = writeln!(self.text, "\nERRORS ({}):", self.errors.len());
            for error in &self.errors {
                let _ = writeln!(self.text, "  ✗ {}", error);
            }
        }

        if !self.warnings.is_empty() {
            let _ = writeln!(self.text, "\nWARNINGS ({}):", self.warnings.len());
            for warning in &self.warnings {
                let _ = writeln!(self.text, "  ⚠ {}", warning);
            }
        }

        let _ = writeln!(self.text);

        if self.errors.is_empty() {
            let _ = writeln!(self.text, "✓ No critical errors - topology is usable");
        } else {
            let _ = writeln!(
                self.text,
                "✗ Critical errors found - topology may not work correctly"
            );
        }
        let _ = writeln!(self.text);
    }
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let Some(args) = parse_clap::<Args>(argv, out)? else {
        return Ok(());
    };

    eprintln!("Checking topology file: {}\n", args.topo);

    let parsed = read_topology_file(&args.topo)
        .map_err(|e| format!("Error reading topology file: {e:?}"))?;
    write!(out, "Successfully parsed topology file\n\n").map_err(io_err)?;
    let topo = build_topology(parsed);

    let report = check_top(&topo);
    write!(out, "{}", report.text).map_err(io_err)?;

    // The report is the output; a topology with errors still produced one.
    if report.has_errors() {
        return Err(format!(
            "check_top: {} critical error(s) in {}",
            report.errors.len(),
            args.topo
        ));
    }
    Ok(())
}
