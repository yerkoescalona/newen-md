//! Perturbation topology (.pttopo) reader and writer.
//!
//! The GROMOS `.pttopo` format defines dual-topology (state A / state B)
//! parameters for atoms that change as λ goes from 0 to 1 (FEP / TI).
//!
//! Key blocks parsed by the reader:
//! - `PERTATOMPARAM`  — per-atom IAC, mass, charge for states A and B; soft-core alphas
//! - `PERTATOMPAIR`   — special LJ between two perturbed atoms (excluded from pairlist)
//! - `PERTBONDSTRETCH` / `PERTBONDSTRETCHH` — perturbed bond types
//! - `PERTBONDANGLE`  / `PERTBONDANGLEH`   — perturbed angle types
//! - `PERTIMPROPERDIH` / `PERTIMPROPERDIHH` — perturbed improper dihedral types
//! - `PERTPROPERDIH`  / `PERTPROPERDIHH`   — perturbed proper dihedral types
//!
//! The soft-variant blocks (PERTBONDSOFT, PERTANGLESOFT, PERTIMPROPERDIHSOFT) are
//! skipped for now — they describe bonds/angles/impropers that are absent in one
//! state (force constant = 0 in state A or B) and require soft-core treatment.

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use gromos_core::topology::{
    PerturbedAngle, PerturbedAtom, PerturbedAtomPair, PerturbedBond, PerturbedDihedral,
    PerturbedSolute, SoftAngle, SoftBond, SoftImproper,
};

use crate::IoError;

/// Parsed perturbation topology — pass to `apply_to_topology()` to wire into
/// a `Topology` that has already been built from the regular `.topo` file.
#[derive(Debug, Default)]
pub struct PerturbedTopology {
    pub perturbed_solute: PerturbedSolute,
    /// is_perturbed[atom_seq] = true for atoms in PERTATOMPARAM (0-indexed)
    pub is_perturbed: Vec<bool>,
    /// n_solute_atoms from the PERTATOMPARAM header (used to size is_perturbed)
    pub n_atoms_hint: usize,
}

/// Read a GROMOS perturbation topology file (`.pttopo`).
///
/// Atoms are 1-indexed in the file; returned structs use 0-indexed atoms.
/// IAC values are 1-indexed in the file; returned structs use 0-indexed IAC.
pub fn read_pttopo<P: AsRef<Path>>(path: P) -> Result<PerturbedTopology, IoError> {
    let content = fs::read_to_string(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;

    // Split into named blocks: block_name → Vec<non-comment lines>
    let blocks = parse_blocks(&content);

    let mut result = PerturbedTopology::default();
    let ps = &mut result.perturbed_solute;

    // ── PERTATOMPAIR ─────────────────────────────────────────────────────────
    // Must be read before PERTATOMPARAM (GROMOS convention — exclusion fixup).
    // Format: count \n i j A_type B_type \n ...   (all 1-indexed)
    for bname in &["PERTATOMPAIR", "PERTATOMPAIR03"] {
        if let Some(lines) = blocks.get(*bname) {
            let (count, data) = split_count(lines)?;
            for line in &data {
                let v = tokenize(line);
                if v.len() < 4 {
                    continue;
                }
                let i = parse_usize(v[0])? - 1; // 1-indexed in file → 0-indexed
                let j = parse_usize(v[1])? - 1;
                // Interaction type: 1-indexed in file → 0-indexed (file 1→0: full LJ, file 2→1: 1-4 LJ)
                let a_type = parse_usize(v[2])? - 1;
                // b_type=0 in file is a sentinel: "atom pair absent in state B" → None
                let b_type = {
                    let raw = parse_usize(v[3])?;
                    if raw == 0 {
                        None
                    } else {
                        Some(raw - 1)
                    }
                };
                let (i, j) = if i <= j { (i, j) } else { (j, i) };
                ps.atom_pairs.push(PerturbedAtomPair {
                    i,
                    j,
                    a_type,
                    b_type,
                });
            }
            if ps.atom_pairs.len() != count {
                log::warn!(
                    "PERTATOMPAIR: expected {count} pairs, got {}",
                    ps.atom_pairs.len()
                );
            }
            break;
        }
    }

    // ── PERTATOMPARAM ────────────────────────────────────────────────────────
    // Format: count \n seq res name a_iac a_mass a_charge b_iac b_mass b_charge lj_soft crf_soft
    // All IAC values 1-indexed; seq is 1-indexed atom number.
    for bname in &["PERTATOMPARAM", "PERTATOM03"] {
        if let Some(lines) = blocks.get(*bname) {
            let (count, data) = split_count(lines)?;
            let mut max_seq = 0usize;
            for line in &data {
                let v = tokenize(line);
                if v.len() < 11 {
                    continue;
                }
                let seq = parse_usize(v[0])? - 1; // 0-indexed
                                                  // v[1] = residue, v[2] = name — skip
                let a_iac = parse_usize(v[3])? - 1; // 0-indexed
                let a_mass = parse_f64(v[4])?;
                let a_charge = parse_f64(v[5])?;
                let b_iac = parse_usize(v[6])? - 1;
                let b_mass = parse_f64(v[7])?;
                let b_charge = parse_f64(v[8])?;
                let lj_soft = parse_f64(v[9])?;
                let crf_soft = parse_f64(v[10])?;
                max_seq = max_seq.max(seq);
                ps.atoms.push(PerturbedAtom {
                    seq,
                    a_iac,
                    a_mass,
                    a_charge,
                    b_iac,
                    b_mass,
                    b_charge,
                    lj_soft,
                    crf_soft,
                });
            }
            if ps.atoms.len() != count {
                log::warn!(
                    "PERTATOMPARAM: expected {count} atoms, got {}",
                    ps.atoms.len()
                );
            }
            // Build is_perturbed: size = max_seq + 1 (will be resized by caller if needed)
            let n = max_seq + 1;
            result.n_atoms_hint = n;
            result.is_perturbed = vec![false; n];
            for a in &ps.atoms {
                if a.seq < result.is_perturbed.len() {
                    result.is_perturbed[a.seq] = true;
                }
            }
            break;
        }
    }

    // ── PERTBONDSTRETCHH / PERTBONDSTRETCH ───────────────────────────────────
    // Format: count \n i j A_type B_type   (1-indexed). The H block first, as gromosXX reads
    // them (`in_perturbation.cc:172-175`): with NTC = 3 this is the perturbed SHAKE order.
    for bname in &["PERTBONDSTRETCHH", "PERTBONDSTRETCH", "PERTBOND03"] {
        if let Some(lines) = blocks.get(*bname) {
            let (count, data) = split_count(lines)?;
            let before = ps.bonds.len();
            for line in &data {
                let v = tokenize(line);
                if v.len() < 4 {
                    continue;
                }
                let i = parse_usize(v[0])? - 1;
                let j = parse_usize(v[1])? - 1;
                let a_type = parse_usize(v[2])? - 1; // file is 1-indexed
                let b_type = parse_usize(v[3])? - 1;
                ps.bonds.push(PerturbedBond {
                    i,
                    j,
                    a_type,
                    b_type,
                    hydrogen: *bname == "PERTBONDSTRETCHH",
                });
            }
            if ps.bonds.len() - before != count {
                log::warn!(
                    "{bname}: expected {count} bonds, got {}",
                    ps.bonds.len() - before
                );
            }
        }
    }

    // ── PERTBONDANGLEH / PERTBONDANGLE ───────────────────────────────────────
    // Format: count \n i j k A_type B_type (H block first, `in_perturbation.cc:590-591`)
    for bname in &["PERTBONDANGLEH", "PERTBONDANGLE", "PERTBANGLE03"] {
        if let Some(lines) = blocks.get(*bname) {
            let (count, data) = split_count(lines)?;
            let before = ps.angles.len();
            for line in &data {
                let v = tokenize(line);
                if v.len() < 5 {
                    continue;
                }
                let i = parse_usize(v[0])? - 1;
                let j = parse_usize(v[1])? - 1;
                let k = parse_usize(v[2])? - 1;
                let a_type = parse_usize(v[3])? - 1;
                let b_type = parse_usize(v[4])? - 1;
                ps.angles.push(PerturbedAngle {
                    i,
                    j,
                    k,
                    a_type,
                    b_type,
                });
            }
            if ps.angles.len() - before != count {
                log::warn!(
                    "{bname}: expected {count} angles, got {}",
                    ps.angles.len() - before
                );
            }
        }
    }

    // ── PERTIMPROPERDIHH / PERTIMPROPERDIH ───────────────────────────────────
    // Format: count \n i j k l A_type B_type (H block first, `in_perturbation.cc:827-828`)
    for bname in &["PERTIMPROPERDIHH", "PERTIMPROPERDIH", "PERTIMPDIHEDRAL03"] {
        if let Some(lines) = blocks.get(*bname) {
            let (count, data) = split_count(lines)?;
            let before = ps.improper_dihedrals.len();
            for line in &data {
                let v = tokenize(line);
                if v.len() < 6 {
                    continue;
                }
                let i = parse_usize(v[0])? - 1;
                let j = parse_usize(v[1])? - 1;
                let k = parse_usize(v[2])? - 1;
                let l = parse_usize(v[3])? - 1;
                let a_type = parse_usize(v[4])? - 1;
                let b_type = parse_usize(v[5])? - 1;
                ps.improper_dihedrals.push(PerturbedDihedral {
                    i,
                    j,
                    k,
                    l,
                    a_type,
                    b_type,
                });
            }
            if ps.improper_dihedrals.len() - before != count {
                log::warn!(
                    "{bname}: expected {count} impropers, got {}",
                    ps.improper_dihedrals.len() - before
                );
            }
        }
    }

    // ── PERTPROPERDIHH / PERTPROPERDIH ───────────────────────────────────────
    // Format: count \n i j k l A_type B_type (H block first, `in_perturbation.cc:1058-1059`)
    for bname in &["PERTPROPERDIHH", "PERTPROPERDIH", "PERTDIHEDRAL03"] {
        if let Some(lines) = blocks.get(*bname) {
            let (count, data) = split_count(lines)?;
            let before = ps.proper_dihedrals.len();
            for line in &data {
                let v = tokenize(line);
                if v.len() < 6 {
                    continue;
                }
                let i = parse_usize(v[0])? - 1;
                let j = parse_usize(v[1])? - 1;
                let k = parse_usize(v[2])? - 1;
                let l = parse_usize(v[3])? - 1;
                let a_type = parse_usize(v[4])? - 1;
                let b_type = parse_usize(v[5])? - 1;
                ps.proper_dihedrals.push(PerturbedDihedral {
                    i,
                    j,
                    k,
                    l,
                    a_type,
                    b_type,
                });
            }
            if ps.proper_dihedrals.len() - before != count {
                log::warn!(
                    "{bname}: expected {count} dihedrals, got {}",
                    ps.proper_dihedrals.len() - before
                );
            }
        }
    }

    // ── PERTBONDSOFT ─────────────────────────────────────────────────────────
    // Format: count \n i j A_type B_type alpha   (1-indexed; B_type=0 → absent in state B)
    // Uses harmonic bond parameters (k_harmonic). Alpha is the softness parameter.
    if let Some(lines) = blocks.get("PERTBONDSOFT") {
        let (count, data) = split_count(lines)?;
        for line in &data {
            let v = tokenize(line);
            if v.len() < 5 {
                continue;
            }
            let i = parse_usize(v[0])? - 1;
            let j = parse_usize(v[1])? - 1;
            let a_type = parse_usize(v[2])? - 1;
            let b_raw = parse_usize(v[3])?;
            let b_type = if b_raw == 0 { None } else { Some(b_raw - 1) };
            let alpha = parse_f64(v[4])?;
            ps.soft_bonds.push(SoftBond {
                i,
                j,
                a_type,
                b_type,
                alpha,
            });
        }
        if ps.soft_bonds.len() != count {
            log::warn!(
                "PERTBONDSOFT: expected {count} bonds, got {}",
                ps.soft_bonds.len()
            );
        }
    }

    // ── PERTANGLESOFT ────────────────────────────────────────────────────────
    // Format: count \n i j k A_type B_type alpha   (1-indexed; B_type=0 → absent in state B)
    // Uses cos-harmonic angle parameters.
    if let Some(lines) = blocks.get("PERTANGLESOFT") {
        let (count, data) = split_count(lines)?;
        for line in &data {
            let v = tokenize(line);
            if v.len() < 6 {
                continue;
            }
            let i = parse_usize(v[0])? - 1;
            let j = parse_usize(v[1])? - 1;
            let k = parse_usize(v[2])? - 1;
            let a_type = parse_usize(v[3])? - 1;
            let b_raw = parse_usize(v[4])?;
            let b_type = if b_raw == 0 { None } else { Some(b_raw - 1) };
            let alpha = parse_f64(v[5])?;
            ps.soft_angles.push(SoftAngle {
                i,
                j,
                k,
                a_type,
                b_type,
                alpha,
            });
        }
        if ps.soft_angles.len() != count {
            log::warn!(
                "PERTANGLESOFT: expected {count} angles, got {}",
                ps.soft_angles.len()
            );
        }
    }

    // ── PERTIMPROPERDIHSOFT ──────────────────────────────────────────────────
    // Format: count \n i j k l A_type B_type alpha  (1-indexed; B_type=0 → absent in state B)
    if let Some(lines) = blocks.get("PERTIMPROPERDIHSOFT") {
        let (count, data) = split_count(lines)?;
        for line in &data {
            let v = tokenize(line);
            if v.len() < 7 {
                continue;
            }
            let i = parse_usize(v[0])? - 1;
            let j = parse_usize(v[1])? - 1;
            let k = parse_usize(v[2])? - 1;
            let l = parse_usize(v[3])? - 1;
            let a_type = parse_usize(v[4])? - 1;
            let b_raw = parse_usize(v[5])?;
            let b_type = if b_raw == 0 { None } else { Some(b_raw - 1) };
            let alpha = parse_f64(v[6])?;
            ps.soft_impropers.push(SoftImproper {
                i,
                j,
                k,
                l,
                a_type,
                b_type,
                alpha,
            });
        }
        if ps.soft_impropers.len() != count {
            log::warn!(
                "PERTIMPROPERDIHSOFT: expected {count} impropers, got {}",
                ps.soft_impropers.len()
            );
        }
    }

    Ok(result)
}

// ─── Helpers ────────────────────────────────────────────────────────────────

/// Parse the file into named blocks: block_name → data lines (comments stripped).
fn parse_blocks(content: &str) -> HashMap<String, Vec<String>> {
    let mut blocks: HashMap<String, Vec<String>> = HashMap::new();
    let mut current: Option<String> = None;

    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END" {
            current = None;
            continue;
        }
        if trimmed == "TITLE" {
            current = Some("TITLE".to_string());
            continue;
        }

        // Block header: starts with uppercase letter, contains only alnum/underscore, no spaces.
        // This distinguishes "PERTATOMPARAM" (header) from "3" (count line, no uppercase start).
        let is_header = !trimmed.contains(' ')
            && trimmed.starts_with(|c: char| c.is_ascii_uppercase())
            && trimmed
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || c == '_');
        if is_header {
            current = Some(trimmed.to_string());
            blocks.entry(trimmed.to_string()).or_default();
            continue;
        }

        if let Some(ref name) = current {
            if name != "TITLE" {
                blocks
                    .entry(name.clone())
                    .or_default()
                    .push(trimmed.to_string());
            }
        }
    }
    blocks
}

/// Split block lines into (count, data_lines): first line is count, rest are data.
fn split_count(lines: &[String]) -> Result<(usize, Vec<String>), IoError> {
    if lines.is_empty() {
        return Ok((0, Vec::new()));
    }
    let count: usize = lines[0]
        .trim()
        .parse()
        .map_err(|_| IoError::ParseError(format!("expected count, got: {}", lines[0])))?;
    Ok((count, lines[1..].to_vec()))
}

fn tokenize(line: &str) -> Vec<&str> {
    line.split_whitespace().collect()
}

fn parse_usize(s: &str) -> Result<usize, IoError> {
    s.parse::<usize>()
        .map_err(|_| IoError::ParseError(format!("expected usize: {s}")))
}

fn parse_f64(s: &str) -> Result<f64, IoError> {
    s.parse::<f64>()
        .map_err(|_| IoError::ParseError(format!("expected f64: {s}")))
}

// ─── Legacy writer (kept for py-gromos tooling) ─────────────────────────────

use std::fs::File;
use std::io::{self, Write};

use gromos_core::topology::Topology;

/// Writer for a simple perturbation topology file.
/// Note: this writes a non-standard format (PERTURBEDATOM block);
/// the reader above handles the real GROMOS PERTATOMPARAM format.
pub struct PtpWriter {
    file: File,
}

impl PtpWriter {
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            file: File::create(path)?,
        })
    }

    pub fn write(
        &mut self,
        topology_a: &Topology,
        topology_b: &Topology,
        alpha_lj: f64,
        alpha_crf: f64,
        title: &str,
    ) -> io::Result<()> {
        if topology_a.num_atoms() != topology_b.num_atoms() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "atom count mismatch",
            ));
        }
        writeln!(self.file, "TITLE\n{title}\nEND")?;
        self.write_perturbed_atoms(topology_a, topology_b, alpha_lj, alpha_crf)?;
        // The other blocks gromos++ `make_pt_top` writes, empty: this writer perturbs atom
        // parameters only (: bonded perturbations come from gromos++ or by hand).
        for (block, what, cols) in [
            (
                "PERTATOMPAIR",
                "perturbed atom pairs",
                "#    i     j  i(A)  i(B)",
            ),
            (
                "PERTBONDSTRETCH",
                "perturbed bonds",
                "#    i     j t(A) t(B)",
            ),
            (
                "PERTBONDSTRETCHH",
                "perturbed bonds",
                "#    i     j t(A) t(B)",
            ),
            (
                "PERTBONDANGLE",
                "perturbed bond angles",
                "#    i     j     k t(A) t(B)",
            ),
            (
                "PERTBONDANGLEH",
                "perturbed bond angles",
                "#    i     j     k t(A) t(B)",
            ),
            (
                "PERTIMPROPERDIH",
                "perturbed improper dihedrals",
                "#    i     j     k     l t(A) t(B)",
            ),
            (
                "PERTIMPROPERDIHH",
                "perturbed improper dihedrals",
                "#    i     j     k     l t(A) t(B)",
            ),
            (
                "PERTPROPERDIH",
                "perturbed dihedrals",
                "#    i     j     k     l t(A) t(B)",
            ),
            (
                "PERTPROPERDIHH",
                "perturbed dihedrals",
                "#    i     j     k     l t(A) t(B)",
            ),
        ] {
            writeln!(self.file, "{block}\n# number of {what}\n0\n{cols}\nEND")?;
        }
        Ok(())
    }

    fn write_perturbed_atoms(
        &mut self,
        a: &Topology,
        b: &Topology,
        alpha_lj: f64,
        alpha_crf: f64,
    ) -> io::Result<()> {
        // gromosXX / gromos++ layout (`PERTATOMPARAM`): NR RES NAME IAC(A) MASS(A) CHARGE(A)
        // IAC(B) MASS(B) CHARGE(B) ALJ ACRF, IACs 1-based in the file.
        let mut lines = Vec::new();
        for i in 0..a.num_atoms() {
            let (iac_a, iac_b) = (a.iac[i], b.iac[i]);
            let (ma, mb) = (a.mass[i], b.mass[i]);
            let (ca, cb) = (a.charge[i], b.charge[i]);
            if iac_a != iac_b || (ma - mb).abs() > 1e-6 || (ca - cb).abs() > 1e-6 {
                lines.push(format!(
                    "{:5}{:5}{:>6}{:4}{:11.5}{:11.5}{:4}{:11.5}{:11.5}{:11.2}{:11.2}",
                    i + 1,
                    a.residue_nr(i).unwrap_or(1),
                    a.atom_name(i).unwrap_or(""),
                    iac_a + 1,
                    ma,
                    ca,
                    iac_b + 1,
                    mb,
                    cb,
                    alpha_lj,
                    alpha_crf
                ));
            }
        }
        writeln!(self.file, "PERTATOMPARAM")?;
        writeln!(self.file, "# number of perturbed atoms")?;
        writeln!(self.file, "{:5}", lines.len())?;
        writeln!(
            self.file,
            "#  NR  RES  NAME IAC(A)  MASS(A) CHARGE(A) IAC(B)  MASS(B) CHARGE(B)       ALJ       ACRF"
        )?;
        for l in lines {
            writeln!(self.file, "{l}")?;
        }
        writeln!(self.file, "END")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(content: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join("gromos_ptp_test.tmp");
        std::fs::File::create(&path)
            .unwrap()
            .write_all(content.as_bytes())
            .unwrap();
        path
    }

    /// Parse the actual ch4_spc_dummy.ptp used in the ch4_water_fep reference test.
    /// This asserts lj_soft and crf_soft are read correctly — the exact values that
    /// feed into the effective soft-core alpha (alpha_lj = lj_soft × global_ALPHLJ).
    #[test]
    fn test_ch4_dummy_ptp() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../crates/gromos-md/tests/ch4_spc_dummy.ptp"
        );
        if !std::path::Path::new(path).exists() {
            eprintln!("SKIP: ch4_spc_dummy.ptp not found");
            return;
        }
        let pt = read_pttopo(path).expect("read_pttopo failed");
        let atoms = &pt.perturbed_solute.atoms;
        assert_eq!(atoms.len(), 1, "CH4→dummy has 1 perturbed atom");
        let cm = &atoms[0];
        assert_eq!(cm.seq, 0, "CM is atom 0 (0-indexed)");
        assert_eq!(cm.a_iac, 16, "CM state-A IAC=17 (0-indexed: 16)");
        assert_eq!(cm.b_iac, 21, "dummy state-B IAC=22 (0-indexed: 21)");
        assert!((cm.a_charge - 0.0).abs() < 1e-6, "CM charge A = 0");
        assert!((cm.b_charge - 0.0).abs() < 1e-6, "dummy charge B = 0");
        // These are the per-atom scaling factors — multiplied by global ALPHLJ/ALPHC
        // in Forcefield::build_pert_info. If these are wrong the soft-core distances
        // will be wrong and the FEP energies will be off.
        assert!((cm.lj_soft - 0.5).abs() < 1e-6, "lj_soft = 0.5");
        assert!((cm.crf_soft - 0.5).abs() < 1e-6, "crf_soft = 0.5");
    }

    /// Minimal inline ptp with the same PERTATOMPARAM structure as ch4_spc_dummy.ptp.
    /// Does not depend on any external file — always runs.
    #[test]
    fn test_pertatomparam_alpha_columns() {
        let content = "\
TITLE
test
END
PERTATOMPARAM
   1
#  NR  RES  NAME IAC(A) MASS(A) CHARGE(A) IAC(B)  MASS(B) CHARGE(B)  ALJ  ACRF
    1    1     CM  17    16.04300  0.00000   22   16.04300  0.00000   0.5   0.5
END
";
        let path = write_tmp(content);
        let pt = read_pttopo(&path).expect("read_pttopo");
        let a = &pt.perturbed_solute.atoms[0];
        assert!(
            (a.lj_soft - 0.5).abs() < 1e-9,
            "lj_soft must be read from column 9 (ALJ)"
        );
        assert!(
            (a.crf_soft - 0.5).abs() < 1e-9,
            "crf_soft must be read from column 10 (ACRF)"
        );
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn hydrogen_blocks_come_first_and_are_marked() {
        // gromosXX reads PERTBONDSTRETCHH before PERTBONDSTRETCH (in_perturbation.cc:172-175)
        let content = "\
TITLE
order
END
PERTBONDSTRETCH
   2
    2    4   26   20
    8   10    9   26
END
PERTBONDSTRETCHH
   1
    4    5    2   36
END
PERTBONDANGLE
   1
    1    2    3    1    2
END
PERTBONDANGLEH
   1
    3    4    5    3    4
END
";
        let path =
            std::env::temp_dir().join(format!("gromos_ptp_order_{}.tmp", std::process::id()));
        std::fs::write(&path, content).unwrap();
        let ps = read_pttopo(&path).expect("read_pttopo").perturbed_solute;
        let bonds: Vec<(usize, usize, bool)> =
            ps.bonds.iter().map(|b| (b.i, b.j, b.hydrogen)).collect();
        assert_eq!(bonds, vec![(3, 4, true), (1, 3, false), (7, 9, false)]);
        let angles: Vec<(usize, usize, usize)> =
            ps.angles.iter().map(|a| (a.i, a.j, a.k)).collect();
        assert_eq!(angles, vec![(2, 3, 4), (0, 1, 2)]);
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_read_pttopo_aladip() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../.local/GROMOS/md++/src/check/data/aladip.pttopo"
        );
        if !std::path::Path::new(path).exists() {
            eprintln!("SKIP: aladip.pttopo not found at {path}");
            return;
        }
        let pt = read_pttopo(path).expect("read_pttopo failed");
        let ps = &pt.perturbed_solute;

        // aladip.pttopo has 3 perturbed atoms, 1 pair, 2 bonds, 2 angles, 2 impropers, 2 dihedrals
        assert_eq!(ps.atoms.len(), 3, "expected 3 perturbed atoms");
        assert_eq!(ps.atom_pairs.len(), 1, "expected 1 perturbed pair");
        assert_eq!(ps.bonds.len(), 2, "expected 2 perturbed bonds");
        assert_eq!(ps.angles.len(), 2, "expected 2 perturbed angles");
        assert_eq!(
            ps.improper_dihedrals.len(),
            2,
            "expected 2 perturbed impropers"
        );
        assert_eq!(
            ps.proper_dihedrals.len(),
            2,
            "expected 2 perturbed dihedrals"
        );

        // Atom 0 (seq=0, 1-indexed seq=1 in file): CB, IAC 14→11 (0-indexed: 13→10)
        let a0 = &ps.atoms[0];
        assert_eq!(a0.seq, 0);
        assert_eq!(a0.a_iac, 13); // 14-1
        assert_eq!(a0.b_iac, 10); // 11-1
        assert!((a0.a_mass - 15.035).abs() < 0.01);
        assert!((a0.a_charge - 0.0).abs() < 1e-6);
        assert!((a0.b_charge - 1.0).abs() < 1e-6);
        assert!((a0.lj_soft - 0.3).abs() < 1e-6);
        assert!((a0.crf_soft - 0.6).abs() < 1e-6);

        // is_perturbed: atoms 0, 3, 10 are perturbed (0-indexed seq from file: 1,4,11 → 0,3,10)
        assert!(pt.is_perturbed[0]);
        assert!(pt.is_perturbed[3]);
        assert!(pt.is_perturbed[10]);
    }
}
