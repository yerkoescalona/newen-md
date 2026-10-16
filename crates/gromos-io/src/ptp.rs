//! Perturbation topology (.pttopo) reader and writer.
//!
//! The GROMOS `.pttopo` format defines dual-topology (state A / state B)
//! parameters for atoms that change as λ goes from 0 to 1 (FEP / TI).
//!
//! Key blocks parsed by the reader:
//! - `PERTATOMPARAM`  — per-atom IAC, mass, charge for states A and B; soft-core alphas
//! - `PERTATOMPAIR`   — atom pairs whose interaction type (excluded / normal / 1-4) is perturbed
//! - `PERTBONDSTRETCH` / `PERTBONDSTRETCHH` — perturbed bond types
//! - `PERTBONDANGLE`  / `PERTBONDANGLEH`   — perturbed angle types
//! - `PERTIMPROPERDIH` / `PERTIMPROPERDIHH` — perturbed improper dihedral types
//! - `PERTPROPERDIH`  / `PERTPROPERDIHH`   — perturbed proper dihedral types
//! - `MPERTATOM`      — the end states of an EDS run (type and charge of each atom per state)
//! - `PERTBONDSOFT` / `PERTANGLESOFT` / `PERTIMPROPERDIHSOFT` — soft-core terms, type 0 for a
//!   term absent in that state (made or broken)

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use gromos_core::topology::{
    EdsPerturbedAtom, PairInteraction, PerturbedAngle, PerturbedAtom, PerturbedAtomPair,
    PerturbedBond, PerturbedDihedral, PerturbedSolute, SoftAngle, SoftBond, SoftImproper,
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
    /// The `MPERTATOM` block (EDS end states), if the file has one.
    pub eds: Option<EdsPerturbation>,
    /// The names of every block in the file except `TITLE`, sorted.
    pub blocks: Vec<String>,
}

/// An `MPERTATOM` block (gromosXX `in_perturbation.cc:1625-1801`): the type and charge of each
/// listed atom in every EDS end state, in file order (the order decides which of two EDS atoms
/// owns their excluded or 1-4 pair, `Topology::apply_eds_perturbation`).
#[derive(Debug, Clone, PartialEq)]
pub struct EdsPerturbation {
    /// NPTB, the number of end states.
    pub num_states: usize,
    /// The end-state names.
    pub state_names: Vec<String>,
    /// One per atom line, 0-indexed seq and types, partner lists empty.
    pub atoms: Vec<EdsPerturbedAtom>,
}

/// Parse the data lines of an `MPERTATOM` block: `NJLA NPTB`, NPTB state names, then NJLA
/// rows `seq name IAC[0] CG[0] … IAC[NPTB-1] CG[NPTB-1]` (1-based seq and IAC). Anything after
/// the last charge (the softness columns some files carry) is ignored, as gromosXX does.
pub fn parse_mpertatom(lines: &[String]) -> Result<EdsPerturbation, IoError> {
    let err = |m: String| IoError::ParseError(format!("MPERTATOM: {m}"));
    let mut it = lines.iter();
    let header = it
        .next()
        .ok_or_else(|| err("missing the NJLA NPTB line".to_string()))?;
    let h = tokenize(header);
    let (Some(njla), Some(nptb)) = (
        h.first().and_then(|v| v.parse::<usize>().ok()),
        h.get(1).and_then(|v| v.parse::<usize>().ok()),
    ) else {
        return Err(err(format!("bad line (NJLA NPTB): {header:?}")));
    };
    let names_line = it
        .next()
        .ok_or_else(|| err("missing the state-name line".to_string()))?;
    let state_names: Vec<String> = tokenize(names_line).iter().map(|s| s.to_string()).collect();
    if state_names.len() < nptb {
        return Err(err(format!(
            "{nptb} states but {} state names: {names_line:?}",
            state_names.len()
        )));
    }
    let state_names = state_names[..nptb].to_vec();
    let mut atoms = Vec::with_capacity(njla);
    for line in it {
        let v = tokenize(line);
        if v.len() < 2 + 2 * nptb {
            return Err(err(format!(
                "bad line (seq name IAC CG per state): {line:?}"
            )));
        }
        let seq = v[0]
            .parse::<usize>()
            .ok()
            .filter(|&s| s >= 1)
            .ok_or_else(|| err(format!("atom sequence number wrong: {line:?}")))?;
        let mut iac = Vec::with_capacity(nptb);
        let mut charge = Vec::with_capacity(nptb);
        for s in 0..nptb {
            let t = v[2 + 2 * s]
                .parse::<usize>()
                .ok()
                .filter(|&t| t >= 1)
                .ok_or_else(|| err(format!("IAC must be a positive integer: {line:?}")))?;
            iac.push(t - 1);
            charge.push(
                v[3 + 2 * s]
                    .parse::<f64>()
                    .map_err(|_| err(format!("bad charge: {line:?}")))?,
            );
        }
        atoms.push(EdsPerturbedAtom {
            seq: seq - 1,
            iac,
            charge,
            exclusions: Vec::new(),
            one_four_pairs: Vec::new(),
        });
    }
    if atoms.len() != njla {
        return Err(err(format!(
            "wrong number of perturbed atoms: NJLA = {njla}, {} rows",
            atoms.len()
        )));
    }
    Ok(EdsPerturbation {
        num_states: nptb,
        state_names,
        atoms,
    })
}

/// The blocks [`parse_pttopo`] applies.
const READ_BLOCKS: [&str; 14] = [
    "PERTATOMPAIR",
    "PERTATOMPARAM",
    "PERTBONDSTRETCHH",
    "PERTBONDSTRETCH",
    "PERTBONDANGLEH",
    "PERTBONDANGLE",
    "PERTIMPROPERDIHH",
    "PERTIMPROPERDIH",
    "PERTPROPERDIHH",
    "PERTPROPERDIH",
    "PERTBONDSOFT",
    "PERTANGLESOFT",
    "PERTIMPROPERDIHSOFT",
    "MPERTATOM",
];

/// Read a GROMOS perturbation topology file (`.pttopo`).
///
/// Atoms are 1-indexed in the file; returned structs use 0-indexed atoms.
/// IAC values are 1-indexed in the file; returned structs use 0-indexed IAC.
pub fn read_pttopo<P: AsRef<Path>>(path: P) -> Result<PerturbedTopology, IoError> {
    let content = fs::read_to_string(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    parse_pttopo(&content)
}

/// [`read_pttopo`] on the text of a perturbation topology.
pub fn parse_pttopo(content: &str) -> Result<PerturbedTopology, IoError> {
    // Split into named blocks: block_name → Vec<non-comment lines>
    let blocks = parse_blocks(content);

    // The pre-2005 block names are errors in gromosXX (`in_perturbation.cc:164-170`, 580-586,
    // 817-823, 1149-1156, 1245-1252). PERTDIHEDRAL03 has no rename message there: it is never
    // read, so the unread-block check (`in_perturbation.cc:1809-1813`) refuses it.
    for (old, new) in [
        ("PERTBOND03", "PERTBONDSTRETCH"),
        ("PERTBANGLE03", "PERTBONDANGLE"),
        ("PERTIMPDIHEDRAL03", "PERTIMPROPERDIH"),
        ("PERTATOMPAIR03", "PERTATOMPAIR"),
        ("PERTATOM03", "PERTATOMPARAM"),
    ] {
        if blocks.contains_key(old) {
            return Err(IoError::ParseError(format!(
                "The {old} block was renamed to {new}."
            )));
        }
    }
    if blocks.contains_key("LAMBDADEP") {
        return Err(IoError::ParseError(
            "LAMBDADEP block in perturbation topology is deprecated. Use LAMBDAS block in the \
             input file"
                .to_string(),
        ));
    }
    // Every other block this reader does not apply is refused, as gromosXX refuses the blocks
    // it does not read (`in_perturbation.cc:1803-1820`). gromosXX does read PERTCONSTRAINT03
    // (NTC ≠ 1) and PERTPOLPARAM, which gromos-rs does not model. SCALEDINTERACTIONS passes:
    // the run refuses NSCALE ≠ 0 instead (`check_perturbation`, gromos-run).
    let mut unread: Vec<&String> = blocks
        .keys()
        .filter(|b| !READ_BLOCKS.contains(&b.as_str()) && b.as_str() != "SCALEDINTERACTIONS")
        .collect();
    unread.sort();
    if let Some(b) = unread.first() {
        return Err(IoError::ParseError(format!("block {b} not supported!")));
    }

    let mut names: Vec<String> = blocks.keys().cloned().collect();
    names.sort();
    let mut result = PerturbedTopology {
        eds: blocks
            .get("MPERTATOM")
            .map(|lines| parse_mpertatom(lines))
            .transpose()?,
        blocks: names,
        ..Default::default()
    };
    let ps = &mut result.perturbed_solute;

    // ── PERTATOMPAIR ─────────────────────────────────────────────────────────
    // Must be read before PERTATOMPARAM (GROMOS convention — exclusion fixup).
    // Format: count \n i j A_type B_type \n ...   (all 1-indexed)
    if let Some(lines) = blocks.get("PERTATOMPAIR") {
        let (count, data) = split_count(lines)?;
        for line in &data {
            let v = tokenize(line);
            if v.len() < 4 {
                continue;
            }
            let i = parse_usize(v[0])? - 1; // 1-indexed in file → 0-indexed
            let j = parse_usize(v[1])? - 1;
            // Interaction types as in the file: 0 excluded, 1 normal, 2 1-4
            let pair_type = |s: &str| {
                s.parse::<i64>()
                    .ok()
                    .and_then(PairInteraction::from_code)
                    .ok_or_else(|| {
                        IoError::ParseError(format!(
                            "PERTATOMPAIR: interaction type must be 0, 1 or 2, got {s}"
                        ))
                    })
            };
            let a_type = pair_type(v[2])?;
            let b_type = pair_type(v[3])?;
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
    }

    // ── PERTATOMPARAM ────────────────────────────────────────────────────────
    // Format: count \n seq res name a_iac a_mass a_charge b_iac b_mass b_charge lj_soft crf_soft
    // All IAC values 1-indexed; seq is 1-indexed atom number.
    if let Some(lines) = blocks.get("PERTATOMPARAM") {
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
    }

    // ── PERTBONDSTRETCHH / PERTBONDSTRETCH ───────────────────────────────────
    // Format: count \n i j A_type B_type   (1-indexed). The H block first, as gromosXX reads
    // them (`in_perturbation.cc:172-175`): with NTC = 3 this is the perturbed SHAKE order.
    for bname in &["PERTBONDSTRETCHH", "PERTBONDSTRETCH"] {
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
    for bname in &["PERTBONDANGLEH", "PERTBONDANGLE"] {
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
    for bname in &["PERTIMPROPERDIHH", "PERTIMPROPERDIH"] {
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
    for bname in &["PERTPROPERDIHH", "PERTPROPERDIH"] {
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
    // Format: count \n i j A_type B_type alpha   (1-indexed; type 0 → absent in that state)
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
            let (a_type, b_type) = soft_types(v[2], v[3], "bond", "PERTBONDSOFT")?;
            // gromosXX's warnings for a bond made or broken (`in_perturbation.cc:444, 448`)
            if a_type.is_none() {
                log::warn!(
                    "You are making a bond, make sure you adapt the excluded atoms to this change."
                );
            } else if b_type.is_none() {
                log::warn!(
                    "You are breaking a bond, make sure you adapt the excluded atoms to this change."
                );
            }
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
    // Format: count \n i j k A_type B_type alpha   (1-indexed; type 0 → absent in that state)
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
            let (a_type, b_type) = soft_types(v[3], v[4], "angle", "PERTANGLESOFT")?;
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
    // Format: count \n i j k l A_type B_type alpha  (1-indexed; type 0 → absent in that state)
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
            let (a_type, b_type) = soft_types(v[4], v[5], "improper", "PERTIMPROPERDIHSOFT")?;
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

        // Block header: outside a block (a TITLE line or a one-word state-name line is data),
        // starts with uppercase letter, contains only alnum/underscore, no spaces.
        // This distinguishes "PERTATOMPARAM" (header) from "3" (count line, no uppercase start).
        let is_header = current.is_none()
            && !trimmed.contains(' ')
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

/// The state A and B types of a soft-bonded line (1-based in the file, 0 = absent in that
/// state, as gromosXX reads them, `in_perturbation.cc:435-449, 765-777, 1003-1015`). Both 0 is
/// gromosXX's "Conversion from bond type 0 to bond type 0" error.
fn soft_types(
    a: &str,
    b: &str,
    term: &str,
    block: &str,
) -> Result<(Option<usize>, Option<usize>), IoError> {
    let a = parse_usize(a)?.checked_sub(1);
    let b = parse_usize(b)?.checked_sub(1);
    if a.is_none() && b.is_none() {
        return Err(IoError::ParseError(format!(
            "Conversion from {term} type 0 to {term} type 0 in {block} block."
        )));
    }
    Ok((a, b))
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

/// Writer for a perturbation topology file (gromos++ `make_pt_top` layout): the atoms whose
/// IAC, mass or charge differ between two topologies in `PERTATOMPARAM`, the bonded blocks
/// empty. [`write_ptp`] writes the same text to any writer.
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
        write_ptp(
            &mut self.file,
            topology_a,
            topology_b,
            alpha_lj,
            alpha_crf,
            title,
        )
    }
}

/// The perturbation topology from state A (`topology_a`) to state B (`topology_b`), every
/// perturbed atom with the soft-core parameters `alpha_lj`, `alpha_crf`.
pub fn write_ptp(
    out: &mut dyn Write,
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
    writeln!(out, "TITLE\n{title}\nEND")?;
    write_perturbed_atoms(out, topology_a, topology_b, alpha_lj, alpha_crf)?;
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
        writeln!(out, "{block}\n# number of {what}\n0\n{cols}\nEND")?;
    }
    Ok(())
}

fn write_perturbed_atoms(
    out: &mut dyn Write,
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
    writeln!(out, "PERTATOMPARAM")?;
    writeln!(out, "# number of perturbed atoms")?;
    writeln!(out, "{:5}", lines.len())?;
    writeln!(
        out,
        "#  NR  RES  NAME IAC(A)  MASS(A) CHARGE(A) IAC(B)  MASS(B) CHARGE(B)       ALJ       ACRF"
    )?;
    for l in lines {
        writeln!(out, "{l}")?;
    }
    writeln!(out, "END")
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

    /// PERTATOMPAIR interaction types are read as they are (0 excluded, 1 normal, 2 1-4), the
    /// atoms ordered i < j; any other type is an error.
    #[test]
    fn pertatompair_types_are_read_raw() {
        let content = "\
TITLE
pairs
END
PERTATOMPAIR
   3
    7   10    2    0
   11    1    0    2
    4    7    1    1
END
";
        // own file names: `write_tmp`'s single path races with the other tests
        let tmp = |tag: &str, content: &str| {
            let p = std::env::temp_dir()
                .join(format!("gromos_ptp_pairs_{tag}_{}.tmp", std::process::id()));
            std::fs::write(&p, content).unwrap();
            p
        };
        let path = tmp("ok", content);
        let pt = read_pttopo(&path).expect("read_pttopo");
        std::fs::remove_file(path).ok();
        let pairs: Vec<_> = pt
            .perturbed_solute
            .atom_pairs
            .iter()
            .map(|p| (p.i, p.j, p.a_type, p.b_type))
            .collect();
        use PairInteraction::{Excluded, Normal, OneFour};
        assert_eq!(
            pairs,
            vec![
                (6, 9, OneFour, Excluded),
                (0, 10, Excluded, OneFour),
                (3, 6, Normal, Normal),
            ]
        );

        let bad = tmp(
            "bad",
            "TITLE\nbad\nEND\nPERTATOMPAIR\n   1\n    1    2    3    0\nEND\n",
        );
        assert!(read_pttopo(&bad).is_err(), "type 3 must be rejected");
        std::fs::remove_file(bad).ok();
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
    fn mpertatom_reads_states_and_ignores_trailing_columns() {
        // the LiveCoMS t_05 layout: tabs, upper-case state names, softness columns at the end
        let content = "TITLE\nEDS\nEND\nMPERTATOM\n 2  3\n# names\n   STATEA\tSTATEB   STATEC\n\
            1  OW\t5  -0.82   22  0.0   5 -0.82   1.00  1.00\n\
            3  HW2 21  0.41   22  0.0  21  0.41\nEND\n";
        let path =
            std::env::temp_dir().join(format!("gromos_ptp_mpert_{}.tmp", std::process::id()));
        std::fs::write(&path, content).unwrap();
        let pt = read_pttopo(&path).expect("read_pttopo");
        std::fs::remove_file(&path).ok();
        assert_eq!(pt.blocks, ["MPERTATOM"]);
        assert!(pt.perturbed_solute.is_empty());
        let eds = pt.eds.expect("MPERTATOM");
        assert_eq!(eds.num_states, 3);
        assert_eq!(eds.state_names, ["STATEA", "STATEB", "STATEC"]);
        assert_eq!(eds.atoms.len(), 2);
        assert_eq!(eds.atoms[0].seq, 0);
        assert_eq!(eds.atoms[0].iac, [4, 21, 4]);
        assert_eq!(eds.atoms[0].charge, [-0.82, 0.0, -0.82]);
        assert_eq!(eds.atoms[1].seq, 2);
        assert_eq!(eds.atoms[1].iac, [20, 21, 20]);
    }

    #[test]
    fn mpertatom_errors_name_the_problem() {
        let parse = |body: &str| {
            let lines: Vec<String> = body.lines().map(str::to_string).collect();
            parse_mpertatom(&lines).unwrap_err().to_string()
        };
        // NJLA says 2, one row
        assert!(parse("2 2\nA B\n1 CB 14 0.0 11 1.0").contains("wrong number"));
        // a row with too few columns for NPTB states
        assert!(parse("1 2\nA B\n1 CB 14 0.0 11").contains("bad line"));
        assert!(parse("x 2\nA B\n").contains("NJLA NPTB"));
        assert!(parse("1 2\nA B\n1 CB 0 0.0 11 1.0").contains("IAC"));
        assert!(parse("1 2\nA\n").contains("state names"));
    }

    /// The pre-2005 block names were read as their successors; gromosXX refuses them.
    #[test]
    fn renamed_03_blocks_are_refused() {
        for (block, body, msg) in [
            (
                "PERTBOND03",
                "1\n1 2 1 2",
                "The PERTBOND03 block was renamed to PERTBONDSTRETCH.",
            ),
            (
                "PERTBANGLE03",
                "1\n1 2 3 1 2",
                "The PERTBANGLE03 block was renamed to PERTBONDANGLE.",
            ),
            (
                "PERTIMPDIHEDRAL03",
                "1\n1 2 3 4 1 2",
                "The PERTIMPDIHEDRAL03 block was renamed to PERTIMPROPERDIH.",
            ),
            (
                "PERTATOMPAIR03",
                "1\n1 4 1 0",
                "The PERTATOMPAIR03 block was renamed to PERTATOMPAIR.",
            ),
            (
                "PERTATOM03",
                "1\n1 1 CB 1 15.0 0.0 2 15.0 0.0 1.0 1.0",
                "The PERTATOM03 block was renamed to PERTATOMPARAM.",
            ),
            (
                "PERTDIHEDRAL03",
                "1\n1 2 3 4 1 2",
                "block PERTDIHEDRAL03 not supported!",
            ),
        ] {
            let text = format!("TITLE\nold\nEND\n{block}\n{body}\nEND\n");
            let e = parse_pttopo(&text).unwrap_err().to_string();
            assert!(e.contains(msg), "{block}: {e}");
        }
    }

    /// A block the reader does not apply used to be dropped silently; it is refused now, as
    /// gromosXX refuses the blocks it does not read. SCALEDINTERACTIONS is left to the NSCALE
    /// check of the run.
    #[test]
    fn unread_blocks_are_refused() {
        let with = |block: &str| {
            parse_pttopo(&format!(
                "TITLE\nt\nEND\nPERTBONDSTRETCH\n0\nEND\n{block}\n1\n1 2 1 2\nEND\n"
            ))
        };
        for block in ["PERTCONSTRAINT03", "PERTPOLPARAM", "PERTATOM", "PERTFOO"] {
            let e = with(block).unwrap_err().to_string();
            assert!(e.contains(&format!("block {block} not supported!")), "{e}");
        }
        assert!(with("LAMBDADEP")
            .unwrap_err()
            .to_string()
            .contains("LAMBDADEP block in perturbation topology is deprecated."));
        assert!(with("SCALEDINTERACTIONS").is_ok());
    }

    /// Type 0 in a soft block is a term absent in that state (made or broken); it used to
    /// underflow in state A. Both 0 is refused, as gromosXX does.
    #[test]
    fn soft_blocks_read_type_0_as_an_absent_state() {
        let read = |block: &str, atoms: &str, types: &str| {
            parse_pttopo(&format!("{block}\n1\n{atoms} {types} 1.5\nEND\n"))
        };
        for (a, b, want) in [("0", "3", (None, Some(2))), ("3", "0", (Some(2), None))] {
            let types = format!("{a} {b}");
            let p = read("PERTBONDSOFT", "1 2", &types)
                .unwrap()
                .perturbed_solute;
            assert_eq!((p.soft_bonds[0].a_type, p.soft_bonds[0].b_type), want);
            let p = read("PERTANGLESOFT", "1 2 3", &types)
                .unwrap()
                .perturbed_solute;
            assert_eq!((p.soft_angles[0].a_type, p.soft_angles[0].b_type), want);
            let p = read("PERTIMPROPERDIHSOFT", "1 2 3 4", &types)
                .unwrap()
                .perturbed_solute;
            assert_eq!(
                (p.soft_impropers[0].a_type, p.soft_impropers[0].b_type),
                want
            );
        }
        for (block, atoms, term) in [
            ("PERTBONDSOFT", "1 2", "bond"),
            ("PERTANGLESOFT", "1 2 3", "angle"),
            ("PERTIMPROPERDIHSOFT", "1 2 3 4", "improper"),
        ] {
            let e = read(block, atoms, "0 0").unwrap_err().to_string();
            assert!(
                e.contains(&format!(
                    "Conversion from {term} type 0 to {term} type 0 in {block} block."
                )),
                "{e}"
            );
        }
    }

    /// A perturbed or soft bond must be a bond of the topology, matched on the ordered pair
    /// as gromosXX's `bondMatcher` does.
    #[test]
    fn perturbation_of_a_non_existing_bond_is_refused() {
        use gromos_core::topology::{Bond, MolTypeAtom};

        let topology = || {
            let mut t = Topology::new();
            for k in 0..3 {
                t.moltypes[0].atoms.push(MolTypeAtom {
                    name: format!("A{}", k + 1),
                    residue_nr: 1,
                    residue_name: "TEST".to_string(),
                    iac: 0,
                    mass: 1.0,
                    charge: 0.0,
                    is_perturbed: false,
                    is_polarisable: false,
                    is_coarse_grained: false,
                });
            }
            t.moltypes[0].bonds.push(Bond {
                i: 0,
                j: 1,
                bond_type: 0,
                hydrogen: false,
            });
            t.resize_atom_arrays();
            t
        };
        let apply = |block: &str, line: &str| {
            let pt = parse_pttopo(&format!("{block}\n1\n{line}\nEND\n")).unwrap();
            topology().apply_perturbation(pt.perturbed_solute, pt.is_perturbed)
        };
        assert!(apply("PERTBONDSTRETCH", "1 2 1 2").is_ok());
        assert!(apply("PERTBONDSOFT", "1 2 1 2 1.0").is_ok());
        assert!(
            apply("PERTBONDSOFT", "2 3 0 1 1.0").is_ok(),
            "a bond being made"
        );
        for (block, line) in [
            ("PERTBONDSTRETCH", "1 3 1 2"),
            ("PERTBONDSTRETCH", "2 1 1 2"),
            ("PERTBONDSTRETCHH", "2 3 1 2"),
            ("PERTBONDSOFT", "2 3 1 2 1.0"),
        ] {
            let e = apply(block, line).unwrap_err();
            assert!(
                e.starts_with(&format!(
                    "Perturbation of a non-existing bond in {block} block."
                )),
                "{block} {line}: {e}"
            );
        }
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
