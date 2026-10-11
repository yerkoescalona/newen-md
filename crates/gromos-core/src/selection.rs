//! Atom selection — gromos-rs AtomSpecifier grammar over a flat index space.
//!
//! Parses the gromos-rs specifier syntax and resolves it to a sorted,
//! deduplicated `Vec<usize>` of 0-based atom indices.
//!
//! # Grammar (subset implemented)
//!
//! ```text
//! spec     ::= term (';' term)*
//! term     ::= 'all' | 'no'
//!            | 'not' '(' spec ')'
//!            | 'minus' '(' spec ')'
//!            | mol_spec
//!            | global_range
//!
//! mol_spec ::= mol_sel ':' atom_sel
//! mol_sel  ::= 'a'                       -- all solute molecules (never solvent)
//!            | 'm'                       -- all solute molecules (role == Solute)
//!            | 's'                       -- all solvent molecules (role == Solvent)
//!            | N                         -- molecule N (1-based)
//!            | N '-' M                   -- molecules N through M
//!
//! atom_sel ::= 'a'                       -- all atoms in the molecule(s)
//!            | 'res' '(' res_spec ')'   -- residue directive
//!            | name_list                 -- atom names, e.g. "CA,N,C"
//!            | range_list               -- local 1-based atom numbers
//!
//! res_spec ::= residue_sel ':' name_or_range
//! residue_sel ::= N | N '-' M | NAME   -- residue numbers or names
//! name_or_range ::= 'a' | name_list | range_list
//!
//! global_range ::= range_list           -- 1-based global atom numbers
//! range_list ::= range_item (',' range_item)*
//! range_item ::= N | N '-' M
//! name_list  ::= NAME (',' NAME)*       -- pure alpha tokens
//! ```
//!
//! Multiple terms are joined by `;` (union). `not()` defines a permanent exclusion
//! mask applied after all union terms. `minus()` removes atoms from the current set.
//!
//! As in gromos++ `AtomSpecifier`, `a:` is every *solute* molecule: on a solvated system
//! `a:OW` selects no water oxygens, `s:OW` does.
//!
//! The implementation deliberately avoids `num_solute_atoms()` as a role boundary.
//! All molecule-level operations route through `topology.molecules` ranges, filtered by
//! the role of the instance that owns them.

use crate::topology::{Role, Topology};
use std::collections::HashSet;

/// Atom selection error with syntax hints.
#[derive(Debug, Clone)]
pub enum SelectionError {
    /// Generic parse failure. `input` is what was given; `hint` suggests the fix.
    ParseError(String),
    /// A number range like `3-1` is backwards.
    InvalidRange(String),
    /// Molecule index is out of range for the topology.
    InvalidMolecule(usize),
    /// Residue number not found.
    InvalidResidue(usize),
    /// Atom name not found in the selected atoms.
    InvalidAtomName(String),
    /// The spec parsed successfully but matched zero atoms.
    EmptySelection,
    /// `s:` was used on a topology without solvent molecules.
    NoSolvent,
}

const SYNTAX_HINT: &str = "\
  Valid syntax examples:
    all              → all atoms
    1:a              → all atoms in molecule 1
    1:1-10           → atoms 1–10 in molecule 1
    1:CA,N,C         → atoms named CA, N or C in molecule 1
    a:CA             → CA atoms in every solute molecule
    m:a              → all solute atoms  (role == Solute)
    m:CA             → CA atoms in solute
    s:a              → all solvent atoms (role == Solvent)
    s:1              → first solvent molecule
    1:res(2:CA)      → CA of residue 2 in molecule 1
    1:res(ALA:a)     → all atoms of ALA residues in molecule 1
    not(a:H) 1:a     → all solute atoms except hydrogens
    1:1-5;1:10-12    → union of two ranges";

impl std::fmt::Display for SelectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SelectionError::ParseError(m) => {
                write!(f, "selection parse error: {m}\n{SYNTAX_HINT}")
            },
            SelectionError::InvalidRange(m) => {
                write!(
                    f,
                    "invalid range '{m}' — ranges must be N-M with N ≤ M (e.g. '1-10')"
                )
            },
            SelectionError::InvalidMolecule(n) => {
                write!(
                    f,
                    "molecule {n} not found in topology\n  \
                    Hint: use '1:a' for the first molecule, 'm:a' for all solute, \
                    's:a' for all solvent"
                )
            },
            SelectionError::InvalidResidue(r) => {
                write!(
                    f,
                    "residue {r} not found\n  \
                    Hint: use '1:res({r}:a)' to select all atoms of residue {r}, \
                    or 'a:res(NAME:a)' to select by residue name"
                )
            },
            SelectionError::InvalidAtomName(n) => {
                write!(
                    f,
                    "atom name '{n}' not found in selected atoms\n  \
                    Hint: check the name in your topology; \
                    names are case-sensitive (e.g. 'CA' not 'ca')"
                )
            },
            SelectionError::EmptySelection => {
                write!(f, "selection matched no atoms\n{SYNTAX_HINT}")
            },
            SelectionError::NoSolvent => {
                write!(
                    f,
                    "'s:' selects solvent, and this topology has no solvent molecules \
                     (its SOLVENTATOM block is empty, or it is not solvated)\n  \
                     Hint: molecules of the solute are selected by number ('1-216:a') or \
                     with 'a:'"
                )
            },
        }
    }
}
impl std::error::Error for SelectionError {}

/// Sorted, deduplicated set of 0-based atom indices.
#[derive(Debug, Clone)]
pub struct AtomSelection {
    indices: Vec<usize>,
    // TODO: use for bounds-checking once from_indices validates lazily
    #[allow(dead_code)]
    n_atoms: usize,
}

impl AtomSelection {
    /// All atoms in the system.
    pub fn all(n_atoms: usize) -> Self {
        Self {
            indices: (0..n_atoms).collect(),
            n_atoms,
        }
    }

    /// From explicit 0-based indices (validated).
    pub fn from_indices(mut indices: Vec<usize>, n_atoms: usize) -> Result<Self, SelectionError> {
        for &i in &indices {
            if i >= n_atoms {
                return Err(SelectionError::InvalidRange(format!(
                    "atom {i} >= {n_atoms}"
                )));
            }
        }
        if indices.is_empty() {
            return Err(SelectionError::EmptySelection);
        }
        indices.sort_unstable();
        indices.dedup();
        Ok(Self { indices, n_atoms })
    }

    /// Parse an AtomSpecifier string against `topology`.
    pub fn from_string(spec: &str, topology: &Topology) -> Result<Self, SelectionError> {
        let n = topology.num_atoms();
        let spec = spec.trim();

        // --- extract not(...) clauses (permanent exclusion mask) ---
        let (not_mask, spec) = extract_parens_clause(spec, "not");
        let (minus_set, spec) = extract_parens_clause(&spec, "minus");

        let not_indices: HashSet<usize> = if not_mask.is_empty() {
            HashSet::new()
        } else {
            parse_union(&not_mask, topology)?.into_iter().collect()
        };
        let minus_indices: HashSet<usize> = if minus_set.is_empty() {
            HashSet::new()
        } else {
            parse_union(&minus_set, topology)?.into_iter().collect()
        };

        let mut result = parse_union(spec.trim(), topology)?;

        // Apply exclusions
        if !not_indices.is_empty() || !minus_indices.is_empty() {
            result.retain(|i| !not_indices.contains(i) && !minus_indices.contains(i));
        }

        if result.is_empty() {
            return Err(SelectionError::EmptySelection);
        }
        result.sort_unstable();
        result.dedup();
        Ok(Self {
            indices: result,
            n_atoms: n,
        })
    }

    /// Return the sorted slice of 0-based atom indices.
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }
    /// Number of atoms in this selection.
    pub fn len(&self) -> usize {
        self.indices.len()
    }
    /// Returns `true` if no atoms are selected.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
    /// Iterate over the 0-based atom indices in sorted order.
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.indices.iter().copied()
    }
}

// ─── Internal parsers ────────────────────────────────────────────────────────

/// Extract `keyword(...)` clauses from `spec`, returning (clause_contents, remaining_spec).
/// Handles only the first occurrence; multiple `not()` are rare in practice.
fn extract_parens_clause(spec: &str, keyword: &str) -> (String, String) {
    let pat = format!("{keyword}(");
    if let Some(start) = spec.find(&pat) {
        let after = &spec[start + pat.len()..];
        if let Some(end) = find_matching_paren(after) {
            let clause = after[..end].to_string();
            let remaining = format!("{}{}", &spec[..start], &after[end + 1..]);
            return (clause, remaining.trim().to_string());
        }
    }
    (String::new(), spec.to_string())
}

/// Find the index of the closing `)` matching the opening `(` that was already consumed.
fn find_matching_paren(s: &str) -> Option<usize> {
    let mut depth = 1usize;
    for (i, c) in s.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth == 0 {
                    return Some(i);
                }
            },
            _ => {},
        }
    }
    None
}

/// Parse a semicolon-separated union of terms.
fn parse_union(spec: &str, topology: &Topology) -> Result<Vec<usize>, SelectionError> {
    let mut result: HashSet<usize> = HashSet::new();
    for term in spec.split(';') {
        let term = term.trim();
        if term.is_empty() {
            continue;
        }
        let indices = parse_term(term, topology)?;
        result.extend(indices);
    }
    Ok(result.into_iter().collect())
}

/// Parse a single non-semicolon term.
fn parse_term(spec: &str, topology: &Topology) -> Result<Vec<usize>, SelectionError> {
    let n = topology.num_atoms();

    match spec.trim() {
        "all" => return Ok((0..n).collect()),
        "no" => return Ok(Vec::new()),
        _ => {},
    }

    // not(...) / minus(...) nested inside a term — treat as complement within this term
    if spec.starts_with("not(") || spec.starts_with("minus(") {
        let kw = if spec.starts_with("not(") {
            "not"
        } else {
            "minus"
        };
        let inner = &spec[kw.len() + 1..];
        if let Some(end) = find_matching_paren(inner) {
            let inner_spec = &inner[..end];
            let inner_ids = parse_union(inner_spec, topology)?;
            let inner_set: HashSet<usize> = inner_ids.into_iter().collect();
            return Ok((0..n).filter(|i| !inner_set.contains(i)).collect());
        }
    }

    if spec.contains(':') {
        parse_prefixed(spec, topology)
    } else {
        parse_global_range_list(spec, n)
    }
}

/// Parse `PREFIX:SELECTOR` form.
fn parse_prefixed(spec: &str, topology: &Topology) -> Result<Vec<usize>, SelectionError> {
    let (prefix, selector) = split_first_colon(spec)?;
    let mol_ranges = resolve_mol_prefix(prefix.trim(), topology)?;
    resolve_atom_selector(selector.trim(), &mol_ranges, topology)
}

/// Resolve the molecule prefix to a list of atom ranges.
/// Returns `Vec<(global_start, global_end_exclusive)>`.
fn resolve_mol_prefix(
    prefix: &str,
    topology: &Topology,
) -> Result<Vec<std::ops::Range<usize>>, SelectionError> {
    let n = topology.num_atoms();
    match prefix {
        "a" => {
            // Every solute molecule, never solvent (gromos++ `a` = molecules 1..numMolecules).
            // No `molecules` ranges means no solvent either: `solvate()` adds them.
            if topology.molecules.is_empty() {
                Ok(vec![0..n])
            } else if topology.instances.is_empty() {
                Ok(topology.molecules.clone())
            } else {
                Ok(molecules_with_role(topology, Role::Solute))
            }
        },
        "s" => {
            // s:N — Nth solvent molecule (role == Solvent), 1-based.
            // s:a — all solvent molecules.
            // Routes through the instancing role attribute; falls back to molecules[1..].
            let solvent = if topology.instances.is_empty() {
                topology.molecules.get(1..).unwrap_or_default().to_vec()
            } else {
                molecules_with_role(topology, Role::Solvent)
            };
            if solvent.is_empty() {
                return Err(SelectionError::NoSolvent);
            }
            Ok(solvent)
        },
        "m" => {
            // m:a — every solute molecule (role == Solute); `m:CA` the CA atoms in them.
            // Falls back to molecules[0] when instances are not populated.
            let solute = if topology.instances.is_empty() {
                topology.molecules.first().cloned().into_iter().collect()
            } else {
                molecules_with_role(topology, Role::Solute)
            };
            if solute.is_empty() {
                return Err(SelectionError::EmptySelection);
            }
            Ok(solute)
        },
        // ── Common mistake detection ──────────────────────────────────────────
        "r" => Err(SelectionError::ParseError(
            "'r:' is not a valid prefix in gromos-rs AtomSpecifier\n  \
                 To select by residue number use: 1:res(N:a) or a:res(N:a)\n  \
                 To select by residue name use:   1:res(NAME:a)"
                .into(),
        )),
        "mol" | "molecule" => Err(SelectionError::ParseError(
            "'mol:' / 'molecule:' is not valid — use the molecule number directly, e.g. '1:a'"
                .into(),
        )),
        "res" | "residue" => Err(SelectionError::ParseError(
            "'res:' / 'residue:' is not a valid prefix — use '1:res(N:a)' or 'a:res(NAME:a)'"
                .into(),
        )),
        other => {
            // Number or number range: "1", "1-3"
            // Detect bare atom names used as prefix (e.g. "CA:a" instead of "a:CA")
            if other.chars().all(|c| c.is_alphabetic()) && other.len() >= 2 && other.len() <= 4 {
                return Err(SelectionError::ParseError(format!(
                    "'{other}' looks like an atom name used as a prefix — \
                        did you mean 'a:{other}' (all molecules) or '1:{other}' (molecule 1)?"
                )));
            }
            let mol_nrs = parse_range_or_list(other)?; // 1-based
            let mut ranges = Vec::new();
            for mn in mol_nrs {
                if mn < 1 {
                    return Err(SelectionError::InvalidMolecule(mn));
                }
                let range = if topology.molecules.is_empty() {
                    if mn == 1 {
                        0..n
                    } else {
                        return Err(SelectionError::InvalidMolecule(mn));
                    }
                } else {
                    topology
                        .molecules
                        .get(mn - 1)
                        .cloned()
                        .ok_or(SelectionError::InvalidMolecule(mn))?
                };
                ranges.push(range);
            }
            Ok(ranges)
        },
    }
}

/// The `molecules` ranges whose atoms belong to an instance of `role`. Matched by atom, not by
/// index: the solute is one instance spanning all its molecules, so `instances[k]` is not
/// `molecules[k]` once the solute has more than one.
fn molecules_with_role(topology: &Topology, role: Role) -> Vec<std::ops::Range<usize>> {
    topology
        .molecules
        .iter()
        .filter(|r| topology.role_of_atom(r.start) == Some(role))
        .cloned()
        .collect()
}

/// Resolve an atom selector string against a set of molecule ranges.
fn resolve_atom_selector(
    selector: &str,
    mol_ranges: &[std::ops::Range<usize>],
    topology: &Topology,
) -> Result<Vec<usize>, SelectionError> {
    if selector == "a" {
        // All atoms in these molecules
        let indices: Vec<usize> = mol_ranges.iter().flat_map(|r| r.clone()).collect();
        return Ok(indices);
    }

    // res(residue_sel:atom_sel) — residue directive
    if let Some(inner) = selector.strip_prefix("res(") {
        if let Some(end) = find_matching_paren(inner) {
            let res_spec = &inner[..end];
            return resolve_res_directive(res_spec, mol_ranges, topology);
        }
        return Err(SelectionError::ParseError(format!(
            "unclosed res(): {selector}"
        )));
    }

    // Determine if selector is names (letters) or numbers (digits/dash/comma)
    let is_names = selector.chars().any(|c| c.is_alphabetic());

    if is_names {
        // Atom name list: "CA,N,C" or a mix like "CA"
        let names: HashSet<&str> = selector.split(',').map(|s| s.trim()).collect();
        let mut indices = Vec::new();
        for r in mol_ranges {
            for i in r.clone() {
                if let Some(name) = topology.atom_name(i) {
                    if names.contains(name) {
                        indices.push(i);
                    }
                }
            }
        }
        Ok(indices)
    } else {
        // Local 1-based atom number range within each molecule
        let local_nrs = parse_range_or_list(selector)?;
        let mut indices = Vec::new();
        for r in mol_ranges {
            let mol_size = r.len();
            for &ln in &local_nrs {
                if ln >= 1 && ln <= mol_size {
                    indices.push(r.start + ln - 1);
                }
            }
        }
        Ok(indices)
    }
}

/// Parse `residue_sel:atom_sel` inside `res(...)`.
fn resolve_res_directive(
    spec: &str,
    mol_ranges: &[std::ops::Range<usize>],
    topology: &Topology,
) -> Result<Vec<usize>, SelectionError> {
    let (res_sel, atom_sel) = split_first_colon(spec)?;
    let res_sel = res_sel.trim();
    let atom_sel = atom_sel.trim();

    // Build residue filter
    let res_by_name = res_sel.chars().any(|c| c.is_alphabetic());
    let residue_filter: Box<dyn Fn(usize) -> bool> = if res_by_name {
        let _names: HashSet<String> = res_sel.split(',').map(|s| s.trim().to_string()).collect();
        Box::new(move |i: usize| {
            // We'll check residue name per-atom below
            let _ = i;
            false // placeholder — filled below per-atom
        })
    } else {
        let residue_nrs: HashSet<usize> = parse_range_or_list(res_sel)?.into_iter().collect();
        Box::new(move |rn: usize| residue_nrs.contains(&rn))
    };

    // Atom filter
    let atom_all = atom_sel == "a";
    let atom_names: HashSet<&str> = if atom_all {
        HashSet::new()
    } else {
        atom_sel.split(',').map(|s| s.trim()).collect()
    };

    let mut indices = Vec::new();
    for r in mol_ranges {
        for i in r.clone() {
            // Residue check
            let residue_ok = if res_by_name {
                topology.residue_name(i).is_some_and(|rn| {
                    res_sel
                        .split(',')
                        .any(|s| s.trim().eq_ignore_ascii_case(rn))
                })
            } else {
                topology.residue_nr(i).is_some_and(&residue_filter)
            };
            if !residue_ok {
                continue;
            }

            // Atom check
            let atom_ok = atom_all
                || topology
                    .atom_name(i)
                    .is_some_and(|name| atom_names.contains(name));
            if atom_ok {
                indices.push(i);
            }
        }
    }
    Ok(indices)
}

/// Parse global 1-based atom number range list: "1-10,15,20-25".
fn parse_global_range_list(spec: &str, n_atoms: usize) -> Result<Vec<usize>, SelectionError> {
    let nrs = parse_range_or_list(spec)?;
    let mut indices = Vec::new();
    for nr in nrs {
        if nr < 1 {
            return Err(SelectionError::InvalidRange(format!(
                "atom number must be ≥1, got {nr}"
            )));
        }
        let idx = nr - 1;
        if idx < n_atoms {
            indices.push(idx);
        }
    }
    Ok(indices)
}

/// Split on the first `:` returning `(before, after)`.
fn split_first_colon(s: &str) -> Result<(&str, &str), SelectionError> {
    s.find(':')
        .map(|i| (&s[..i], &s[i + 1..]))
        .ok_or_else(|| SelectionError::ParseError(format!("expected ':' in '{s}'")))
}

/// Parse a comma-separated list of numbers and ranges ("1", "1-5", "1,3,5-8").
/// Returns the expanded set of values (1-based as given).
fn parse_range_or_list(spec: &str) -> Result<Vec<usize>, SelectionError> {
    let mut values = Vec::new();
    for part in spec.split(',') {
        let part = part.trim();
        if let Some(dash) = part.find('-') {
            let a: usize = part[..dash].trim().parse().map_err(|_| {
                SelectionError::ParseError(format!("not a number: '{}'", &part[..dash]))
            })?;
            let b: usize = part[dash + 1..].trim().parse().map_err(|_| {
                SelectionError::ParseError(format!("not a number: '{}'", &part[dash + 1..]))
            })?;
            if a > b {
                return Err(SelectionError::InvalidRange(format!("{a}-{b}")));
            }
            values.extend(a..=b);
        } else {
            let v: usize = part
                .parse()
                .map_err(|_| SelectionError::ParseError(format!("not a number: '{part}'")))?;
            values.push(v);
        }
    }
    Ok(values)
}

// ─── Tests ───────────────────────────────────────────────────────────────────
//
// All expected indices are verified against the gromos-rs `atominfo` binary
// (gromosPlsPls BUILD/programs/atominfo) on the aladip topology.
//
// Aladip has 12 atoms, 1 molecule, 3 residues:
//   Res 1 GLY:  atom 1=CB  2=C   3=O         (0-based: 0,1,2)
//   Res 2 ALA:  atom 4=N   5=H   6=CA 7=CB 8=C 9=O  (0-based: 3,4,5,6,7,8)
//   Res 3 GLY:  atom 10=N  11=H  12=CB        (0-based: 9,10,11)

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::{Atom, Topology};

    /// Construct a topology matching the aladip system used by gromos-rs atominfo reference runs.
    fn aladip_topo() -> Topology {
        let atoms: &[(&str, usize, &str)] = &[
            // (name, residue_nr, residue_name)
            ("CB", 1, "GLY"),
            ("C", 1, "GLY"),
            ("O", 1, "GLY"),
            ("N", 2, "ALA"),
            ("H", 2, "ALA"),
            ("CA", 2, "ALA"),
            ("CB", 2, "ALA"),
            ("C", 2, "ALA"),
            ("O", 2, "ALA"),
            ("N", 3, "GLY"),
            ("H", 3, "GLY"),
            ("CB", 3, "GLY"),
        ];
        let mut topo = Topology::new();
        for &(name, res_nr, res_name) in atoms {
            topo.moltypes[0].atoms.push(Atom {
                name: name.into(),
                residue_nr: res_nr,
                residue_name: res_name.into(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
            topo.iac.push(0);
            topo.mass.push(12.0);
            topo.charge.push(0.0);
        }
        // Populate molecules[0] = 0..12 so molecule_nr() and atom_name() work.
        topo.init_solute_moltype();
        topo
    }

    // ── m: facade (role-based solute selection) ───────────────────────────────

    #[test] // m:1 = first solute molecule (role == Solute) = all 12 aladip atoms
    fn test_m_solute_all() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("m:a", &t).unwrap();
        assert_eq!(s.len(), 12);
    }

    #[test] // m:CA = CA atoms in solute molecules (same as a:CA for single-solute system)
    fn test_m_solute_name() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("m:CA", &t).unwrap();
        assert_eq!(s.indices(), &[5]); // CA is atom 6 (0-based 5) in ALA residue
    }

    // ── basic ────────────────────────────────────────────────────────────────

    #[test] // atominfo: 1:a → atoms 1-12
    fn test_all_atoms() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:a", &t).unwrap();
        assert_eq!(s.indices(), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
    }

    #[test] // atominfo: keyword "all"
    fn test_all_keyword() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("all", &t).unwrap();
        assert_eq!(s.len(), 12);
    }

    // ── global atom number ranges ─────────────────────────────────────────────

    #[test] // atominfo: 1:1-5 → atoms 1-5
    fn test_mol_local_range() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:1-5", &t).unwrap();
        assert_eq!(s.indices(), &[0, 1, 2, 3, 4]);
    }

    #[test] // atominfo: 1:1 → atom 1
    fn test_mol_single_atom() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:1", &t).unwrap();
        assert_eq!(s.indices(), &[0]);
    }

    #[test] // atominfo: 1:5-8 → atoms 5-8 (0-based: 4,5,6,7)
    fn test_mol_mid_range() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:5-8", &t).unwrap();
        assert_eq!(s.indices(), &[4, 5, 6, 7]);
    }

    #[test] // atominfo: 1:2,5,8,11 → 0-based: 1,4,7,10
    fn test_mol_comma_list() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:2,5,8,11", &t).unwrap();
        assert_eq!(s.indices(), &[1, 4, 7, 10]);
    }

    // ── atom name selections ──────────────────────────────────────────────────

    #[test] // atominfo: a:CA → 1:6 (0-based: 5)
    fn test_all_mol_ca() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("a:CA", &t).unwrap();
        assert_eq!(s.indices(), &[5]);
    }

    #[test] // atominfo: a:N → 1:4,10 (0-based: 3,9)
    fn test_all_mol_n() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("a:N", &t).unwrap();
        assert_eq!(s.indices(), &[3, 9]);
    }

    #[test] // atominfo: a:H → 1:5,11 (0-based: 4,10)
    fn test_all_mol_h() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("a:H", &t).unwrap();
        assert_eq!(s.indices(), &[4, 10]);
    }

    #[test] // atominfo: 1:CA → 1:6 (0-based: 5)
    fn test_mol_single_name() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:CA", &t).unwrap();
        assert_eq!(s.indices(), &[5]);
    }

    #[test] // atominfo: 1:CA,N,C,O → 0-based: 1,2,3,5,7,8,9
    fn test_mol_name_list() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:CA,N,C,O", &t).unwrap();
        assert_eq!(s.indices(), &[1, 2, 3, 5, 7, 8, 9]);
    }

    #[test] // atominfo: 1:N,CA,C,O → same 7 atoms
    fn test_mol_name_list_reorder() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:N,CA,C,O", &t).unwrap();
        assert_eq!(s.indices(), &[1, 2, 3, 5, 7, 8, 9]);
    }

    // ── res() directive by residue number ────────────────────────────────────

    #[test] // atominfo: 1:res(1:a) → atoms 1-3 (0-based: 0,1,2)
    fn test_res_nr_all() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(1:a)", &t).unwrap();
        assert_eq!(s.indices(), &[0, 1, 2]);
    }

    #[test] // atominfo: 1:res(2:a) → atoms 4-9 (0-based: 3..=8)
    fn test_res_nr2_all() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(2:a)", &t).unwrap();
        assert_eq!(s.indices(), &[3, 4, 5, 6, 7, 8]);
    }

    #[test] // atominfo: 1:res(1-2:CA) → only ALA (res 2) has CA → 0-based: 5
    fn test_res_nr_range_ca() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(1-2:CA)", &t).unwrap();
        assert_eq!(s.indices(), &[5]);
    }

    #[test] // atominfo: 1:res(1,3:a) → res 1 + res 3 (0-based: 0,1,2,9,10,11)
    fn test_res_nr_comma_list() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(1,3:a)", &t).unwrap();
        assert_eq!(s.indices(), &[0, 1, 2, 9, 10, 11]);
    }

    #[test] // atominfo: 1:res(1,3:CB) → GLY CB → 0-based: 0,11
    fn test_res_nr_comma_name() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(1,3:CB)", &t).unwrap();
        assert_eq!(s.indices(), &[0, 11]);
    }

    #[test] // atominfo: 1:res(2:N,CA,C) → 0-based: 3,5,7
    fn test_res_nr_multi_name() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(2:N,CA,C)", &t).unwrap();
        assert_eq!(s.indices(), &[3, 5, 7]);
    }

    // ── res() directive by residue name ──────────────────────────────────────

    #[test] // atominfo: 1:res(ALA:a) → atoms 4-9 (0-based: 3..=8)
    fn test_res_name_all() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(ALA:a)", &t).unwrap();
        assert_eq!(s.indices(), &[3, 4, 5, 6, 7, 8]);
    }

    #[test] // atominfo: 1:res(GLY:CB) → 0-based: 0,11
    fn test_res_name_cb() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(GLY:CB)", &t).unwrap();
        assert_eq!(s.indices(), &[0, 11]);
    }

    #[test] // atominfo: 1:res(GLY:N,H,CB) → 0-based: 0,9,10,11
    fn test_res_name_multi_atom() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(GLY:N,H,CB)", &t).unwrap();
        assert_eq!(s.indices(), &[0, 9, 10, 11]);
    }

    #[test] // atominfo: 1:res(ILE,ALA:N) → only ALA has N → 0-based: 3
    fn test_res_name_multi_residue() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(ILE,ALA:N)", &t).unwrap();
        assert_eq!(s.indices(), &[3]);
    }

    #[test] // atominfo: 1:res(GLY,ALA:N) → 0-based: 3,9
    fn test_res_name_gly_ala_n() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:res(GLY,ALA:N)", &t).unwrap();
        assert_eq!(s.indices(), &[3, 9]);
    }

    // ── semicolon union ───────────────────────────────────────────────────────

    #[test] // atominfo: 1:1-3;1:10-12 → 0-based: 0,1,2,9,10,11
    fn test_semicolon_union() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:1-3;1:10-12", &t).unwrap();
        assert_eq!(s.indices(), &[0, 1, 2, 9, 10, 11]);
    }

    #[test] // atominfo: 1:1;1:3;1:5 → 0-based: 0,2,4
    fn test_semicolon_singles() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:1;1:3;1:5", &t).unwrap();
        assert_eq!(s.indices(), &[0, 2, 4]);
    }

    #[test] // atominfo: 1:1-6;1:7-12 → all 12
    fn test_semicolon_full_coverage() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:1-6;1:7-12", &t).unwrap();
        assert_eq!(s.len(), 12);
    }

    // ── not() / minus() exclusions ───────────────────────────────────────────

    #[test] // atominfo: not(1:1-6) 1:a → atoms 7-12 (0-based: 6..=11)
    fn test_not_range() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("not(1:1-6) 1:a", &t).unwrap();
        assert_eq!(s.indices(), &[6, 7, 8, 9, 10, 11]);
    }

    #[test] // atominfo: 1:a minus(1:1-3) → atoms 4-12 (0-based: 3..=11)
    fn test_minus_range() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:a minus(1:1-3)", &t).unwrap();
        assert_eq!(s.indices(), &[3, 4, 5, 6, 7, 8, 9, 10, 11]);
    }

    #[test] // atominfo: not(a:H) 1:a → all except H atoms (0-based: 4,10 excluded)
    fn test_not_name() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("not(a:H) 1:a", &t).unwrap();
        assert_eq!(s.indices(), &[0, 1, 2, 3, 5, 6, 7, 8, 9, 11]);
    }

    #[test] // atominfo: 1:a minus(a:H) → same result as not(a:H) 1:a
    fn test_minus_name() {
        let t = aladip_topo();
        let s = AtomSelection::from_string("1:a minus(a:H)", &t).unwrap();
        assert_eq!(s.indices(), &[0, 1, 2, 3, 5, 6, 7, 8, 9, 11]);
    }

    // ── multi-molecule solute ─────────────────────────────────────────────────
    //
    // water_216_box.topo (gromos-md reference): 216 three-atom SOL molecules (OW, HW1, HW2)
    // all in the solute — SOLUTEMOLECULES lists 216 ranges, SOLVENTATOM is empty. Built here
    // the way `gromos_io::topology::build_topology` builds it: one solute moltype, one
    // `molecules` range per SOLUTEMOLECULES entry.

    fn water_216_solute_topo() -> Topology {
        let mut topo = Topology::new();
        for mol in 0..216 {
            for (name, charge) in [("OW", -0.82), ("HW1", 0.41), ("HW2", 0.41)] {
                topo.moltypes[0].atoms.push(Atom {
                    name: name.into(),
                    residue_nr: mol + 1,
                    residue_name: "SOL".into(),
                    iac: 0,
                    mass: 1.0,
                    charge,
                    is_perturbed: false,
                    is_polarisable: false,
                    is_coarse_grained: false,
                });
                topo.iac.push(0);
                topo.mass.push(1.0);
                topo.charge.push(charge);
            }
            topo.molecules.push(3 * mol..3 * mol + 3);
        }
        topo.init_solute_moltype();
        topo
    }

    #[test] // a:OW = the OW of every molecule, not only of the first
    fn test_multi_molecule_solute_name_in_all_molecules() {
        let t = water_216_solute_topo();
        let s = AtomSelection::from_string("a:OW", &t).unwrap();
        assert_eq!(s.len(), 216);
        assert!(s.iter().all(|i| i % 3 == 0));
    }

    #[test] // s: on a topology without solvent is an error that says so
    fn test_multi_molecule_solute_has_no_solvent() {
        let t = water_216_solute_topo();
        let err = AtomSelection::from_string("s:OW", &t).unwrap_err();
        assert!(matches!(err, SelectionError::NoSolvent), "{err}");
        assert!(err.to_string().contains("no solvent"), "{err}");
    }

    #[test] // 1-216:1 = atom 1 of each molecule
    fn test_multi_molecule_solute_local_atom_in_range() {
        let t = water_216_solute_topo();
        let s = AtomSelection::from_string("1-216:1", &t).unwrap();
        assert_eq!(s.len(), 216);
        assert_eq!(s.indices()[..3], [0, 3, 6]);
    }

    #[test] // 2:1 = atom 1 of molecule 2, and it is named
    fn test_multi_molecule_solute_single_atom() {
        let t = water_216_solute_topo();
        let s = AtomSelection::from_string("2:1", &t).unwrap();
        assert_eq!(s.indices(), &[3]);
        assert_eq!(t.atom_name(3), Some("OW"));
        assert_eq!(t.residue_nr(3), Some(2));
    }

    /// `t` with `nsm` SPC waters (OW, HW1, HW2) appended as solvent.
    fn solvated(mut t: Topology, nsm: usize) -> Topology {
        use crate::topology::SolventAtomTemplate;
        t.solvent_atom_template = ["OW", "HW1", "HW2"]
            .iter()
            .map(|n| SolventAtomTemplate {
                iac: 0,
                name: (*n).into(),
                mass: 1.0,
                charge: 0.0,
            })
            .collect();
        t.solvate(nsm);
        t
    }

    #[test] // m:a = every solute molecule; s:a = the solvent once the topology is solvated
    fn test_multi_molecule_solute_with_solvent() {
        let t = solvated(water_216_solute_topo(), 2);
        assert_eq!(AtomSelection::from_string("m:a", &t).unwrap().len(), 648);
        let solvent = AtomSelection::from_string("s:a", &t).unwrap();
        assert_eq!(solvent.indices(), &[648, 649, 650, 651, 652, 653]);
        assert_eq!(AtomSelection::from_string("s:OW", &t).unwrap().len(), 2);
    }

    // ── a: is solute only (gromos++ `a` = molecules 1..numMolecules) ─────────

    #[test] // aladip + 20 waters: OW exists only in the solvent, so a:OW is empty
    fn test_a_name_only_in_solvent_selects_nothing() {
        let t = solvated(aladip_topo(), 20);
        assert_eq!(t.num_atoms(), 72);
        assert!(matches!(
            AtomSelection::from_string("a:OW", &t),
            Err(SelectionError::EmptySelection)
        ));
        let ow = AtomSelection::from_string("s:OW", &t).unwrap();
        assert_eq!(ow.len(), 20);
        assert!(ow.iter().all(|i| i >= 12 && (i - 12) % 3 == 0));
    }

    #[test] // OW in solute and solvent: a:OW takes the 216 solute oxygens only
    fn test_a_name_in_both_roles_selects_solute_only() {
        let t = solvated(water_216_solute_topo(), 2);
        let s = AtomSelection::from_string("a:OW", &t).unwrap();
        assert_eq!(s.len(), 216);
        assert!(s.iter().all(|i| i < 648 && i % 3 == 0));
    }

    #[test] // a:a = every atom of every solute molecule, no solvent
    fn test_a_all_excludes_solvent() {
        let t = solvated(aladip_topo(), 20);
        let s = AtomSelection::from_string("a:a", &t).unwrap();
        assert_eq!(s.indices(), (0..12).collect::<Vec<_>>());
        let t = solvated(water_216_solute_topo(), 2);
        assert_eq!(AtomSelection::from_string("a:a", &t).unwrap().len(), 648);
    }

    #[test] // local numbers and res() under a: stay inside the solute too
    fn test_a_local_range_and_res_exclude_solvent() {
        let t = solvated(aladip_topo(), 20);
        let s = AtomSelection::from_string("a:1-3", &t).unwrap();
        assert_eq!(s.indices(), &[0, 1, 2]);
        let s = AtomSelection::from_string("a:res(1:a)", &t).unwrap();
        assert_eq!(s.indices(), &[0, 1, 2]);
        assert_eq!(
            AtomSelection::from_string("a:N", &t).unwrap().indices(),
            &[3, 9]
        );
    }
}
