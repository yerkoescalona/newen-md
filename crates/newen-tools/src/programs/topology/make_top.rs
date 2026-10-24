//! make_top - Build GROMOS topology from building blocks
//!
//! Usage: make_top @build <mtb> @param <ifp> @seq <residues...> [@solv <solvent>]
//!    or: make_top @f <argfile>
//!    or: make_top --build <mtb> --param <ifp> --seq <residues...> [--solv <solvent>]
//!
//! Assembles a complete molecular topology from MTB building blocks
//! and IFP force field parameters, following the GROMOS convention.
//! Supports both GROMOS @-prefix and standard --prefix argument styles.

use crate::programs::{io_err, parse_clap};
use clap::Parser;
use newen_io::args::{Program, Status};
use newen_io::ifp::{self, ForceFieldParameters};
use newen_io::mtb::{self, BbAtom, BbEnd, BbSolute, BlockRef, BuildingBlocks};
use newen_io::topology::{parse_topology, ParsedTopology};
use std::collections::BTreeSet;
use std::fmt::Write as FmtWrite;
use std::io::Write;

pub const USAGE: &str = "make_top - Build a GROMOS topology from building blocks\n\nUsage: make_top @build <mtb> [<mtb> ...] @param <ifp> @seq <residue> [<residue> ...] [@solv <name>]\n\nArguments:\n  @build  Molecular building-block (MTB) file(s)\n  @param  Force-field parameter (IFP) file\n  @seq    Residue sequence to build\n  @solv   Solvent building-block name";

pub const KNOWN: &[&str] = &["build", "param", "seq", "solv"];

pub const PROGRAM: Program = Program {
    name: "make_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// A built topology: the file text and the same topology as data.
#[derive(Debug, Clone)]
pub struct MakeTopResult {
    /// The topology, read back from `text` — what `Topology` wraps in Python.
    pub topology: ParsedTopology,
    /// The topology file, in gromos++ `OutTopology` layout.
    pub text: String,
}

/// Build a GROMOS topology from building blocks and force-field parameters.
///
/// Supports GROMOS @-prefix args: make_top @build 54a7.mtb @param 54a7.ifp @seq NH3+ ALA COO- @solv H2O
#[derive(Parser)]
#[command(name = "make_top", version, about)]
struct Args {
    /// Building block (MTB) files; later files add blocks to the first
    #[arg(long, num_args = 1..)]
    build: Vec<String>,

    /// Force field parameter (IFP) file
    #[arg(long)]
    param: String,

    /// Residue sequence (end groups + building blocks)
    #[arg(long, num_args = 1..)]
    seq: Vec<String>,

    /// Solvent building block name (e.g. H2O)
    #[arg(long)]
    solv: Option<String>,
}

// ── The linear topology gromos++ assembles (`gcore::LinearTopology`) ──────────────────────
//
// Atom indices are 0-based, as in gromos++ after `InBuildingBlock` subtracted one from every
// index in the MTB file (`--i`); negative values are the "previous residue" references of the
// building blocks. IAC −1 marks an atom that a begin group replaces (removed afterwards), IAC −2
// an end-group atom whose charge is only added to the last atom of that name.

#[derive(Clone, Debug)]
struct LtAtom {
    name: String,
    iac: i32,
    mass: f64,
    charge: f64,
    chargegroup: usize,
    exclusions: BTreeSet<i64>,
    one_four: BTreeSet<usize>,
    res: usize,
}

/// gromos++ `std::set<Bond>` order: (i, j, type).
type Bond = (i64, i64, i64);
/// `std::set<Angle>` order: (j, i, k, type) — stored in that key order.
type AngleKey = (i64, i64, i64, i64);
/// `std::set<Improper>` order: (i, j, k, l, type).
type Improper = (i64, i64, i64, i64, i64);
/// `std::set<Dihedral>` order: (j, k, i, l, type) — stored in that key order.
type DihedralKey = (i64, i64, i64, i64, i64);

fn angle_key(i: i64, j: i64, k: i64, t: i64) -> AngleKey {
    (j, i, k, t)
}
fn angle_atoms(a: &AngleKey) -> (i64, i64, i64, i64) {
    (a.1, a.0, a.2, a.3)
}
fn dihedral_key(i: i64, j: i64, k: i64, l: i64, t: i64) -> DihedralKey {
    (j, k, i, l, t)
}
fn dihedral_atoms(d: &DihedralKey) -> (i64, i64, i64, i64, i64) {
    (d.2, d.0, d.1, d.3, d.4)
}

#[derive(Default)]
struct Lt {
    atoms: Vec<LtAtom>,
    res_names: Vec<String>,
    bonds: BTreeSet<Bond>,
    angles: BTreeSet<AngleKey>,
    impropers: BTreeSet<Improper>,
    dihedrals: BTreeSet<DihedralKey>,
}

/// A building block as gromos++ sees it: 0-based indices, resolved masses.
struct Block {
    name: String,
    rep: i32,
    atoms: Vec<LtAtom>,
    /// per-atom exclusions, 0-based relative to the block
    excl: Vec<Vec<i64>>,
    /// preceding exclusions: (atoms before the block, 0-based relative) with their lists
    pexcl: Vec<Vec<i64>>,
    bonds: Vec<(i64, i64, i64)>,
    angles: Vec<(i64, i64, i64, i64)>,
    impropers: Vec<(i64, i64, i64, i64, i64)>,
    dihedrals: Vec<(i64, i64, i64, i64, i64)>,
}

fn bb_atom(a: &BbAtom, ffp: &ForceFieldParameters) -> LtAtom {
    LtAtom {
        name: a.name.clone(),
        iac: a.iac - 1,
        mass: ffp.get_mass(a.mass_code).unwrap_or(0.0),
        charge: a.charge,
        chargegroup: a.chargegroup,
        exclusions: BTreeSet::new(),
        one_four: BTreeSet::new(),
        res: 0,
    }
}

fn from_solute(bb: &BbSolute, ffp: &ForceFieldParameters) -> Block {
    Block {
        name: bb.name.clone(),
        rep: 0,
        atoms: bb.atoms.iter().map(|a| bb_atom(a, ffp)).collect(),
        excl: bb
            .atoms
            .iter()
            .map(|a| a.exclusions.iter().map(|&e| e as i64 - 1).collect())
            .collect(),
        pexcl: bb
            .preceding_exclusions
            .iter()
            .map(|p| p.excluded_atoms.iter().map(|&e| e as i64 - 1).collect())
            .collect(),
        bonds: bb
            .bonds
            .iter()
            .map(|&(i, j, t)| (i as i64 - 1, j as i64 - 1, t as i64))
            .collect(),
        angles: bb
            .angles
            .iter()
            .map(|&(i, j, k, t)| (i as i64 - 1, j as i64 - 1, k as i64 - 1, t as i64))
            .collect(),
        impropers: bb
            .improper_dihedrals
            .iter()
            .map(|&(i, j, k, l, t)| {
                (
                    i as i64 - 1,
                    j as i64 - 1,
                    k as i64 - 1,
                    l as i64 - 1,
                    t as i64,
                )
            })
            .collect(),
        dihedrals: bb
            .proper_dihedrals
            .iter()
            .map(|&(i, j, k, l, t)| {
                (
                    i as i64 - 1,
                    j as i64 - 1,
                    k as i64 - 1,
                    l as i64 - 1,
                    t as i64,
                )
            })
            .collect(),
    }
}

fn from_end(bb: &BbEnd, ffp: &ForceFieldParameters) -> Block {
    Block {
        name: bb.name.clone(),
        rep: bb.n_replace,
        atoms: bb.atoms.iter().map(|a| bb_atom(a, ffp)).collect(),
        excl: bb
            .atoms
            .iter()
            .map(|a| a.exclusions.iter().map(|&e| e as i64 - 1).collect())
            .collect(),
        pexcl: Vec::new(),
        bonds: bb
            .bonds
            .iter()
            .map(|&(i, j, t)| (i as i64 - 1, j as i64 - 1, t as i64))
            .collect(),
        angles: bb
            .angles
            .iter()
            .map(|&(i, j, k, t)| (i as i64 - 1, j as i64 - 1, k as i64 - 1, t as i64))
            .collect(),
        impropers: bb
            .improper_dihedrals
            .iter()
            .map(|&(i, j, k, l, t)| {
                (
                    i as i64 - 1,
                    j as i64 - 1,
                    k as i64 - 1,
                    l as i64 - 1,
                    t as i64,
                )
            })
            .collect(),
        dihedrals: bb
            .proper_dihedrals
            .iter()
            .map(|&(i, j, k, l, t)| {
                (
                    i as i64 - 1,
                    j as i64 - 1,
                    k as i64 - 1,
                    l as i64 - 1,
                    t as i64,
                )
            })
            .collect(),
    }
}

impl Lt {
    fn set_res_name(&mut self, resnum: usize, name: &str) {
        if self.res_names.len() <= resnum {
            self.res_names.resize(resnum + 1, String::new());
        }
        self.res_names[resnum] = name.to_string();
    }

    fn has_bond(&self, i: i64, j: i64) -> bool {
        self.bonds.iter().any(|b| b.0 == i && b.1 == j)
    }
    fn has_angle(&self, i: i64, j: i64, k: i64) -> bool {
        self.angles
            .iter()
            .any(|a| angle_atoms(a).0 == i && a.0 == j && a.2 == k)
    }
    fn has_improper(&self, i: i64, j: i64, k: i64, l: i64) -> bool {
        self.impropers
            .iter()
            .any(|p| p.0 == i && p.1 == j && p.2 == k && p.3 == l)
    }

    /// gromos++ `bondedAtoms`: atoms bonded (as the first partner) to one of `atoms`, at an
    /// index not above `offset` and not in `atoms` — the previous residue's link candidates.
    fn bonded_atoms(&self, atoms: &BTreeSet<i64>, offset: i64) -> BTreeSet<i64> {
        let min = atoms.iter().copied().fold(offset, i64::min);
        let mut out = BTreeSet::new();
        for &a in atoms {
            for b in &self.bonds {
                if b.1 == a && b.0 <= min && !atoms.contains(&b.0) {
                    out.insert(b.0);
                }
            }
        }
        out
    }

    /// gromos++ `LinearTopology::removeAtoms`: drop atoms with a negative IAC (replaced by a
    /// begin group) and renumber everything.
    fn remove_atoms(&mut self) {
        let removed: BTreeSet<usize> = self
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.iac < 0)
            .map(|(i, _)| i)
            .collect();
        if removed.is_empty() {
            return;
        }
        let mut ren = vec![-1i64; self.atoms.len()];
        let mut corr = 0usize;
        for (i, r) in ren.iter_mut().enumerate() {
            if removed.contains(&i) {
                corr += 1;
            } else {
                *r = (i - corr) as i64;
            }
        }
        let map = |x: i64| -> Option<i64> {
            if x < 0 || x as usize >= ren.len() {
                Some(x) // out of range: keep (later checks drop it)
            } else if ren[x as usize] < 0 {
                None
            } else {
                Some(ren[x as usize])
            }
        };
        let mut atoms = Vec::new();
        for (i, a) in self.atoms.iter().enumerate() {
            if removed.contains(&i) {
                continue;
            }
            let mut a = a.clone();
            a.exclusions = a.exclusions.iter().filter_map(|&e| map(e)).collect();
            atoms.push(a);
        }
        self.atoms = atoms;
        self.bonds = self
            .bonds
            .iter()
            .filter_map(|b| Some((map(b.0)?, map(b.1)?, b.2)))
            .collect();
        self.angles = self
            .angles
            .iter()
            .filter_map(|a| Some((map(a.0)?, map(a.1)?, map(a.2)?, a.3)))
            .collect();
        self.impropers = self
            .impropers
            .iter()
            .filter_map(|p| Some((map(p.0)?, map(p.1)?, map(p.2)?, map(p.3)?, p.4)))
            .collect();
        self.dihedrals = self
            .dihedrals
            .iter()
            .filter_map(|d| Some((map(d.0)?, map(d.1)?, map(d.2)?, map(d.3)?, d.4)))
            .collect();
    }

    /// gromos++ `LinearTopology::get14s`: atoms exactly three bonds apart (and not closer)
    /// that are not excluded are 1-4 partners.
    fn get14s(&mut self) {
        let n = self.atoms.len();
        let mut neigh: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
        for b in &self.bonds {
            if b.0 >= 0 && b.1 >= 0 && (b.0 as usize) < n && (b.1 as usize) < n {
                neigh[b.0 as usize].insert(b.1 as usize);
                neigh[b.1 as usize].insert(b.0 as usize);
            }
        }
        for i in 0..n {
            let first = neigh[i].clone();
            let mut second = BTreeSet::new();
            for &f in &first {
                second.extend(neigh[f].iter().copied());
            }
            let mut third = BTreeSet::new();
            for &s in &second {
                third.extend(neigh[s].iter().copied());
            }
            let mut e = BTreeSet::new();
            for &t in &third {
                if i < t
                    && !first.contains(&t)
                    && !second.contains(&t)
                    && !self.atoms[i].exclusions.contains(&(t as i64))
                {
                    e.insert(t);
                }
            }
            self.atoms[i].one_four = e;
        }
    }
}

/// gromos++ `addSolute`.
fn add_solute(
    lt: &mut Lt,
    bb: &Block,
    resnum: usize,
    resname: &str,
    rep: i64,
    nn: i64,
) -> Result<(), String> {
    let na = lt.atoms.len() as i64;
    let strt = na - rep;
    let beg = if strt < 0 { -strt } else { 0 };
    for i in beg..rep {
        let e: BTreeSet<i64> = bb.excl[i as usize].iter().map(|x| x + strt).collect();
        lt.atoms[(strt + i) as usize].exclusions = e;
    }
    if rep == 0 {
        let pexl = strt - bb.pexcl.len() as i64;
        if pexl < 0 {
            return Err(format!(
                "{resname}: preceding exclusions, but no preceding atoms"
            ));
        }
        for (i, list) in bb.pexcl.iter().enumerate() {
            let e: BTreeSet<i64> = list.iter().map(|x| x + strt).collect();
            lt.atoms[pexl as usize + i].exclusions = e;
        }
    }
    for i in rep as usize..bb.atoms.len() {
        let mut a = bb.atoms[i].clone();
        a.exclusions = bb.excl[i].iter().map(|x| x + strt).collect();
        a.res = resnum;
        lt.atoms.push(a);
    }
    lt.set_res_name(resnum, resname);
    let offset = strt;
    for &(i, j, t) in &bb.bonds {
        let b = (i + offset, j + offset, t);
        if !(rep > 0 && lt.has_bond(b.0, b.1)) {
            lt.bonds.insert(b);
        }
    }
    for &(i, j, k, t) in &bb.angles {
        let (i, j, k) = (i + offset, j + offset, k + offset);
        let found = rep > 0 && lt.has_angle(i, j, k);
        if (i >= na || j >= na || k >= na) && !found {
            lt.angles.insert(angle_key(i, j, k, t));
        }
    }
    for &(i, j, k, l, t) in &bb.impropers {
        let (i, j, k, l) = (i + offset, j + offset, k + offset, l + offset);
        let found = rep > 0 && lt.has_improper(i, j, k, l);
        if (i >= na || j >= na || k >= na || l >= na) && !found {
            lt.impropers.insert((i, j, k, l, t));
        }
    }
    if rep == 0 {
        for &(d0, d1, d2, d3, t) in &bb.dihedrals {
            let (mut corr0, corr1, corr2, mut corr3) = (offset, offset, offset, offset);
            if d3 == -2 {
                corr3 = 0;
            }
            if d0 == -3 {
                corr0 = 0;
                for b in &lt.bonds {
                    if b.1 == d1 + offset {
                        corr0 = b.0 + 3;
                    }
                }
            }
            lt.dihedrals.insert(dihedral_key(
                d0 + corr0,
                d1 + corr1,
                d2 + corr2,
                d3 + corr3,
                t,
            ));
        }
    } else if rep > 0 {
        let mut old: Vec<(i64, i64, i64, i64, i64)> =
            lt.dihedrals.iter().map(dihedral_atoms).collect();
        lt.dihedrals.clear();
        for &(d0, d1, d2, d3, t) in &bb.dihedrals {
            let mut corr = offset;
            if d0 == -3 {
                for b in &lt.bonds {
                    if b.1 == d1 + offset {
                        corr = b.0 + 3;
                    }
                }
            }
            if d0 == -2 {
                for b in &lt.bonds {
                    if b.1 == d1 + offset {
                        corr = b.0 + 2;
                    }
                }
            }
            let b = (d0 + corr, d1 + offset, d2 + offset, d3 + offset, t);
            let mut add = true;
            for o in old.iter_mut() {
                if o.0 == b.0 && o.1 == b.1 && o.2 == b.2 {
                    o.3 = b.3;
                    add = false;
                }
            }
            if b.0 < nn || b.1 < nn || b.2 < nn || b.3 < nn {
                add = false;
            }
            if b.1 == -6 + offset && b.2 == -5 + offset && b.3 == -4 + offset {
                add = true;
            }
            if b.0 == -5 + offset && b.1 == -6 + offset {
                add = true;
            }
            if b.3 == -6 + offset {
                add = true;
            }
            if add {
                lt.dihedrals.insert(dihedral_key(b.0, b.1, b.2, b.3, b.4));
            }
        }
        for o in old {
            lt.dihedrals.insert(dihedral_key(o.0, o.1, o.2, o.3, o.4));
        }
    }
    Ok(())
}

/// gromos++ `addBegin`: the atoms of a begin group, exclusions shifted by the atoms so far.
fn add_begin(lt: &mut Lt, bb: &Block, resnum: usize) -> i64 {
    let na = lt.atoms.len() as i64;
    for (i, a) in bb.atoms.iter().enumerate() {
        let mut a = a.clone();
        a.exclusions = bb.excl[i].iter().map(|x| x + na).collect();
        a.res = resnum;
        lt.atoms.push(a);
    }
    bb.rep as i64
}

/// gromos++ `addEnd`: replace the last |rep| atoms of the chain by the end group.
fn add_end(lt: &mut Lt, bb: &Block, resnum: usize) -> Result<(), String> {
    let strt = lt.atoms.len() as i64 + bb.rep as i64;
    let mut search = Vec::new();
    for (i, a) in bb.atoms.iter().enumerate() {
        if a.iac < -2 {
            return Err(format!("IAC < -1 found in end group {}", bb.name));
        }
        if a.iac == -2 {
            search.push(i);
            eprintln!(
                "# WARNING: for atom {} in MTBUILDBLEND group {} only the CHARGE is transferred to the last atom with this name in the chain.",
                a.name, bb.name
            );
        }
    }
    let n_pop = (-bb.rep) as usize - search.len();
    for _ in 0..n_pop {
        lt.atoms.pop();
    }
    for &s in &search {
        let mut j = lt.atoms.len() - 1;
        while lt.atoms[j].name != bb.atoms[s].name {
            j -= 1;
        }
        lt.atoms[j].charge += bb.atoms[s].charge;
    }
    for (i, a) in bb.atoms.iter().enumerate() {
        if a.iac >= -1 {
            let mut a = a.clone();
            a.exclusions = bb.excl[i].iter().map(|x| x + strt).collect();
            a.res = resnum;
            lt.atoms.push(a);
        }
    }
    Ok(())
}

/// gromos++ `addCovEnd`: the covalent terms of an end group; for a terminal group the
/// references to the previous residue are resolved through the bond graph.
fn add_cov_end(lt: &mut Lt, bb: &Block, offset: i64) -> Result<(), String> {
    let name = &bb.name;
    for &(i, j, t) in &bb.bonds {
        let mut b = (i + offset, j + offset, t);
        if bb.rep < 0 {
            if i < 0 || bb.atoms[i as usize].iac == -2 {
                let atoms: BTreeSet<i64> = [j + offset].into_iter().collect();
                let cand = lt.bonded_atoms(&atoms, offset);
                if cand.is_empty() {
                    return Err(format!(
                        "bond to the previous residue in {name} is not found"
                    ));
                }
                if cand.len() != 1 {
                    return Err(format!(
                        "bond to the previous residue in {name} is ambiguous"
                    ));
                }
                b.0 = *cand.iter().next().unwrap();
            }
            lt.bonds.retain(|x| !(x.0 == b.0 && x.1 == b.1));
            lt.bonds.insert(b);
        } else {
            lt.bonds.insert(b);
        }
    }
    for &(i, j, k, t) in &bb.angles {
        let mut a = (i + offset, j + offset, k + offset, t);
        if bb.rep < 0 && (i < 0 || bb.atoms[i as usize].iac == -2) {
            let atoms: BTreeSet<i64> = [j + offset].into_iter().collect();
            let cand = lt.bonded_atoms(&atoms, offset);
            if cand.is_empty() {
                return Err(format!(
                    "angle to the previous residue in {name} is not found"
                ));
            }
            if cand.len() != 1 {
                return Err(format!(
                    "angle to the previous residue in {name} is ambiguous"
                ));
            }
            a.0 = *cand.iter().next().unwrap();
        }
        lt.angles.retain(|x| {
            let (xi, xj, xk, _) = angle_atoms(x);
            !(xi == a.0 && xj == a.1 && xk == a.2)
        });
        lt.angles.insert(angle_key(a.0, a.1, a.2, a.3));
    }
    for &(i, j, k, l, t) in &bb.impropers {
        let mut list = [i + offset, j + offset, k + offset, l + offset];
        if bb.rep < 0 {
            let raw = [i, j, k, l];
            let mut atoms = BTreeSet::new();
            let mut negs = Vec::new();
            for (n, &r) in raw.iter().enumerate() {
                if r >= 0 && bb.atoms[r as usize].iac >= -1 {
                    atoms.insert(r + offset);
                } else {
                    negs.push(n);
                }
            }
            if atoms.len() < 4 {
                let cand = lt.bonded_atoms(&atoms, offset);
                if cand.is_empty() {
                    return Err(format!(
                        "improper to the previous residue in {name} is not found"
                    ));
                }
                if cand.len() != 1 {
                    return Err(format!(
                        "improper to the previous residue in {name} is ambiguous"
                    ));
                }
                if negs.len() > 1 {
                    return Err(format!("multiple negative values in {name} impropers"));
                }
                list[negs[0]] = *cand.iter().next().unwrap();
            }
        }
        let set: BTreeSet<i64> = list.iter().copied().collect();
        lt.impropers.retain(|x| {
            !(set.contains(&x.0) && set.contains(&x.1) && set.contains(&x.2) && set.contains(&x.3))
        });
        lt.impropers.insert((list[0], list[1], list[2], list[3], t));
    }
    let mut added: BTreeSet<DihedralKey> = BTreeSet::new();
    for &(i, j, k, l, t) in &bb.dihedrals {
        let mut list = [i + offset, j + offset, k + offset, l + offset];
        if bb.rep < 0 {
            let raw = [i, j, k, l];
            let mut atoms = BTreeSet::new();
            let mut negs = Vec::new();
            for (n, &r) in raw.iter().enumerate() {
                if r >= 0 && bb.atoms[r as usize].iac >= -1 {
                    atoms.insert(r + offset);
                } else {
                    negs.push(n);
                }
            }
            if atoms.len() < 4 {
                let cand = lt.bonded_atoms(&atoms, offset);
                if cand.is_empty() {
                    return Err(format!(
                        "dihedral to the previous residue in {name} is not found"
                    ));
                }
                if cand.len() != 1 {
                    eprintln!("Warning: dihedral to the previous residue in {name} is ambiguous");
                }
                if negs.len() > 1 {
                    return Err(format!("multiple negative values in {name} dihedrals"));
                }
                list[negs[0]] = *cand.iter().next().unwrap();
            }
        }
        let key = dihedral_key(list[0], list[1], list[2], list[3], t);
        lt.dihedrals.retain(|x| {
            let (xi, xj, xk, xl, _) = dihedral_atoms(x);
            !(xi == list[0] && xj == list[1] && xk == list[2] && xl == list[3]) || added.contains(x)
        });
        lt.dihedrals.insert(key);
        added.insert(key);
    }
    Ok(())
}

fn is_h(mass: f64) -> bool {
    (mass - 1.008).abs() < 1e-9
}

/// Build a topology from building blocks, gromos++ `make_top`'s residue loop.
///
/// `seq` is the residue sequence; `solvent_name` names a solvent building block to append. The
/// result carries both the file text and the topology parsed back from it.
pub fn make_top(
    bb: &BuildingBlocks,
    ffp: &ForceFieldParameters,
    seq: &[String],
    solvent_name: Option<&str>,
) -> Result<MakeTopResult, String> {
    // A missing solvent block is fatal, and gromos++ only notices it half-way through writing.
    // Check before a single line of output exists, so a failure never yields half a topology.
    if let Some(sname) = solvent_name {
        if bb.find_solvent(sname).is_none() {
            return Err(format!("Error: solvent building block '{sname}' not found"));
        }
    }

    // gromos++ make_top.cc: the residue loop with its status machine
    //   0 solute, 1 begin group, 2 first solute after a begin group, 3 end group
    let mut lt = Lt::default();
    let mut status = 0;
    let mut resnum = 0usize;
    let mut repforward = 0i64;
    let mut first_atom = 0i64;
    for (seq_idx, res_name) in seq.iter().enumerate() {
        let block = match bb.find_block(res_name) {
            Some(BlockRef::End(e)) => {
                status = if e.n_replace < 0 { 3 } else { 1 };
                from_end(e, ffp)
            },
            Some(BlockRef::Solute(s)) => {
                status = if status == 1 { 2 } else { 0 };
                from_solute(s, ffp)
            },
            None => {
                return Err(format!(
                "Error: building block '{}' not found in the MTB file(s) (sequence position {})",
                res_name,
                seq_idx + 1
            ))
            },
        };
        match status {
            0 => {
                add_solute(&mut lt, &block, resnum, res_name, 0, first_atom)
                    .map_err(|e| format!("Error: {e}"))?;
                resnum += 1;
            },
            1 => {
                repforward = add_begin(&mut lt, &block, resnum);
                first_atom = lt.atoms.len() as i64 - block.atoms.len() as i64;
                add_cov_end(&mut lt, &block, first_atom).map_err(|e| format!("Error: {e}"))?;
            },
            2 => {
                add_solute(&mut lt, &block, resnum, res_name, repforward, first_atom)
                    .map_err(|e| format!("Error: {e}"))?;
                lt.remove_atoms();
                resnum += 1;
            },
            _ => {
                if resnum == 0 {
                    return Err(format!(
                        "Error: end group {res_name} cannot be at the start of the sequence"
                    ));
                }
                resnum -= 1;
                add_end(&mut lt, &block, resnum).map_err(|e| format!("Error: {e}"))?;
                let off = lt.atoms.len() as i64 - block.atoms.len() as i64;
                add_cov_end(&mut lt, &block, off).map_err(|e| format!("Error: {e}"))?;
                resnum += 1;
            },
        }
    }
    lt.get14s();
    let n = lt.atoms.len() as i64;
    for (i, a) in lt.atoms.iter_mut().enumerate() {
        let before = a.exclusions.len();
        a.exclusions.retain(|&e| e >= 0 && e < n);
        if a.exclusions.len() != before {
            eprintln!(
                "WARNING: exclusion skipped since it is not within the solute (atom {})",
                i + 1
            );
        }
    }
    let text = write_topology(&lt, bb, ffp, solvent_name)?;
    let topology = parse_topology(text.as_bytes())
        .map_err(|e| format!("make_top produced a topology it cannot read back: {e}"))?;
    Ok(MakeTopResult { topology, text })
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let Some(args) = parse_clap::<Args>(argv, out)? else {
        return Ok(());
    };

    let mut bb: Option<BuildingBlocks> = None;
    for path in &args.build {
        let b = mtb::read_mtb_file(path)
            .map_err(|e| format!("Error reading MTB file '{path}': {e:?}"))?;
        match bb.as_mut() {
            None => bb = Some(b),
            Some(first) => {
                first.solute_blocks.extend(b.solute_blocks);
                first.end_blocks.extend(b.end_blocks);
                first.solvent_blocks.extend(b.solvent_blocks);
            },
        }
    }
    let bb = bb.ok_or("Error: at least one @build file is required")?;
    let ffp = ifp::read_ifp_file(&args.param)
        .map_err(|e| format!("Error reading IFP file '{}': {e:?}", args.param))?;

    eprintln!("Building topology for {} residues", args.seq.len());
    eprintln!("Force field: {}", ffp.force_field);

    let result = make_top(&bb, &ffp, &args.seq, args.solv.as_deref())?;
    write!(out, "{}", result.text).map_err(io_err)?;
    Ok(())
}

/// Molecules from bond connectivity (union–find), as gromos++ `parse_molecules` does.
fn molecules(lt: &Lt) -> Vec<usize> {
    let n = lt.atoms.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for b in &lt.bonds {
        if b.0 >= 0 && b.1 >= 0 && (b.0 as usize) < n && (b.1 as usize) < n {
            let (a, c) = (
                find(&mut parent, b.0 as usize),
                find(&mut parent, b.1 as usize),
            );
            if a != c {
                parent[a.max(c)] = a.min(c);
            }
        }
    }
    // last atom of each molecule, molecules being contiguous in a linear topology
    let mut last_atoms = Vec::new();
    for i in 0..n {
        let next_same = i + 1 < n && find(&mut parent, i + 1) == find(&mut parent, i);
        if !next_same {
            last_atoms.push(i + 1);
        }
    }
    last_atoms
}

fn print_list(s: &mut String, prefix_width: usize, items: &[i64]) {
    // gromos++ OutTopology: six per line, continuation lines indented
    for (k, x) in items.iter().enumerate() {
        if k % 6 == 0 && k != 0 {
            let _ = writeln!(s);
            let _ = write!(s, "{:width$}", "", width = prefix_width);
        }
        let _ = write!(s, " {:5}", x);
    }
    let _ = writeln!(s);
}

/// The topology file in gromos++ `OutTopology` layout.
fn write_topology(
    lt: &Lt,
    bb: &BuildingBlocks,
    ffp: &ForceFieldParameters,
    solvent_name: Option<&str>,
) -> Result<String, String> {
    let mut s = String::new();
    let _ = writeln!(s, "TITLE");
    let _ = writeln!(s, "MAKE_TOP topology (newen), using:");
    let _ = writeln!(s, "Force-field code: {}", ffp.force_field);
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "PHYSICALCONSTANTS");
    let _ = writeln!(
        s,
        "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)"
    );
    let _ = writeln!(s, "{:15.4}", bb.fpepsi);
    let _ = writeln!(s, "# HBAR: Planck's constant HBAR = H/(2* PI)");
    let _ = writeln!(s, "{:15.7}", bb.hbar);
    let _ = writeln!(s, "# SPDL: Speed of light (nm/ps)");
    let _ = writeln!(s, "{:15.3}", bb.spdl);
    let _ = writeln!(s, "# BOLTZ: Boltzmann's constant kB");
    let _ = writeln!(s, "{:15.8}", bb.boltz);
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "TOPVERSION");
    let _ = writeln!(s, "2.0");
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "ATOMTYPENAME");
    let _ = writeln!(s, "# NRATT: number of van der Waals atom types");
    let _ = writeln!(s, "{:5}", ffp.num_atom_types);
    let _ = writeln!(s, "# TYPE: atom type names");
    for alj in &ffp.atom_lj {
        let _ = writeln!(s, "{}", alj.name);
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "RESNAME");
    let _ = writeln!(s, "# NRAA2: number of residues in a solute molecule");
    let _ = writeln!(s, "{:5}", lt.res_names.len());
    let _ = writeln!(s, "# AANM: residue names");
    for name in &lt.res_names {
        let _ = writeln!(s, "{}", name);
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "SOLUTEATOM");
    let _ = writeln!(s, "#   NRP: number of solute atoms");
    let _ = writeln!(s, "{:5}", lt.atoms.len());
    let _ = writeln!(s, "#  ATNM: atom number");
    let _ = writeln!(s, "#  MRES: residue number");
    let _ = writeln!(s, "#  PANM: atom name of solute atom");
    let _ = writeln!(s, "#   IAC: integer (van der Waals) atom type code");
    let _ = writeln!(s, "#  MASS: mass of solute atom");
    let _ = writeln!(s, "#    CG: charge of solute atom");
    let _ = writeln!(s, "#   CGC: charge group code (0 or 1)");
    let _ = writeln!(s, "#   INE: number of excluded atoms");
    let _ = writeln!(s, "# INE14: number of 1-4 interactions");
    let _ = writeln!(s, "# ATNM MRES PANM IAC     MASS       CG  CGC INE");
    let _ = writeln!(s, "#                                           INE14");
    for (i, a) in lt.atoms.iter().enumerate() {
        let _ = write!(
            s,
            "{:6} {:4} {:4} {:3} {:8.5} {:8.5} {:2} {:5}",
            i + 1,
            a.res + 1,
            a.name,
            a.iac + 1,
            a.mass,
            a.charge,
            a.chargegroup,
            a.exclusions.len()
        );
        let excl: Vec<i64> = a.exclusions.iter().map(|e| e + 1).collect();
        print_list(&mut s, 47, &excl);
        let _ = write!(s, "{:46}{}", "", a.one_four.len());
        let one_four: Vec<i64> = a.one_four.iter().map(|&e| e as i64 + 1).collect();
        print_list(&mut s, 45, &one_four);
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "BONDSTRETCHTYPE");
    let _ = writeln!(s, "#  NBTY: number of covalent bond types");
    let _ = writeln!(s, "{:5}", ffp.bond_types.len());
    let _ = writeln!(s, "#  CB: force constant");
    let _ = writeln!(s, "#  HB: harmonic force constant");
    let _ = writeln!(s, "#  B0: bond length at minimum energy");
    let _ = writeln!(s, "#         CB         HB         B0");
    for bt in &ffp.bond_types {
        let _ = writeln!(
            s,
            "{:15.7e}{:15.7e}{:13.7e}",
            bt.k_quartic, bt.k_harmonic, bt.r0
        );
    }
    let _ = writeln!(s, "END");
    let atoms = &lt.atoms;
    let h = |i: i64| is_h(atoms[i as usize].mass);
    let (bondh, bond): (Vec<&Bond>, Vec<&Bond>) = lt.bonds.iter().partition(|b| h(b.0) || h(b.1));
    let _ = writeln!(s, "BONDH");
    let _ = writeln!(s, "#  NBONH: number of bonds involving H atoms in solute");
    let _ = writeln!(s, "{:5}", bondh.len());
    let _ = writeln!(
        s,
        "#  IBH, JBH: atom sequence numbers of atoms forming a bond"
    );
    let _ = writeln!(s, "#  ICBH: bond type code");
    let _ = writeln!(s, "#   IBH    JBH ICBH");
    for b in &bondh {
        let _ = writeln!(s, "{:7}{:7}{:5}", b.0 + 1, b.1 + 1, b.2);
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "BOND");
    let _ = writeln!(
        s,
        "#  NBON: number of bonds NOT involving H atoms in solute"
    );
    let _ = writeln!(s, "{:5}", bond.len());
    let _ = writeln!(
        s,
        "#  IB, JB: atom sequence numbers of atoms forming a bond"
    );
    let _ = writeln!(s, "#  ICB: bond type code");
    let _ = writeln!(s, "#    IB     JB  ICB");
    for b in &bond {
        let _ = writeln!(s, "{:7}{:7}{:5}", b.0 + 1, b.1 + 1, b.2);
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "BONDANGLEBENDTYPE");
    let _ = writeln!(s, "#  NTTY: number of bond angle types");
    let _ = writeln!(s, "{:5}", ffp.angle_types.len());
    let _ = writeln!(s, "#  CT: force constant (based on potential");
    let _ = writeln!(s, "#      harmonic in the angle cosine)");
    let _ = writeln!(s, "#  CHT: force constant (based on potential");
    let _ = writeln!(s, "#      harmonic in the angle)");
    let _ = writeln!(s, "#  T0: bond angle at minimum energy in degrees");
    let _ = writeln!(s, "#         CT         CHT          T0");
    for at in &ffp.angle_types {
        let _ = writeln!(
            s,
            "{:15.7e}{:15.7e}{:13.7e}",
            at.k_non_harmonic, at.k_harmonic, at.theta0
        );
    }
    let _ = writeln!(s, "END");
    let (angleh, angle): (Vec<&AngleKey>, Vec<&AngleKey>) =
        lt.angles.iter().partition(|a| h(a.0) || h(a.1) || h(a.2));
    let _ = writeln!(s, "BONDANGLEH");
    let _ = writeln!(
        s,
        "#  NTHEH: number of bond angles involving H atoms in solute"
    );
    let _ = writeln!(s, "{:5}", angleh.len());
    let _ = writeln!(s, "#  ITH, JTH, KTH: atom sequence numbers");
    let _ = writeln!(s, "#    of atoms forming a bond angle in solute");
    let _ = writeln!(s, "#  ICTH: bond angle type code");
    let _ = writeln!(s, "#   ITH    JTH    KTH ICTH");
    for a in &angleh {
        let (i, j, k, t) = angle_atoms(a);
        let _ = writeln!(s, "{:7}{:7}{:7}{:5}", i + 1, j + 1, k + 1, t);
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "BONDANGLE");
    let _ = writeln!(s, "#  NTHE: number of bond angles NOT");
    let _ = writeln!(s, "#        involving H atoms in solute");
    let _ = writeln!(s, "{:5}", angle.len());
    let _ = writeln!(s, "#  IT, JT, KT: atom sequence numbers of atoms");
    let _ = writeln!(s, "#     forming a bond angle");
    let _ = writeln!(s, "#  ICT: bond angle type code");
    let _ = writeln!(s, "#    IT     JT     KT  ICT");
    for a in &angle {
        let (i, j, k, t) = angle_atoms(a);
        let _ = writeln!(s, "{:7}{:7}{:7}{:5}", i + 1, j + 1, k + 1, t);
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "IMPDIHEDRALTYPE");
    let _ = writeln!(s, "#  NQTY: number of improper dihedrals");
    let _ = writeln!(s, "{:5}", ffp.improper_types.len());
    let _ = writeln!(
        s,
        "#  CQ: force constant of improper dihedral per degrees square"
    );
    let _ = writeln!(
        s,
        "#  Q0: improper dihedral angle at minimum energy in degrees"
    );
    let _ = writeln!(s, "#         CQ          Q0");
    for it in &ffp.improper_types {
        let _ = writeln!(s, "{:15.7e}{:15.7e}", it.k, it.xi0);
    }
    let _ = writeln!(s, "END");
    let (imph, imp): (Vec<&Improper>, Vec<&Improper>) = lt
        .impropers
        .iter()
        .partition(|p| h(p.0) || h(p.1) || h(p.2) || h(p.3));
    let _ = writeln!(s, "IMPDIHEDRALH");
    let _ = writeln!(s, "#  NQHIH: number of improper dihedrals");
    let _ = writeln!(s, "#         involving H atoms in the solute");
    let _ = writeln!(s, "{:5}", imph.len());
    let _ = writeln!(s, "#  IQH,JQH,KQH,LQH: atom sequence numbers");
    let _ = writeln!(s, "#     of atoms forming an improper dihedral");
    let _ = writeln!(s, "#  ICQH: improper dihedral type code");
    let _ = writeln!(s, "#   IQH    JQH    KQH    LQH ICQH");
    for p in &imph {
        let _ = writeln!(
            s,
            "{:7}{:7}{:7}{:7}{:5}",
            p.0 + 1,
            p.1 + 1,
            p.2 + 1,
            p.3 + 1,
            p.4
        );
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "IMPDIHEDRAL");
    let _ = writeln!(s, "#  NQHI: number of improper dihedrals NOT");
    let _ = writeln!(s, "#    involving H atoms in solute");
    let _ = writeln!(s, "{:5}", imp.len());
    let _ = writeln!(s, "#  IQ,JQ,KQ,LQ: atom sequence numbers of atoms");
    let _ = writeln!(s, "#    forming an improper dihedral");
    let _ = writeln!(s, "#  ICQ: improper dihedral type code");
    let _ = writeln!(s, "#    IQ     JQ     KQ     LQ  ICQ");
    for p in &imp {
        let _ = writeln!(
            s,
            "{:7}{:7}{:7}{:7}{:5}",
            p.0 + 1,
            p.1 + 1,
            p.2 + 1,
            p.3 + 1,
            p.4
        );
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "TORSDIHEDRALTYPE");
    let _ = writeln!(s, "#  NPTY: number of dihedral types");
    let _ = writeln!(s, "{:5}", ffp.torsion_types.len());
    let _ = writeln!(s, "#  CP: force constant");
    let _ = writeln!(s, "#  PD: phase-shift angle");
    let _ = writeln!(s, "#  NP: multiplicity");
    let _ = writeln!(s, "#       CP        PD  NP");
    for tt in &ffp.torsion_types {
        let _ = writeln!(s, "{:10.5}{:10.5}{:4}", tt.k, tt.delta, tt.multiplicity);
    }
    let _ = writeln!(s, "END");
    let (dihh, dih): (Vec<&DihedralKey>, Vec<&DihedralKey>) = lt.dihedrals.iter().partition(|d| {
        let (i, j, k, l, _) = dihedral_atoms(d);
        h(i) || h(j) || h(k) || h(l)
    });
    let _ = writeln!(s, "DIHEDRALH");
    let _ = writeln!(
        s,
        "#  NPHIH: number of dihedrals involving H atoms in solute"
    );
    let _ = writeln!(s, "{:5}", dihh.len());
    let _ = writeln!(s, "#  IPH, JPH, KPH, LPH: atom sequence numbers");
    let _ = writeln!(s, "#    of atoms forming a dihedral");
    let _ = writeln!(s, "#  ICPH: dihedral type code");
    let _ = writeln!(s, "#   IPH    JPH    KPH    LPH ICPH");
    for d in &dihh {
        let (i, j, k, l, t) = dihedral_atoms(d);
        let _ = writeln!(s, "{:7}{:7}{:7}{:7}{:5}", i + 1, j + 1, k + 1, l + 1, t);
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "DIHEDRAL");
    let _ = writeln!(
        s,
        "#  NPHI: number of dihedrals NOT involving H atoms in solute"
    );
    let _ = writeln!(s, "{:5}", dih.len());
    let _ = writeln!(s, "#  IP, JP, KP, LP: atom sequence numbers");
    let _ = writeln!(s, "#     of atoms forming a dihedral");
    let _ = writeln!(s, "#  ICP: dihedral type code");
    let _ = writeln!(s, "#    IP     JP     KP     LP  ICP");
    for d in &dih {
        let (i, j, k, l, t) = dihedral_atoms(d);
        let _ = writeln!(s, "{:7}{:7}{:7}{:7}{:5}", i + 1, j + 1, k + 1, l + 1, t);
    }
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "CROSSDIHEDRALH");
    let _ = writeln!(
        s,
        "#  NPHIH: number of cross dihedrals involving H atoms in solute"
    );
    let _ = writeln!(s, "    0");
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "CROSSDIHEDRAL");
    let _ = writeln!(
        s,
        "#  NPPC: number of cross dihedrals NOT involving H atoms in solute"
    );
    let _ = writeln!(s, "    0");
    let _ = writeln!(s, "END");
    let _ = writeln!(s, "LJPARAMETERS");
    let _ = writeln!(
        s,
        "#  NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2"
    );
    let n_types = ffp.num_atom_types;
    let _ = writeln!(s, "{:5}", n_types * (n_types + 1) / 2);
    let _ = writeln!(s, "#  IAC,JAC: integer (van der Waals) atom type code");
    let _ = writeln!(s, "#  C12: r**(-12) term in nonbonded interactions");
    let _ = writeln!(s, "#   C6: r**(-6) term in nonbonded interactions");
    let _ = writeln!(s, "# CS12: r**(-12) term in 1-4 nonbonded interactions");
    let _ = writeln!(s, "#  CS6: r**(-6) term in 1-4 nonbonded interactions");
    let _ = writeln!(
        s,
        "# IAC  JAC           C12            C6          CS12           CS6"
    );
    for j in 1..=n_types {
        for i in 1..=j {
            let c6 = ffp.compute_c6(i, j);
            let c12 = ffp.compute_c12(i, j);
            let cs6 = ffp.compute_c6_14(i, j);
            let cs12 = ffp.compute_c12_14(i, j);
            let _ = writeln!(
                s,
                "{:5}{:5}{:15.6e}{:15.6e}{:15.6e}{:15.6e}",
                i, j, c12, c6, cs12, cs6
            );
        }
    }
    let _ = writeln!(s, "END");
    let mols = molecules(lt);
    let _ = writeln!(s, "SOLUTEMOLECULES");
    let _ = writeln!(s, "# NSPM: number of separate molecules in solute block");
    let _ = writeln!(s, "# NSP[1...NSPM]: atom sequence number of last atom");
    let _ = writeln!(s, "#                of the successive submolecules");
    let _ = writeln!(s, "#      NSPM  NSP[1...NSPM]");
    let _ = writeln!(s, "{:5}", mols.len());
    let mols_i: Vec<i64> = mols.iter().map(|&m| m as i64).collect();
    for chunk in mols_i.chunks(10) {
        let _ = writeln!(
            s,
            "{}",
            chunk.iter().map(|m| format!("{:6}", m)).collect::<String>()
        );
    }
    let _ = writeln!(s, "END");
    for (block, comment) in [
        ("TEMPERATUREGROUPS", "temperature"),
        ("PRESSUREGROUPS", "pressure"),
    ] {
        let _ = writeln!(s, "{block}");
        let _ = writeln!(s, "# NSTM: number of {comment} atom groups");
        let _ = writeln!(s, "# NST[1...NSTM]: atom sequence number of last atom");
        let _ = writeln!(
            s,
            "#                of the successive {comment} atom groups"
        );
        let _ = writeln!(s, "#      NSTM  NST[1...NSTM]");
        let _ = writeln!(s, "    1");
        let _ = writeln!(s, "{:6}", lt.atoms.len());
        let _ = writeln!(s, "END");
    }
    let _ = writeln!(s, "LJEXCEPTIONS");
    let _ = writeln!(
        s,
        "# This block defines special LJ-interactions based on atom numbers"
    );
    let _ = writeln!(
        s,
        "# This overrules the normal LJ-parameters (including 1-4s)"
    );
    let _ = writeln!(s, "# NEX: number of exceptions");
    let _ = writeln!(s, "    0");
    let _ = writeln!(s, "END");
    if let Some(sname) = solvent_name {
        let sol = bb
            .find_solvent(sname)
            .ok_or_else(|| format!("Error: solvent building block '{sname}' not found"))?;
        let _ = writeln!(s, "SOLVENTATOM");
        let _ = writeln!(s, "#  NRAM: number of atoms per solvent molecule");
        let _ = writeln!(s, "{:5}", sol.atoms.len());
        let _ = writeln!(s, "#     I: solvent atom sequence number");
        let _ = writeln!(s, "#  IACS: integer (van der Waals) atom type code");
        let _ = writeln!(s, "#  ANMS: atom name of solvent atom");
        let _ = writeln!(s, "#  MASS: mass of solvent atom");
        let _ = writeln!(s, "#   CGS: charge of solvent atom");
        let _ = writeln!(s, "#  I  ANMS IACS      MASS        CGS");
        for (idx, (aname, iac, mass_code, charge)) in sol.atoms.iter().enumerate() {
            let mass = ffp.get_mass(*mass_code).unwrap_or(0.0);
            let _ = writeln!(
                s,
                "{:5} {:>5}{:5}{:11.5}{:11.5}",
                idx + 1,
                aname,
                iac,
                mass,
                charge
            );
        }
        let _ = writeln!(s, "END");
        let _ = writeln!(s, "SOLVENTCONSTR");
        let _ = writeln!(s, "#  NCONS: number of constraints");
        let _ = writeln!(s, "{:5}", sol.constraints.len());
        let _ = writeln!(
            s,
            "#  ICONS, JCONS: atom sequence numbers forming constraint"
        );
        let _ = writeln!(s, "#   CONS constraint length");
        let _ = writeln!(s, "#ICONS JCONS         CONS");
        for (ic, jc, dist) in &sol.constraints {
            let _ = writeln!(s, "{:5}{:5}{:15.7}", ic, jc, dist);
        }
        let _ = writeln!(s, "END");
    }
    Ok(s)
}
