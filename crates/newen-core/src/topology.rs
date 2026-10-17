//! Topology module - molecular structure and force field parameters
//!
//! This is a direct Rust translation of GROMOS topology structures from:
//! - md++/src/topology/topology.h
//! - md++/src/topology/solute.h
//! - md++/src/topology/solvent.h
//!
//! # The instancing model (Phase 1)
//!
//! Why it is shaped this way: `.claude/CONTEXT.md` §"The instancing model".
//!
//! The legacy `Solute` / `Vec<Solvent>` fields remain the authoritative source
//! for bonded interactions during the transition.  The new fields alongside them
//! implement the *representation* half of it:
//!
//! - `MoleculeType` stores per-atom parameters **once** per unique molecule type.
//! - `MoleculeInstance` places one copy at an `atom_offset` in the flat arrays
//!   and carries a `Role` attribute (Solute / Solvent).
//! - `Topology::moltypes` + `Topology::instances` are populated by
//!   `init_solute_moltype()` and `solvate()`.
//!
//! **Phase 1 invariants:**
//! - `instances` are in system order (non-decreasing `atom_offset`). They match `molecules`
//!   1-to-1 only when the solute is one molecule: the solute is a single instance spanning
//!   every SOLUTEMOLECULES range, so look an atom's instance up by offset
//!   (`instance_of_atom`), never by `molecule_nr`.
//! - `role_of_atom(i)` → `Some(Role::Solvent)` iff atom `i` belongs to a
//!   solvent instance.  Returns `None` when the instance model is not yet
//!   initialised (pre-`solvate()` test topologies).
//! - `s:` in `AtomSelection` routes through `role_of_atom` when instances are
//!   available; falls back to `molecules[1..]` otherwise.
//!
//! Phase 2 (future): move bonds/angles/dihedrals into `MoleculeType`; make the
//! flat `iac`/`mass`/`charge` arrays into rebuildable caches; remove `Solute`
//! and `Vec<Solvent>`; add `promote(instance)` with CG/exclusion renumbering.

use crate::lambdas::{IndividualLambdas, InteractionLambda, PrecalcLambdas};
use crate::math::Vec3;

// ─── instancing model ───────────────────────────────────────────────────────

/// Semantic role of a molecule instance in the simulation.
///
/// Routes SHAKE/SETTLE dispatch, the solvent fast-path in the nonbonded
/// innerloop, and the `s:`/`m:` AtomSpecifier grammar.
/// In Phase 1 this is set at `solvate()` time; Phase 2 adds `promote()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Role {
    /// Molecule is part of the solute (non-solvent) component.
    Solute,
    /// Molecule is a solvent molecule.
    Solvent,
}

/// Per-atom parameters stored *once* inside a `MoleculeType`.
///
/// Phase 2d: merged with `Atom` — carries the same fields including the
/// boolean attribute flags.  `Atom` is now an alias kept for backwards compat
/// during the remaining migration.
#[derive(Debug, Clone, Default)]
pub struct MolTypeAtom {
    /// Atom name (e.g. `"CA"`, `"OW"`).
    pub name: String,
    /// 1-based residue number within the molecule type.
    pub residue_nr: usize,
    /// Residue name (e.g. `"ALA"`, `"SOL"`).
    pub residue_name: String,
    /// Integer atom code (0-based atom type index into the LJ matrix).
    pub iac: usize,
    /// Atomic mass (u).
    pub mass: f64,
    /// Partial charge (elementary charge units).
    pub charge: f64,
    /// `true` if this atom appears in the PERTATOMPARAM block (FEP).
    pub is_perturbed: bool,
    /// `true` if this atom uses the COS polarisation model.
    pub is_polarisable: bool,
    /// `true` if this atom is a coarse-grained bead.
    pub is_coarse_grained: bool,
}

/// Legacy alias — `Atom` is now the same type as `MolTypeAtom`.
/// Kept so existing code compiles unchanged during Phase 2d migration;
/// will be removed once all sites use `MolTypeAtom` directly.
pub type Atom = MolTypeAtom;

/// Topology stored once for a unique molecule type.
///
/// Phase 1: atom parameters.
/// Phase 2: bonds, angles, dihedrals — migrated here from `Solute`.
/// Indices within bonded terms are *local* (0-based within this moltype).
/// For the solute at atom_offset=0, local == global, so no conversion needed
/// during the migration.
#[derive(Debug, Clone)]
pub struct MoleculeType {
    /// Human-readable name ("SOLUTE", "SPC", …).
    pub name: String,
    /// One entry per atom in the template, in order.
    pub atoms: Vec<MolTypeAtom>,
    /// Covalent bonds within this molecule type (local atom indices).
    pub bonds: Vec<Bond>,
    /// Bond-angle terms within this molecule type (local atom indices).
    pub angles: Vec<Angle>,
    /// Proper dihedral terms within this molecule type (local atom indices).
    pub proper_dihedrals: Vec<Dihedral>,
    /// Improper dihedral terms within this molecule type (local atom indices).
    pub improper_dihedrals: Vec<Dihedral>,
    /// Cross-dihedral coupling terms within this molecule type (local atom indices).
    pub cross_dihedrals: Vec<CrossDihedral>,
}

impl MoleculeType {
    /// Number of atoms in one copy of this molecule type.
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
}

/// One copy of a `MoleculeType` placed in the system.
#[derive(Debug, Clone)]
pub struct MoleculeInstance {
    /// Index into `Topology::moltypes`.
    pub moltype_id: usize,
    /// First atom's 0-based global index in the flat position/force arrays.
    pub atom_offset: usize,
    /// Semantic role — drives dispatch and selection grammar.
    pub role: Role,
}

/// Two-body bonded term (bonds, harmonic constraints)
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Index into `Topology::bond_parameters`.
    pub bond_type: usize,
    /// The bond is in the topology's BONDH block (a bond to hydrogen): the set NTC = 2
    /// constrains (gromosXX `in_topology.cc:2143-2178`).
    pub hydrogen: bool,
}

/// Three-body bonded term (angles)
#[derive(Debug, Clone, Copy)]
pub struct Angle {
    /// First atom index.
    pub i: usize,
    /// Central atom index.
    pub j: usize,
    /// Third atom index.
    pub k: usize,
    /// Index into `Topology::angle_parameters`.
    pub angle_type: usize,
}

/// Four-body bonded term (proper and improper dihedrals)
#[derive(Debug, Clone, Copy)]
pub struct Dihedral {
    /// First atom index (i–j–k–l sequence).
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Third atom index.
    pub k: usize,
    /// Fourth atom index.
    pub l: usize,
    /// Index into `Topology::dihedral_parameters` or `improper_dihedral_parameters`.
    pub dihedral_type: usize,
}

/// Cross-dihedral term (8 atoms)
#[derive(Debug, Clone, Copy)]
pub struct CrossDihedral {
    /// Atom a (first dihedral plane).
    pub a: usize,
    /// Atom b.
    pub b: usize,
    /// Atom c.
    pub c: usize,
    /// Atom d (shared with second plane).
    pub d: usize,
    /// Atom e.
    pub e: usize,
    /// Atom f.
    pub f: usize,
    /// Atom g.
    pub g: usize,
    /// Atom h (second dihedral plane).
    pub h: usize,
    /// Index into the cross-dihedral parameter table.
    pub cross_dihedral_type: usize,
}

/// Perturbed atom — state A and B LJ+charge parameters (PERTATOMPARAM block).
/// `lj_soft` and `crf_soft` are the per-atom soft-core alphas from the file;
/// they are scaled by ALPHLJ/ALPHC from the imd PERTURBATION block at run time.
#[derive(Debug, Clone)]
pub struct PerturbedAtom {
    /// 0-based atom index in the topology.
    pub seq: usize,
    /// State A integer atom code (0-indexed).
    pub a_iac: usize,
    /// State A atomic mass (u).
    pub a_mass: f64,
    /// State A partial charge.
    pub a_charge: f64,
    /// State B integer atom code (0-indexed).
    pub b_iac: usize,
    /// State B atomic mass (u).
    pub b_mass: f64,
    /// State B partial charge.
    pub b_charge: f64,
    /// Per-atom soft-core α for LJ (scaled by global ALPHLJ at runtime).
    pub lj_soft: f64,
    /// Per-atom soft-core α for CRF (scaled by global ALPHC at runtime).
    pub crf_soft: f64,
}

/// Net charge (e) of a set of atoms in the two end states of a perturbation (1.8).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NetCharge {
    /// State A.
    pub a: f64,
    /// State B.
    pub b: f64,
}

impl NetCharge {
    /// Charges summed in f64 differ from an exact sum by far less than this.
    pub const TOLERANCE: f64 = 1e-6;

    /// ΔQ = Q_B − Q_A, exactly 0 within [`NetCharge::TOLERANCE`].
    pub fn delta(&self) -> f64 {
        let d = self.b - self.a;
        if d.abs() < Self::TOLERANCE {
            0.0
        } else {
            d
        }
    }

    /// The net charge of each molecule in `molecules` that holds an atom marked in
    /// `perturbed`, with its index: `charge_a` and `charge_b` are the per-atom charges of the
    /// two states. An empty `molecules` is one molecule over every atom.
    pub fn of_perturbed_molecules(
        molecules: &[std::ops::Range<usize>],
        charge_a: &[f64],
        charge_b: &[f64],
        perturbed: &[bool],
    ) -> Vec<(usize, NetCharge)> {
        let whole = 0..charge_a.len();
        let molecules = if molecules.is_empty() {
            std::slice::from_ref(&whole)
        } else {
            molecules
        };
        molecules
            .iter()
            .enumerate()
            .filter(|(_, m)| (m.start..m.end).any(|i| perturbed.get(i).copied().unwrap_or(false)))
            .map(|(k, m)| {
                (
                    k,
                    NetCharge {
                        a: m.clone().filter_map(|i| charge_a.get(i)).sum(),
                        b: m.clone().filter_map(|i| charge_b.get(i)).sum(),
                    },
                )
            })
            .collect()
    }
}

/// `A 1.000 B 0.000 (ΔQ = -1.000)`: three decimals, a rounded −0 printed as 0.
impl std::fmt::Display for NetCharge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let r = |x: f64| (x * 1e3).round() / 1e3 + 0.0;
        write!(
            f,
            "A {:.3} B {:.3} (ΔQ = {:.3})",
            r(self.a),
            r(self.b),
            r(self.delta())
        )
    }
}

/// How a perturbed atom pair interacts in one state (`PERTATOMPAIR` columns `i(A)`, `i(B)`).
///
/// The file values are read as they are, as gromosXX does (`in_perturbation.cc:1186`,
/// `perturbed_nonbonded_pair.cc:328`): `0` excluded, `1` normal, `2` 1-4.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairInteraction {
    /// File `0`: excluded — only the reaction-field term of an excluded pair (with NSLFEXCL).
    Excluded,
    /// File `1`: a normal pair — LJ with c6/c12 and the full reaction field.
    Normal,
    /// File `2`: a 1-4 pair — LJ with cs6/cs12 and the full reaction field.
    OneFour,
}

impl PairInteraction {
    /// The interaction type of a `PERTATOMPAIR` file value; `None` for anything but 0, 1, 2.
    pub fn from_code(code: i64) -> Option<Self> {
        match code {
            0 => Some(Self::Excluded),
            1 => Some(Self::Normal),
            2 => Some(Self::OneFour),
            _ => None,
        }
    }
}

/// Perturbed atom pair (`PERTATOMPAIR` block): a pair whose interaction type differs between
/// states A and B. Like gromosXX (`in_perturbation.cc:1207-1224`),
/// [`Topology::apply_perturbation`] takes such a pair out of the exclusion and 1-4 lists and
/// keeps it off the pairlist, so the perturbed atom-pair term is the only place it is computed.
///
/// Atom indices are 0-indexed (the file's are 1-indexed).
#[derive(Debug, Clone, Copy)]
pub struct PerturbedAtomPair {
    /// First atom index (0-indexed, always i < j).
    pub i: usize,
    /// Second atom index (0-indexed).
    pub j: usize,
    /// State A interaction type.
    pub a_type: PairInteraction,
    /// State B interaction type.
    pub b_type: PairInteraction,
}

/// One atom of an EDS multi-state solute (`MPERTATOM`, gromosXX `EDS_Perturbed_Atom`): its
/// atom type and charge in every end state, and the excluded and 1-4 partners whose
/// state-dependent terms it owns. `iac`/`charge` are indexed by state; `topo.iac`/`topo.charge`
/// keep the atom's topology values.
#[derive(Debug, Clone, PartialEq)]
pub struct EdsPerturbedAtom {
    /// 0-based atom index.
    pub seq: usize,
    /// Integer atom code per end state (0-indexed).
    pub iac: Vec<usize>,
    /// Partial charge per end state.
    pub charge: Vec<f64>,
    /// Excluded partners (sorted): their reaction-field term is computed per state for this
    /// atom and nowhere else ([`Topology::apply_eds_perturbation`]).
    pub exclusions: Vec<usize>,
    /// 1-4 partners (sorted), likewise per state.
    pub one_four_pairs: Vec<usize>,
}

/// The EDS multi-state solute (`MPERTATOM` block of the perturbation topology).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct EdsSolute {
    /// Number of end states (`NPTB`, the `EDS` block's NUMSTATES).
    pub num_states: usize,
    /// The end-state names of the block's second line.
    pub state_names: Vec<String>,
    /// The EDS atoms in ascending `seq` (gromosXX keeps them in a `std::map`).
    pub atoms: Vec<EdsPerturbedAtom>,
}

impl EdsSolute {
    /// No EDS atom.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }
}

/// Perturbed bond (for FEP calculations)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedBond {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// State A bond type index.
    pub a_type: usize,
    /// State B bond type index.
    pub b_type: usize,
    /// A bond to hydrogen: the topology bond it replaces is in BONDH (or, with no such bond,
    /// it was read from PERTBONDSTRETCHH).
    pub hydrogen: bool,
}

/// Perturbed angle (for FEP calculations)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedAngle {
    /// First atom index.
    pub i: usize,
    /// Central atom index.
    pub j: usize,
    /// Third atom index.
    pub k: usize,
    /// State A angle type index.
    pub a_type: usize,
    /// State B angle type index.
    pub b_type: usize,
}

/// Perturbed dihedral (for FEP calculations)
#[derive(Debug, Clone, Copy)]
pub struct PerturbedDihedral {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Third atom index.
    pub k: usize,
    /// Fourth atom index.
    pub l: usize,
    /// State A dihedral type index.
    pub a_type: usize,
    /// State B dihedral type index.
    pub b_type: usize,
}

/// Bond force field parameters (GROMOS format)
#[derive(Debug, Clone, Copy)]
pub struct BondParameters {
    /// Quartic force constant (kJ mol⁻¹ nm⁻⁴).
    pub k_quartic: f64,
    /// Harmonic force constant (kJ mol⁻¹ nm⁻²).
    pub k_harmonic: f64,
    /// Equilibrium bond length (nm).
    pub r0: f64,
}

/// Angle force field parameters (GROMOS format)
#[derive(Debug, Clone, Copy)]
pub struct AngleParameters {
    /// Force constant for the cos-harmonic form (kJ mol⁻¹).
    pub k_cosine: f64,
    /// Force constant for the harmonic-in-angle form (kJ mol⁻¹ rad⁻²).
    pub k_harmonic: f64,
    /// Equilibrium bond angle (rad).
    pub theta0: f64,
}

/// Dihedral force field parameters
#[derive(Debug, Clone, Copy)]
pub struct DihedralParameters {
    /// Torsion force constant (kJ mol⁻¹).
    pub k: f64,
    /// Phase shift δ (rad); typically 0 or π.
    pub pd: f64,
    /// Precomputed cos(δ) for efficiency.
    pub cospd: f64,
    /// Torsion multiplicity n.
    pub m: i32,
}

/// Improper dihedral force field parameters
#[derive(Debug, Clone, Copy)]
pub struct ImproperDihedralParameters {
    /// Force constant (kJ mol⁻¹ rad⁻²).
    pub k: f64,
    /// Equilibrium improper dihedral angle (rad).
    pub q0: f64,
}

/// Lennard-Jones parameters
#[derive(Debug, Clone, Copy)]
pub struct LJParameters {
    /// C6 attractive coefficient (kJ mol⁻¹ nm⁶); contributes −C6/r⁶.
    pub c6: f64,
    /// C12 repulsive coefficient (kJ mol⁻¹ nm¹²); contributes +C12/r¹².
    pub c12: f64,
    /// 1-4 (soft-core) C6 used in free energy calculations.
    pub cs6: f64,
    /// 1-4 (soft-core) C12 used in free energy calculations.
    pub cs12: f64,
}

impl LJParameters {
    /// Create standard LJ parameters (without softcore)
    pub fn new(c6: f64, c12: f64) -> Self {
        Self {
            c6,
            c12,
            cs6: c6,
            cs12: c12,
        }
    }

    /// Create LJ parameters with explicit 1-4 (cs6, cs12) values
    pub fn new_with_14(c6: f64, c12: f64, cs6: f64, cs12: f64) -> Self {
        Self { c6, c12, cs6, cs12 }
    }

    /// Calculate sigma (size parameter) from C6 and C12
    /// sigma = (C12/C6)^(1/6)
    pub fn sigma(&self) -> f64 {
        if self.c6 > 0.0 {
            (self.c12 / self.c6).powf(1.0 / 6.0)
        } else {
            0.0
        }
    }

    /// Calculate epsilon (well depth) from C6 and C12
    /// epsilon = C6^2 / (4 * C12)
    pub fn epsilon(&self) -> f64 {
        if self.c12 > 0.0 {
            self.c6 * self.c6 / (4.0 * self.c12)
        } else {
            0.0
        }
    }
}

impl Default for LJParameters {
    fn default() -> Self {
        Self {
            c6: 0.0,
            c12: 0.0,
            cs6: 0.0,
            cs12: 0.0,
        }
    }
}

/// Exclusion list for an atom (atoms that don't interact via nonbonded forces)
///
/// Uses a sorted Vec for fast binary-search lookups. Matches the GROMOS C++
/// convention (sorted vector with binary search). For typical exclusion counts
/// (2-6 per atom), this is much faster than HashSet due to cache locality.
pub type Exclusions = Vec<usize>;

/// Membership test on one atom's sorted exclusion list.
///
/// Relies on two invariants every construction path maintains (the topology
/// reader, [`Topology::build_exclusions`], and the solvent expansion):
/// lists are **sorted** and **symmetric** — `j ∈ exclusions[i]` iff
/// `i ∈ exclusions[j]` — which is what lets callers query one side only.
/// Both are checked by `exclusion_invariants` in the newen-io tests.
///
/// The `last < j` early-out is gromosXX's (`topology/exclusions.h:79`,
/// `Exclusions::is_excluded`): exclusions are intramolecular and the lists
/// are short, so for a partner in another molecule the answer is usually
/// settled by looking at the last element, with no search at all.
#[inline]
fn exclusion_list_contains(list: &[usize], j: usize) -> bool {
    match list.last() {
        None => false,
        Some(&last) if last < j => false,
        _ => list.binary_search(&j).is_ok(),
    }
}

/// Chargegroup - atoms whose charge sum is zero (or small)
///
/// Chargegroups are used for cutoff optimizations:
/// - Distance calculated between chargegroup centers
/// - All atoms in both groups included if distance < cutoff
#[derive(Debug, Clone)]
pub struct ChargeGroup {
    /// 0-based global atom indices in this charge group.
    pub atoms: Vec<usize>,
}

impl ChargeGroup {
    /// Calculate center of geometry for this chargegroup
    pub fn center_of_geometry(&self, positions: &[Vec3]) -> Vec3 {
        if self.atoms.is_empty() {
            return Vec3::ZERO;
        }

        let sum: Vec3 = self.atoms.iter().map(|&i| positions[i]).sum();

        sum / self.atoms.len() as f64
    }
}

/// Solute atom data.
///
/// Phase 2c: bonded terms (bonds/angles/dihedrals) have been moved to
/// Solute — Phase 2d: atoms moved to MoleculeType. This struct is now empty
/// and kept only as a named placeholder until the `topo.solute` field
/// is removed in Phase 2e (when `Vec<Solvent>` is also removed).
/// All data is now in `topo.moltypes[0]`.
#[derive(Debug, Clone, Default)]
pub struct Solute;

impl Solute {
    /// Create the empty solute placeholder.
    pub fn new() -> Self {
        Self
    }
    /// Number of solute atoms — delegates to Topology::num_solute_atoms().
    /// Call `topo.num_solute_atoms()` directly instead.
    pub fn num_atoms(&self) -> usize {
        0
    }
}

/// Soft-core perturbed harmonic bond (PERTBONDSOFT).
/// Used when a bond is absent in one state (K=0), requiring soft-core treatment.
/// A type `None` means the bond is absent in that state: K = 0 and the other state's r0
/// (gromosXX `perturbed_soft_bond_interaction.cc:249-263`). At most one of the two is `None`.
#[derive(Debug, Clone)]
pub struct SoftBond {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// State A bond type index (uses `k_harmonic`); `None` means K_A=0, r0_A=r0_B.
    pub a_type: Option<usize>,
    /// State B bond type index; `None` means K_B=0, r0_B=r0_A.
    pub b_type: Option<usize>,
    /// Soft-core coupling strength α.
    pub alpha: f64,
}

/// Soft-core perturbed cos-harmonic angle (PERTANGLESOFT).
/// A type `None` means the angle is absent in that state: K = 0 and the other state's cos0.
#[derive(Debug, Clone)]
pub struct SoftAngle {
    /// First atom index.
    pub i: usize,
    /// Central atom index.
    pub j: usize,
    /// Third atom index.
    pub k: usize,
    /// State A angle type index; `None` means K_A=0, cos0_A=cos0_B.
    pub a_type: Option<usize>,
    /// State B angle type index; `None` means K_B=0, cos0_B=cos0_A.
    pub b_type: Option<usize>,
    /// Soft-core coupling strength α.
    pub alpha: f64,
}

/// Soft-core perturbed improper dihedral (PERTIMPROPERDIHSOFT).
/// A type `None` means the improper is absent in that state: K = 0 and the other state's q0.
#[derive(Debug, Clone)]
pub struct SoftImproper {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Third atom index.
    pub k: usize,
    /// Fourth atom index.
    pub l: usize,
    /// State A improper dihedral type index; `None` means K_A=0, q0_A=q0_B.
    pub a_type: Option<usize>,
    /// State B improper dihedral type index; `None` means K_B=0, q0_B=q0_A.
    pub b_type: Option<usize>,
    /// Soft-core coupling strength α.
    pub alpha: f64,
}

/// Perturbed solute for FEP calculations
///
/// Contains dual-topology (A/B state) bonded terms for free energy perturbation
#[derive(Debug, Clone)]
pub struct PerturbedSolute {
    /// Perturbed atoms: state A/B IAC, mass, charge, soft-core alphas.
    /// Indexed sparsely — only atoms listed in PERTATOMPARAM.
    pub atoms: Vec<PerturbedAtom>,
    /// Perturbed atom pairs (PERTATOMPAIR): special LJ between two perturbed atoms.
    pub atom_pairs: Vec<PerturbedAtomPair>,
    /// Perturbed bond terms (PERTBOND block).
    pub bonds: Vec<PerturbedBond>,
    /// Perturbed angle terms (PERTANGLE block).
    pub angles: Vec<PerturbedAngle>,
    /// Perturbed proper dihedral terms (PERTPROPERDIH block).
    pub proper_dihedrals: Vec<PerturbedDihedral>,
    /// Perturbed improper dihedral terms (PERTIMPROPERDIH block).
    pub improper_dihedrals: Vec<PerturbedDihedral>,
    /// Soft-core perturbed bonds (PERTBONDSOFT): bond absent in one state.
    pub soft_bonds: Vec<SoftBond>,
    /// Soft-core perturbed angles (PERTANGLESOFT): angle absent in one state.
    pub soft_angles: Vec<SoftAngle>,
    /// Soft-core perturbed impropers (PERTIMPROPERDIHSOFT): improper absent in one state.
    pub soft_impropers: Vec<SoftImproper>,
    /// Perturbed bonds a constraint algorithm owns (gromosXX
    /// `perturbed_solute().distance_constraints()`): with NTC = 3 every perturbed bond moves
    /// here, out of `bonds`, and perturbed SHAKE holds it at the λ-mixed length.
    pub distance_constraints: Vec<PerturbedBond>,
}

impl PerturbedSolute {
    /// Create an empty perturbed-solute container.
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            atom_pairs: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            proper_dihedrals: Vec::new(),
            improper_dihedrals: Vec::new(),
            soft_bonds: Vec::new(),
            soft_angles: Vec::new(),
            soft_impropers: Vec::new(),
            distance_constraints: Vec::new(),
        }
    }

    /// Returns `true` if no perturbed atoms or bonded terms have been added.
    /// No perturbation of any kind — every list, including the soft-core bonded terms and the
    /// perturbed atom pairs (a `.ptp` with only `PERTBONDSOFT` used to count as empty, so its
    /// soft bonds were never evaluated: the regular bond was removed and nothing replaced it).
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
            && self.atom_pairs.is_empty()
            && self.bonds.is_empty()
            && self.angles.is_empty()
            && self.proper_dihedrals.is_empty()
            && self.improper_dihedrals.is_empty()
            && self.soft_bonds.is_empty()
            && self.soft_angles.is_empty()
            && self.soft_impropers.is_empty()
            && self.distance_constraints.is_empty()
    }
}

impl Default for PerturbedSolute {
    fn default() -> Self {
        Self::new()
    }
}

/// Solvent molecule structure (typically water)
/// Solvent atom template (one molecule's worth of atom properties)
#[derive(Debug, Clone)]
pub struct SolventAtomTemplate {
    /// Integer atom code (0-based atom type).
    pub iac: usize,
    /// Atom name (e.g. `"OW"`, `"HW1"`).
    pub name: String,
    /// Atomic mass (u).
    pub mass: f64,
    /// Partial charge.
    pub charge: f64,
}

/// Solvent constraint template (within one molecule)
#[derive(Debug, Clone)]
pub struct SolventConstraintTemplate {
    /// First atom index (0-based within the solvent molecule template).
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Constraint distance (nm).
    pub length: f64,
}

/// Distance restraint specification (virtual atom type 0 only)
#[derive(Debug, Clone)]
pub struct DistanceRestraintSpec {
    /// First restrained atom (0-based).
    pub atom1: usize,
    /// Second restrained atom (0-based).
    pub atom2: usize,
    /// Reference distance (nm).
    pub r0: f64,
    /// Force constant / weight.
    pub w0: f64,
    /// Restraint type flag (RAH).
    pub rah: i32,
}

/// Perturbed distance restraint specification (virtual atom type 0 only)
#[derive(Debug, Clone)]
pub struct PerturbedDistanceRestraintSpec {
    /// First restrained atom (0-based).
    pub atom1: usize,
    /// Second restrained atom (0-based).
    pub atom2: usize,
    /// State A exponent n.
    pub n: i32,
    /// State B exponent m.
    pub m: i32,
    /// State A reference distance (nm).
    pub a_r0: f64,
    /// State B reference distance (nm).
    pub b_r0: f64,
    /// State A force constant / weight.
    pub a_w0: f64,
    /// State B force constant / weight.
    pub b_w0: f64,
    /// Restraint type flag (RAH).
    pub rah: i32,
}

/// Main topology structure containing all molecular information
#[derive(Debug, Clone)]
pub struct Topology {
    // Solute is now a shell (atoms in moltypes[0]).
    // Phase 2e: Vec<Solvent> removed — solvent info lives in moltypes + instances.
    /// Legacy solute placeholder (Phase 2e: will be removed; use `moltypes[0]` directly).
    pub solute: Solute,
    /// FEP dual-topology perturbed terms (atoms, bonds, angles, dihedrals).
    pub perturbed_solute: PerturbedSolute,
    /// EDS end states of the solute atoms listed in `MPERTATOM` (empty without EDS).
    pub eds_solute: EdsSolute,

    // Per-atom properties (flat arrays for all atoms)
    /// Integer atom codes (atom types).
    ///
    /// **INVARIANT: 0-indexed** (0..N\_types-1), matching LJ matrix rows/cols.
    /// GROMOS files store IAC as 1-indexed; `newen-io/topology.rs` subtracts 1
    /// at the parse boundary.  `ptp.rs` does the same for `PerturbedAtom.a_iac`
    /// and `b_iac`.  Never store a raw 1-indexed file value here.
    pub iac: Vec<usize>,
    /// Flat atomic masses array (u); indexed by global atom index.
    pub mass: Vec<f64>,
    /// Precomputed 1/mass for each atom; avoids division in integrators.
    pub inverse_mass: Vec<f64>,
    /// Flat atomic charges array; indexed by global atom index.
    pub charge: Vec<f64>,

    // Exclusions (atoms that don't interact via nonbonded)
    /// Per-atom exclusion lists; `exclusions[i]` lists atoms that do not interact with i.
    pub exclusions: Vec<Exclusions>,
    /// 1-4 atom pairs (receive special LJ/CRF scaling).
    pub one_four_pairs: Vec<Vec<usize>>,
    /// Pairs in neither list above but kept off the pairlist (gromosXX's `all_exclusion` is
    /// built before the perturbation topology is read), per atom, sorted and symmetric: the
    /// `PERTATOMPAIR` pairs and the excluded and 1-4 pairs `MPERTATOM` hands to an EDS atom.
    /// Empty without such pairs.
    pub atom_pair_exclusions: Vec<Vec<usize>>,

    // Chargegroups
    /// All charge groups in the system (solute first, then solvent).
    pub chargegroups: Vec<ChargeGroup>,
    /// Maps each atom to its charge group index.
    pub atom_to_chargegroup: Vec<usize>,
    /// Number of solute charge groups; CGs `[0..n)` are solute, `[n..)` are solvent.
    pub num_solute_chargegroups: usize,

    // Temperature and pressure coupling groups
    /// Atoms in each temperature coupling group (MULTIBATH/TCOUPLE).
    pub temperature_groups: Vec<Vec<usize>>,
    /// Pressure coupling groups as atom index ranges (GROMOS PRESSUREGROUPS block).
    pub pressure_groups: Vec<std::ops::Range<usize>>,
    /// Temperature coupling group index for each atom.
    pub atom_to_temperature_group: Vec<usize>,
    /// Pressure coupling group index for each atom.
    pub atom_to_pressure_group: Vec<usize>,

    // Energy groups (for energy monitoring)
    /// Atoms in each energy monitoring group.
    pub energy_groups: Vec<Vec<usize>>,
    /// Energy group index for each atom.
    pub atom_to_energy_group: Vec<usize>,

    // Molecule boundaries (for molecule-based operations)
    /// Per-molecule atom ranges; `molecules[m]` is `start..end` for molecule m.
    pub molecules: Vec<std::ops::Range<usize>>,

    // ── instancing model ─────────────────────────────────────────────────────
    // Populated by init_solute_moltype() + solvate().
    // One solute instance (however many molecules the solute has), then one per solvent
    // molecule — so instances[k] is molecules[k] only for a one-molecule solute.
    // Empty until those methods are called (e.g. in raw test topologies).
    /// Registry: one entry per unique molecule type (parameters stored once).
    pub moltypes: Vec<MoleculeType>,
    /// All instances in system order (non-decreasing `atom_offset`).
    pub instances: Vec<MoleculeInstance>,

    // Force field parameters
    /// Bond force field parameters indexed by bond type.
    pub bond_parameters: Vec<BondParameters>,
    /// Angle force field parameters indexed by angle type.
    pub angle_parameters: Vec<AngleParameters>,
    /// Proper and improper dihedral parameters indexed by dihedral type.
    pub dihedral_parameters: Vec<DihedralParameters>,
    /// Improper dihedral parameters indexed by improper dihedral type.
    pub improper_dihedral_parameters: Vec<ImproperDihedralParameters>,
    /// LJ parameter matrix indexed `[type_i][type_j]`.
    pub lj_parameters: Vec<Vec<LJParameters>>,

    /// Atom parameters for one solvent molecule (from topology SOLUTEATOM/SOLVENTATOM).
    pub solvent_atom_template: Vec<SolventAtomTemplate>,
    /// Constraint template within one solvent molecule.
    pub solvent_constraint_template: Vec<SolventConstraintTemplate>,

    /// Solute bonds that a constraint algorithm owns (`CONSTRAINT` NTC ≥ 2), as ordered global
    /// atom pairs. gromosXX moves these out of the bond *force* list when it creates the
    /// constraints (`create_constraints.cc`), so a constrained bond contributes no potential
    /// energy and no force — its length is imposed instead. Empty for NTC = 1.
    pub constrained_bonds: std::collections::HashSet<(usize, usize)>,

    /// Charge-group codes (CGC) for solute atoms, from the topology CGSOLUTE block.
    pub chargegroup_codes: Vec<usize>,

    /// Distance restraint specifications (loaded from the `.distres` file).
    pub distance_restraints: Vec<DistanceRestraintSpec>,
    /// Perturbed distance restraint specifications (FEP/TI).
    pub perturbed_distance_restraints: Vec<PerturbedDistanceRestraintSpec>,

    /// Per-solute-atom FEP flag; `true` if the atom appears in PERTATOMPARAM.
    pub is_perturbed: Vec<bool>,

    /// The global λ of a perturbed run (gromosXX `topo.lambda()`, RLAM).
    pub lambda: f64,
    /// NLAM, the exponent of the nonbonded λ dependence (gromosXX `topo.lambda_exp()`).
    pub lambda_exp: i32,
    /// The individual λ of every interaction kind and energy-group pair, evaluated at
    /// [`Self::lambda`] by [`Self::update_for_lambda`].
    pub individual_lambdas: IndividualLambdas,
    /// PRECALCLAM: the extended-TI λ grid and the steps it is evaluated on; `None` = off.
    pub precalc_lambdas: Option<PrecalcLambdas>,
}

impl Topology {
    /// Create an empty topology.
    pub fn new() -> Self {
        Self {
            solute: Solute::new(),
            perturbed_solute: PerturbedSolute::new(),
            eds_solute: EdsSolute::default(),
            iac: Vec::new(),
            mass: Vec::new(),
            inverse_mass: Vec::new(),
            charge: Vec::new(),
            exclusions: Vec::new(),
            one_four_pairs: Vec::new(),
            atom_pair_exclusions: Vec::new(),
            chargegroups: Vec::new(),
            atom_to_chargegroup: Vec::new(),
            num_solute_chargegroups: 0,
            temperature_groups: Vec::new(),
            pressure_groups: Vec::new(), // populated from PRESSUREGROUPS topology block
            atom_to_temperature_group: Vec::new(),
            atom_to_pressure_group: Vec::new(),
            energy_groups: Vec::new(),
            atom_to_energy_group: Vec::new(),
            molecules: Vec::new(),
            // moltypes[0] is always the solute MoleculeType (GROMACS moltype[0] convention).
            // Initialised empty here so callers can write `topo.moltypes[0].bonds.push(..)`
            // without any guard or ensure_* call.
            moltypes: vec![MoleculeType {
                name: "SOLUTE".to_string(),
                atoms: Vec::new(),
                bonds: Vec::new(),
                angles: Vec::new(),
                proper_dihedrals: Vec::new(),
                improper_dihedrals: Vec::new(),
                cross_dihedrals: Vec::new(),
            }],
            instances: vec![MoleculeInstance {
                moltype_id: 0,
                atom_offset: 0,
                role: Role::Solute,
            }],
            bond_parameters: Vec::new(),
            angle_parameters: Vec::new(),
            dihedral_parameters: Vec::new(),
            improper_dihedral_parameters: Vec::new(),
            lj_parameters: Vec::new(),
            solvent_atom_template: Vec::new(),
            solvent_constraint_template: Vec::new(),
            constrained_bonds: std::collections::HashSet::new(),
            chargegroup_codes: Vec::new(),
            distance_restraints: Vec::new(),
            perturbed_distance_restraints: Vec::new(),
            is_perturbed: Vec::new(),
            lambda: 0.0,
            lambda_exp: 1,
            individual_lambdas: IndividualLambdas::identity(1),
            precalc_lambdas: None,
        }
    }

    /// Total number of atoms in the system — derived from instances.
    pub fn num_atoms(&self) -> usize {
        if !self.instances.is_empty() {
            self.instances
                .iter()
                .map(|inst| {
                    self.moltypes
                        .get(inst.moltype_id)
                        .map(|mt| mt.num_atoms())
                        .unwrap_or(0)
                })
                .sum()
        } else {
            self.num_solute_atoms()
        }
    }

    // ── Phase 2e helpers — replace .solvents[0].X accesses ──────────────────

    /// Number of solvent molecule instances.
    pub fn num_solvent_molecules(&self) -> usize {
        self.instances
            .iter()
            .filter(|i| i.role == Role::Solvent)
            .count()
    }

    /// Number of atoms per solvent molecule (from the solvent MoleculeType).
    /// Falls back to `solvent_atom_template.len()` for pre-solvate topologies.
    pub fn atoms_per_solvent(&self) -> usize {
        self.instances
            .iter()
            .find(|i| i.role == Role::Solvent)
            .and_then(|i| self.moltypes.get(i.moltype_id))
            .map(|mt| mt.num_atoms())
            .unwrap_or(self.solvent_atom_template.len())
    }

    /// Name of the solvent molecule type.
    pub fn solvent_name(&self) -> Option<&str> {
        self.instances
            .iter()
            .find(|i| i.role == Role::Solvent)
            .and_then(|i| self.moltypes.get(i.moltype_id))
            .map(|mt| mt.name.as_str())
    }

    /// Number of solute atoms — derived from instances when available.
    pub fn num_solute_atoms(&self) -> usize {
        if !self.instances.is_empty() {
            self.instances
                .iter()
                .filter(|inst| inst.role == Role::Solute)
                .map(|inst| {
                    self.moltypes
                        .get(inst.moltype_id)
                        .map(|mt| mt.num_atoms())
                        .unwrap_or(0)
                })
                .sum()
        } else {
            // no instances: the solute is `moltypes[0]`
            self.moltypes.first().map_or(0, |mt| mt.num_atoms())
        }
    }

    /// Number of atom types
    pub fn num_atom_types(&self) -> usize {
        self.lj_parameters.len()
    }

    /// Get LJ parameters for atom pair (i, j)
    #[inline]
    pub fn lj_parameter(&self, type_i: usize, type_j: usize) -> &LJParameters {
        &self.lj_parameters[type_i][type_j]
    }

    /// Check if atoms i and j are excluded from nonbonded interactions
    /// (1-2 and 1-3 bonded pairs only — used for RF excluded interactions)
    ///
    /// One lookup suffices: exclusion lists are symmetric (see
    /// [`exclusion_list_contains`]), so `j ∈ exclusions[i]` ⇔ `i ∈ exclusions[j]`.
    #[inline]
    pub fn is_excluded(&self, i: usize, j: usize) -> bool {
        exclusion_list_contains(&self.exclusions[i], j)
    }

    /// Check if atoms i and j are excluded, a 1-4 pair or a `PERTATOMPAIR` pair.
    /// Used for pairlist exclusion (GROMOS: all_exclusion).
    #[inline]
    pub fn is_excluded_or_14(&self, i: usize, j: usize) -> bool {
        // Both lists are symmetric, so each needs only one direction. This is
        // the pairlist's innermost test — one call per candidate atom pair on
        // every rebuild — which is why the second lookup was worth removing.
        exclusion_list_contains(&self.exclusions[i], j)
            || (!self.one_four_pairs.is_empty() && self.one_four_pairs[i].contains(&j))
            || self
                .atom_pair_exclusions
                .get(i)
                .is_some_and(|l| l.contains(&j))
    }

    /// Get chargegroup index for atom i
    #[inline]
    pub fn chargegroup(&self, i: usize) -> usize {
        self.atom_to_chargegroup[i]
    }

    /// Get temperature group index for atom i
    #[inline]
    pub fn temperature_group(&self, i: usize) -> usize {
        self.atom_to_temperature_group[i]
    }

    /// Get energy group index for atom i
    #[inline]
    pub fn energy_group(&self, i: usize) -> usize {
        self.atom_to_energy_group[i]
    }

    /// gromosXX `Topology::update_for_lambda` (`topology.cc:1419-1478`): evaluate the
    /// individual λ at [`Self::lambda`] over the energy groups, then mix each perturbed atom's
    /// mass with the mass λ of its own group, `m = (1 − L)·m_A + L·m_B` (no NLAM exponent).
    /// Call it after the energy groups are assigned and whenever λ changes.
    pub fn update_for_lambda(&mut self) {
        let n_groups = self.energy_groups.len().max(1);
        self.individual_lambdas.resize(n_groups);
        self.individual_lambdas.update(self.lambda);
        for a in &self.perturbed_solute.atoms {
            let g = self.atom_to_energy_group.get(a.seq).copied().unwrap_or(0);
            let l = self
                .individual_lambdas
                .lambda(InteractionLambda::Mass, g, g);
            let mass = (1.0 - l) * a.a_mass + l * a.b_mass;
            if a.seq < self.mass.len() && mass > 0.0 {
                self.mass[a.seq] = mass;
                self.inverse_mass[a.seq] = 1.0 / mass;
            }
        }
    }

    // ── instancing accessors ─────────────────────────────────────────────────

    /// Index into `instances` of the instance that owns global atom `i`.
    ///
    /// Found by atom offset, not by `molecule_nr`: the solute is one instance however many
    /// molecules it has, so instance and molecule numbers agree only for a one-molecule
    /// solute. Instances are in system order; an empty instance (the solute of a pure-solvent
    /// system shares offset 0 with the first solvent molecule) owns no atom.
    pub fn instance_of_atom(&self, i: usize) -> Option<usize> {
        let k = self
            .instances
            .partition_point(|inst| inst.atom_offset <= i)
            .checked_sub(1)?;
        let inst = &self.instances[k];
        let n = self.moltypes.get(inst.moltype_id)?.atoms.len();
        (i - inst.atom_offset < n).then_some(k)
    }

    /// Look up the MolTypeAtom for global atom index `i` via instances.
    pub fn moltype_atom(&self, i: usize) -> Option<&MolTypeAtom> {
        let inst = &self.instances[self.instance_of_atom(i)?];
        self.moltypes[inst.moltype_id]
            .atoms
            .get(i - inst.atom_offset)
    }

    /// Role of the molecule instance that owns atom `i`.
    /// Returns `None` when the instance model hasn't been initialised yet.
    pub fn role_of_atom(&self, i: usize) -> Option<Role> {
        self.instance_of_atom(i).map(|k| self.instances[k].role)
    }

    /// `true` iff atom `i` belongs to a `Role::Solvent` instance.
    pub fn is_solvent_atom(&self, i: usize) -> bool {
        self.role_of_atom(i) == Some(Role::Solvent)
    }

    /// Change the role of molecule instance `mol_idx` to `Role::Solute`.
    ///
    /// This is the "promote one solvent molecule to solute" operation.
    /// Phase 2 will add CG/exclusion renumbering; for now it's a label flip.
    pub fn promote(&mut self, mol_idx: usize) -> Result<(), String> {
        match self.instances.get_mut(mol_idx) {
            Some(inst) => {
                inst.role = Role::Solute;
                Ok(())
            },
            None => Err(format!(
                "molecule index {mol_idx} out of range (have {})",
                self.instances.len()
            )),
        }
    }

    /// Build the solute `MoleculeType` + `MoleculeInstance` from the existing
    /// `solute.atoms` and flat parameter arrays.
    ///
    /// Safe to call repeatedly (no-op if already initialised).
    /// Called by `solvate()` and by `newen-io::build_topology()`.
    /// Ensure a solute MoleculeType + instance exist at index 0.
    ///
    /// Creates an EMPTY moltype (no atoms, no bonds) if none exists yet.
    /// Sync flat iac/mass/charge into `moltypes[0].atoms` after the flat arrays
    /// are populated by `build_topology()`.  `moltypes[0]` always exists (created
    /// in `Topology::new()`) so no guard is needed.
    pub fn init_solute_moltype(&mut self) {
        let n_sol = self.moltypes[0].atoms.len();
        if n_sol == 0 {
            return;
        }

        // Overwrite iac/mass/charge from the flat arrays (source of truth at load time).
        let mt_atoms: Vec<MolTypeAtom> = self.moltypes[0]
            .atoms
            .iter()
            .enumerate()
            .map(|(i, a)| MolTypeAtom {
                name: a.name.clone(),
                residue_nr: a.residue_nr,
                residue_name: a.residue_name.clone(),
                iac: self.iac.get(i).copied().unwrap_or(a.iac),
                mass: self.mass.get(i).copied().unwrap_or(a.mass),
                charge: self.charge.get(i).copied().unwrap_or(a.charge),
                is_perturbed: a.is_perturbed,
                is_polarisable: a.is_polarisable,
                is_coarse_grained: a.is_coarse_grained,
            })
            .collect();
        self.moltypes[0].atoms = mt_atoms;

        // Ensure molecules[0] is set
        if self.molecules.is_empty() && n_sol > 0 {
            self.molecules.push(0..n_sol);
        }
    }

    /// Atom name for atom `i` (covers solute + all expanded solvent). Returns `None` if out of range.
    ///
    /// Prefers the moltype path when instances are populated; falls back to
    /// the legacy `solute.atoms` + `solvent_atom_template` approach.
    pub fn atom_name(&self, i: usize) -> Option<&str> {
        self.moltype_atom(i).map(|a| a.name.as_str())
    }

    /// 1-based residue number for atom `i`.
    pub fn residue_nr(&self, i: usize) -> Option<usize> {
        self.moltype_atom(i).map(|a| a.residue_nr)
    }

    /// Residue name for atom `i` (e.g. `"ALA"`, `"SOL"`).
    pub fn residue_name(&self, i: usize) -> Option<&str> {
        self.moltype_atom(i).map(|a| a.residue_name.as_str())
    }

    /// 0-based molecule index for atom `i` (uses the `molecules` ranges populated by `solvate()`).
    /// Returns `None` if `molecules` is empty or atom is out of all molecule ranges.
    pub fn molecule_nr(&self, i: usize) -> Option<usize> {
        // molecules is sorted by start; binary-search on start to find candidate
        let idx = self.molecules.partition_point(|r| r.start <= i);
        if idx > 0 {
            let mol = idx - 1;
            if self.molecules[mol].contains(&i) {
                return Some(mol);
            }
        }
        None
    }

    // ── Phase 3: instance-iterating bonded force iterators ───────────────────
    // Each iterator walks all MoleculeType instances, translates local bond indices
    // (0-based within the moltype) to global indices using inst.atom_offset, and
    // yields values with global i/j/k/l. This handles flexible solute, flexible
    // solvent, and any future repeated molecule type in one unified loop.

    /// Iterate over all bond terms with global atom indices.
    pub fn all_bonds_global(&self) -> impl Iterator<Item = Bond> + '_ {
        self.instances.iter().flat_map(move |inst| {
            let mt = &self.moltypes[inst.moltype_id];
            let off = inst.atom_offset;
            mt.bonds.iter().map(move |b| Bond {
                i: off + b.i,
                j: off + b.j,
                ..*b
            })
        })
    }

    /// Iterate over all angle terms with global atom indices.
    pub fn all_angles_global(&self) -> impl Iterator<Item = Angle> + '_ {
        self.instances.iter().flat_map(move |inst| {
            let mt = &self.moltypes[inst.moltype_id];
            let off = inst.atom_offset;
            mt.angles.iter().map(move |a| Angle {
                i: off + a.i,
                j: off + a.j,
                k: off + a.k,
                angle_type: a.angle_type,
            })
        })
    }

    /// Iterate over all proper dihedral terms with global atom indices.
    pub fn all_proper_dihedrals_global(&self) -> impl Iterator<Item = Dihedral> + '_ {
        self.instances.iter().flat_map(move |inst| {
            let mt = &self.moltypes[inst.moltype_id];
            let off = inst.atom_offset;
            mt.proper_dihedrals.iter().map(move |d| Dihedral {
                i: off + d.i,
                j: off + d.j,
                k: off + d.k,
                l: off + d.l,
                dihedral_type: d.dihedral_type,
            })
        })
    }

    /// Iterate over all improper dihedral terms with global atom indices.
    pub fn all_improper_dihedrals_global(&self) -> impl Iterator<Item = Dihedral> + '_ {
        self.instances.iter().flat_map(move |inst| {
            let mt = &self.moltypes[inst.moltype_id];
            let off = inst.atom_offset;
            mt.improper_dihedrals.iter().map(move |d| Dihedral {
                i: off + d.i,
                j: off + d.j,
                k: off + d.k,
                l: off + d.l,
                dihedral_type: d.dihedral_type,
            })
        })
    }

    /// Iterate over all cross-dihedral coupling terms with global atom indices.
    pub fn all_cross_dihedrals_global(&self) -> impl Iterator<Item = CrossDihedral> + '_ {
        self.instances.iter().flat_map(move |inst| {
            let mt = &self.moltypes[inst.moltype_id];
            let off = inst.atom_offset;
            mt.cross_dihedrals.iter().map(move |c| CrossDihedral {
                a: off + c.a,
                b: off + c.b,
                c: off + c.c,
                d: off + c.d,
                e: off + c.e,
                f: off + c.f,
                g: off + c.g,
                h: off + c.h,
                cross_dihedral_type: c.cross_dihedral_type,
            })
        })
    }

    /// Build exclusions from bonded topology
    ///
    /// Excludes:
    /// - 1-2 bonded neighbors
    /// - 1-3 angle neighbors
    /// - Optionally 1-4 dihedral neighbors (or use special 1-4 scaling)
    pub fn build_exclusions(&mut self, exclude_14: bool) {
        let n_atoms = self.num_atoms();
        self.exclusions = vec![Vec::new(); n_atoms];

        // Collect into locals first — moltypes and exclusions are different struct fields
        // so Rust allows simultaneous borrows, but collecting avoids lifetime complexity.
        let bonds: Vec<Bond> = self.moltypes[0].bonds.clone();
        let angles: Vec<Angle> = self.moltypes[0].angles.clone();
        let dihedrals: Vec<Dihedral> = self.moltypes[0].proper_dihedrals.clone();

        // Exclude bonded neighbors (1-2)
        for bond in &bonds {
            self.exclusions[bond.i].push(bond.j);
            self.exclusions[bond.j].push(bond.i);
        }

        // Exclude angle neighbors (1-3)
        for angle in &angles {
            self.exclusions[angle.i].push(angle.k);
            self.exclusions[angle.k].push(angle.i);
        }

        // Optionally exclude or mark 1-4 pairs
        if exclude_14 {
            for dihedral in &dihedrals {
                self.exclusions[dihedral.i].push(dihedral.l);
                self.exclusions[dihedral.l].push(dihedral.i);
            }
        } else {
            // Store as special 1-4 pairs
            self.one_four_pairs = vec![Vec::new(); n_atoms];
            for dihedral in &dihedrals {
                self.one_four_pairs[dihedral.i].push(dihedral.l);
                self.one_four_pairs[dihedral.l].push(dihedral.i);
            }
        }

        // Sort and deduplicate for fast binary_search lookups
        for excl in &mut self.exclusions {
            excl.sort_unstable();
            excl.dedup();
        }
    }

    /// Initialize LJ parameter matrix from combining rules
    ///
    /// GROMOS uses geometric mean combining rules:
    /// - C6_ij = sqrt(C6_ii * C6_jj)
    /// - C12_ij = sqrt(C12_ii * C12_jj)
    pub fn build_lj_matrix(&mut self, lj_types: &[LJParameters]) {
        let n_types = lj_types.len();
        self.lj_parameters = vec![vec![LJParameters::new(0.0, 0.0); n_types]; n_types];

        for i in 0..n_types {
            for j in 0..n_types {
                let c6 = (lj_types[i].c6 * lj_types[j].c6).sqrt();
                let c12 = (lj_types[i].c12 * lj_types[j].c12).sqrt();
                self.lj_parameters[i][j] = LJParameters::new(c6, c12);
            }
        }
    }

    /// Resize all per-atom arrays to match number of atoms
    pub fn resize_atom_arrays(&mut self) {
        let n_atoms = self.num_atoms();
        self.iac.resize(n_atoms, 0);
        self.mass.resize(n_atoms, 0.0);
        self.inverse_mass.resize(n_atoms, 0.0);
        self.charge.resize(n_atoms, 0.0);
        self.exclusions.resize(n_atoms, Vec::new());
        self.one_four_pairs.resize(n_atoms, Vec::new());
        self.atom_to_chargegroup.resize(n_atoms, 0);
        self.atom_to_temperature_group.resize(n_atoms, 0);
        self.atom_to_pressure_group.resize(n_atoms, 0);
        self.atom_to_energy_group.resize(n_atoms, 0);
    }

    /// Rebuild `iac`, `mass`, `charge` flat caches from `moltypes` + `instances`.
    ///
    /// This is the Phase 2a step: the flat arrays are no longer expanded
    /// manually in `solvate()` — they are derived from the instance registry.
    /// Calling this after adding/removing instances keeps the caches coherent.
    pub fn rebuild_flat_arrays(&mut self) {
        self.iac.clear();
        self.mass.clear();
        self.charge.clear();
        for inst in &self.instances {
            if let Some(mt) = self.moltypes.get(inst.moltype_id) {
                for a in &mt.atoms {
                    self.iac.push(a.iac);
                    self.mass.push(a.mass);
                    self.charge.push(a.charge);
                }
            }
        }
        self.compute_inverse_masses();
    }

    /// Compute inverse masses for all atoms (for efficient integration)
    pub fn compute_inverse_masses(&mut self) {
        self.inverse_mass = self
            .mass
            .iter()
            .map(|&m| if m > 0.0 { 1.0 / m } else { 0.0 })
            .collect();
    }

    /// Merge a perturbation topology (`.ptp`) into this topology, as gromosXX does when
    /// `NTG != 0`: the perturbed solute terms are stored in `perturbed_solute`, and the bonds,
    /// angles, improper and proper dihedrals they perturb are removed from the regular solute
    /// so the two sets stay disjoint. `is_perturbed` is resized to the solute atom count.
    /// Each `PERTATOMPAIR` pair leaves the exclusion and 1-4 lists and goes to
    /// `atom_pair_exclusions` (gromosXX `in_perturbation.cc:1207-1224`).
    ///
    /// Call after `build_topology()` and before `solvate()`; applying twice replaces the
    /// previous perturbation's terms but cannot restore bonds or exclusions already removed.
    /// Returns the solute's net charge in the two states ([`Topology::perturbation_net_charge`]),
    /// or an error, leaving the topology unchanged, when a perturbed or soft bond is not a bond
    /// of the topology.
    pub fn apply_perturbation(
        &mut self,
        perturbed_solute: PerturbedSolute,
        is_perturbed: Vec<bool>,
    ) -> Result<NetCharge, String> {
        self.check_perturbed_bonds_exist(&perturbed_solute)?;
        self.merge_perturbation(perturbed_solute, is_perturbed);
        Ok(self.perturbation_net_charge())
    }

    /// gromosXX's "Perturbation of a non-existing bond" errors (`in_perturbation.cc:232-247`,
    /// 271-286, 420-425). Its match is on the ordered pair, so `j i` does not find bond `i j`.
    /// gromosXX looks in the distance constraints instead when NTC = 3 (and, for PERTBONDSOFT,
    /// always in both); every BONDH/BONDS bond of the topology counts here.
    fn check_perturbed_bonds_exist(&self, ps: &PerturbedSolute) -> Result<(), String> {
        use std::collections::HashSet;

        let bonds: HashSet<(usize, usize)> = self
            .moltypes
            .first()
            .map(|m| m.bonds.iter().map(|b| (b.i, b.j)).collect())
            .unwrap_or_default();
        let missing = |i: usize, j: usize, block: &str| {
            Err(format!(
                "Perturbation of a non-existing bond in {block} block. (atoms {} {})",
                i + 1,
                j + 1
            ))
        };
        for b in &ps.bonds {
            if !bonds.contains(&(b.i, b.j)) {
                let block = if b.hydrogen {
                    "PERTBONDSTRETCHH"
                } else {
                    "PERTBONDSTRETCH"
                };
                return missing(b.i, b.j, block);
            }
        }
        // A soft bond absent in state A is a bond being made: it need not exist.
        for b in ps.soft_bonds.iter().filter(|b| b.a_type.is_some()) {
            if !bonds.contains(&(b.i, b.j)) {
                return missing(b.i, b.j, "PERTBONDSOFT");
            }
        }
        Ok(())
    }

    /// The solute's per-atom charges in state A and state B: the topology's charges with each
    /// perturbed atom's `PERTATOMPARAM` charge of that state in its place.
    pub fn end_state_charges(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.num_solute_atoms().min(self.charge.len());
        let mut a = self.charge[..n].to_vec();
        let mut b = a.clone();
        for p in &self.perturbed_solute.atoms {
            if p.seq < n {
                a[p.seq] = p.a_charge;
                b[p.seq] = p.b_charge;
            }
        }
        (a, b)
    }

    /// Net charge of the whole solute in state A and state B.
    pub fn perturbation_net_charge(&self) -> NetCharge {
        let (a, b) = self.end_state_charges();
        NetCharge {
            a: a.iter().sum(),
            b: b.iter().sum(),
        }
    }

    /// Net charge in state A and B of each solute molecule holding a perturbed atom, with the
    /// molecule's index.
    pub fn perturbed_molecule_charges(&self) -> Vec<(usize, NetCharge)> {
        let (a, b) = self.end_state_charges();
        let mut perturbed = vec![false; a.len()];
        for p in &self.perturbed_solute.atoms {
            if let Some(x) = perturbed.get_mut(p.seq) {
                *x = true;
            }
        }
        let solute: Vec<_> = self
            .molecules
            .iter()
            .filter(|m| m.end <= a.len())
            .cloned()
            .collect();
        NetCharge::of_perturbed_molecules(&solute, &a, &b, &perturbed)
    }

    fn merge_perturbation(
        &mut self,
        perturbed_solute: PerturbedSolute,
        mut is_perturbed: Vec<bool>,
    ) {
        use std::collections::HashSet;

        let n = self.num_solute_atoms();
        is_perturbed.resize(n, false);
        self.is_perturbed = is_perturbed;
        self.perturbed_solute = perturbed_solute;

        for ap in &self.perturbed_solute.atom_pairs {
            for (a, b) in [(ap.i, ap.j), (ap.j, ap.i)] {
                if let Some(l) = self.exclusions.get_mut(a) {
                    l.retain(|&x| x != b);
                }
                if let Some(l) = self.one_four_pairs.get_mut(a) {
                    l.retain(|&x| x != b);
                }
                if self.atom_pair_exclusions.len() <= a {
                    self.atom_pair_exclusions.resize(a + 1, Vec::new());
                }
                let l = &mut self.atom_pair_exclusions[a];
                if let Err(pos) = l.binary_search(&b) {
                    l.insert(pos, b);
                }
            }
        }

        let bondh: HashSet<(usize, usize)> = self.moltypes[0]
            .bonds
            .iter()
            .filter(|b| b.hydrogen)
            .map(|b| (b.i.min(b.j), b.i.max(b.j)))
            .collect();
        let regular: HashSet<(usize, usize)> = self.moltypes[0]
            .bonds
            .iter()
            .map(|b| (b.i.min(b.j), b.i.max(b.j)))
            .collect();
        for b in &mut self.perturbed_solute.bonds {
            let key = (b.i.min(b.j), b.i.max(b.j));
            if regular.contains(&key) {
                b.hydrogen = bondh.contains(&key);
            }
        }
        let pb: HashSet<(usize, usize)> = self
            .perturbed_solute
            .bonds
            .iter()
            .map(|b| (b.i.min(b.j), b.i.max(b.j)))
            .chain(
                self.perturbed_solute
                    .soft_bonds
                    .iter()
                    .map(|b| (b.i.min(b.j), b.i.max(b.j))),
            )
            .collect();
        self.moltypes[0]
            .bonds
            .retain(|b| !pb.contains(&(b.i.min(b.j), b.i.max(b.j))));

        let pa: HashSet<(usize, usize, usize)> = self
            .perturbed_solute
            .angles
            .iter()
            .map(|a| (a.i, a.j, a.k))
            .chain(
                self.perturbed_solute
                    .soft_angles
                    .iter()
                    .map(|a| (a.i, a.j, a.k)),
            )
            .collect();
        self.moltypes[0]
            .angles
            .retain(|a| !pa.contains(&(a.i, a.j, a.k)));

        let pi: HashSet<(usize, usize, usize, usize)> = self
            .perturbed_solute
            .improper_dihedrals
            .iter()
            .map(|d| (d.i, d.j, d.k, d.l))
            .chain(
                self.perturbed_solute
                    .soft_impropers
                    .iter()
                    .map(|d| (d.i, d.j, d.k, d.l)),
            )
            .collect();
        self.moltypes[0]
            .improper_dihedrals
            .retain(|d| !pi.contains(&(d.i, d.j, d.k, d.l)));

        let pd: HashSet<(usize, usize, usize, usize)> = self
            .perturbed_solute
            .proper_dihedrals
            .iter()
            .map(|d| (d.i, d.j, d.k, d.l))
            .collect();
        self.moltypes[0]
            .proper_dihedrals
            .retain(|d| !pd.contains(&(d.i, d.j, d.k, d.l)));
    }

    /// Merge the EDS end states of an `MPERTATOM` block, as gromosXX's `In_Perturbation` does
    /// (`in_perturbation.cc:1757-1777`). `atoms` come in file order with empty partner lists.
    ///
    /// Each atom takes over its remaining excluded and 1-4 partners: they leave
    /// `exclusions`/`one_four_pairs` (so the standard reaction-field and 1-4 loops no longer
    /// see the pair) and stay off the pairlist through `atom_pair_exclusions` (gromosXX's
    /// `all_exclusion` is not updated). A pair of two EDS atoms belongs to the one listed
    /// first. The atoms are stored in ascending `seq`, the order gromosXX's map visits them.
    ///
    /// Call after `build_topology()` (and `apply_perturbation`) and before `solvate()`.
    pub fn apply_eds_perturbation(
        &mut self,
        num_states: usize,
        state_names: Vec<String>,
        atoms: Vec<EdsPerturbedAtom>,
    ) -> Result<(), String> {
        let n_solute = self.num_solute_atoms();
        let n_types = self.lj_parameters.len();
        let mut seen = vec![false; n_solute];
        for a in &atoms {
            if a.seq >= n_solute {
                return Err(format!(
                    "MPERTATOM: atom sequence number {} is not a solute atom (the solute has {n_solute})",
                    a.seq + 1
                ));
            }
            if std::mem::replace(&mut seen[a.seq], true) {
                return Err(format!("MPERTATOM: atom {} listed twice", a.seq + 1));
            }
            if a.iac.len() != num_states || a.charge.len() != num_states {
                return Err(format!(
                    "MPERTATOM: atom {} has {} types and {} charges for {num_states} states",
                    a.seq + 1,
                    a.iac.len(),
                    a.charge.len()
                ));
            }
            if let Some(&t) = a.iac.iter().find(|&&t| t >= n_types) {
                return Err(format!(
                    "MPERTATOM: atom {} has atom type {} but the topology has {n_types}",
                    a.seq + 1,
                    t + 1
                ));
            }
        }
        let mut atoms = atoms;
        for atom in &mut atoms {
            let e = atom.seq;
            for (lists, owned) in [
                (&mut self.exclusions, &mut atom.exclusions),
                (&mut self.one_four_pairs, &mut atom.one_four_pairs),
            ] {
                let Some(list) = lists.get_mut(e) else {
                    continue;
                };
                *owned = std::mem::take(list);
                for &j in owned.iter() {
                    if let Some(l) = lists.get_mut(j) {
                        l.retain(|&x| x != e);
                    }
                }
            }
            for &j in atom.exclusions.iter().chain(&atom.one_four_pairs) {
                for (a, b) in [(e, j), (j, e)] {
                    if self.atom_pair_exclusions.len() <= a {
                        self.atom_pair_exclusions.resize(a + 1, Vec::new());
                    }
                    let l = &mut self.atom_pair_exclusions[a];
                    if let Err(pos) = l.binary_search(&b) {
                        l.insert(pos, b);
                    }
                }
            }
        }
        atoms.sort_by_key(|a| a.seq);
        self.eds_solute = EdsSolute {
            num_states,
            state_names,
            atoms,
        };
        Ok(())
    }

    /// Expand solvent molecules (GROMOS: topo.solvate(0, nsm))
    ///
    /// Uses the solvent template stored on the topology to create `nsm` copies.
    /// Expands flat arrays (mass, charge, iac), chargegroups, and exclusions.
    /// Must be called after build_topology() and before reading coordinates.
    pub fn solvate(&mut self, nsm: usize) {
        if nsm == 0 || self.solvent_atom_template.is_empty() {
            return;
        }

        let atoms_per_solvent = self.solvent_atom_template.len();
        let n_solute = self.num_solute_atoms();

        // Ensure the solute moltype + instance are initialised before we add solvent.
        self.init_solute_moltype();

        // Add solute as one molecule (if it has atoms)
        if n_solute > 0 && self.molecules.is_empty() {
            self.molecules.push(0..n_solute);
        }

        let n_solvent_atoms = atoms_per_solvent * nsm;
        let n_total = n_solute + n_solvent_atoms;

        // Phase 2e: Solvent struct removed — solvent data lives in moltypes + instances.

        // Each solvent molecule is its own chargegroup
        for mol in 0..nsm {
            let base = n_solute + mol * atoms_per_solvent;
            let atoms: Vec<usize> = (base..base + atoms_per_solvent).collect();
            self.chargegroups.push(ChargeGroup { atoms });
        }

        // Rebuild atom_to_chargegroup mapping
        self.atom_to_chargegroup = vec![0; n_total];
        for (cg_idx, cg) in self.chargegroups.iter().enumerate() {
            for &atom in &cg.atoms {
                if atom < n_total {
                    self.atom_to_chargegroup[atom] = cg_idx;
                }
            }
        }

        // Extend exclusions for solvent atoms
        self.exclusions.resize(n_total, Vec::new());
        // `one_four_pairs` is either empty (no 1-4 data at all) or indexed by
        // every atom; solvent atoms have no 1-4 partners but must still be
        // addressable, or an atomic-cutoff pairlist on a solvated system
        // indexes past the end.
        if !self.one_four_pairs.is_empty() {
            self.one_four_pairs.resize(n_total, Vec::new());
        }

        // Solvent intra-molecular exclusions: all atoms within each molecule exclude each other
        for mol in 0..nsm {
            let base = n_solute + mol * atoms_per_solvent;
            // Populate molecules list
            self.molecules.push(base..base + atoms_per_solvent);
            // Each solvent molecule is its own pressure and temperature group (GROMOS convention)
            self.pressure_groups.push(base..base + atoms_per_solvent);
            self.temperature_groups
                .push((base..base + atoms_per_solvent).collect());
            for a in 0..atoms_per_solvent {
                for b in (a + 1)..atoms_per_solvent {
                    self.exclusions[base + a].push(base + b);
                    self.exclusions[base + b].push(base + a);
                }
            }
        }

        // Sort solvent exclusions for binary_search
        for i in n_solute..n_total {
            self.exclusions[i].sort_unstable();
        }

        // Add the solvent MoleculeType (stored once) + one Instance per molecule.
        {
            let sv_name = "SOLV".to_string();

            let sv_moltype = MoleculeType {
                name: sv_name,
                atoms: self
                    .solvent_atom_template
                    .iter()
                    .enumerate()
                    .map(|(idx, sa)| MolTypeAtom {
                        name: sa.name.clone(),
                        residue_nr: idx + 1,
                        residue_name: "SOLV".to_string(),
                        iac: sa.iac,
                        mass: sa.mass,
                        charge: sa.charge,
                        ..MolTypeAtom::default()
                    })
                    .collect(),
                // Solvent bonded terms live in solvent_constraint_template (constraints only).
                // No proper bonds/angles/dihedrals for rigid SPC-like solvents.
                bonds: Vec::new(),
                angles: Vec::new(),
                proper_dihedrals: Vec::new(),
                improper_dihedrals: Vec::new(),
                cross_dihedrals: Vec::new(),
            };
            let solvent_moltype_id = self.moltypes.len();
            self.moltypes.push(sv_moltype);

            for mol in 0..nsm {
                let atom_offset = n_solute + mol * atoms_per_solvent;
                self.instances.push(MoleculeInstance {
                    moltype_id: solvent_moltype_id,
                    atom_offset,
                    role: Role::Solvent,
                });
            }
        }

        // Phase 2a: derive flat arrays from the now-complete instance registry.
        // Must run before the LJ matrix resize which reads self.iac.
        self.rebuild_flat_arrays();

        // Rebuild LJ matrix if solvent introduces new IAC types
        let max_iac = self.iac.iter().max().copied().unwrap_or(0);
        let n_types = max_iac + 1;
        if n_types > self.lj_parameters.len() {
            let old_len = self.lj_parameters.len();
            self.lj_parameters
                .resize(n_types, vec![LJParameters::default(); n_types]);
            for row in self.lj_parameters.iter_mut() {
                row.resize(n_types, LJParameters::default());
            }
            log::debug!(
                "LJ matrix expanded from {}x{} to {}x{} for solvent types",
                old_len,
                old_len,
                n_types,
                n_types
            );
        }

        log::debug!("Solvated: {} solute + {} solvent ({} molecules × {} atoms) = {} total, {} chargegroups",
            n_solute, n_solvent_atoms, nsm, atoms_per_solvent, n_total,
            self.chargegroups.len());
    }
}

impl Default for Topology {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solute_atoms_without_instances_come_from_the_solute_moltype() {
        let mut topo = Topology::new();
        topo.instances.clear();
        assert_eq!(topo.num_solute_atoms(), 0);
        assert_eq!(topo.num_atoms(), 0);
        for name in ["C", "O"] {
            topo.moltypes[0].atoms.push(MolTypeAtom {
                name: name.into(),
                residue_nr: 1,
                residue_name: "RES".into(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        assert_eq!(topo.num_solute_atoms(), 2);
        assert_eq!(topo.num_atoms(), 2);
    }

    /// Five solute atoms with symmetric exclusion and 1-4 lists and two atom types.
    fn five_atoms() -> Topology {
        let mut topo = Topology::new();
        for i in 0..5 {
            topo.moltypes[0].atoms.push(MolTypeAtom {
                name: format!("A{}", i + 1),
                residue_nr: 1,
                residue_name: "RES".into(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.lj_parameters = vec![vec![LJParameters::new(0.0, 0.0); 2]; 2];
        // excluded: (0,1) (0,2) (1,2) (2,3); 1-4: (0,3) (1,4)
        topo.exclusions = vec![vec![1, 2], vec![0, 2], vec![0, 1, 3], vec![2], vec![]];
        topo.one_four_pairs = vec![vec![3], vec![4], vec![], vec![0], vec![1]];
        topo
    }

    fn eds_atom(seq: usize) -> EdsPerturbedAtom {
        EdsPerturbedAtom {
            seq,
            iac: vec![0, 1],
            charge: vec![0.1, -0.1],
            exclusions: Vec::new(),
            one_four_pairs: Vec::new(),
        }
    }

    #[test]
    fn eds_atoms_take_their_excluded_and_one_four_pairs() {
        let mut topo = five_atoms();
        // atom 3 (seq 2) is listed before atom 1 (seq 0)
        topo.apply_eds_perturbation(
            2,
            vec!["A".into(), "B".into()],
            vec![eds_atom(2), eds_atom(0)],
        )
        .unwrap();
        let eds = &topo.eds_solute;
        assert_eq!(eds.atoms.iter().map(|a| a.seq).collect::<Vec<_>>(), [0, 2]);
        // atom 2 (listed first) owns all of its partners, the EDS-EDS pair (0,2) included
        assert_eq!(eds.atoms[1].exclusions, [0, 1, 3]);
        assert_eq!(eds.atoms[1].one_four_pairs, Vec::<usize>::new());
        // atom 0 keeps what atom 2 left: (0,1) excluded, (0,3) 1-4
        assert_eq!(eds.atoms[0].exclusions, [1]);
        assert_eq!(eds.atoms[0].one_four_pairs, [3]);
        // the standard lists no longer hold any pair of an EDS atom, and stay symmetric
        assert_eq!(
            topo.exclusions,
            vec![vec![], vec![], vec![], vec![], vec![]]
        );
        assert_eq!(
            topo.one_four_pairs,
            vec![vec![], vec![4], vec![], vec![], vec![1]]
        );
        // but the pairlist keeps every one of them out (gromosXX `all_exclusion`)
        for (i, j) in [(0, 1), (0, 2), (1, 2), (2, 3), (0, 3), (1, 4)] {
            assert!(
                topo.is_excluded_or_14(i, j) && topo.is_excluded_or_14(j, i),
                "({i},{j})"
            );
        }
        assert!(!topo.is_excluded_or_14(0, 4));
    }

    #[test]
    fn the_first_listed_eds_atom_owns_an_eds_eds_pair() {
        let mut topo = five_atoms();
        topo.apply_eds_perturbation(2, vec![], vec![eds_atom(0), eds_atom(2)])
            .unwrap();
        let eds = &topo.eds_solute;
        assert_eq!(eds.atoms[0].exclusions, [1, 2]);
        assert_eq!(eds.atoms[1].exclusions, [1, 3]);
    }

    #[test]
    fn mpertatom_rows_are_checked_against_the_topology() {
        let err = |atoms: Vec<EdsPerturbedAtom>| {
            five_atoms()
                .apply_eds_perturbation(2, vec![], atoms)
                .unwrap_err()
        };
        assert!(err(vec![eds_atom(5)]).contains("not a solute atom"));
        assert!(err(vec![eds_atom(1), eds_atom(1)]).contains("listed twice"));
        let mut a = eds_atom(1);
        a.iac = vec![0, 2];
        assert!(err(vec![a]).contains("atom type 3"));
        let mut a = eds_atom(1);
        a.charge.pop();
        assert!(err(vec![a]).contains("for 2 states"));
    }

    #[test]
    fn test_lj_parameters() {
        // For sigma = 0.34 nm and epsilon = 0.65 kJ/mol (typical for methane)
        // c6 = 4 * epsilon * sigma^6 = 4 * 0.65 * 0.34^6 = 0.00398
        // c12 = 4 * epsilon * sigma^12 = 4 * 0.65 * 0.34^12 = 4.6e-7
        // Or simpler: sigma = (c12/c6)^(1/6), epsilon = c6^2/(4*c12)

        // Use values that give known sigma/epsilon
        let sigma_expected: f64 = 0.34;
        let epsilon_expected: f64 = 0.65;
        let c6 = 4.0 * epsilon_expected * sigma_expected.powi(6);
        let c12 = 4.0 * epsilon_expected * sigma_expected.powi(12);

        let lj = LJParameters::new(c6, c12);

        // Check sigma
        let sigma = lj.sigma();
        assert!(
            (sigma - sigma_expected).abs() < 0.01,
            "Sigma mismatch: got {}",
            sigma
        );

        // Check epsilon
        let epsilon = lj.epsilon();
        assert!(
            (epsilon - epsilon_expected).abs() < 0.01,
            "Epsilon mismatch: got {}",
            epsilon
        );
    }

    #[test]
    fn test_topology_creation() {
        let mut topo = Topology::new();
        topo.solvent_atom_template = vec![SolventAtomTemplate {
            iac: 0,
            name: "OW".into(),
            mass: 15.9994,
            charge: -0.82,
        }];
        topo.solvate(100);
        assert_eq!(topo.num_atoms(), 100);
    }

    #[test]
    fn test_exclusions() {
        let mut topo = Topology::new();

        // Add 3 atoms first so num_atoms() returns 3
        for i in 0..3 {
            topo.moltypes[0].atoms.push(MolTypeAtom {
                name: format!("A{}", i + 1),
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

        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.moltypes[0].bonds.push(Bond {
            i: 1,
            j: 2,
            bond_type: 0,
            hydrogen: false,
        });
        topo.moltypes[0].angles.push(Angle {
            i: 0,
            j: 1,
            k: 2,
            angle_type: 0,
        });

        topo.resize_atom_arrays();
        topo.build_exclusions(false);

        // 0-1 should be excluded (bonded)
        assert!(topo.is_excluded(0, 1));

        // 0-2 should be excluded (1-3 neighbors)
        assert!(topo.is_excluded(0, 2));
    }

    /// A PERTATOMPAIR pair leaves the exclusion and 1-4 lists but stays off the pairlist,
    /// whatever it was in the topology (gromosXX `in_perturbation.cc:1207-1224`).
    #[test]
    fn perturbed_atom_pairs_leave_exclusion_and_one_four_lists() {
        let mut topo = Topology::new();
        for i in 0..4 {
            topo.moltypes[0].atoms.push(MolTypeAtom {
                name: format!("A{}", i + 1),
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
        topo.resize_atom_arrays();
        topo.exclusions = vec![vec![1, 2], vec![0, 2], vec![0, 1], vec![]];
        topo.one_four_pairs = vec![vec![3], vec![], vec![], vec![0]];
        let mut ps = PerturbedSolute::new();
        for (i, j) in [(0, 2), (0, 3), (1, 3)] {
            ps.atom_pairs.push(PerturbedAtomPair {
                i,
                j,
                a_type: PairInteraction::Normal,
                b_type: PairInteraction::Excluded,
            });
        }
        topo.apply_perturbation(ps, vec![false; 4]).unwrap();

        assert_eq!(topo.exclusions, vec![vec![1], vec![0, 2], vec![1], vec![]]);
        assert!(topo.one_four_pairs.iter().all(Vec::is_empty));
        for (i, j) in [(0, 2), (0, 3), (1, 3), (0, 1)] {
            assert!(topo.is_excluded_or_14(i, j), "({i}, {j}) on the pairlist");
            assert!(topo.is_excluded_or_14(j, i), "({j}, {i}) on the pairlist");
        }
        assert!(!topo.is_excluded_or_14(2, 3));
    }

    #[test]
    fn test_lj_matrix() {
        let mut topo = Topology::new();

        // Two atom types
        let lj_types = vec![
            LJParameters::new(0.001, 0.0001), // Type 0
            LJParameters::new(0.002, 0.0002), // Type 1
        ];

        topo.build_lj_matrix(&lj_types);

        // Check diagonal elements
        assert!((topo.lj_parameter(0, 0).c6 - 0.001).abs() < 1e-10);
        assert!((topo.lj_parameter(1, 1).c6 - 0.002).abs() < 1e-10);

        // Check off-diagonal (geometric mean)
        let c6_01 = (0.001_f64 * 0.002_f64).sqrt();
        assert!((topo.lj_parameter(0, 1).c6 - c6_01).abs() < 1e-10);
    }

    #[test]
    fn net_charge_of_the_perturbed_molecules() {
        // two molecules; atom 3 (the second's first) goes from +1 to 0
        let molecules = [0..2, 2..4];
        let a = [-0.5, 0.5, 1.0, 0.0];
        let b = [-0.5, 0.5, 0.0, 0.0];
        let perturbed = [false, false, true, false];
        let q = NetCharge::of_perturbed_molecules(&molecules, &a, &b, &perturbed);
        assert_eq!(q, vec![(1, NetCharge { a: 1.0, b: 0.0 })]);
        assert_eq!(q[0].1.delta(), -1.0);
        // no molecule boundaries: one molecule over every atom
        let q = NetCharge::of_perturbed_molecules(&[], &a, &b, &perturbed);
        assert_eq!(q, vec![(0, NetCharge { a: 1.0, b: 0.0 })]);
        // a sum that is zero up to rounding is exactly zero
        let q = NetCharge {
            a: 0.1 + 0.2,
            b: 0.3,
        };
        assert_eq!(q.delta().to_bits(), 0.0_f64.to_bits());
        let q = NetCharge { a: -1e-17, b: 1.0 };
        assert_eq!(q.to_string(), "A 0.000 B 1.000 (ΔQ = 1.000)");
    }
}
