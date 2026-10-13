//! EDS end-state nonbonded terms (gromosXX `Perturbed_Nonbonded_Innerloop`, the EDS cases).
//!
//! An atom listed in `MPERTATOM` has an atom type and a charge per end state. Every pair it is
//! part of (pairlist pairs, its 1-4 pairs, its excluded pairs and its reaction-field self term)
//! is computed once per state into an [`EdsStorage`] — per-state forces, virials and energies —
//! and never into the standard force, energy or group matrices; the EDS algorithm mixes the
//! states afterwards. Pairs without an EDS atom stay with the standard kernels.
//!
//! Only the lj_crf interaction function is ported (cases 0 EDS–normal and 1 EDS–EDS of
//! `perturbed_nonbonded_innerloop.cc:84-165`); the EDS–perturbed cases are refused when the run
//! is built.

use super::{CRFParameters, EnergyGroups, ForceStorage, LJParamMatrix, LJParameters};
use gromos_core::math::{BoundaryCondition, Vec3};
use gromos_core::topology::Topology;

/// Marks an atom without end states in [`EdsAtoms::slot`].
const NOT_EDS: u32 = u32::MAX;

/// The EDS atoms of a topology in the layout the kernels read.
#[derive(Debug, Clone)]
pub struct EdsAtoms {
    /// Number of end states.
    pub num_states: usize,
    /// Per atom of the system: its index among the EDS atoms, or `u32::MAX`.
    slot: Vec<u32>,
    /// Atom index of every EDS atom, ascending.
    pub seq: Vec<usize>,
    /// `[k * num_states + s]`: type of EDS atom k in state s (0-indexed).
    iac: Vec<u32>,
    /// `[k * num_states + s]`: charge of EDS atom k in state s.
    charge: Vec<f64>,
    /// The excluded partners each EDS atom owns.
    exclusions: Vec<Vec<usize>>,
    /// The 1-4 partners each EDS atom owns.
    one_four_pairs: Vec<Vec<usize>>,
}

impl EdsAtoms {
    /// The EDS atoms of `topo` (`None` without any), for a system of `n_atoms` atoms.
    pub fn from_topology(topo: &Topology, n_atoms: usize) -> Option<Self> {
        let eds = &topo.eds_solute;
        if eds.is_empty() {
            return None;
        }
        let n = eds.num_states;
        let mut slot = vec![NOT_EDS; n_atoms];
        let mut out = Self {
            num_states: n,
            slot: Vec::new(),
            seq: Vec::with_capacity(eds.atoms.len()),
            iac: Vec::with_capacity(eds.atoms.len() * n),
            charge: Vec::with_capacity(eds.atoms.len() * n),
            exclusions: Vec::with_capacity(eds.atoms.len()),
            one_four_pairs: Vec::with_capacity(eds.atoms.len()),
        };
        for (k, a) in eds.atoms.iter().enumerate() {
            slot[a.seq] = k as u32;
            out.seq.push(a.seq);
            out.iac.extend(a.iac.iter().map(|&t| t as u32));
            out.charge.extend_from_slice(&a.charge);
            out.exclusions.push(a.exclusions.clone());
            out.one_four_pairs.push(a.one_four_pairs.clone());
        }
        out.slot = slot;
        Some(out)
    }

    /// Whether atom `i` has end states.
    #[inline]
    pub fn is_eds(&self, i: usize) -> bool {
        self.slot[i] != NOT_EDS
    }

    /// Types and charges of atom `i` in every state; `None` for an atom without end states.
    #[inline]
    fn states(&self, i: usize) -> Option<(&[u32], &[f64])> {
        let k = self.slot[i];
        (k != NOT_EDS).then(|| {
            let r = k as usize * self.num_states..(k as usize + 1) * self.num_states;
            (&self.iac[r.clone()], &self.charge[r])
        })
    }
}

/// Per-state forces, virials and energies of the EDS terms, allocated once.
#[derive(Debug, Clone, Default)]
pub struct EdsStorage {
    num_states: usize,
    n_atoms: usize,
    /// `[s * n_atoms + i]`: force on atom i in state s.
    pub forces: Vec<Vec3>,
    /// `virial[s][a][b] = Σ r_a f_b` (the `ForceStorage::virial` convention).
    pub virial: Vec<[[f64; 3]; 3]>,
    /// Energy of state s.
    pub energy: Vec<f64>,
}

impl EdsStorage {
    pub fn new(num_states: usize, n_atoms: usize) -> Self {
        Self {
            num_states,
            n_atoms,
            forces: vec![Vec3::ZERO; num_states * n_atoms],
            virial: vec![[[0.0; 3]; 3]; num_states],
            energy: vec![0.0; num_states],
        }
    }

    pub fn clear(&mut self) {
        self.forces.fill(Vec3::ZERO);
        self.virial.fill([[0.0; 3]; 3]);
        self.energy.fill(0.0);
    }

    /// Add `other` (same layout) to this storage.
    pub fn merge(&mut self, other: &EdsStorage) {
        for (d, s) in self.forces.iter_mut().zip(&other.forces) {
            *d += *s;
        }
        for (d, s) in self.virial.iter_mut().zip(&other.virial) {
            for (dr, sr) in d.iter_mut().zip(s) {
                for (x, y) in dr.iter_mut().zip(sr) {
                    *x += y;
                }
            }
        }
        for (d, s) in self.energy.iter_mut().zip(&other.energy) {
            *d += s;
        }
    }

    pub fn num_states(&self) -> usize {
        self.num_states
    }

    /// The forces of state `s`, one per atom.
    pub fn state_forces(&self, s: usize) -> &[Vec3] {
        &self.forces[s * self.n_atoms..(s + 1) * self.n_atoms]
    }
}

/// Move every pair with an EDS atom from `list` to `eds`, the EDS atom first — gromosXX's
/// `insert_pair` tries `(a1, a2)` and then `(a2, a1)` (`standard_pairlist_algorithm.cc:499-505,
/// 680-698`), for the short- and the long-range lists alike. Run after every pairlist update.
pub fn split_eds_pairs(list: &mut Vec<(u32, u32)>, atoms: &EdsAtoms, eds: &mut Vec<(u32, u32)>) {
    list.retain(|&(a, b)| {
        if atoms.is_eds(a as usize) {
            eds.push((a, b));
            false
        } else if atoms.is_eds(b as usize) {
            eds.push((b, a));
            false
        } else {
            true
        }
    });
}

/// gromosXX `Eds_Nonbonded_Term::eds_lj_crf_interaction` (`perturbed_nonbonded_term.cc:359-381`):
/// the force factor (times r) and the energy of one state.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn eds_lj_crf(
    dist2: f64,
    dist6: f64,
    disti: f64,
    c6: f64,
    c12: f64,
    q: f64,
    coulomb_scaling: f64,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
) -> (f64, f64) {
    let dist6i = 1.0 / dist6;
    let q_eps = q * four_pi_eps_i;
    let c12_dist6i = c12 * dist6i;
    let dist2i = 1.0 / dist2;
    let elj = (c12_dist6i - c6) * dist6i;
    let ecrf = q_eps * (disti * coulomb_scaling - crf.crf_2cut3i * dist2 - crf.crf_cut);
    let flj = (c12_dist6i + c12_dist6i - c6) * 6.0 * dist6i * dist2i;
    let fcrf = q_eps * (disti * coulomb_scaling * dist2i + 2.0 * crf.crf_2cut3i);
    (flj + fcrf, elj + ecrf)
}

/// One pair `(i, j)` with EDS atom `i`, in every state: case 0 (j without end states, its
/// topology type and charge) or case 1 (both EDS atoms, state by state). `one_four` selects the
/// cs6/cs12 parameters.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn eds_pair<BC: BoundaryCondition, const VIRIAL: bool, const ONE_FOUR: bool>(
    i: usize,
    j: usize,
    positions: &[Vec3],
    atoms: &EdsAtoms,
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    boundary: &BC,
    four_pi_eps_i: f64,
    coulomb_scaling: f64,
    out: &mut EdsStorage,
) {
    let Some((iac_i, q_i)) = atoms.states(i) else {
        return;
    };
    let r = boundary.nearest_image(positions[i], positions[j]);
    let dist2 = r.length_squared();
    let disti = 1.0 / dist2.sqrt();
    let dist6 = dist2 * dist2 * dist2;
    let state_j = atoms.states(j);
    let n_atoms = out.n_atoms;
    for s in 0..atoms.num_states {
        let (t_j, q_j) = match state_j {
            Some((t, q)) => (t[s], q[s]),
            None => (iac[j], charges[j]),
        };
        let lj: LJParameters = lj_params.get(iac_i[s] as usize, t_j as usize);
        let (c6, c12) = if ONE_FOUR {
            (lj.cs6, lj.cs12)
        } else {
            (lj.c6, lj.c12)
        };
        let (f, e) = eds_lj_crf(
            dist2,
            dist6,
            disti,
            c6,
            c12,
            q_i[s] * q_j,
            coulomb_scaling,
            crf,
            four_pi_eps_i,
        );
        let term = r * f;
        out.forces[s * n_atoms + i] += term;
        out.forces[s * n_atoms + j] -= term;
        if VIRIAL {
            let (rv, tv) = ([r.x, r.y, r.z], [term.x, term.y, term.z]);
            let v = &mut out.virial[s];
            for (b, rb) in rv.iter().enumerate() {
                for (a, ta) in tv.iter().enumerate() {
                    v[b][a] += rb * ta;
                }
            }
        }
        out.energy[s] += e;
    }
}

/// The pairlist pairs of the EDS atoms (each `(i, j)` with EDS atom `i`, as [`split_eds_pairs`]
/// leaves them), every state, into `out`; the per-state virial only with `VIRIAL` (gromosXX adds
/// the pairs' end-state virials only with a virial, `perturbed_nonbonded_set.cc:306-318, 367-376`).
#[allow(clippy::too_many_arguments)]
pub fn eds_pair_interactions<BC: BoundaryCondition, const VIRIAL: bool>(
    pairs: &[(u32, u32)],
    positions: &[Vec3],
    atoms: &EdsAtoms,
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    boundary: &BC,
    four_pi_eps_i: f64,
    out: &mut EdsStorage,
) {
    for &(i, j) in pairs {
        eds_pair::<BC, VIRIAL, false>(
            i as usize,
            j as usize,
            positions,
            atoms,
            charges,
            iac,
            lj_params,
            crf,
            boundary,
            four_pi_eps_i,
            1.0,
            out,
        );
    }
}

/// The 1-4 pairs the EDS atoms own, in ascending EDS atom order (gromosXX
/// `_perturbed_one_four_outerloop`, `perturbed_nonbonded_outerloop.cc:151-205`): cs6/cs12 and
/// the Coulomb scaling of the standard 1-4 term; the per-state virial always.
#[allow(clippy::too_many_arguments)]
pub fn eds_one_four_interactions<BC: BoundaryCondition>(
    positions: &[Vec3],
    atoms: &EdsAtoms,
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    boundary: &BC,
    four_pi_eps_i: f64,
    coulomb_scaling: f64,
    out: &mut EdsStorage,
) {
    for (k, &i) in atoms.seq.iter().enumerate() {
        for &j in &atoms.one_four_pairs[k] {
            eds_pair::<BC, true, true>(
                i,
                j,
                positions,
                atoms,
                charges,
                iac,
                lj_params,
                crf,
                boundary,
                four_pi_eps_i,
                coulomb_scaling,
                out,
            );
        }
    }
}

/// gromosXX `eds_rf_interaction` (`perturbed_nonbonded_term.cc:547-559`): the force and energy
/// of the distance-dependent reaction-field term of an excluded pair with charge product `q`.
#[inline(always)]
fn eds_rf(r: Vec3, q: f64, crf: &CRFParameters, four_pi_eps_i: f64) -> (Vec3, f64) {
    let dist2 = r.length_squared();
    let f = r * (q * four_pi_eps_i * 2.0 * crf.crf_2cut3i);
    let e = q * four_pi_eps_i * (-crf.crf_2cut3i * dist2 - crf.crf_cut);
    (f, e)
}

/// The reaction-field self and excluded-pair terms of the EDS atoms (gromosXX
/// `eds_RF_excluded_interaction_innerloop` and `eds_perturbed_RF_exclusions_loop`,
/// `perturbed_nonbonded_innerloop.cc:3243-3337, 3470-3600`), with NSLFEXCL.
///
/// The standard self term the reaction-field loop booked with the atom's topology charge is
/// taken back out of `standard` (its CRF total and group cell `[eg(i)][eg(i)]`); every state gets
/// `½·e_rf(q_s²)` and, for each owned excluded partner, the term with the state's charges. The
/// per-state virial is always accumulated.
#[allow(clippy::too_many_arguments)]
pub fn eds_rf_excluded_interactions<BC: BoundaryCondition>(
    positions: &[Vec3],
    atoms: &EdsAtoms,
    charges: &[f64],
    crf: &CRFParameters,
    boundary: &BC,
    four_pi_eps_i: f64,
    out: &mut EdsStorage,
    standard: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    let n_atoms = out.n_atoms;
    for (k, &i) in atoms.seq.iter().enumerate() {
        let (_, e_self) = eds_rf(Vec3::ZERO, charges[i] * charges[i], crf, four_pi_eps_i);
        standard.e_crf -= 0.5 * e_self;
        if egroups.is_grouped() {
            standard.groups.crf[egroups.cell(i, i)] -= 0.5 * e_self;
        }
        let q_i = &atoms.charge[k * atoms.num_states..(k + 1) * atoms.num_states];
        for (s, &q) in q_i.iter().enumerate() {
            let (_, e) = eds_rf(Vec3::ZERO, q * q, crf, four_pi_eps_i);
            out.energy[s] += 0.5 * e;
        }
        for &j in &atoms.exclusions[k] {
            let r = boundary.nearest_image(positions[i], positions[j]);
            let q_j = atoms.states(j).map(|(_, q)| q);
            for (s, &qi) in q_i.iter().enumerate() {
                let qj = q_j.map_or(charges[j], |q| q[s]);
                let (f, e) = eds_rf(r, qi * qj, crf, four_pi_eps_i);
                out.energy[s] += e;
                out.forces[s * n_atoms + i] += f;
                out.forces[s * n_atoms + j] -= f;
                let (rv, fv) = ([r.x, r.y, r.z], [f.x, f.y, f.z]);
                let v = &mut out.virial[s];
                for (a, ra) in rv.iter().enumerate() {
                    for (b, fb) in fv.iter().enumerate() {
                        v[a][b] += ra * fb;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::math::Vacuum;
    use gromos_core::topology::{EdsPerturbedAtom, EdsSolute};

    fn lj() -> LJParamMatrix {
        let p = |c6: f64, c12: f64| LJParameters {
            c6,
            c12,
            cs6: 0.5 * c6,
            cs12: 0.5 * c12,
        };
        LJParamMatrix::from_nested(&[
            vec![p(2.0e-3, 2.0e-6), p(3.0e-3, 4.0e-6)],
            vec![p(3.0e-3, 4.0e-6), p(5.0e-3, 9.0e-6)],
        ])
    }

    /// Atoms 0 and 2 have two end states; atom 1 has none. Atom 0 owns its exclusion with 2.
    fn atoms() -> EdsAtoms {
        let mut topo = Topology::new();
        topo.eds_solute = EdsSolute {
            num_states: 2,
            state_names: vec![],
            atoms: vec![
                EdsPerturbedAtom {
                    seq: 0,
                    iac: vec![0, 1],
                    charge: vec![0.4, -0.2],
                    exclusions: vec![2],
                    one_four_pairs: vec![1],
                },
                EdsPerturbedAtom {
                    seq: 2,
                    iac: vec![1, 0],
                    charge: vec![-0.3, 0.1],
                    exclusions: vec![],
                    one_four_pairs: vec![],
                },
            ],
        };
        EdsAtoms::from_topology(&topo, 3).unwrap()
    }

    #[test]
    fn pairs_with_an_eds_atom_move_with_the_eds_atom_first() {
        let a = atoms();
        let mut list = vec![(0, 1), (1, 2), (1, 1), (2, 0)];
        let mut eds = Vec::new();
        split_eds_pairs(&mut list, &a, &mut eds);
        assert_eq!(list, [(1, 1)]);
        assert_eq!(eds, [(0, 1), (2, 1), (2, 0)]);
    }

    /// Each state of a pair is the plain LJ + reaction-field pair with that state's parameters.
    #[test]
    fn a_state_is_the_pair_term_with_the_states_type_and_charge() {
        let a = atoms();
        let crf = CRFParameters::new(1.4, 1.0, 61.0, 0.0);
        let pos = [
            Vec3::ZERO,
            Vec3::new(0.31, 0.05, -0.02),
            Vec3::new(-0.1, 0.22, 0.07),
        ];
        let charges = [0.9, 0.25, -0.7]; // topology charges: unused for EDS atoms
        let iac = [1, 1, 0];
        let fpe = 138.9354;
        let mut out = EdsStorage::new(2, 3);
        eds_pair_interactions::<_, true>(
            &[(0, 1), (0, 2)],
            &pos,
            &a,
            &charges,
            &iac,
            &lj(),
            &crf,
            &Vacuum,
            fpe,
            &mut out,
        );
        let pair = |ri: Vec3, rj: Vec3, ti: usize, tj: usize, q: f64| {
            let p = lj().get(ti, tj);
            let r = (ri - rj).length();
            let e_lj = p.c12 / r.powi(12) - p.c6 / r.powi(6);
            let e_crf = q * fpe * (1.0 / r - crf.crf_2cut3i * r * r - crf.crf_cut);
            e_lj + e_crf
        };
        // state 0: atom 0 type 0 charge 0.4; atom 2 type 1 charge -0.3; atom 1 its topology
        let e0 = pair(pos[0], pos[1], 0, 1, 0.4 * 0.25) + pair(pos[0], pos[2], 0, 1, 0.4 * -0.3);
        let e1 = pair(pos[0], pos[1], 1, 1, -0.2 * 0.25) + pair(pos[0], pos[2], 1, 0, -0.2 * 0.1);
        assert!((out.energy[0] - e0).abs() < 1e-12, "{} {e0}", out.energy[0]);
        assert!((out.energy[1] - e1).abs() < 1e-12, "{} {e1}", out.energy[1]);
        // Newton's third law per state, and the virial trace is Σ r·f
        for s in 0..2 {
            let f = out.state_forces(s);
            assert!((f[0] + f[1] + f[2]).length() < 1e-12);
        }
    }

    #[test]
    fn the_rf_self_term_moves_from_the_standard_total_to_the_states() {
        let a = atoms();
        let crf = CRFParameters::new(1.4, 1.0, 61.0, 0.0);
        let pos = [
            Vec3::ZERO,
            Vec3::new(0.3, 0.0, 0.0),
            Vec3::new(0.0, 0.2, 0.0),
        ];
        let charges = [0.9, 0.25, -0.7];
        let fpe = 138.9354;
        let mut out = EdsStorage::new(2, 3);
        let mut standard = ForceStorage::new(3);
        eds_rf_excluded_interactions(
            &pos,
            &a,
            &charges,
            &crf,
            &Vacuum,
            fpe,
            &mut out,
            &mut standard,
            EnergyGroups::single(),
        );
        let self_term = |q: f64| -0.5 * q * q * fpe * crf.crf_cut;
        let excl = |q: f64| q * fpe * (-crf.crf_2cut3i * 0.2 * 0.2 - crf.crf_cut);
        let removed = self_term(0.9) + self_term(-0.7);
        assert!((standard.e_crf + removed).abs() < 1e-12);
        let s0 = self_term(0.4) + self_term(-0.3) + excl(0.4 * -0.3);
        let s1 = self_term(-0.2) + self_term(0.1) + excl(-0.2 * 0.1);
        assert!((out.energy[0] - s0).abs() < 1e-12, "{} {s0}", out.energy[0]);
        assert!((out.energy[1] - s1).abs() < 1e-12, "{} {s1}", out.energy[1]);
    }
}
