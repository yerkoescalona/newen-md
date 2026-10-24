//! `ElectrostaticEmbedding` — the QM↔MM Coulomb coupling term, with forces on **both** sides.
//!
//! ** P2.7 Step 3.** This is the concrete answer to "how do forces get back onto MM
//! atoms": the QM zone carries charges (derived by whatever engine owns it), the environment
//! carries the force field's partial charges, and the Coulomb interaction between them produces
//! energy plus equal-and-opposite forces on both sets. P5 called for a provider that
//! contributes to an *arbitrary* atom subset rather than only its own region — this is the first
//! provider that actually does so, putting forces on MM atoms outside `region`.
//!
//! **Which of GROMOS's three EE force paths this is.** Poliak et al. 2025 (*J. Comput. Chem.*
//! 46:e70053) describe three ways forces reach MM atoms, depending on what the QM program
//! returns: (a) the program hands back MM forces directly (MNDO, DFTB, ORCA, Turbomole, xtb),
//! (b) it returns the electric field at each MM site and GROMOS computes `f_i = q_i·E_i`
//! (Gaussian), or (c) it returns QM partial charges and **GROMOS computes the pairwise Coulomb
//! forces itself** (MOPAC). This implements **(c)**.
//!
//! **This is an alternative to path (a), not something to combine with it.**
//! [`crate::nonbonded::XtbInteraction`]'s `Embedding::Electrostatic` mode now implements path (a)
//! directly (xtb's own `pcharge`/`pcgrad` files — real SCF polarization, real MM forces from the
//! engine itself, confirmed against a real gromosXX-equivalent run). Registering *both* this
//! provider and an xtb-electrostatic `XtbInteraction` on the same QM/MM boundary computes the
//! coupling twice — the exact double-counting hazard `zones.rs` (assumption A5) exists to guard
//! against. Use this provider only for engines that don't return MM gradients (confirmed against
//! gromosXX's real source: only `xtb_worker.cc`, `orca_worker.cc`, `turbomole_worker.cc`
//! implement `parse_mm_gradients` — MOPAC/MNDO/DFTB do not, so they need this path).
//!
//! **Deliberate scope limits**, so this isn't mistaken for full QM/MM:
//! - **Charges are an input, not a prediction.** `region_charges` is supplied by the caller, not
//!   refreshed from an engine each step. **Correction, found by reading gromosXX's real source
//!   (`interaction/qmmm/qm_zone.cc`) rather than trusting secondhand paper summaries:** there is
//!   no separate "extra force" correction term for fluctuating charges anywhere in gromosXX
//!   (`grep -rn "MEDC"` over the whole source returns nothing — an earlier version of this note
//!   cited one that doesn't exist in the implementation). What gromosXX actually does is simpler
//!   and narrower: **dynamic per-step charges exist only for *mechanical* embedding**
//!   (`qm_zone.cc:177-196`, gated on `qmmm_mechanical && qm_ch_dynamic`) — the engine's charges
//!   are written straight into `topo.charge()` and the ordinary classical loop picks them up next
//!   step, with the `dq/dR` derivative simply neglected (an accepted approximation of that
//!   scheme, not a solved problem). For *electrostatic* embedding specifically, gromosXX doesn't
//!   use this path's static-QM-charge approach at all — it uses path (a) instead, which needs no
//!   charge-refresh logic because the QM program differentiates the whole coupled energy itself.
//! - **Plain Coulomb, not reaction-field.** GROMOS's classical MM–MM electrostatics use RF
//!   (`CRFParameters`); the QM↔MM embedding term here is bare cutoff Coulomb, matching the
//!   point-charge picture. Whether GROMOS applies RF to the QM/MM cross term is a separate
//!   question this does not answer or assume.
//! - **No mutual polarization.** Poliak: EE "includes electronic polarization of the QM zone by
//!   the MM atoms, but not vice versa."
//! - **No link atoms.** The tutorial system has `QMLI=0` throughout (no bond crosses the
//!   boundary), and Poliak states boundary charge redistribution is unimplemented in gromosXX
//! too — see P2.7 A4.

use newen_core::configuration::Configuration;
use newen_core::math::{Mat3, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::SpatialIndex;
use newen_core::topology::Topology;

use crate::provider::{Contribution, Embedding, PotentialProvider, ProviderError, ProviderExtra};

/// QM↔MM electrostatic embedding via explicit pairwise Coulomb (Poliak path (c)).
#[derive(Debug, Clone)]
pub struct ElectrostaticEmbedding {
    /// Cutoff for gathering MM point charges around the region (GROMOS `RCUTQM`), in nm.
    pub cutoff: f64,
    /// QM-derived partial charge per **global** atom index. Only entries for atoms in `region`
    /// are read; sized like `Topology::charge` so indices need no remapping.
    region_charges: Vec<f64>,
}

impl ElectrostaticEmbedding {
    /// `region_charges` is indexed by global atom index (same convention as `Topology::charge`).
    pub fn new(cutoff: f64, region_charges: Vec<f64>) -> Self {
        Self {
            cutoff,
            region_charges,
        }
    }
}

impl PotentialProvider for ElectrostaticEmbedding {
    fn contribute(
        &mut self,
        region: &AtomSelection,
        topo: &Topology,
        conf: &Configuration,
        neigh: &dyn SpatialIndex,
    ) -> Result<Contribution, ProviderError> {
        if region.is_empty() {
            return Err(ProviderError::InvalidRegion("empty region".to_string()));
        }
        let n_atoms = topo.num_atoms();
        if self.region_charges.len() < n_atoms {
            return Err(ProviderError::InvalidRegion(format!(
                "region_charges has {} entries, need {} (indexed by global atom index)",
                self.region_charges.len(),
                n_atoms
            )));
        }

        let mut in_region = vec![false; n_atoms];
        for i in region.iter() {
            in_region[i] = true;
        }

        let positions = &conf.current().pos;
        let four_pi_eps_i = newen_core::units::four_pi_eps_i;

        // Sparse force accumulator over *all* touched atoms — QM and MM alike (P5).
        let mut force_acc: Vec<Vec3> = vec![Vec3::ZERO; n_atoms];
        let mut touched = vec![false; n_atoms];
        let mut energy = 0.0;
        // Atomic pairwise virial, same convention as `ForceStorage::virial` in
        // `nonbonded/innerloops.rs`: `W[a][b] = sum_pairs r[a] * f[b]`, `r = r_qm - r_mm`,
        // `f = f_qm` (force on the QM atom). P2.8-3.
        let mut virial = Mat3::ZERO;

        for (lo, hi, shift) in neigh.neighbor_pairs(region, self.cutoff) {
            // Only cross-boundary pairs are embedding pairs. Intra-region pairs are the QM
            // engine's business (double-counting them here is exactly the P2.7 A5 hazard);
            // `neighbor_pairs` never returns pairs with both endpoints outside `region`.
            let (qm, mm) = match (in_region[lo], in_region[hi]) {
                (true, false) => (lo, hi),
                (false, true) => (hi, lo),
                _ => continue,
            };

            let q_qm = self.region_charges[qm];
            let q_mm = topo.charge[mm];
            if q_qm == 0.0 || q_mm == 0.0 {
                continue;
            }

            // `shift` is defined so that `pos[lo] - (pos[hi] + shift)` is the minimum-image
            // vector; recover it in QM→MM orientation regardless of which endpoint was `lo`.
            let r_lo_hi = (positions[lo] - positions[hi]) - shift;
            let r_vec = if qm == lo { r_lo_hi } else { -r_lo_hi };

            let r2 = r_vec.length_squared();
            if r2 < 1e-12 {
                return Err(ProviderError::ComputationFailed(format!(
                    "QM atom {qm} and MM atom {mm} are coincident; Coulomb term is singular"
                )));
            }
            let r = r2.sqrt();

            // E = q_qm q_mm / (4πε₀ r);  F_qm = q_qm q_mm r̂ / (4πε₀ r²) = ... * r_vec / r³
            let qq = q_qm * q_mm * four_pi_eps_i;
            energy += qq / r;

            let f_qm = r_vec * (qq / (r2 * r));
            force_acc[qm] += f_qm;
            force_acc[mm] -= f_qm; // Newton's third law: the MM atom feels the reaction
            touched[qm] = true;
            touched[mm] = true;

            // Outer product r_vec ⊗ f_qm, matching `ForceStorage`'s `virial[a][b] += r[a]*f[b]`:
            // glam's `Mat3::from_cols` takes columns, and column `b` of that matrix is
            // `r_vec * f_qm[b]` (row `a` of column `b` is `r_vec[a] * f_qm[b]`).
            virial += Mat3::from_cols(r_vec * f_qm.x, r_vec * f_qm.y, r_vec * f_qm.z);
        }

        let forces = (0..n_atoms)
            .filter(|&i| touched[i])
            .map(|i| (i, force_acc[i]))
            .collect();

        Ok(Contribution {
            energy,
            forces,
            virial,
            extra: ProviderExtra::default(),
        })
    }

    fn name(&self) -> &str {
        "QM/MM electrostatic embedding"
    }

    fn embedding(&self) -> Embedding {
        Embedding::Electrostatic
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::{Periodicity, Vacuum};
    use newen_core::spatial_index::ConfigurationSpatialIndex;
    use newen_core::topology::MolTypeAtom;

    /// Two atoms: index 0 is the "QM" region, index 1 is an "MM" point charge.
    fn two_charge_system(q_qm: f64, q_mm: f64, separation: f64) -> (Topology, Configuration) {
        let mut topo = Topology::new();
        topo.charge = vec![0.0, q_mm]; // MM charge lives in the topology, QM charge does not
        topo.iac = vec![0, 0];
        topo.mass = vec![12.0, 16.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 16.0];
        topo.exclusions = vec![Vec::new(), Vec::new()];
        topo.moltypes[0].atoms = (0..2)
            .map(|i| MolTypeAtom {
                name: format!("A{i}"),
                residue_nr: 1,
                residue_name: "MOL".to_string(),
                iac: 0,
                mass: 12.0,
                charge: topo.charge[i],
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            })
            .collect();
        let _ = q_qm;

        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos = vec![Vec3::ZERO, Vec3::new(separation, 0.0, 0.0)];
        (topo, conf)
    }

    /// Strongest possible check: two point charges have an exactly known Coulomb energy and
    /// force, so this is a closed-form oracle, not a self-consistency check.
    #[test]
    fn matches_analytic_coulomb_for_two_point_charges() {
        let (q_qm, q_mm, r) = (0.5, -0.8, 0.4);
        let (topo, conf) = two_charge_system(q_qm, q_mm, r);

        let mut provider = ElectrostaticEmbedding::new(1.4, vec![q_qm, 0.0]);
        let region = AtomSelection::from_indices(vec![0], 2).unwrap();
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let c = provider
            .contribute(&region, &topo, &conf, &index)
            .expect("embedding should succeed");

        let expected_energy = q_qm * q_mm * newen_core::units::four_pi_eps_i / r;
        assert!(
            (c.energy - expected_energy).abs() < 1e-10,
            "energy {} vs analytic {}",
            c.energy,
            expected_energy
        );

        // Opposite charges attract: the QM atom at the origin is pulled toward +x.
        let expected_force_mag = (q_qm * q_mm * newen_core::units::four_pi_eps_i / (r * r)).abs();
        let f_qm = c.forces.iter().find(|(i, _)| *i == 0).unwrap().1;
        let f_mm = c.forces.iter().find(|(i, _)| *i == 1).unwrap().1;

        assert!(
            (f_qm.x - expected_force_mag).abs() < 1e-10,
            "f_qm = {f_qm:?}"
        );
        assert!(
            (f_qm + f_mm).length() < 1e-12,
            "forces must be equal and opposite, got {f_qm:?} and {f_mm:?}"
        );
    }

    /// ** P2.8-3, closed-form half.** For a pure `q_qm*q_mm/r` pair (Coulomb, no cutoff
    /// switching), Euler's theorem for a homogeneous degree-(-1) potential gives
    /// `r · dE/dr = -E`, and the atomic pairwise virial's trace is exactly `r_vec · f_qm`
    /// (`f_qm = -dE/dr_vec`, so `r_vec · f_qm = -r·dE/dr = E`). So `trace(virial) == energy`
    /// holds *exactly* here, not just approximately — a stronger, cheaper check than finite
    /// differences, done first for the same reason P2.7 Step 3 did the analytic force check
    /// before the FD one.
    #[test]
    fn virial_trace_equals_energy_for_a_pure_coulomb_pair() {
        let (q_qm, q_mm, r) = (0.5, -0.8, 0.4);
        let (topo, conf) = two_charge_system(q_qm, q_mm, r);

        let mut provider = ElectrostaticEmbedding::new(1.4, vec![q_qm, 0.0]);
        let region = AtomSelection::from_indices(vec![0], 2).unwrap();
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let c = provider
            .contribute(&region, &topo, &conf, &index)
            .expect("embedding should succeed");

        let trace = c.virial.x_axis.x + c.virial.y_axis.y + c.virial.z_axis.z;
        assert!(
            (trace - c.energy).abs() < 1e-10,
            "trace(virial) = {trace} vs energy = {} (should be exactly equal for a 1/r pair)",
            c.energy
        );
    }

    /// ** P2.8-3, finite-difference half.** Scaling both atoms' positions away from a
    /// fixed center by `s` is the vacuum analogue of isotropic volume scaling (`dV/V = 3 ds/s`
    /// at `s=1`) — the standard virial-pressure relation `P ~ -dE/dV` reduces, for this one pair,
    /// to `trace(virial) = r_vec · f_qm = -dE/ds|_{s=1}` (verified analytically above; checked
    /// here against a genuine numerical derivative instead of trusting the derivation). Distinct
    /// from the force FD check in `force_on_mm_atom_matches_finite_difference`, which perturbs
    /// one atom at fixed separation; this perturbs the *separation itself*.
    #[test]
    fn virial_trace_matches_finite_difference_scaling_derivative() {
        let (q_qm, q_mm, r0) = (0.5, -0.8, 0.4);
        let region = AtomSelection::from_indices(vec![0], 2).unwrap();
        let periodicity = Periodicity::Vacuum(Vacuum);

        let energy_at_scale = |s: f64| -> f64 {
            let (topo, conf) = two_charge_system(q_qm, q_mm, s * r0);
            let mut provider = ElectrostaticEmbedding::new(1.4, vec![q_qm, 0.0]);
            let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
            provider
                .contribute(&region, &topo, &conf, &index)
                .unwrap()
                .energy
        };

        let h = 1e-6;
        let fd_denergy_ds = (energy_at_scale(1.0 + h) - energy_at_scale(1.0 - h)) / (2.0 * h);

        let (topo, conf) = two_charge_system(q_qm, q_mm, r0);
        let mut provider = ElectrostaticEmbedding::new(1.4, vec![q_qm, 0.0]);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let c = provider.contribute(&region, &topo, &conf, &index).unwrap();
        let trace = c.virial.x_axis.x + c.virial.y_axis.y + c.virial.z_axis.z;

        assert!(
            (trace - (-fd_denergy_ds)).abs() < 1e-6,
            "trace(virial) = {trace} vs -dE/ds (finite diff) = {}",
            -fd_denergy_ds
        );
    }

    /// The Step 3 exit criterion: forces land on **MM atoms outside the region**, and they are
    /// the true gradient of the energy with respect to *that MM atom's* position.
    #[test]
    fn force_on_mm_atom_matches_finite_difference() {
        let (q_qm, q_mm, r) = (0.5, -0.8, 0.4);
        let (topo, conf) = two_charge_system(q_qm, q_mm, r);

        let mut provider = ElectrostaticEmbedding::new(1.4, vec![q_qm, 0.0]);
        let region = AtomSelection::from_indices(vec![0], 2).unwrap();
        let periodicity = Periodicity::Vacuum(Vacuum);

        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let c = provider
            .contribute(&region, &topo, &conf, &index)
            .expect("embedding should succeed");
        let f_mm = c.forces.iter().find(|(i, _)| *i == 1).unwrap().1;
        assert!(
            c.forces.iter().any(|&(i, _)| i == 1),
            "an MM atom outside the region must receive a force"
        );

        let mut energy_at = |positions: &[Vec3]| -> f64 {
            let mut c2 = Configuration::new(2, 1, 1);
            c2.current_mut().pos = positions.to_vec();
            let idx = ConfigurationSpatialIndex::new(&c2, &periodicity);
            provider
                .contribute(&region, &topo, &c2, &idx)
                .unwrap()
                .energy
        };

        let h = 1e-6;
        let base = conf.current().pos.clone();
        for axis in 0..3 {
            let delta = match axis {
                0 => Vec3::new(h, 0.0, 0.0),
                1 => Vec3::new(0.0, h, 0.0),
                _ => Vec3::new(0.0, 0.0, h),
            };
            let mut plus = base.clone();
            let mut minus = base.clone();
            plus[1] += delta; // perturb the MM atom, not the QM atom
            minus[1] -= delta;

            let fd = -(energy_at(&plus) - energy_at(&minus)) / (2.0 * h);
            let analytic = match axis {
                0 => f_mm.x,
                1 => f_mm.y,
                _ => f_mm.z,
            };
            assert!(
                (fd - analytic).abs() < 1e-4,
                "MM-atom force axis {axis}: finite-diff {fd} vs provider {analytic}"
            );
        }
    }

    /// Intra-region pairs must not contribute — they belong to the QM engine, and counting
    /// them here is the P2.7 A5 double-counting hazard in miniature.
    #[test]
    fn intra_region_pairs_are_not_counted() {
        let (topo, conf) = two_charge_system(0.5, -0.8, 0.4);
        let mut provider = ElectrostaticEmbedding::new(1.4, vec![0.5, -0.8]);

        // Both atoms in the region: there is no QM↔MM pair left, so no embedding energy.
        let region = AtomSelection::all(2);
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let c = provider.contribute(&region, &topo, &conf, &index).unwrap();
        assert_eq!(c.energy, 0.0, "intra-region pairs must not be embedded");
        assert!(c.forces.is_empty());
    }
}
