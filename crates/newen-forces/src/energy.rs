//! Single-point energy calculation.
//!
//! Computes the full potential energy (bonded + nonbonded) for a given set of
//! positions without running an MD step.  Designed to be called from analysis
//! tools (e.g. `ener`) without depending on `newen-integrators`.
//!
//! # Usage
//! ```rust,ignore
//! use newen_forces::energy::{single_point_energy, EnergyParams};
//!
//! let params = EnergyParams::from_imd(&imd);
//! let e = single_point_energy(&topo, &positions, box_dims, &params);
//! println!("E_pot = {} kJ/mol", e.potential_total);
//! ```

use newen_core::{
    configuration::{Box as SimBox, Configuration},
    math::{Periodicity, Rectangular, Vacuum, Vec3},
    pairlist::{PairlistContainer, StandardPairlistAlgorithm},
    topology::Topology,
};

use crate::{
    bonded::calculate_bonded_forces_ntf,
    nonbonded::{
        lj_crf_innerloop_novirial, one_four_interaction_loop, rf_excluded_interactions,
        solvent_innerloop_novirial, CRFParameters, EnergyGroups, ForceStorage, LJParamMatrix,
        LJParameters,
    },
};

/// Parameters needed for single-point energy evaluation.
#[derive(Debug, Clone)]
pub struct EnergyParams {
    /// Short-range (and long-range) cutoff in nm.
    pub cutoff: f64,
    /// Interior dielectric constant (always 1.0 in GROMOS).
    pub epsilon: f64,
    /// Reaction-field dielectric constant.
    pub rf_epsilon: f64,
    /// Reaction-field inverse Debye screening length (nm⁻¹).
    pub rf_kappa: f64,
    /// Pairlist update every N steps (use 1 for single-point).
    pub pairlist_freq: usize,
    /// NTF flags: [bond, angle, improper, dihedral]
    pub ntf: [bool; 4],
    /// Atoms per solvent molecule (3 for SPC/TIP3P).
    pub atoms_per_solvent: usize,
    /// Use quartic bond potential (true) or harmonic (false).
    pub quartic_bonds: bool,
}

impl Default for EnergyParams {
    fn default() -> Self {
        Self {
            cutoff: 1.4,
            epsilon: 1.0,
            rf_epsilon: 78.5,
            rf_kappa: 0.0,
            pairlist_freq: 1,
            ntf: [true, true, true, true],
            atoms_per_solvent: 3,
            quartic_bonds: true,
        }
    }
}

/// Energy components from a single-point calculation.
#[derive(Debug, Clone, Default)]
pub struct SinglePointEnergy {
    pub bond: f64,
    pub angle: f64,
    pub dihedral: f64,
    pub improper: f64,
    pub lj: f64,
    pub crf: f64,
    pub special: f64,
    pub potential: f64,
}

impl SinglePointEnergy {
    fn update_potential(&mut self) {
        self.potential = self.bond
            + self.angle
            + self.dihedral
            + self.improper
            + self.lj
            + self.crf
            + self.special;
    }
}

/// Compute the full potential energy for a single configuration.
///
/// Builds a fresh pairlist, evaluates all bonded and nonbonded terms,
/// and returns decomposed energy components.
pub fn single_point_energy(
    topo: &Topology,
    positions: &[Vec3],
    box_dims: Vec3,
    params: &EnergyParams,
) -> SinglePointEnergy {
    let n_atoms = topo.num_atoms();

    // Build a minimal Configuration wrapping the positions
    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = positions.to_vec();
    conf.current_mut().vel = vec![Vec3::ZERO; n_atoms];
    conf.current_mut().force = vec![Vec3::ZERO; n_atoms];
    if box_dims.x > 0.0 {
        conf.current_mut().box_config = SimBox::rectangular(box_dims.x, box_dims.y, box_dims.z);
    }

    // Periodicity
    let periodicity = if box_dims.x > 0.0 {
        Periodicity::Rectangular(Rectangular::new(box_dims))
    } else {
        Periodicity::Vacuum(Vacuum)
    };

    // Build pairlist
    let mut pairlist = PairlistContainer::new(params.cutoff, params.cutoff, 0.0);
    pairlist.update_frequency = 1;
    let pl_algo = StandardPairlistAlgorithm::new(!topo.chargegroups.is_empty());
    pl_algo.update(topo, &conf, &mut pairlist, &periodicity);
    pairlist.step();

    let charges = topo.charge.clone();
    let iac_u32: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let lj_nested: Vec<Vec<LJParameters>> = topo
        .lj_parameters
        .iter()
        .map(|r| r.iter().map(LJParameters::from).collect())
        .collect();
    let lj_mat = LJParamMatrix::from_nested(&lj_nested);
    let crf = CRFParameters::new(
        params.cutoff,
        params.epsilon,
        params.rf_epsilon,
        params.rf_kappa,
    );

    let pairlist_u32: Vec<(u32, u32)> = pairlist
        .solute_short
        .iter()
        .chain(pairlist.solute_long.iter())
        .copied()
        .collect();
    let solvent_u32: Vec<(u32, u32)> = pairlist
        .solvent_short
        .iter()
        .chain(pairlist.solvent_long.iter())
        .copied()
        .collect();

    // ── Bonded ──────────────────────────────────────────────────────────────
    let bonded = calculate_bonded_forces_ntf(
        topo,
        &conf,
        &periodicity,
        params.quartic_bonds,
        crate::bonded::CovalentForm::default(),
        params.ntf[0],
        params.ntf[1],
        params.ntf[3],
        params.ntf[2],
        &mut crate::bonded::BondedGroupEnergies::default(),
    );

    // ── Nonbonded ────────────────────────────────────────────────────────────
    let mut nb = ForceStorage::new(n_atoms);

    if !pairlist_u32.is_empty() {
        lj_crf_innerloop_novirial(
            &conf.current().pos,
            &charges,
            &iac_u32,
            &pairlist_u32,
            &lj_mat,
            &crf,
            &periodicity,
            newen_core::units::four_pi_eps_i,
            &mut nb,
            EnergyGroups::single(),
        );
    }
    if !solvent_u32.is_empty() {
        solvent_innerloop_novirial(
            &conf.current().pos,
            &charges,
            &iac_u32,
            &solvent_u32,
            &lj_mat,
            &crf,
            &periodicity,
            params.atoms_per_solvent,
            newen_core::units::four_pi_eps_i,
            &mut nb,
            EnergyGroups::single(),
        );
    }

    rf_excluded_interactions(
        &charges,
        &topo.exclusions,
        &conf.current().pos,
        &crf,
        &periodicity,
        newen_core::units::four_pi_eps_i,
        &mut nb,
        topo.num_solute_atoms(),
        EnergyGroups::single(),
    );

    if !topo.one_four_pairs.is_empty() {
        one_four_interaction_loop(
            &topo.one_four_pairs,
            &conf.current().pos,
            &charges,
            &iac_u32,
            &lj_mat,
            &crf,
            &periodicity,
            newen_core::units::four_pi_eps_i,
            &mut nb,
            1.0,
            EnergyGroups::single(),
        );
    }

    // ── Assemble ─────────────────────────────────────────────────────────────
    let mut e = SinglePointEnergy {
        bond: bonded.energy,
        angle: 0.0, // split not yet available from ForceEnergy
        dihedral: 0.0,
        improper: 0.0,
        lj: nb.e_lj,
        crf: nb.e_crf,
        special: 0.0,
        potential: 0.0,
    };
    e.update_potential();
    // Bond total already includes angle+dihedral+improper from calculate_bonded_forces_ntf
    // Adjust labelling: report everything under bond for now
    e.potential = bonded.energy + nb.e_lj + nb.e_crf;
    e
}
