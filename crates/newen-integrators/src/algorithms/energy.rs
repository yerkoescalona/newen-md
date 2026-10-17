//! Energy calculation algorithm.
//!
//! Equivalent to GROMOS `algorithm::Energy_Calculation`.
//! Finalizes total energy = E_kin + E_pot.

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::topology::Topology;

/// Finalizes total energy computation.
///
/// This is the last algorithm in the standard sequence.
/// It ensures E_total = E_kinetic + E_potential is up to date.
#[derive(Debug, Clone)]
pub struct EnergyCalculation;

impl EnergyCalculation {
    pub fn new() -> Self {
        Self
    }
}

impl Default for EnergyCalculation {
    fn default() -> Self {
        Self::new()
    }
}

impl Algorithm for EnergyCalculation {
    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        // GROMOS convention: energies are stored in old() after exchange_state.
        // Forcefield wrote potential energies to what is now old(),
        // TemperatureCalculation wrote kinetic energy to old().
        // Finalize the potential total here.
        let state = conf.old_mut();
        state.energies.update_potential_total();
        state.energies.book_totals_into_single_group();
        // gromosXX `calculate_totals` (energy.cc:397-423) on both records, PRECALCLAM part
        state.energies.precalc.calculate_totals();
        let d = &mut state.perturbed_energy_derivatives;
        d.book_totals_into_single_group();
        d.precalc.calculate_totals();
        // gromosXX energy_calculation.cc:87-89
        d.entropy_term = d.potential_total * state.energies.potential_total;
        log::debug!(
            "  E_pot={:.10e}  E_kin={:.10e}  E_tot={:.10e}",
            state.energies.potential_total,
            state.energies.kinetic_total,
            state.energies.total()
        );
        Ok(())
    }

    fn name(&self) -> &str {
        "Energy_Calculation"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn potential_total_is_the_sum_of_the_potential_terms() {
        let topo = Topology::new();
        let mut conf = Configuration::new(1, 1, 1);
        let e = &mut conf.old_mut().energies;
        e.bond_total = 1.0;
        e.angle_total = 2.0;
        e.dihedral_total = 3.0;
        e.improper_total = 4.0;
        e.lj_total = -10.0;
        e.crf_total = -20.0;
        e.kinetic_total = 100.0; // not part of the potential
        e.distanceres_total = 7.0; // special, not part of the potential
        let sim = SimulationState::new(0.002, 1);
        EnergyCalculation::new()
            .apply(&topo, &mut conf, &sim)
            .unwrap();
        assert_eq!(conf.old().energies.potential_total, -20.0);
        assert_eq!(conf.old().energies.total(), 100.0 - 20.0 + 7.0);
    }
}
