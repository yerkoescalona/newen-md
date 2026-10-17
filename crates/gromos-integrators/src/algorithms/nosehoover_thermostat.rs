//! Nosé-Hoover (chain) thermostat algorithm.
//!
//! Implements both single NHC (nhc=1) and NHC chain (nhc≥2), equivalent to
//! GROMOS `algorithm::NoseHoover_Thermostat`.
//!
//! Placed between LeapFrogVelocity and LeapFrogPosition in the sequence.
//! Uses the "new" kinetic energy from the previous step's TemperatureCalculation
//! (stored in energies.kinetic_energy_new, equivalent to GROMOS multibath.bath.ekin).
//!
//! # Single NHC (nhc = 1)
//! ```text
//!   T_free = 2·E_kin / (dof·k_B)
//!   ζ[0]  += dt/τ² · (T_free/T₀ - 1)
//!   scale  = 1 - ζ[0]·dt
//! ```
//!
//! # NHC chain (nhc = N ≥ 2)  — chain updated from tail to head
//! ```text
//!   τ_vec[0]   = τ²
//!   τ_vec[i>0] = τ²/dof
//!
//!   ζ[N-1] += (τ_vec[N-2]·ζ[N-2]² - 1/dof) / τ_vec[N-1] · dt
//!   for i in (N-2)..1:
//!     ζ[i] += ((τ_vec[i-1]·ζ[i-1]² - 1/dof)/τ_vec[i] - ζ[i]·ζ[i+1]) · dt
//!   ζ[0] += ((T_free/T₀ - 1)/τ_vec[0] - ζ[1]·ζ[0]) · dt
//!   scale = 1 - ζ[0]·dt
//! ```
//!
//! Source: md++/src/algorithm/temperature/nosehoover_thermostat.cc
//!         md++/src/algorithm/temperature/thermostat.cc (scale method)

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::topology::Topology;
use gromos_core::units::kB;

use super::multibath::BathRange;

/// Nosé-Hoover thermostat parameters for a single temperature bath.
#[derive(Debug, Clone)]
pub struct NoseHooverThermostatParams {
    /// Reference temperature T₀ (K)
    pub temperature: f64,
    /// Coupling time constant τ (ps). τ < 0 → no coupling
    pub tau: f64,
    /// Degrees of freedom for this bath
    pub dof: f64,
    /// Chain length: 1 = single NHC, N ≥ 2 = NHC chain
    pub nhc: usize,
}

/// Nosé-Hoover (chain) thermostat.
///
/// GROMOS sequence position: after Leap_Frog_Velocity, before Leap_Frog_Position.
///
/// The ζ chain state is persisted across steps in `zeta`. On construction all
/// ζ values are initialised to zero (GROMOS convention: see `init()` in the
/// C++ source).
#[derive(Debug, Clone)]
pub struct NoseHooverThermostat {
    pub params: Vec<NoseHooverThermostatParams>,
    /// The DOFSET ranges that map atoms to `params` (gromosXX `bath_index`), in file order.
    pub ranges: Vec<BathRange>,
    /// Scaling factor per bath of the current step.
    scale: Vec<f64>,
    /// Per-bath ζ chain state — length = nhc per bath
    zeta: Vec<Vec<f64>>,
}

impl NoseHooverThermostat {
    /// The scaling factor of every bath at the last step (gromosXX `bath_struct::scale`).
    pub fn scales(&self) -> &[f64] {
        &self.scale
    }

    /// Create a single-bath Nosé-Hoover thermostat with chain length 1.
    ///
    /// The bath covers all atoms `[0, n_atoms)`.
    pub fn new_single_bath(temperature: f64, tau: f64, dof: f64, n_atoms: usize) -> Self {
        Self {
            params: vec![NoseHooverThermostatParams {
                temperature,
                tau,
                dof,
                nhc: 1,
            }],
            ranges: whole_system(n_atoms),
            scale: vec![1.0],
            zeta: vec![vec![0.0; 1]],
        }
    }

    /// Create a single-bath Nosé-Hoover chain thermostat with chain length `nhc ≥ 2`.
    ///
    /// The bath covers all atoms `[0, n_atoms)`.
    ///
    /// # Panics
    /// Panics in debug mode if `nhc < 2`. Use [`new_single_bath`] for the single-NHC case.
    ///
    /// [`new_single_bath`]: NoseHooverThermostat::new_single_bath
    pub fn new_chain_bath(
        temperature: f64,
        tau: f64,
        dof: f64,
        n_atoms: usize,
        nhc: usize,
    ) -> Self {
        debug_assert!(
            nhc >= 2,
            "new_chain_bath requires nhc >= 2; use new_single_bath for nhc=1"
        );
        Self {
            params: vec![NoseHooverThermostatParams {
                temperature,
                tau,
                dof,
                nhc,
            }],
            ranges: whole_system(n_atoms),
            scale: vec![1.0],
            zeta: vec![vec![0.0; nhc]],
        }
    }
}

impl Algorithm for NoseHooverThermostat {
    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        let dt = sim.dt;

        self.scale.resize(self.params.len(), 1.0);
        for bath_idx in 0..self.params.len() {
            self.scale[bath_idx] = 1.0;
            // Copy all params as plain scalars (f64/usize are Copy) so that
            // the borrow of self.params is released before we mutate self.zeta.
            let temperature = self.params[bath_idx].temperature;
            let tau = self.params[bath_idx].tau;
            let dof = self.params[bath_idx].dof;
            let nhc = self.params[bath_idx].nhc;

            // τ < 0 → no coupling for this bath (GROMOS: bath.tau == -1 check)
            if tau < 0.0 {
                continue;
            }

            if nhc == 0 {
                return Err(format!("NoseHoover bath {bath_idx}: nhc must be >= 1"));
            }

            // E_kin from previous step (set by TemperatureCalculation): per bath with several
            // baths, the total with one
            let ekin = conf
                .current()
                .energies
                .kinetic_energy_new_bath
                .get(bath_idx)
                .copied()
                .unwrap_or(conf.current().energies.kinetic_energy_new);

            // T_free = 2·E_kin / (dof·k_B)
            let mut free_temp = if dof > 0.0 {
                2.0 * ekin / (dof * kB)
            } else {
                0.0
            };

            // Guard against division-by-zero / zero temperature (GROMOS: epsilon check)
            if free_temp < f64::EPSILON {
                free_temp = temperature;
            }

            let scale = if nhc == 1 {
                // ── Single NHC ──────────────────────────────────────────────
                // ζ[0] += dt/τ² · (T_free/T₀ - 1)
                self.zeta[bath_idx][0] += dt / (tau * tau) * (free_temp / temperature - 1.0);
                1.0 - self.zeta[bath_idx][0] * dt
            } else {
                // ── NHC chain ───────────────────────────────────────────────
                // τ_vec[0] = τ²,  τ_vec[i>0] = τ²/dof
                let tau2 = tau * tau;
                // Inline closure — captures tau2 and dof by value (both Copy).
                let tv = |i: usize| -> f64 {
                    if i == 0 {
                        tau2
                    } else {
                        tau2 / dof
                    }
                };

                // 1. Tail: ζ[N-1] += (τ_vec[N-2]·ζ[N-2]² - 1/dof) / τ_vec[N-1] · dt
                {
                    let z_nm2 = self.zeta[bath_idx][nhc - 2];
                    let delta = (tv(nhc - 2) * z_nm2 * z_nm2 - 1.0 / dof) / tv(nhc - 1) * dt;
                    self.zeta[bath_idx][nhc - 1] += delta;
                }

                // 2. Middle links: i from N-2 down to 1
                //    ζ[i] += ((τ_vec[i-1]·ζ[i-1]² - 1/dof)/τ_vec[i] - ζ[i]·ζ[i+1]) · dt
                //    Note: ζ[i+1] is already updated (tail-to-head order, mirrors GROMOS).
                for i in (1..nhc - 1).rev() {
                    let z_prev = self.zeta[bath_idx][i - 1]; // old ζ[i-1]
                    let z_curr = self.zeta[bath_idx][i]; // old ζ[i]
                    let z_next = self.zeta[bath_idx][i + 1]; // already updated ζ[i+1]
                    let delta =
                        ((tv(i - 1) * z_prev * z_prev - 1.0 / dof) / tv(i) - z_curr * z_next) * dt;
                    self.zeta[bath_idx][i] += delta;
                }

                // 3. Head: ζ[0] += ((T_free/T₀ - 1)/τ_vec[0] - ζ[1]·ζ[0]) · dt
                //    Uses freshly updated ζ[1] (mirrors GROMOS).
                {
                    let z0 = self.zeta[bath_idx][0]; // old ζ[0]
                    let z1 = self.zeta[bath_idx][1]; // already updated ζ[1]
                    let delta = ((free_temp / temperature - 1.0) / tv(0) - z1 * z0) * dt;
                    self.zeta[bath_idx][0] += delta;
                }

                1.0 - self.zeta[bath_idx][0] * dt
            };

            log::debug!(
                "  NoseHoover bath {}: E_kin={:.6e}, T_free={:.2}, T0={:.2}, \
                 scale={:.10}, zeta[0]={:.10}",
                bath_idx,
                ekin,
                free_temp,
                temperature,
                scale,
                self.zeta[bath_idx][0],
            );

            self.scale[bath_idx] = scale;
        }

        // gromosXX `Thermostat::scale`: every DOFSET range is scaled by its bath
        let mut first = 0;
        for r in &self.ranges {
            let scale = self.scale[r.com_bath];
            if (scale - 1.0).abs() >= 1e-15 {
                for v in &mut conf.current_mut().vel[first..=r.last_atom] {
                    *v *= scale;
                }
            }
            first = r.last_atom + 1;
        }

        Ok(())
    }

    fn name(&self) -> &str {
        "NoseHoover_Thermostat"
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}

/// One DOFSET range over all `n_atoms`, coupled to bath 0.
fn whole_system(n_atoms: usize) -> Vec<BathRange> {
    vec![BathRange {
        last_atom: n_atoms - 1,
        com_bath: 0,
        ir_bath: 0,
    }]
}

// ─────────────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::algorithm::SimulationState;
    use gromos_core::configuration::Configuration;
    use gromos_core::math::Vec3;
    use gromos_core::topology::Topology;

    // ── helpers ──────────────────────────────────────────────────────────────

    fn make_sim(dt: f64) -> SimulationState {
        SimulationState::new(dt, 1)
    }

    /// Build a Configuration with `n_atoms` atoms all carrying velocity (1,0,0),
    /// and set `kinetic_energy_new` to the given value.
    fn make_conf(n_atoms: usize, ekin: f64) -> (Topology, Configuration) {
        let topo = Topology::new();
        let mut conf = Configuration::new(n_atoms, 1, 1);
        for i in 0..n_atoms {
            conf.current_mut().vel[i] = Vec3::new(1.0, 0.0, 0.0);
        }
        conf.current_mut().energies.kinetic_energy_new = ekin;
        (topo, conf)
    }

    // ── Single NHC (nhc = 1) ─────────────────────────────────────────────────

    /// Setup chosen so that the arithmetic is exact in f64:
    ///
    ///   T₀   = 300 K, τ = 0.1 ps, dof = 3, dt = 0.002 ps
    ///   ekin = T₀ · 2 · dof · k_B / 2  (free_temp = 2·T₀ = 600 K)
    ///
    ///   Δζ[0] = dt/τ² · (T_free/T₀ - 1) = 0.002/0.01 · 1 = 0.2
    ///   scale = 1 - 0.2·0.002 = 0.9996
    #[test]
    fn test_single_nhc_zeta_and_scale() {
        const T0: f64 = 300.0;
        const TAU: f64 = 0.1;
        const DOF: f64 = 3.0;
        const DT: f64 = 0.002;
        const N: usize = 2;

        // ekin that gives free_temp = 2·T₀ = 600 K
        let ekin = 2.0 * T0 * DOF * kB / 2.0; // = 600·3·k_B/2

        let (topo, mut conf) = make_conf(N, ekin);
        let sim = make_sim(DT);

        let mut nhc = NoseHooverThermostat::new_single_bath(T0, TAU, DOF, N);
        nhc.apply(&topo, &mut conf, &sim).unwrap();

        // Expected: Δζ = 0.002/0.01 · 1.0 = 0.2
        let expected_zeta = 0.2_f64;
        let expected_scale = 1.0 - expected_zeta * DT; // = 0.9996

        assert!(
            (nhc.zeta[0][0] - expected_zeta).abs() < 1e-12,
            "zeta[0] = {}, expected {expected_zeta}",
            nhc.zeta[0][0]
        );

        // All velocities should be scaled by `scale`
        for i in 0..N {
            let vx = conf.current().vel[i].x;
            assert!(
                (vx - expected_scale).abs() < 1e-12,
                "vel[{i}].x = {vx}, expected {expected_scale}"
            );
        }
    }

    /// When the system is cooler than T₀ the thermostat should add energy
    /// (scale > 1, ζ[0] decreases / goes negative).
    #[test]
    fn test_single_nhc_cold_bath_increases_scale() {
        const T0: f64 = 300.0;
        const TAU: f64 = 0.1;
        const DOF: f64 = 3.0;
        const DT: f64 = 0.002;
        const N: usize = 2;

        // ekin gives free_temp = T₀/2 = 150 K (system is cold)
        let ekin = (T0 / 2.0) * DOF * kB / 2.0;

        let (topo, mut conf) = make_conf(N, ekin);
        let sim = make_sim(DT);

        let mut nhc = NoseHooverThermostat::new_single_bath(T0, TAU, DOF, N);
        nhc.apply(&topo, &mut conf, &sim).unwrap();

        // Δζ = dt/τ² · (T_free/T₀ - 1) = 0.008 · (0.5 - 1) = −0.004  → scale > 1
        assert!(
            nhc.zeta[0][0] < 0.0,
            "ζ[0] should be negative for a cold bath, got {}",
            nhc.zeta[0][0]
        );
        let scale = 1.0 - nhc.zeta[0][0] * DT;
        assert!(
            scale > 1.0,
            "scale should be > 1 for a cold bath, got {scale}"
        );
    }

    /// τ < 0 must skip coupling entirely — velocities must not change.
    #[test]
    fn test_single_nhc_no_coupling_when_tau_negative() {
        const N: usize = 3;
        let ekin = 10.0;

        let (topo, mut conf) = make_conf(N, ekin);
        let sim = make_sim(0.002);

        let mut nhc = NoseHooverThermostat::new_single_bath(300.0, -1.0, 3.0, N);
        nhc.apply(&topo, &mut conf, &sim).unwrap();

        // ζ and velocities must remain unchanged
        assert_eq!(nhc.zeta[0][0], 0.0);
        for i in 0..N {
            assert_eq!(conf.current().vel[i].x, 1.0);
        }
    }

    // ── NHC chain (nhc = 3) ───────────────────────────────────────────────────

    /// Chain of length 3, all ζ = 0, free_temp = 2·T₀.
    ///
    /// First step from all-zero ζ with T_free = 2·T₀:
    ///
    ///   τ_vec = [τ², τ²/dof, τ²/dof] = [0.01, 0.01/3, 0.01/3]
    ///
    ///   Tail  (i=2): Δζ = (τ_vec[1]·0² - 1/dof) / τ_vec[2] · dt
    ///                   = (0 - 1/3) / (1/3) · 0.002 / 100 · 300
    ///                   = -1 · 0.2 = -0.2
    ///
    ///   Mid   (i=1): Δζ = ((τ_vec[0]·0² - 1/dof)/τ_vec[1] - 0·(-0.2)) · dt
    ///                   = (-1/3)/(1/300) · 0.002  [same arithmetic] = -0.2
    ///
    ///   Head  (i=0): Δζ = ((T_free/T₀ - 1)/τ_vec[0] - ζ[1]·0) · dt
    ///                   = (1/0.01) · 0.002 = 0.2
    ///
    ///   scale = 1 - 0.2·0.002 = 0.9996
    #[test]
    fn test_chain_nhc_first_step_from_zeros() {
        const T0: f64 = 300.0;
        const TAU: f64 = 0.1;
        const DOF: f64 = 3.0;
        const DT: f64 = 0.002;
        const NHC: usize = 3;
        const N: usize = 2;

        let ekin = 2.0 * T0 * DOF * kB / 2.0; // free_temp = 2·T₀

        let (topo, mut conf) = make_conf(N, ekin);
        let sim = make_sim(DT);

        let mut nhc = NoseHooverThermostat::new_chain_bath(T0, TAU, DOF, N, NHC);
        nhc.apply(&topo, &mut conf, &sim).unwrap();

        let z = &nhc.zeta[0];
        assert!((z[0] - 0.2).abs() < 1e-12, "ζ[0] = {}, expected 0.2", z[0]);
        assert!(
            (z[1] - (-0.2)).abs() < 1e-12,
            "ζ[1] = {}, expected -0.2",
            z[1]
        );
        assert!(
            (z[2] - (-0.2)).abs() < 1e-12,
            "ζ[2] = {}, expected -0.2",
            z[2]
        );

        let expected_scale = 1.0 - 0.2 * DT; // 0.9996
        for i in 0..N {
            let vx = conf.current().vel[i].x;
            assert!(
                (vx - expected_scale).abs() < 1e-12,
                "vel[{i}].x = {vx}, expected {expected_scale}"
            );
        }
    }

    /// Second step — verifies that chain coupling produces different ζ[0] than
    /// the single-NHC would after two steps of identical forcing.
    ///
    /// We apply two steps with the same ekin and check that ζ[1] and ζ[2]
    /// remain non-trivially updated (the chain term −ζ[i]·ζ[i+1] is non-zero
    /// from step 2 onward), confirming propagation through the chain.
    #[test]
    fn test_chain_nhc_propagates_across_two_steps() {
        const T0: f64 = 300.0;
        const TAU: f64 = 0.1;
        const DOF: f64 = 3.0;
        const DT: f64 = 0.002;
        const NHC: usize = 3;
        const N: usize = 2;

        let ekin = 2.0 * T0 * DOF * kB / 2.0; // free_temp = 2·T₀

        let (topo, mut conf) = make_conf(N, ekin);
        // Keep ekin fixed for both steps (we're testing the chain, not the
        // velocity update that would follow in a real MD loop).
        conf.current_mut().energies.kinetic_energy_new = ekin;

        let sim = make_sim(DT);

        let mut nhc_chain = NoseHooverThermostat::new_chain_bath(T0, TAU, DOF, N, NHC);
        nhc_chain.apply(&topo, &mut conf, &sim).unwrap();

        // Reset ekin for second step (velocities were scaled, but we want to
        // test chain arithmetic independently).
        conf.current_mut().energies.kinetic_energy_new = ekin;
        // Reset velocities so the scaling test is clean.
        for i in 0..N {
            conf.current_mut().vel[i] = Vec3::new(1.0, 0.0, 0.0);
        }

        nhc_chain.apply(&topo, &mut conf, &sim).unwrap();

        // After step 1: ζ = [0.2, −0.2, −0.2]
        // Step 2 tail uses ζ[1]² > 0, so the tail term is not just −1/dof
        // anymore — ζ[2] must differ from −0.4 (what two identical tail steps
        // with zero coupling would give).
        //
        // Exact value from the recurrence:
        //   Δζ₂ = (τ_vec[1]·(−0.2)² − 1/3) / τ_vec[2] · dt
        //        = ((0.01/3)·0.04 − 1/3) / (0.01/3) · 0.002
        //        = (0.04 − 100) · 0.002        [after dividing by 0.01/3]
        //
        // (0.01/3)·0.04 = 0.0004/3; divided by 0.01/3 = *300 → 0.04
        // −1/3 divided by 0.01/3 = −100
        // total: (0.04 − 100)·0.002 = −99.96·0.002 = −0.19992
        // ζ[2] = −0.2 + (−0.19992) = −0.39992

        let z = &nhc_chain.zeta[0];

        assert!(
            (z[2] - (-0.39992)).abs() < 1e-9,
            "ζ[2] after step 2 = {}, expected ≈ -0.39992",
            z[2]
        );

        // The chain cross-term (−ζ[i]·ζ[i+1]) is non-zero in the middle update,
        // so ζ[1] after step 2 must differ from what two identical single-NHC
        // steps (no chain coupling) would produce.
        // For a sanity check: ζ[1] must be negative and more negative than −0.2.
        assert!(
            z[1] < -0.2,
            "ζ[1] after step 2 should be < −0.2 (chain coupling active), got {}",
            z[1]
        );

        // ζ[0] must be positive and larger than after step 1 (system still hot).
        assert!(
            z[0] > 0.2,
            "ζ[0] after step 2 should be > 0.2 (continued driving), got {}",
            z[0]
        );
    }
}
