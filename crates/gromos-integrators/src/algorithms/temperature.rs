//! Temperature / kinetic energy calculation algorithm.
//!
//! Equivalent to GROMOS `algorithm::Temperature_Calculation`.
//! Computes kinetic energy using the GROMOS convention:
//! E_kin = 0.5 * sum_i( m_i * (|v_new_i|^2 + |v_old_i|^2) / 2 )
//! where v_new = conf.current().vel, v_old = conf.old().vel
//!
//! **Init behavior**: `init()` calls `apply()` once before the MD loop to
//! pre-compute `kinetic_energy_new` from the initial velocities. This is
//! required so that the Berendsen thermostat has a valid E_kin for its very
//! first scaling step. This matches GROMOS `Temperature_Calculation::init()`
//! (see md++/src/algorithm/temperature/temperature_calculation.cc line ~248-260).
//!
//! Note: this is NOT a GROMOS bug — it is physically correct. The thermostat
//! needs to know the current kinetic energy to compute the first scaling factor.
//! Without this init, E_kin_new=0 at step 0, the thermostat falls back to
//! T_free=T0 and scale=1.0, producing an unscaled first step that cascades
//! into diverging trajectories.
//!
//! Source: md++/src/configuration/state_properties.cc `molecular_translational_ekin()`

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::lambdas::InteractionLambda;
use gromos_core::math::Vec3;
use gromos_core::topology::Topology;

use super::multibath::{in_bath, BathRange};

/// Calculates kinetic energy and temperature using GROMOS's averaged formula.
///
/// In GROMOS, E_kin = 0.5 * sum_i( m_i * (|v_new|^2 + |v_old|^2) / 2 )
/// This averages the kinetic energy between the current and old velocities.
/// The result is stored in conf.old().energies (GROMOS convention).
///
/// Per bath, the kinetic energy is split as gromosXX splits it
/// (`temperature_calculation.cc:60-113`, `state_properties.cc:153-186`): every temperature
/// group contributes the translational kinetic energy of its centre of mass (`com`) to the COM
/// bath and the rest (`ir`, internal and rotational) to the IR bath of its first atom
/// ([`in_bath`]). On a perturbed run the mass derivative of the kinetic energy goes to
/// `perturbed_energy_derivatives`, booked to the COM bath of each atom's range
/// (`temperature_calculation.cc:116-237`).
#[derive(Debug, Clone, Default)]
pub struct TemperatureCalculation {
    /// Number of temperature baths.
    n_baths: usize,
    /// The DOFSET ranges, in file order. Empty = one bath over the whole system, where only the
    /// totals feed the thermostat.
    ranges: Vec<BathRange>,
    /// Temperature groups as (first, last) atoms, built from the topology on first use.
    groups: Vec<(usize, usize)>,
    /// Perturbed solute atoms and their `(m_A, m_B)`, in atom order.
    perturbed_mass: Vec<(usize, f64, f64)>,
    /// Number of atoms `groups` and `perturbed_mass` were built for.
    prepared_for: Option<usize>,
    com: Vec<f64>,
    ir: Vec<f64>,
    new_bath: Vec<f64>,
    de_kin: Vec<f64>,
    /// `(|v_old|² + |v_new|²)/2` of each perturbed atom, this step.
    avg_v2: Vec<f64>,
}

impl TemperatureCalculation {
    pub fn new() -> Self {
        Self::default()
    }

    /// The kinetic energy is accumulated per bath as well as in total, as gromosXX does
    /// (`Multibath::calculate_kinetic_energy`): each bath's thermostat scales on its own energy.
    /// `ranges` are the DOFSET lines that map atoms to the `n_baths` baths.
    pub fn with_baths(n_baths: usize, ranges: &[BathRange]) -> Self {
        Self {
            n_baths,
            ranges: ranges.to_vec(),
            ..Self::default()
        }
    }

    fn prepare(&mut self, topo: &Topology, n_atoms: usize) {
        if self.prepared_for == Some(n_atoms) {
            return;
        }
        self.groups = temperature_group_ranges(topo, n_atoms);
        let n_solute = topo.num_solute_atoms().min(n_atoms);
        let mut pert: Vec<(usize, f64, f64)> = topo
            .perturbed_solute
            .atoms
            .iter()
            .filter(|a| a.seq < n_solute)
            .map(|a| (a.seq, a.a_mass, a.b_mass))
            .collect();
        pert.sort_by_key(|&(i, _, _)| i);
        self.perturbed_mass = pert;
        self.prepared_for = Some(n_atoms);
    }
}

/// The temperature groups of the topology as (first, last) atoms. Atoms the topology does not
/// group follow gromosXX's defaults (`topology.cc:591-594`, `solvate`): the rest of the solute
/// is one group, every solvent molecule another.
fn temperature_group_ranges(topo: &Topology, n_atoms: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut next = 0usize;
    for g in &topo.temperature_groups {
        if let (Some(&first), Some(&last)) = (g.first(), g.last()) {
            if first == next && last < n_atoms {
                out.push((first, last));
                next = last + 1;
            }
        }
    }
    let n_solute = topo.num_solute_atoms().min(n_atoms);
    if next < n_solute {
        out.push((next, n_solute - 1));
        next = n_solute;
    }
    let per_molecule = topo.solvent_atom_template.len();
    while per_molecule > 0 && next + per_molecule <= n_atoms {
        out.push((next, next + per_molecule - 1));
        next += per_molecule;
    }
    if next < n_atoms {
        out.push((next, n_atoms - 1));
    }
    out
}

impl Algorithm for TemperatureCalculation {
    /// Initialize by running apply() once to pre-compute kinetic_energy_new.
    ///
    /// GROMOS: Temperature_Calculation::init() calls apply() so that
    /// multibath.bath.ekin is populated before the first MD step. Without
    /// this, the Berendsen thermostat would see E_kin=0 at step 0 and
    /// skip scaling entirely.
    fn init(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        log::debug!("TemperatureCalculation::init() — pre-computing E_kin from initial velocities");
        self.apply(topo, conf, sim)
    }

    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        // GROMOS convention: E_kin = 0.5 * sum_i( m_i * (|v_new|^2 + |v_old|^2) / 2 )
        // After Leap_Frog_Velocity's exchange_state():
        //   current().vel = new velocities v(t+dt/2)
        //   old().vel     = previous velocities v(t-dt/2)
        // The average gives the kinetic energy at time t.
        let n_atoms = topo.inverse_mass.len();
        let mut e_kin = 0.0;
        let mut e_kin_new = 0.0;
        for i in 0..n_atoms {
            let v_new = conf.current().vel[i];
            let v_old = conf.old().vel[i];
            let m = topo.mass[i];
            e_kin += 0.5 * m * (v_new.length_squared() + v_old.length_squared()) / 2.0;
            e_kin_new += 0.5 * m * v_new.length_squared();
        }

        self.prepare(topo, n_atoms);
        let n_baths = if self.ranges.is_empty() {
            1
        } else {
            self.n_baths.max(1)
        };
        for v in [&mut self.com, &mut self.ir, &mut self.new_bath] {
            v.clear();
            v.resize(n_baths, 0.0);
        }
        {
            let vel = &conf.current().vel;
            let old_vel = &conf.old().vel;
            for &(first, last) in &self.groups {
                let mut tot_mass = 0.0;
                let mut com_v = Vec3::ZERO;
                let mut ekin = 0.0;
                let mut new_com_v = Vec3::ZERO;
                let mut new_ekin = 0.0;
                for a in first..=last {
                    let m = topo.mass[a];
                    tot_mass += m;
                    let (new_v, old_v) = (vel[a], old_vel[a]);
                    let v = 0.5 * (new_v + old_v);
                    com_v += m * v;
                    ekin += m * (new_v.length_squared() + old_v.length_squared()) / 2.0;
                    new_com_v += m * new_v;
                    new_ekin += m * new_v.length_squared();
                }
                com_v /= tot_mass;
                new_com_v /= tot_mass;
                let com_ekin = 0.5 * tot_mass * com_v.length_squared();
                ekin *= 0.5;
                let new_com_ekin = 0.5 * tot_mass * new_com_v.length_squared();
                new_ekin *= 0.5;

                let (com_bath, ir_bath) = in_bath(&self.ranges, first);
                self.new_bath[com_bath] += new_com_ekin;
                self.new_bath[ir_bath] += new_ekin - new_com_ekin;
                self.com[com_bath] += com_ekin;
                self.ir[ir_bath] += ekin - com_ekin;
            }
        }

        // dE_kin/dλ over the perturbed solute atoms, with the mass L′ of each atom's energy
        // group (gromosXX `temperature_calculation.cc:163-178`)
        let perturbed = !self.perturbed_mass.is_empty();
        if perturbed {
            self.de_kin.clear();
            self.de_kin.resize(n_baths, 0.0);
            let vel = &conf.current().vel;
            let old_vel = &conf.old().vel;
            // PRECALCLAM: the kinetic energy of the perturbed atoms with their masses at each
            // λ_k, over the current mass λ (temperature_calculation.cc:200-224); zeroed every
            // step (gromosXX `zero(false, true)` and the force field's zeroing)
            let precalc_on = topo.precalc_lambdas.is_some_and(|p| p.active(sim.step));
            self.avg_v2.clear();
            for &(i, _, _) in &self.perturbed_mass {
                self.avg_v2
                    .push((old_vel[i].length_squared() + vel[i].length_squared()) * 0.5);
            }
            let old = conf.old_mut();
            old.energies.precalc.zero_kinetic();
            old.perturbed_energy_derivatives.precalc.zero_kinetic();
            let n_lambdas = old.energies.precalc.nr_lambdas();
            for (&(i, ma, mb), &avg_v2) in self.perturbed_mass.iter().zip(&self.avg_v2) {
                let g = topo.atom_to_energy_group.get(i).copied().unwrap_or(0);
                let l_deriv = topo
                    .individual_lambdas
                    .derivative(InteractionLambda::Mass, g, g);
                let (bath, _) = in_bath(&self.ranges, i);
                self.de_kin[bath] -= l_deriv * (mb - ma) * avg_v2;
                if precalc_on {
                    let grid = topo.precalc_lambdas.expect("precalc_on");
                    let slam = topo
                        .individual_lambdas
                        .lambda(InteractionLambda::Mass, g, g);
                    let mass_s = (1.0 - slam) * ma + slam * mb;
                    let mass_s2 = mass_s * mass_s;
                    for k in 0..n_lambdas {
                        let lam = grid.lambda(k);
                        let mass_p = (1.0 - lam) * ma + lam * mb;
                        let mass_p2 = mass_p * mass_p;
                        old.energies.precalc.ab_kinetic[k] += 0.5 * (mass_s2 / mass_p) * avg_v2;
                        old.perturbed_energy_derivatives.precalc.ab_kinetic[k] +=
                            -0.5 * (mb - ma) * (mass_s2 / mass_p2) * avg_v2;
                    }
                }
            }
            for e in self.de_kin.iter_mut() {
                *e *= 0.5;
            }
        }

        let old = conf.old_mut();
        let e = &mut old.energies;
        e.kinetic_total = e_kin;
        // Store "new" E_kin for thermostat scaling (GROMOS: multibath.bath.ekin)
        e.kinetic_energy_new = e_kin_new;
        for v in [
            &mut e.kinetic_energy,
            &mut e.com_kinetic_energy,
            &mut e.ir_kinetic_energy,
        ] {
            v.resize(n_baths, 0.0);
        }
        for b in 0..n_baths {
            e.com_kinetic_energy[b] = self.com[b];
            e.ir_kinetic_energy[b] = self.ir[b];
            e.kinetic_energy[b] = self.com[b] + self.ir[b];
        }
        if !self.ranges.is_empty() {
            e.kinetic_energy_new_bath.clear();
            e.kinetic_energy_new_bath.extend_from_slice(&self.new_bath);
        }
        if perturbed {
            let d = &mut old.perturbed_energy_derivatives;
            d.kinetic_energy.resize(n_baths, 0.0);
            d.kinetic_energy.copy_from_slice(&self.de_kin);
            d.kinetic_total = 0.0;
            for &k in &self.de_kin {
                d.kinetic_total += k;
            }
        }
        log::debug!("  E_kin={:.10e}  E_kin_new={:.10e}", e_kin, e_kin_new);
        Ok(())
    }

    fn name(&self) -> &str {
        "Temperature_Calculation"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::math::Vec3;

    #[test]
    fn kinetic_energy_is_the_average_of_the_old_and_new_half_steps() {
        // gromosXX: E_kin = ½ m (v_old² + v_new²) / 2 — the leap-frog average
        let mut topo = Topology::new();
        topo.mass = vec![2.0];
        topo.inverse_mass = vec![0.5];
        let mut conf = Configuration::new(1, 1, 1);
        conf.old_mut().vel[0] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().vel[0] = Vec3::new(3.0, 0.0, 0.0);
        let sim = SimulationState::new(0.002, 1);
        TemperatureCalculation::new()
            .apply(&topo, &mut conf, &sim)
            .unwrap();
        // ½·2·(1 + 9)/2 = 5 ; new-only: ½·2·9 = 9
        assert_eq!(conf.old().energies.kinetic_total, 5.0);
        assert_eq!(conf.old().energies.kinetic_energy_new, 9.0);
    }

    #[test]
    fn each_temperature_group_splits_into_com_and_internal_energy() {
        // two atoms moving apart in one group: no centre-of-mass motion, all internal
        let mut topo = Topology::new();
        topo.mass = vec![1.0, 1.0];
        topo.inverse_mass = vec![1.0, 1.0];
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().vel = vec![Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0)];
        conf.copy_current_to_old();
        let sim = SimulationState::new(0.002, 1);
        TemperatureCalculation::new()
            .apply(&topo, &mut conf, &sim)
            .unwrap();
        let e = &conf.old().energies;
        assert_eq!(e.com_kinetic_energy, vec![0.0]);
        assert_eq!(e.ir_kinetic_energy, vec![1.0]);
        assert_eq!(e.kinetic_energy, vec![1.0]);

        // one group per atom, one bath per group: all centre-of-mass motion
        topo.temperature_groups = vec![vec![0], vec![1]];
        let ranges = [
            BathRange {
                last_atom: 0,
                com_bath: 0,
                ir_bath: 0,
            },
            BathRange {
                last_atom: 1,
                com_bath: 1,
                ir_bath: 1,
            },
        ];
        let mut tc = TemperatureCalculation::with_baths(2, &ranges);
        tc.apply(&topo, &mut conf, &sim).unwrap();
        let e = &conf.old().energies;
        assert_eq!(e.com_kinetic_energy, vec![0.5, 0.5]);
        assert_eq!(e.ir_kinetic_energy, vec![0.0, 0.0]);
        assert_eq!(e.kinetic_energy_new_bath, vec![0.5, 0.5]);
        assert_eq!(e.kinetic_total, 1.0);
    }

    #[test]
    fn a_perturbed_mass_gives_a_kinetic_lambda_derivative() {
        // gromosXX: dE_kin/dλ = −½ (m_B − m_A) (|v_old|² + |v_new|²)/2
        let mut topo = Topology::new();
        topo.mass = vec![2.0];
        topo.inverse_mass = vec![0.5];
        topo.moltypes[0]
            .atoms
            .push(gromos_core::topology::MolTypeAtom {
                name: "A".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: 2.0,
                charge: 0.0,
                is_perturbed: true,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        topo.perturbed_solute
            .atoms
            .push(gromos_core::topology::PerturbedAtom {
                seq: 0,
                a_iac: 0,
                a_mass: 1.0,
                a_charge: 0.0,
                b_iac: 0,
                b_mass: 3.0,
                b_charge: 0.0,
                lj_soft: 0.0,
                crf_soft: 0.0,
            });
        let mut conf = Configuration::new(1, 1, 1);
        conf.old_mut().vel[0] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().vel[0] = Vec3::new(3.0, 0.0, 0.0);
        let sim = SimulationState::new(0.002, 1);
        TemperatureCalculation::new()
            .apply(&topo, &mut conf, &sim)
            .unwrap();
        // avg v² = 5, (m_B − m_A) = 2 → −0.5 · 2 · 5
        let d = &conf.old().perturbed_energy_derivatives;
        assert_eq!(d.kinetic_energy, vec![-5.0]);
        assert_eq!(d.kinetic_total, -5.0);
    }

    #[test]
    fn init_computes_the_same_thing_as_apply() {
        let mut topo = Topology::new();
        topo.mass = vec![1.0, 1.0];
        topo.inverse_mass = vec![1.0, 1.0];
        let mut conf = Configuration::new(2, 1, 1);
        conf.old_mut().vel[1] = Vec3::new(0.0, 2.0, 0.0);
        conf.current_mut().vel[1] = Vec3::new(0.0, 2.0, 0.0);
        let sim = SimulationState::new(0.002, 1);
        TemperatureCalculation::new()
            .init(&topo, &mut conf, &sim)
            .unwrap();
        assert_eq!(conf.old().energies.kinetic_total, 2.0);
    }
}
