//! Bonded interaction calculations.
//!
//! Submodules by term type; this module re-exports everything and provides
//! the shared types and top-level combiners.

pub mod angles;
pub mod bonds;
pub mod dihedrals;
pub mod improper;
pub mod perturbed;

pub use angles::*;
pub use bonds::*;
pub use dihedrals::*;
pub use improper::*;
pub use perturbed::{
    calculate_perturbed_angle_forces, calculate_perturbed_bond_forces,
    calculate_perturbed_bonded_forces, calculate_perturbed_dihedral_forces,
    calculate_perturbed_improper_dihedral_forces, BondedPrecalc, TermLambda,
};

use gromos_core::configuration::Configuration;
use gromos_core::math::Periodicity;
use gromos_core::math::Vec3;
use gromos_core::topology::Topology;

// ─── Shared result types ─────────────────────────────────────────────────────

/// Energy + per-atom forces + virial tensor from a bonded term.
///
/// `energy` is always the combined total. The `*_energy` fields are the
/// per-term breakdown (bond/angle/dihedral/improper); they're only populated
/// by combiners like [`calculate_bonded_forces_ntf`] that know which term is
/// which — a single-term calculator (e.g. `calculate_bond_forces_quartic`)
/// leaves them at 0 and relies on `energy` alone.
#[derive(Debug, Clone)]
pub struct ForceEnergy {
    pub energy: f64,
    pub forces: Vec<Vec3>,
    /// virial[a][b] = Σ r[b] * f[a]  (GROMOS convention)
    pub virial: [[f64; 3]; 3],
    pub bond_energy: f64,
    pub angle_energy: f64,
    pub dihedral_energy: f64,
    pub improper_energy: f64,
}

/// Energy + per-atom forces + dH/dλ from a perturbed (FEP/TI) bonded term.
#[derive(Debug, Clone)]
pub struct ForceEnergyLambda {
    pub energy: f64,
    pub forces: Vec<Vec3>,
    pub lambda_derivative: f64,
    /// `energy` split by term kind, so the perturbed contributions land in the same `.tre`
    /// columns as their unperturbed counterparts (bond/angle/improper/dihedral) — gromosXX
    /// books them there; soft bonds/angles/impropers count as their kind.
    pub bond_energy: f64,
    pub angle_energy: f64,
    pub improper_energy: f64,
    pub dihedral_energy: f64,
    /// `lambda_derivative` split the same way (gromosXX `perturbed_energy_derivatives`
    /// `bond_energy`, … — the per-term columns of FREEENERDERIVS03).
    pub bond_dhdl: f64,
    pub angle_dhdl: f64,
    pub improper_dhdl: f64,
    pub dihedral_dhdl: f64,
}

/// Bonded energies per energy group, each term booked to the group of its first atom
/// (gromosXX `bond_energy[atom_energy_group(i)]`, `angle_energy`, `improper_energy`,
/// `dihedral_energy`). Empty with one energy group, where nothing is booked per group.
#[derive(Debug, Clone, Default)]
pub struct BondedGroupEnergies {
    pub bond: Vec<f64>,
    pub angle: Vec<f64>,
    pub improper: Vec<f64>,
    pub dihedral: Vec<f64>,
}

impl BondedGroupEnergies {
    /// Zeroed arrays for `n_groups` groups; empty for one group.
    pub fn new(n_groups: usize) -> Self {
        let n = if n_groups > 1 { n_groups } else { 0 };
        Self {
            bond: vec![0.0; n],
            angle: vec![0.0; n],
            improper: vec![0.0; n],
            dihedral: vec![0.0; n],
        }
    }

    /// Whether the arrays are sized for more than one group.
    #[inline]
    pub fn is_grouped(&self) -> bool {
        !self.bond.is_empty()
    }

    pub fn clear(&mut self) {
        self.bond.fill(0.0);
        self.angle.fill(0.0);
        self.improper.fill(0.0);
        self.dihedral.fill(0.0);
    }
}

impl ForceEnergy {
    pub fn new(num_atoms: usize) -> Self {
        Self {
            energy: 0.0,
            forces: vec![Vec3::ZERO; num_atoms],
            virial: [[0.0; 3]; 3],
            bond_energy: 0.0,
            angle_energy: 0.0,
            dihedral_energy: 0.0,
            improper_energy: 0.0,
        }
    }

    pub fn add(&mut self, other: &ForceEnergy) {
        self.energy += other.energy;
        for i in 0..self.forces.len().min(other.forces.len()) {
            self.forces[i] += other.forces[i];
        }
        for a in 0..3 {
            for b in 0..3 {
                self.virial[a][b] += other.virial[a][b];
            }
        }
        self.bond_energy += other.bond_energy;
        self.angle_energy += other.angle_energy;
        self.dihedral_energy += other.dihedral_energy;
        self.improper_energy += other.improper_energy;
    }
}

impl ForceEnergyLambda {
    pub fn new(num_atoms: usize) -> Self {
        Self {
            energy: 0.0,
            forces: vec![Vec3::ZERO; num_atoms],
            lambda_derivative: 0.0,
            bond_energy: 0.0,
            angle_energy: 0.0,
            improper_energy: 0.0,
            dihedral_energy: 0.0,
            bond_dhdl: 0.0,
            angle_dhdl: 0.0,
            improper_dhdl: 0.0,
            dihedral_dhdl: 0.0,
        }
    }

    pub fn add(&mut self, other: &ForceEnergyLambda) {
        self.energy += other.energy;
        self.lambda_derivative += other.lambda_derivative;
        self.bond_energy += other.bond_energy;
        self.angle_energy += other.angle_energy;
        self.improper_energy += other.improper_energy;
        self.dihedral_energy += other.dihedral_energy;
        self.bond_dhdl += other.bond_dhdl;
        self.angle_dhdl += other.angle_dhdl;
        self.improper_dhdl += other.improper_dhdl;
        self.dihedral_dhdl += other.dihedral_dhdl;
        for i in 0..self.forces.len().min(other.forces.len()) {
            self.forces[i] += other.forces[i];
        }
    }
}

// ─── Top-level combiners ─────────────────────────────────────────────────────

/// Calculate all bonded forces.
pub fn calculate_bonded_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    use_quartic_bonds: bool,
) -> ForceEnergy {
    calculate_bonded_forces_ntf(
        topo,
        conf,
        pbc,
        use_quartic_bonds,
        CovalentForm::default(),
        true,
        true,
        true,
        true,
        &mut BondedGroupEnergies::default(),
    )
}

/// Calculate bonded forces gated by NTF flags (FORCE block in GROMOS).
/// Which functional form each covalent term takes — gromosXX's `COVALENTFORM` block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CovalentForm {
    /// NTBBH = 0: quartic bonds (the GROMOS default); 1: harmonic
    pub quartic_bonds: bool,
    /// NTBAH = 0: cosine-harmonic angles (the GROMOS default); 1: harmonic
    pub cosine_harmonic_angles: bool,
    /// NTBDN = 0: arbitrary dihedral phase shifts (the GROMOS default); 1: 0° and 180° only
    pub arbitrary_phase_shifts: bool,
}

impl Default for CovalentForm {
    fn default() -> Self {
        Self {
            quartic_bonds: true,
            cosine_harmonic_angles: true,
            arbitrary_phase_shifts: true,
        }
    }
}

#[allow(clippy::too_many_arguments)] // the NTF flags of the FORCE block, one per term
pub fn calculate_bonded_forces_ntf(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    use_quartic_bonds: bool,
    form: CovalentForm,
    ntf_bond: bool,
    ntf_angle: bool,
    ntf_dihedral: bool,
    ntf_improper: bool,
    groups: &mut BondedGroupEnergies,
) -> ForceEnergy {
    if groups.is_grouped() {
        bonded_forces_ntf::<true>(
            topo,
            conf,
            pbc,
            use_quartic_bonds,
            form,
            [ntf_bond, ntf_angle, ntf_dihedral, ntf_improper],
            groups,
        )
    } else {
        bonded_forces_ntf::<false>(
            topo,
            conf,
            pbc,
            use_quartic_bonds,
            form,
            [ntf_bond, ntf_angle, ntf_dihedral, ntf_improper],
            groups,
        )
    }
}

/// [`calculate_bonded_forces_ntf`]; with `GROUPS` every term is also booked into `groups`.
fn bonded_forces_ntf<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    use_quartic_bonds: bool,
    form: CovalentForm,
    [ntf_bond, ntf_angle, ntf_dihedral, ntf_improper]: [bool; 4],
    groups: &mut BondedGroupEnergies,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    if ntf_bond {
        let mut bf = if use_quartic_bonds && form.quartic_bonds {
            calculate_bond_forces_quartic_grouped::<GROUPS>(topo, conf, pbc, &mut groups.bond)
        } else {
            calculate_bond_forces_harmonic_grouped::<GROUPS>(topo, conf, pbc, &mut groups.bond)
        };
        log::debug!("  bonded: bond={:.6e}", bf.energy);
        bf.bond_energy = bf.energy;
        result.add(&bf);
    }
    if ntf_angle {
        // COVALENTFORM NTBAH: cosine-harmonic (the GROMOS default) or harmonic
        let mut af = if form.cosine_harmonic_angles {
            calculate_angle_forces_grouped::<GROUPS>(topo, conf, pbc, &mut groups.angle)
        } else {
            calculate_harmonic_angle_forces_grouped::<GROUPS>(topo, conf, pbc, &mut groups.angle)
        };
        log::debug!("  bonded: angle={:.6e}", af.energy);
        af.angle_energy = af.energy;
        result.add(&af);
    }
    if ntf_dihedral {
        // COVALENTFORM NTBDN: arbitrary phase shifts (the GROMOS default) or 0°/180° only
        let mut df = if form.arbitrary_phase_shifts {
            calculate_dihedral_new_forces_grouped::<GROUPS>(topo, conf, pbc, &mut groups.dihedral)
        } else {
            calculate_dihedral_forces_grouped::<GROUPS>(topo, conf, pbc, &mut groups.dihedral)
        };
        log::debug!("  bonded: dihe={:.6e}", df.energy);
        df.dihedral_energy = df.energy;
        result.add(&df);
    }
    if ntf_improper {
        let mut imf = calculate_improper_dihedral_forces_grouped::<GROUPS>(
            topo,
            conf,
            pbc,
            &mut groups.improper,
        );
        log::debug!("  bonded: impr={:.6e}", imf.energy);
        imf.improper_energy = imf.energy;
        result.add(&imf);
    }

    log::debug!(
        "  bonded total={:.6e}  max|f|={:.6e}",
        result.energy,
        result
            .forces
            .iter()
            .map(|f| f.length())
            .fold(0.0_f64, f64::max)
    );
    result
}

#[cfg(test)]
mod tests {
    /// The unit tests below are single isolated molecules: no box, so no minimum image.
    fn vacuum() -> gromos_core::math::Periodicity {
        gromos_core::math::Periodicity::Vacuum(gromos_core::math::Vacuum)
    }

    use super::*;
    use gromos_core::configuration::Configuration;
    use gromos_core::topology::Topology;

    #[test]
    fn test_bonded_forces_complete() {
        let topo = Topology::new();
        let conf = Configuration::new(0, 1, 1);
        let result = calculate_bonded_forces(&topo, &conf, &vacuum(), true);
        assert_eq!(result.energy, 0.0);
        assert!(result.forces.is_empty());
    }
}
