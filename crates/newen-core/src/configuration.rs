//! Configuration module - system state and energy management
//!
//! Direct Rust translation of GROMOS configuration structures from:
//! - md++/src/configuration/configuration.h
//! - md++/src/configuration/energy.h

use crate::math::{Mat3, Vec3};

/// Simulation box representation
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoxType {
    /// No periodic boundary (gas-phase / vacuum simulations).
    Vacuum,
    /// Orthorhombic (rectangular) periodic box.
    Rectangular,
    /// General triclinic periodic box.
    Triclinic,
    /// Truncated-octahedral box stored as a triclinic cell.
    TruncatedOctahedral,
}

/// Simulation box with periodic boundary conditions
#[derive(Debug, Clone)]
pub struct Box {
    /// Periodicity type (vacuum, rectangular, triclinic, truncated-octahedral).
    pub box_type: BoxType,
    /// Box vectors as columns of a 3×3 matrix (nm).
    pub vectors: Mat3,
    /// Inverse of `vectors`; used for wrapping positions into the primary cell.
    pub inv_vectors: Mat3,
}

impl Box {
    /// Create vacuum box (no periodicity)
    pub fn vacuum() -> Self {
        Self {
            box_type: BoxType::Vacuum,
            vectors: Mat3::IDENTITY,
            inv_vectors: Mat3::IDENTITY,
        }
    }

    /// Create rectangular box with dimensions (lx, ly, lz)
    pub fn rectangular(lx: f64, ly: f64, lz: f64) -> Self {
        let vectors = Mat3::from_cols(
            Vec3::new(lx, 0.0, 0.0),
            Vec3::new(0.0, ly, 0.0),
            Vec3::new(0.0, 0.0, lz),
        );

        let inv_vectors = Mat3::from_cols(
            Vec3::new(1.0 / lx, 0.0, 0.0),
            Vec3::new(0.0, 1.0 / ly, 0.0),
            Vec3::new(0.0, 0.0, 1.0 / lz),
        );

        Self {
            box_type: BoxType::Rectangular,
            vectors,
            inv_vectors,
        }
    }

    /// Create triclinic box from box vectors
    pub fn triclinic(vectors: Mat3) -> Self {
        Self {
            box_type: BoxType::Triclinic,
            vectors,
            inv_vectors: vectors.inverse(),
        }
    }

    /// Create a truncated-octahedron box from its (lower-triangular)
    /// triclinic box vectors, as produced by
    /// `math::truncoct_triclinic_box(cubic, true)` from the legacy NTB=-1
    /// "cube edge length" BOX block (P1.4).
    pub fn truncated_octahedral(vectors: Mat3) -> Self {
        Self {
            box_type: BoxType::TruncatedOctahedral,
            vectors,
            inv_vectors: vectors.inverse(),
        }
    }

    /// Get box volume
    pub fn volume(&self) -> f64 {
        // Volume = determinant of box matrix
        self.vectors.determinant()
    }

    /// Get box dimensions (for rectangular box)
    pub fn dimensions(&self) -> Vec3 {
        Vec3::new(
            self.vectors.x_axis.x,
            self.vectors.y_axis.y,
            self.vectors.z_axis.z,
        )
    }
}

/// Energy storage for the system, laid out like gromosXX's `configuration::Energy`
/// (`configuration/energy.h`): totals, per-energy-group arrays, energy-group-pair matrices and
/// per-bath kinetic energies. A `State` carries two of them: the energies and their λ
/// derivatives (`State::perturbed_energy_derivatives`).
///
/// The `*_total` scalars are the engine's totals; the per-group arrays and the matrices are what
/// the `# bonded`, `# nonbonded` and `# special` sections of `ENERGY03` / `FREEENERDERIVS03`
/// print. Neither object is cleared each step: whoever fills a field assigns it every step.
#[derive(Debug, Clone)]
pub struct Energy {
    /// Total kinetic energy (kJ mol⁻¹).
    pub kinetic_total: f64,
    /// Total potential energy (kJ mol⁻¹).
    pub potential_total: f64,

    /// Sum of all covalent bond stretch energies.
    pub bond_total: f64,
    /// Sum of all bond-angle bending energies.
    pub angle_total: f64,
    /// Sum of all proper dihedral torsion energies.
    pub dihedral_total: f64,
    /// Sum of all improper dihedral (planarity / chirality) energies.
    pub improper_total: f64,
    /// Sum of all cross-dihedral coupling energies (GROMOS CMAP-equivalent).
    pub cross_dihedral_total: f64,

    /// Total Lennard-Jones (van der Waals) nonbonded energy.
    pub lj_total: f64,
    /// Total Coulomb reaction-field electrostatic energy.
    pub crf_total: f64,
    /// Total lattice-sum (Ewald / P3M) electrostatic energy: `ls_pair_total + ls_self_total +
    /// ls_surface_total` (gromosXX `energy.cc:584`).
    pub ls_total: f64,
    /// The pair part of the lattice-sum energy: real space + k space + A term (`energy.cc:582`).
    pub ls_pair_total: f64,
    /// Real-space contribution to the lattice-sum energy.
    pub ls_realspace_total: f64,
    /// k-space contribution to the lattice-sum energy.
    pub ls_kspace_total: f64,
    /// The A term of the lattice sum (GROMOS Vol. 2 `e_a_def`).
    pub ls_a_term_total: f64,
    /// The self term of the lattice sum (Vol. 2 `e_slf`).
    pub ls_self_total: f64,
    /// The surface term of the lattice sum (Vol. 2 `e_srf`).
    pub ls_surface_total: f64,

    /// Sum of special interaction energies (position restraints, etc.).
    pub special_total: f64,
    /// Energy of each additive term of the run (`Term("xtb", …)`, `Term("schnet", …)`), keyed by
    /// the term's registry name (`xtb`, or `xtb:1` when a kind occurs more than once), in plan
    /// order. Their sum is booked in `special_total`; this keeps them visible one by one
    /// (3.9 G10) so a compensating error between two terms cannot hide.
    pub term_energies: Vec<(String, f64)>,
    /// Sum of position-restraint energies (also part of `special_total`).
    pub posrest_total: f64,
    /// Sum of distance-restraint energies.
    pub distanceres_total: f64,
    /// Solvent-accessible surface-area (SASA) solvation energy.
    pub sasa_total: f64,
    /// Energy contribution from holonomic constraints (SHAKE/LINCS).
    pub constraint_total: f64,
    /// gromosXX `entropy_term`: on the derivative object, dH/dλ · H of the potential
    /// (`energy_calculation.cc:87-89`); zero on the energies.
    pub entropy_term: f64,

    /// Kinetic energy from the current (new) velocities, used for thermostat scaling.
    pub kinetic_energy_new: f64,

    /// Kinetic energy per temperature bath: `com_kinetic_energy + ir_kinetic_energy`.
    pub kinetic_energy: Vec<f64>,
    /// Centre-of-mass translational kinetic energy of the temperature groups, per bath.
    pub com_kinetic_energy: Vec<f64>,
    /// Internal and rotational kinetic energy of the temperature groups, per bath.
    pub ir_kinetic_energy: Vec<f64>,

    /// Per-bath kinetic energy from the current (new) velocities — what a thermostat scales on.
    /// Empty for a single-bath run, where `kinetic_energy_new` carries the same number.
    pub kinetic_energy_new_bath: Vec<f64>,

    /// Bond energy per energy group (booked to the group of the bond's first atom).
    pub bond_energy: Vec<f64>,
    /// Angle energy per energy group.
    pub angle_energy: Vec<f64>,
    /// Improper-dihedral energy per energy group.
    pub improper_energy: Vec<f64>,
    /// Proper-dihedral energy per energy group.
    pub dihedral_energy: Vec<f64>,
    /// Cross-dihedral energy per energy group.
    pub crossdihedral_energy: Vec<f64>,
    /// Position-restraint energy per energy group.
    pub posrest_energy: Vec<f64>,
    /// Distance-restraint energy per energy group.
    pub distanceres_energy: Vec<f64>,
    /// Distance-field restraint energy per energy group.
    pub disfieldres_energy: Vec<f64>,
    /// Angle-restraint energy per energy group.
    pub angrest_energy: Vec<f64>,
    /// Dihedral-restraint energy per energy group.
    pub dihrest_energy: Vec<f64>,
    /// Constraint energy per energy group.
    pub constraints_energy: Vec<f64>,
    /// SASA energy per energy group.
    pub sasa_energy: Vec<f64>,
    /// SASA volume energy per energy group.
    pub sasa_volume_energy: Vec<f64>,
    /// RDC restraint energy per energy group.
    pub rdc_energy: Vec<f64>,

    /// Per-energy-group-pair LJ energies (energy_groups × energy_groups matrix).
    pub lj_energy: Vec<Vec<f64>>,
    /// Per-energy-group-pair CRF energies (energy_groups × energy_groups matrix).
    pub crf_energy: Vec<Vec<f64>>,
    /// Per-energy-group-pair lattice-sum k-space energies.
    pub ls_k_energy: Vec<Vec<f64>>,
    /// Per-energy-group-pair shifted-RF extra energy, original form.
    pub shift_extra_orig: Vec<Vec<f64>>,
    /// Per-energy-group-pair shifted-RF extra energy, physical form.
    pub shift_extra_phys: Vec<Vec<f64>>,

    /// Scalar virial (trace of the virial tensor), used for pressure coupling.
    pub virial_total: f64,

    /// PRECALCLAM (extended TI): the end-state energies on the λ grid; empty when off.
    pub precalc: PrecalcEnergies,

    /// EDS reference-state energy V_R (gromosXX `eds_vr`; part of `special_total`, not of
    /// `potential_total`). On the derivatives: Σ p_i dV_i/dλ.
    pub eds_vr: f64,
    /// EDS end-state energies V_i (gromosXX `eds_vi`; on the derivatives dV_i/dλ); empty
    /// without EDS.
    pub eds_vi: Vec<f64>,
    /// The special-term part of each V_i (gromosXX `eds_vi_special`).
    pub eds_vi_special: Vec<f64>,
    /// A-EDS: the unboosted mixed energy V_mix (gromosXX `eds_vmix`; in no total).
    pub eds_vmix: f64,
    /// A-EDS: E_max after this step's search update (gromosXX `eds_emax`).
    pub eds_emax: f64,
    /// A-EDS: E_min after this step's search update (gromosXX `eds_emin`).
    pub eds_emin: f64,
    /// A-EDS E_max/E_min search: the target state's mean of V_i − E_i^R (gromosXX `eds_globmin`).
    pub eds_globmin: f64,
    /// A-EDS E_max/E_min search: the target state's standard deviation of V_i
    /// (gromosXX `eds_globminfluc`).
    pub eds_globminfluc: f64,
    /// A-EDS: the offsets E_i^R after this step's search update (gromosXX `eds_eir`, the
    /// `# eds` offset column); sized with the end states, zero in plain EDS.
    pub eds_eir: Vec<f64>,
}

/// The extended-TI part of an [`Energy`] record (gromosXX `A_lj_energy` … `AB_disfld`,
/// `A_dihedral`, `B_dihedral`, `configuration/energy.h`, item 1.8). On the energies it holds the
/// A- and B-state energies at each λ_k of the PRECALCLAM grid; on the derivatives what
/// gromosXX books there (`A_de = +∂A/∂λ`, `B_de = −∂B/∂λ` for the nonbonded soft-core parts,
/// the full ∂E/∂λ at λ_k for the `AB_*` columns). The four group-pair matrices are flat
/// `[k][g1][g2]`; [`Self::calculate_totals`] folds them as gromosXX does.
#[derive(Debug, Clone, Default)]
pub struct PrecalcEnergies {
    n_groups: usize,
    /// A-state Lennard-Jones energy, `[k][g1][g2]`.
    pub a_lj: Vec<f64>,
    /// B-state Lennard-Jones energy, `[k][g1][g2]`.
    pub b_lj: Vec<f64>,
    /// A-state reaction-field energy, `[k][g1][g2]`.
    pub a_crf: Vec<f64>,
    /// B-state reaction-field energy, `[k][g1][g2]`.
    pub b_crf: Vec<f64>,
    /// Total of `a_lj` per λ_k (`calculate_totals`).
    pub a_lj_total: Vec<f64>,
    /// Total of `b_lj` per λ_k.
    pub b_lj_total: Vec<f64>,
    /// Total of `a_crf` per λ_k.
    pub a_crf_total: Vec<f64>,
    /// Total of `b_crf` per λ_k.
    pub b_crf_total: Vec<f64>,
    /// Kinetic energy of the perturbed atoms with their masses at λ_k.
    pub ab_kinetic: Vec<f64>,
    /// Perturbed bonds at λ_k.
    pub ab_bond: Vec<f64>,
    /// Perturbed angles at λ_k.
    pub ab_angle: Vec<f64>,
    /// Perturbed impropers at λ_k.
    pub ab_improper: Vec<f64>,
    /// Perturbed distance restraints at λ_k.
    pub ab_disres: Vec<f64>,
    /// Perturbed angle restraints (not modelled: zero).
    pub ab_angres: Vec<f64>,
    /// Perturbed dihedral restraints (not modelled: zero).
    pub ab_dihres: Vec<f64>,
    /// Perturbed distance-field restraints (not modelled: zero).
    pub ab_disfld: Vec<f64>,
    /// State-A energy of the perturbed proper dihedrals (no λ grid).
    pub a_dihedral: f64,
    /// State-B energy of the perturbed proper dihedrals.
    pub b_dihedral: f64,
}

impl PrecalcEnergies {
    /// Size for `n_lambdas` grid points and `n_groups` energy groups (gromosXX
    /// `Energy::resize`, `energy.cc:309-345`), zeroed.
    pub fn resize(&mut self, n_lambdas: usize, n_groups: usize) {
        let n_groups = n_groups.max(1);
        self.n_groups = n_groups;
        let cells = n_lambdas * n_groups * n_groups;
        for m in [
            &mut self.a_lj,
            &mut self.b_lj,
            &mut self.a_crf,
            &mut self.b_crf,
        ] {
            m.clear();
            m.resize(cells, 0.0);
        }
        for v in self.per_lambda_mut() {
            v.clear();
            v.resize(n_lambdas, 0.0);
        }
        self.a_dihedral = 0.0;
        self.b_dihedral = 0.0;
    }

    /// Number of λ grid points (0 = off).
    pub fn nr_lambdas(&self) -> usize {
        self.ab_kinetic.len()
    }

    /// Energy groups the matrices are sized for.
    pub fn n_groups(&self) -> usize {
        self.n_groups
    }

    fn per_lambda_mut(&mut self) -> [&mut Vec<f64>; 12] {
        [
            &mut self.a_lj_total,
            &mut self.b_lj_total,
            &mut self.a_crf_total,
            &mut self.b_crf_total,
            &mut self.ab_kinetic,
            &mut self.ab_bond,
            &mut self.ab_angle,
            &mut self.ab_improper,
            &mut self.ab_disres,
            &mut self.ab_angres,
            &mut self.ab_dihres,
            &mut self.ab_disfld,
        ]
    }

    /// gromosXX `zero(potential = true, kinetic = false)`: everything but `ab_kinetic`.
    pub fn zero_potential(&mut self) {
        for m in [
            &mut self.a_lj,
            &mut self.b_lj,
            &mut self.a_crf,
            &mut self.b_crf,
        ] {
            m.fill(0.0);
        }
        let [a, b, c, d, _kinetic, e, f, g, h, i, j, k] = self.per_lambda_mut();
        for v in [a, b, c, d, e, f, g, h, i, j, k] {
            v.fill(0.0);
        }
        self.a_dihedral = 0.0;
        self.b_dihedral = 0.0;
    }

    /// gromosXX `zero(potential = false, kinetic = true)`: `ab_kinetic` only.
    pub fn zero_kinetic(&mut self) {
        self.ab_kinetic.fill(0.0);
    }

    /// gromosXX `Energy::calculate_totals` (`energy.cc:397-423`): fold every matrix onto its
    /// lower triangle (`[k][j][i] += [k][i][j]`, `[k][i][j] = 0` for `j < i`) and sum the
    /// triangle in gromosXX's loop order.
    pub fn calculate_totals(&mut self) {
        let n = self.n_groups;
        let cells = n * n;
        for k in 0..self.nr_lambdas() {
            let base = k * cells;
            for (m, total) in [
                (&mut self.a_lj, &mut self.a_lj_total),
                (&mut self.b_lj, &mut self.b_lj_total),
                (&mut self.a_crf, &mut self.a_crf_total),
                (&mut self.b_crf, &mut self.b_crf_total),
            ] {
                total[k] = 0.0;
                for j in 0..n {
                    for i in j..n {
                        if i != j {
                            m[base + i * n + j] += m[base + j * n + i];
                            m[base + j * n + i] = 0.0;
                        }
                        total[k] += m[base + i * n + j];
                    }
                }
            }
        }
    }

    /// The twelve columns of `# precalclam` row `k` (gromosXX `out_configuration.cc:3237-3252`).
    pub fn row(&self, k: usize) -> [f64; 12] {
        [
            self.a_lj_total[k],
            self.b_lj_total[k],
            self.a_crf_total[k],
            self.b_crf_total[k],
            self.ab_kinetic[k],
            self.ab_bond[k],
            self.ab_angle[k],
            self.ab_improper[k],
            self.ab_disres[k],
            self.ab_angres[k],
            self.ab_dihres[k],
            self.ab_disfld[k],
        ]
    }
}

impl Energy {
    /// Create a zeroed energy accumulator for the given number of temperature and energy groups.
    pub fn new(num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        let mut e = Self {
            kinetic_total: 0.0,
            potential_total: 0.0,
            bond_total: 0.0,
            angle_total: 0.0,
            dihedral_total: 0.0,
            improper_total: 0.0,
            cross_dihedral_total: 0.0,
            lj_total: 0.0,
            crf_total: 0.0,
            ls_total: 0.0,
            ls_pair_total: 0.0,
            ls_realspace_total: 0.0,
            ls_kspace_total: 0.0,
            ls_a_term_total: 0.0,
            ls_self_total: 0.0,
            ls_surface_total: 0.0,
            special_total: 0.0,
            term_energies: Vec::new(),
            posrest_total: 0.0,
            distanceres_total: 0.0,
            sasa_total: 0.0,
            constraint_total: 0.0,
            entropy_term: 0.0,
            kinetic_energy_new: 0.0,
            kinetic_energy: Vec::new(),
            com_kinetic_energy: Vec::new(),
            ir_kinetic_energy: Vec::new(),
            kinetic_energy_new_bath: Vec::new(),
            bond_energy: Vec::new(),
            angle_energy: Vec::new(),
            improper_energy: Vec::new(),
            dihedral_energy: Vec::new(),
            crossdihedral_energy: Vec::new(),
            posrest_energy: Vec::new(),
            distanceres_energy: Vec::new(),
            disfieldres_energy: Vec::new(),
            angrest_energy: Vec::new(),
            dihrest_energy: Vec::new(),
            constraints_energy: Vec::new(),
            sasa_energy: Vec::new(),
            sasa_volume_energy: Vec::new(),
            rdc_energy: Vec::new(),
            lj_energy: Vec::new(),
            crf_energy: Vec::new(),
            ls_k_energy: Vec::new(),
            shift_extra_orig: Vec::new(),
            shift_extra_phys: Vec::new(),
            virial_total: 0.0,
            precalc: PrecalcEnergies::default(),
            eds_vr: 0.0,
            eds_vi: Vec::new(),
            eds_vi_special: Vec::new(),
            eds_vmix: 0.0,
            eds_emax: 0.0,
            eds_emin: 0.0,
            eds_globmin: 0.0,
            eds_globminfluc: 0.0,
            eds_eir: Vec::new(),
        };
        e.resize(num_energy_groups, num_temperature_groups);
        e
    }

    /// Size the per-group arrays to `n_groups`, the group-pair matrices to `n_groups × n_groups`
    /// and the per-bath arrays to `n_baths` (gromosXX `Energy::resize(num, numb)`), zeroing them.
    pub fn resize(&mut self, n_groups: usize, n_baths: usize) {
        for v in self.group_arrays_mut() {
            v.clear();
            v.resize(n_groups, 0.0);
        }
        for m in self.group_matrices_mut() {
            m.clear();
            m.resize(n_groups, vec![0.0; n_groups]);
        }
        for v in [
            &mut self.kinetic_energy,
            &mut self.com_kinetic_energy,
            &mut self.ir_kinetic_energy,
        ] {
            v.clear();
            v.resize(n_baths, 0.0);
        }
        let n_lambdas = self.precalc.nr_lambdas();
        self.precalc.resize(n_lambdas, n_groups);
    }

    /// Number of energy groups the per-group arrays are sized for.
    pub fn num_energy_groups(&self) -> usize {
        self.bond_energy.len()
    }

    fn group_arrays_mut(&mut self) -> [&mut Vec<f64>; 14] {
        [
            &mut self.bond_energy,
            &mut self.angle_energy,
            &mut self.improper_energy,
            &mut self.dihedral_energy,
            &mut self.crossdihedral_energy,
            &mut self.posrest_energy,
            &mut self.distanceres_energy,
            &mut self.disfieldres_energy,
            &mut self.angrest_energy,
            &mut self.dihrest_energy,
            &mut self.constraints_energy,
            &mut self.sasa_energy,
            &mut self.sasa_volume_energy,
            &mut self.rdc_energy,
        ]
    }

    fn group_matrices_mut(&mut self) -> [&mut Vec<Vec<f64>>; 5] {
        [
            &mut self.lj_energy,
            &mut self.crf_energy,
            &mut self.ls_k_energy,
            &mut self.shift_extra_orig,
            &mut self.shift_extra_phys,
        ]
    }

    /// With a single energy group every group cell is the total of its term, as gromosXX's
    /// group sums are: book the totals into group 0 (and pair (0, 0)). No-op for any other
    /// group count — those cells are booked by the terms themselves.
    pub fn book_totals_into_single_group(&mut self) {
        if self.bond_energy.len() != 1 || self.lj_energy.len() != 1 {
            return;
        }
        self.bond_energy[0] = self.bond_total;
        self.angle_energy[0] = self.angle_total;
        self.improper_energy[0] = self.improper_total;
        self.dihedral_energy[0] = self.dihedral_total;
        self.crossdihedral_energy[0] = self.cross_dihedral_total;
        self.posrest_energy[0] = self.posrest_total;
        self.distanceres_energy[0] = self.distanceres_total;
        self.constraints_energy[0] = self.constraint_total;
        self.sasa_energy[0] = self.sasa_total;
        self.lj_energy[0][0] = self.lj_total;
        self.crf_energy[0][0] = self.crf_total;
    }

    /// Total energy (kinetic + potential + special), mirroring GROMOS:
    /// total = potential_total + kinetic_total + special_total
    /// where special_total includes distanceres_total, posrest, angrest, etc.
    pub fn total(&self) -> f64 {
        self.kinetic_total + self.potential_total + self.special_total + self.distanceres_total
    }

    /// Update potential total from components (bonded + nonbonded, without special)
    pub fn update_potential_total(&mut self) {
        self.potential_total = self.bond_total
            + self.angle_total
            + self.dihedral_total
            + self.improper_total
            + self.cross_dihedral_total
            + self.lj_total
            + self.crf_total
            + self.ls_total
            + self.sasa_total;
    }

    /// Clear all energies (sizes are kept).
    pub fn clear(&mut self) {
        self.kinetic_total = 0.0;
        self.potential_total = 0.0;
        self.bond_total = 0.0;
        self.angle_total = 0.0;
        self.dihedral_total = 0.0;
        self.improper_total = 0.0;
        self.cross_dihedral_total = 0.0;
        self.lj_total = 0.0;
        self.crf_total = 0.0;
        self.ls_total = 0.0;
        self.ls_pair_total = 0.0;
        self.ls_realspace_total = 0.0;
        self.ls_kspace_total = 0.0;
        self.ls_a_term_total = 0.0;
        self.ls_self_total = 0.0;
        self.ls_surface_total = 0.0;
        self.special_total = 0.0;
        self.term_energies.clear();
        self.posrest_total = 0.0;
        self.distanceres_total = 0.0;
        self.sasa_total = 0.0;
        self.constraint_total = 0.0;
        self.entropy_term = 0.0;
        self.kinetic_energy_new = 0.0;
        self.kinetic_energy_new_bath.fill(0.0);
        self.virial_total = 0.0;

        self.kinetic_energy.fill(0.0);
        self.com_kinetic_energy.fill(0.0);
        self.ir_kinetic_energy.fill(0.0);
        for v in self.group_arrays_mut() {
            v.fill(0.0);
        }
        for m in self.group_matrices_mut() {
            for row in m.iter_mut() {
                row.fill(0.0);
            }
        }
        self.precalc.zero_potential();
        self.precalc.zero_kinetic();
        // gromosXX `Energy::zero` (`energy.cc:151-163`)
        self.eds_vr = 0.0;
        self.eds_vi.fill(0.0);
        self.eds_vi_special.fill(0.0);
        self.eds_vmix = 0.0;
        self.eds_emax = 0.0;
        self.eds_emin = 0.0;
        self.eds_globmin = 0.0;
        self.eds_globminfluc = 0.0;
        self.eds_eir.fill(0.0);
    }

    /// Size the EDS end-state arrays for `num_states` states (gromosXX
    /// `Configuration::init`, `configuration.cc:403-419`), zeroed.
    pub fn resize_eds(&mut self, num_states: usize) {
        self.eds_vr = 0.0;
        for v in [
            &mut self.eds_vi,
            &mut self.eds_vi_special,
            &mut self.eds_eir,
        ] {
            v.clear();
            v.resize(num_states, 0.0);
        }
    }
}

/// State of the system at a single time point
///
/// Contains positions, velocities, forces, and all associated properties
#[derive(Debug, Clone)]
pub struct State {
    /// Atomic positions (nm).
    pub pos: Vec<Vec3>,
    /// Atomic velocities (nm ps⁻¹).
    pub vel: Vec<Vec3>,
    /// Total forces on each atom (kJ mol⁻¹ nm⁻¹).
    pub force: Vec<Vec3>,
    /// Holonomic constraint forces (kJ mol⁻¹ nm⁻¹).
    pub constraint_force: Vec<Vec3>,
    /// Simulation box geometry.
    pub box_config: Box,
    /// Virial tensor (kJ mol⁻¹); used for pressure calculation.
    pub virial_tensor: Mat3,
    /// Kinetic-energy tensor (kJ mol⁻¹); used for pressure calculation.
    pub kinetic_energy_tensor: Mat3,
    /// Pressure tensor (bar); computed by `PressureCalculation`.
    pub pressure_tensor: Mat3,
    /// All energy components for this state.
    pub energies: Energy,
    /// dH/dλ of every energy component, laid out like `energies` (gromosXX
    /// `perturbed_energy_derivatives`); zero unless the run is perturbed.
    pub perturbed_energy_derivatives: Energy,
}

impl State {
    /// Create a zeroed state for `n_atoms` atoms with the given group counts.
    pub fn new(n_atoms: usize, num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            pos: vec![Vec3::ZERO; n_atoms],
            vel: vec![Vec3::ZERO; n_atoms],
            force: vec![Vec3::ZERO; n_atoms],
            constraint_force: vec![Vec3::ZERO; n_atoms],
            box_config: Box::vacuum(),
            virial_tensor: Mat3::ZERO,
            kinetic_energy_tensor: Mat3::ZERO,
            pressure_tensor: Mat3::ZERO,
            energies: Energy::new(num_temperature_groups, num_energy_groups),
            perturbed_energy_derivatives: Energy::new(num_temperature_groups, num_energy_groups),
        }
    }

    /// Clear forces and virial
    pub fn clear_forces(&mut self) {
        self.force.fill(Vec3::ZERO);
        self.virial_tensor = Mat3::ZERO;
        self.kinetic_energy_tensor = Mat3::ZERO;
    }

    /// Clear constraint forces
    pub fn clear_constraint_forces(&mut self) {
        self.constraint_force.fill(Vec3::ZERO);
    }

    /// Calculate kinetic energy
    ///
    /// KE = 0.5 * sum_i(m_i * v_i^2)
    pub fn calculate_kinetic_energy(&mut self, masses: &[f64]) {
        self.energies.kinetic_total = self
            .vel
            .iter()
            .zip(masses.iter())
            .map(|(v, &m)| 0.5 * m * v.length_squared())
            .sum();
    }

    /// Calculate kinetic energy per group
    pub fn calculate_kinetic_energy_groups(
        &mut self,
        masses: &[f64],
        temperature_groups: &[Vec<usize>],
    ) {
        for (group_idx, group) in temperature_groups.iter().enumerate() {
            self.energies.kinetic_energy[group_idx] = group
                .iter()
                .map(|&i| 0.5 * masses[i] * self.vel[i].length_squared())
                .sum();
        }

        // Total kinetic energy is sum of all groups
        self.energies.kinetic_total = self.energies.kinetic_energy.iter().sum();
    }

    /// Calculate temperature from kinetic energy
    ///
    /// T = 2 * KE / (k_B * N_dof)
    /// where N_dof = 3N - N_constraints
    pub fn temperature(&self, num_degrees_of_freedom: usize) -> f64 {
        2.0 * self.energies.kinetic_total / (crate::units::kB * num_degrees_of_freedom as f64)
    }

    /// The pressure, trace / 3 of `pressure_tensor`: what the barostat couples to and what
    /// gromosXX writes (`out_configuration.cc:3322-3323`). Zero unless the run computes it
    /// (pressure coupling).
    pub fn pressure(&self) -> f64 {
        let p = &self.pressure_tensor;
        (p.x_axis.x + p.y_axis.y + p.z_axis.z) / 3.0
    }
}

/// Configuration - manages current and previous states with double buffering
///
/// GROMOS uses double buffering (pointer swapping) for efficient state management.
/// This is translated to Rust using two State structs and an index.
#[derive(Debug, Clone)]
pub struct Configuration {
    state1: State,
    state2: State,
    current_idx: usize, // 0 or 1

    /// Lattice shifts for periodic boundary tracking (FEP)
    /// Stores how many times each atom crossed periodic boundaries
    /// Dimension: [n_atoms][3] for x, y, z shifts
    pub lattice_shifts: Vec<[i32; 3]>,

    /// EDS end-state forces and virials of the current step (gromosXX `conf.special().eds`);
    /// empty without EDS.
    pub eds: EdsEndStates,
    /// The A-EDS parameter-search state the run carries from step to step and from job to job
    /// (the configuration's `AEDSSEARCH` block); `None` without A-EDS.
    pub aeds: Option<AedsSearch>,
    /// The Nosé–Hoover variables ζ, one chain per bath (gromosXX `bath_struct::zeta`, the
    /// configuration's `NHCVARIABLES` block): before the run, the ones read from the start
    /// configuration; once the thermostat runs, its current ones. `None` without them.
    pub nose_hoover: Option<Vec<Vec<f64>>>,
}

/// The mutable state of an A-EDS run: E_max, E_min and the offsets E_i^R with everything their
/// searches accumulate (gromosXX keeps it in `sim.param().eds`, `simulation/parameter.h:3190-3312`,
/// and writes it as the final configuration's `AEDSSEARCH` block, `out_configuration.cc:3985-4004`).
/// Per-state vectors are indexed by end state.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AedsSearch {
    /// E_max (kJ/mol)
    pub emax: f64,
    /// E_min (kJ/mol)
    pub emin: f64,
    /// highest V_mix at a state change since the last round trip
    pub searchemax: f64,
    /// number of round trips (every state visited) so far
    pub emaxcounts: u32,
    /// the state of the previous step
    pub oldstate: u32,
    /// `false` while E_min is held at or above the lowest average end-state energy (RESTREMIN
    /// 1, until the first round trip)
    pub fullemin: bool,
    /// E_i^R
    pub eir: Vec<f64>,
    /// log of the memory-averaged exp(−β(V_i* − V_R*))
    pub lnexpde: Vec<f64>,
    /// −lnexpde/β
    pub statefren: Vec<f64>,
    /// visited since the last round trip
    pub visited: Vec<bool>,
    /// steps spent in state i
    pub visitcounts: Vec<u32>,
    /// running mean of V_i over the steps in state i
    pub avgenergy: Vec<f64>,
    /// running mean of V_i − E_i^R over the steps in state i
    pub eiravgenergy: Vec<f64>,
    /// Welford sum of squares of V_i
    pub bigs: Vec<f64>,
    /// standard deviation of V_i over the steps in state i
    pub stdevenergy: Vec<f64>,
    /// FORMS 5, 6: steps with V_i − V_R* within kT above the state's free energy (written,
    /// never used)
    pub framecounts: Vec<u32>,
    /// The next step starts a fresh search (NTIAEDSS = 1; gromosXX `initaedssearch`). Not part
    /// of the `AEDSSEARCH` block.
    pub initialise: bool,
}

impl AedsSearch {
    /// A fresh state for `num_states` states (gromosXX `read_AEDS`, `in_parameter.cc:4124-4136`).
    pub fn new(emax: f64, emin: f64, eir: Vec<f64>, fullemin: bool) -> Self {
        let n = eir.len();
        Self {
            emax,
            emin,
            searchemax: 0.0,
            emaxcounts: 0,
            oldstate: 0,
            fullemin,
            eir,
            lnexpde: vec![0.0; n],
            statefren: vec![0.0; n],
            visited: vec![false; n],
            visitcounts: vec![0; n],
            avgenergy: vec![0.0; n],
            eiravgenergy: vec![0.0; n],
            bigs: vec![0.0; n],
            stdevenergy: vec![0.0; n],
            framecounts: vec![0; n],
            initialise: true,
        }
    }

    /// Number of end states.
    pub fn num_states(&self) -> usize {
        self.eir.len()
    }
}

/// The state-dependent forces and virials of an EDS run, one entry per end state, written by
/// the force field and mixed into the forces by the EDS algorithm (gromosXX
/// `configuration::Configuration::special_struct::eds_struct`). Not double-buffered.
#[derive(Debug, Clone, Default)]
pub struct EdsEndStates {
    /// `forces[s][i]`: the force on atom i from the terms of end state s.
    pub forces: Vec<Vec<Vec3>>,
    /// `virials[s][a][b] = Σ r_a f_b` of end state s (the convention of `ForceStorage::virial`).
    pub virials: Vec<[[f64; 3]; 3]>,
}

impl EdsEndStates {
    /// Size for `num_states` end states over `n_atoms` atoms, zeroed.
    pub fn resize(&mut self, num_states: usize, n_atoms: usize) {
        self.forces.clear();
        self.forces.resize(num_states, vec![Vec3::ZERO; n_atoms]);
        self.virials.clear();
        self.virials.resize(num_states, [[0.0; 3]; 3]);
    }

    /// Number of end states (0 without EDS).
    pub fn num_states(&self) -> usize {
        self.forces.len()
    }
}

impl Configuration {
    /// Create a new double-buffered configuration for `n_atoms` atoms.
    pub fn new(n_atoms: usize, num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            state1: State::new(n_atoms, num_temperature_groups, num_energy_groups),
            state2: State::new(n_atoms, num_temperature_groups, num_energy_groups),
            current_idx: 0,
            lattice_shifts: vec![[0, 0, 0]; n_atoms],
            eds: EdsEndStates::default(),
            aeds: None,
            nose_hoover: None,
        }
    }

    /// Get reference to current state
    #[inline]
    pub fn current(&self) -> &State {
        if self.current_idx == 0 {
            &self.state1
        } else {
            &self.state2
        }
    }

    /// Get mutable reference to current state
    #[inline]
    pub fn current_mut(&mut self) -> &mut State {
        if self.current_idx == 0 {
            &mut self.state1
        } else {
            &mut self.state2
        }
    }

    /// Get reference to old (previous) state
    #[inline]
    pub fn old(&self) -> &State {
        if self.current_idx == 0 {
            &self.state2
        } else {
            &self.state1
        }
    }

    /// Get mutable reference to old state
    #[inline]
    pub fn old_mut(&mut self) -> &mut State {
        if self.current_idx == 0 {
            &mut self.state2
        } else {
            &mut self.state1
        }
    }

    /// The current state and the EDS end states together (both borrowed mutably).
    pub fn current_and_eds_mut(&mut self) -> (&mut State, &mut EdsEndStates) {
        let state = if self.current_idx == 0 {
            &mut self.state1
        } else {
            &mut self.state2
        };
        (state, &mut self.eds)
    }

    /// Exchange current and old states (zero-cost pointer swap)
    ///
    /// This is the Rust equivalent of GROMOS C++ pointer swapping
    #[inline]
    pub fn exchange_state(&mut self) {
        self.current_idx = 1 - self.current_idx;
    }

    /// Copy current state to old state
    pub fn copy_current_to_old(&mut self) {
        if self.current_idx == 0 {
            self.state2 = self.state1.clone();
        } else {
            self.state1 = self.state2.clone();
        }
    }
}

/// Stochastic dynamics additional variables
#[derive(Debug, Clone)]
pub struct StochasticVariables {
    /// Per-atom friction coefficients γ (ps⁻¹) for Langevin dynamics.
    pub gamma: Vec<f64>,
    /// Per-atom random forces (kJ mol⁻¹ nm⁻¹) drawn each step.
    pub random_force: Vec<Vec3>,
}

impl StochasticVariables {
    /// Create zeroed stochastic variables for `n_atoms` atoms.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            gamma: vec![0.0; n_atoms],
            random_force: vec![Vec3::ZERO; n_atoms],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn precalc_totals_fold_the_group_matrix_as_gromosxx_does() {
        // two λ points, two groups: [k][g1][g2], the upper triangle moves onto the lower one
        let mut p = PrecalcEnergies::default();
        p.resize(2, 2);
        for k in 0..2 {
            for (c, v) in [1.0, 2.0, 4.0, 8.0].iter().enumerate() {
                p.a_lj[k * 4 + c] = v * (k + 1) as f64;
                p.b_crf[k * 4 + c] = -v;
            }
        }
        p.ab_kinetic.fill(3.0);
        p.calculate_totals();
        assert_eq!(p.a_lj_total, vec![15.0, 30.0]);
        assert_eq!(p.b_crf_total, vec![-15.0, -15.0]);
        // [0][1] folded into [1][0]
        assert_eq!(&p.a_lj[0..4], &[1.0, 0.0, 6.0, 8.0]);
        p.zero_potential();
        assert!(p.a_lj.iter().chain(&p.a_lj_total).all(|&v| v == 0.0));
        assert_eq!(
            p.ab_kinetic,
            vec![3.0, 3.0],
            "the kinetic column has its own flag"
        );
        p.zero_kinetic();
        assert_eq!(p.ab_kinetic, vec![0.0, 0.0]);
        assert_eq!(p.row(1), [0.0; 12]);
    }

    #[test]
    fn test_box_rectangular() {
        let box_config = Box::rectangular(3.0, 4.0, 5.0);

        assert_eq!(box_config.box_type, BoxType::Rectangular);
        assert_eq!(box_config.dimensions(), Vec3::new(3.0, 4.0, 5.0));
        assert!((box_config.volume() - 60.0).abs() < 1e-6);
    }

    #[test]
    fn test_energy_total() {
        let mut energy = Energy::new(1, 1);

        energy.kinetic_total = 100.0;
        energy.bond_total = 50.0;
        energy.lj_total = 30.0;
        energy.crf_total = 20.0;

        energy.update_potential_total();

        assert_eq!(energy.potential_total, 100.0);
        assert_eq!(energy.total(), 200.0);
    }

    #[test]
    fn energy_is_laid_out_per_group_pair_and_bath() {
        let mut e = Energy::new(2, 3);
        assert_eq!(e.num_energy_groups(), 3);
        assert_eq!((e.lj_energy.len(), e.ls_k_energy[2].len()), (3, 3));
        assert_eq!((e.kinetic_energy.len(), e.ir_kinetic_energy.len()), (2, 2));
        e.resize(1, 4);
        assert_eq!((e.rdc_energy.len(), e.shift_extra_phys.len()), (1, 1));
        assert_eq!(e.com_kinetic_energy.len(), 4);
        // a single group carries the totals
        e.bond_total = 1.5;
        e.lj_total = -2.0;
        e.posrest_total = 0.25;
        e.book_totals_into_single_group();
        assert_eq!(
            (e.bond_energy[0], e.lj_energy[0][0], e.posrest_energy[0]),
            (1.5, -2.0, 0.25)
        );
        e.clear();
        assert_eq!((e.bond_energy[0], e.lj_energy[0][0]), (0.0, 0.0));
        // more than one group: the cells are the terms' to book
        let mut two = Energy::new(1, 2);
        two.bond_total = 1.0;
        two.book_totals_into_single_group();
        assert_eq!(two.bond_energy, vec![0.0, 0.0]);
    }

    #[test]
    fn test_state_kinetic_energy() {
        let mut state = State::new(2, 1, 1);
        let masses = vec![1.0, 2.0];

        // Set velocities
        state.vel[0] = Vec3::new(1.0, 0.0, 0.0);
        state.vel[1] = Vec3::new(0.0, 2.0, 0.0);

        state.calculate_kinetic_energy(&masses);

        // KE = 0.5 * (1.0 * 1.0^2 + 2.0 * 2.0^2) = 0.5 * (1 + 8) = 4.5
        assert!((state.energies.kinetic_total - 4.5).abs() < 1e-6);
    }

    #[test]
    fn test_temperature_calculation() {
        let mut state = State::new(100, 1, 1);
        let masses = vec![1.0; 100];

        // Set random velocities for testing
        for vel in &mut state.vel {
            *vel = Vec3::new(0.1, 0.1, 0.1);
        }

        state.calculate_kinetic_energy(&masses);

        // 3N degrees of freedom (no constraints)
        let temp = state.temperature(300);

        assert!(temp > 0.0);
        assert!(temp < 100.0); // Sanity check
    }

    #[test]
    fn test_state_exchange() {
        let mut conf = Configuration::new(2, 1, 1);

        // Set some values in current state
        conf.current_mut().pos[0] = Vec3::new(1.0, 2.0, 3.0);

        // Exchange
        conf.exchange_state();

        // Old state should now have the values
        assert_eq!(conf.old().pos[0], Vec3::new(1.0, 2.0, 3.0));

        // Current state should be fresh
        assert_eq!(conf.current().pos[0], Vec3::ZERO);
    }

    #[test]
    fn test_clear_forces() {
        let mut state = State::new(10, 1, 1);

        // Set some forces
        state.force[0] = Vec3::new(1.0, 2.0, 3.0);
        state.virial_tensor = Mat3::from_cols(
            Vec3::new(10.0, 0.0, 0.0),
            Vec3::new(0.0, 10.0, 0.0),
            Vec3::new(0.0, 0.0, 10.0),
        );

        state.clear_forces();

        assert_eq!(state.force[0], Vec3::ZERO);
        assert_eq!(state.virial_tensor, Mat3::ZERO);
    }
}
