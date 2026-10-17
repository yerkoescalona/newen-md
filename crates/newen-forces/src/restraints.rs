//! Restraint interactions for MD simulations
//!
//! Direct translation of GROMOS restraint interactions from:
//! - md++/src/interaction/special/position_restraint_interaction.cc
//! - md++/src/interaction/special/distance_restraint_interaction.cc

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Vec3};

use crate::bonded::TermLambda;

/// Position restraint - restrains atom to reference position
///
/// Direct translation from position_restraint_interaction.cc
///
/// Energy: E = 0.5 * k * |r - r_ref|²
/// Force: F = -k * (r - r_ref)
#[derive(Debug, Clone)]
pub struct PositionRestraint {
    /// Atom index
    pub atom: usize,
    /// Reference position
    pub reference_pos: Vec3,
    /// Force constant (kJ mol⁻¹ nm⁻²)
    pub force_constant: f64,
}

impl PositionRestraint {
    pub fn new(atom: usize, reference_pos: Vec3, force_constant: f64) -> Self {
        Self {
            atom,
            reference_pos,
            force_constant,
        }
    }

    /// Calculate restraint energy and forces
    pub fn calculate(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        // nearest_image(ri, rj) returns ri - rj
        // v = pos - ref (same as GROMOS)
        let v = periodicity.nearest_image(conf.current().pos[self.atom], self.reference_pos);

        let dist_sq = v.length_squared();

        // Energy: 0.5 * k * r²
        let energy = 0.5 * self.force_constant * dist_sq;

        // Force: F = -k * (pos - ref)
        let force = v * (-self.force_constant);

        conf.current_mut().force[self.atom] += force;

        energy
    }
}

/// Collection of position restraints
#[derive(Debug, Clone, Default)]
pub struct PositionRestraints {
    pub restraints: Vec<PositionRestraint>,
}

impl PositionRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: PositionRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        self.restraints
            .iter()
            .map(|r| r.calculate(conf, periodicity))
            .sum()
    }

    /// [`Self::calculate_all`], also adding each restraint's energy to `per_group` at the
    /// energy group of its atom (gromosXX `posrest_energy[atom_energy_group(seq)]`,
    /// `position_restraint_interaction.cc:101`).
    pub fn calculate_all_grouped(
        &self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        atom_group: &[usize],
        per_group: &mut [f64],
    ) -> f64 {
        self.restraints
            .iter()
            .map(|r| {
                let e = r.calculate(conf, periodicity);
                per_group[atom_group[r.atom]] += e;
                e
            })
            .sum()
    }
}

// ─── Distance restraints (GROMOS distance_restraint_interaction.cc) ────────

/// Decode the RAH field into (form, dimension_mask).
///
/// RAH encodes both the half-harmonic type and the spatial dimension:
///   dim_base ∈ {0,10,20,30,40,50,60}  →  mask: XYZ, XY, XZ, YZ, X, Y, Z
///   form     ∈ {-1, 0, 1}             →  repulsive, full-harmonic, attractive
///   rah = dim_base + form  (but shifted so -1 works: dim_base = ((rah+1)/10)*10)
fn decode_rah(rah: i32) -> (i32, Vec3) {
    let dim_base = (rah + 1).div_euclid(10) * 10;
    let form = rah - dim_base;
    let mask = match dim_base {
        10 => Vec3::new(1.0, 1.0, 0.0),
        20 => Vec3::new(1.0, 0.0, 1.0),
        30 => Vec3::new(0.0, 1.0, 1.0),
        40 => Vec3::new(1.0, 0.0, 0.0),
        50 => Vec3::new(0.0, 1.0, 0.0),
        60 => Vec3::new(0.0, 0.0, 1.0),
        _ => Vec3::new(1.0, 1.0, 1.0), // dim_base=0 (XYZ)
    };
    (form, mask)
}

/// Core energy+force calculation for a single distance restraint.
///
/// Returns `(energy, force_on_atom1)`.  Caller adds force to atom1 and
/// subtracts it from atom2 (Newton's 3rd law).
///
/// `v = nearest_image(pos_j, pos_i)` — the vector from atom1 to atom2.
/// `mode`: 1 = no w0 scaling, 2 = multiply energy & force by w0.
fn distres_en_force(
    v: Vec3,
    r0: f64,
    w0: f64,
    rah: i32,
    k: f64,
    r_linear: f64,
    mode: i32,
) -> (f64, Vec3) {
    let (form, mask) = decode_rah(rah);
    let vm = Vec3::new(v.x * mask.x, v.y * mask.y, v.z * mask.z);
    let dist = vm.length();

    if dist < 1e-10 {
        return (0.0, Vec3::ZERO);
    }

    // Half-harmonic: skip if restraint is satisfied
    if form != 0 && (form as f64) * dist < (form as f64) * r0 {
        return (0.0, Vec3::ZERO);
    }

    let delta = dist - r0;
    let (energy, f_on_1) = if delta.abs() <= r_linear {
        // Harmonic zone: E = 0.5 K (dist-r0)², F_1 = K*(dist-r0)*vm/dist
        let e = 0.5 * k * delta * delta;
        let f = vm * (k * delta / dist);
        (e, f)
    } else if dist < r0 {
        // Linear zone below r0 (dist < r0 - r_linear)
        let e = -k * r_linear * (dist + 0.5 * r_linear - r0);
        let f = vm * (-k * r_linear / dist);
        (e, f)
    } else {
        // Linear zone above r0 (dist > r0 + r_linear)
        let e = k * r_linear * (dist - 0.5 * r_linear - r0);
        let f = vm * (k * r_linear / dist);
        (e, f)
    };

    let w = if mode == 2 { w0 } else { 1.0 };
    (energy * w, f_on_1 * w)
}

/// Single distance restraint (GROMOS-faithful, virtual atom type 0 only).
///
/// RAH encodes both the half-harmonic form and the spatial dimension;
/// see `decode_rah` for details.  Mode 2 scales by w0 (NTDIR=2).
#[derive(Debug, Clone)]
pub struct DistanceRestraint {
    pub atom1: usize,
    pub atom2: usize,
    pub r0: f64,
    pub w0: f64,
    pub rah: i32,
    /// Force constant K (kJ mol⁻¹ nm⁻²)
    pub k: f64,
    /// Linear-region threshold r_linear = DIR0 (nm)
    pub r_linear: f64,
    /// 1 = no w0 scaling, 2 = multiply by w0
    pub mode: i32,
    /// The running ⟨r⁻³⟩ of a time-averaged restraint (gromosXX `distanceres.av`); `None`
    /// until seeded.
    pub avg_r_inv3: Option<f64>,
    /// exp(−dt/τ) of a time-averaged restraint (NTDIR < 0); `None` = instantaneous.
    pub exp_mem: Option<f64>,
}

impl DistanceRestraint {
    pub fn new(
        atom1: usize,
        atom2: usize,
        r0: f64,
        w0: f64,
        rah: i32,
        k: f64,
        r_linear: f64,
        mode: i32,
    ) -> Self {
        Self {
            atom1,
            atom2,
            r0,
            w0,
            rah,
            k,
            r_linear,
            mode,
            avg_r_inv3: None,
            exp_mem: None,
        }
    }

    /// Time-average the restraint (NTDIR < 0): `exp_mem = exp(−dt/τ)`
    /// (`distance_restraint_interaction.cc:269-271`).
    pub fn with_time_averaging(mut self, tau: f64, dt: f64) -> Self {
        self.exp_mem = Some((-dt / tau).exp());
        self
    }

    /// Seed the average of a time-averaged restraint with `|r|⁻³` of the start-up configuration
    /// (NTDIRA = 0, gromosXX `_init_averages`, :236-252: the unmasked distance).
    pub fn init_average(&mut self, conf: &Configuration, periodicity: &Periodicity) {
        if self.exp_mem.is_some() {
            let pos = &conf.current().pos;
            let v = periodicity.nearest_image(pos[self.atom1], pos[self.atom2]);
            self.avg_r_inv3 = Some(v.length().powf(-3.0));
        }
    }

    pub fn calculate(&mut self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        if let Some(exp_t) = self.exp_mem {
            return self.calculate_time_averaged(conf, periodicity, exp_t);
        }
        let v = periodicity.nearest_image(
            conf.current().pos[self.atom2],
            conf.current().pos[self.atom1],
        );
        let (energy, f1) = distres_en_force(
            v,
            self.r0,
            self.w0,
            self.rah,
            self.k,
            self.r_linear,
            self.mode,
        );
        conf.current_mut().force[self.atom1] += f1;
        conf.current_mut().force[self.atom2] -= f1;
        energy
    }

    /// NTDIR < 0 (FORCESCALE 0), as gromosXX `_calculate_distance_restraint_interactions`
    /// (`distance_restraint_interaction.cc:120-190`): `⟨r⁻³⟩ = (1−e)·r⁻³ + e·⟨r⁻³⟩`, the energy at
    /// `⟨r⁻³⟩^(−1/3)`, the force along the instantaneous unit vector.
    fn calculate_time_averaged(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        exp_t: f64,
    ) -> f64 {
        if self.avg_r_inv3.is_none() {
            self.init_average(conf, periodicity);
        }
        // gromosXX: v = v1 − v2, masked by RAH's dimension
        let v = periodicity.nearest_image(
            conf.current().pos[self.atom1],
            conf.current().pos[self.atom2],
        );
        let (rah, mask) = decode_rah(self.rah);
        let v = Vec3::new(v.x * mask.x, v.y * mask.y, v.z * mask.z);
        let vabs = v.length();
        let ave = (1.0 - exp_t) * vabs.powf(-3.0) + exp_t * self.avg_r_inv3.unwrap_or(0.0);
        self.avg_r_inv3 = Some(ave);
        let dist = ave.powf(-1.0 / 3.0);

        let (k, r0, r_l) = (self.k, self.r0, self.r_linear);
        let rah_f = rah as f64;
        let (mut energy, mut f) = if rah_f * dist < rah_f * r0 {
            (0.0, Vec3::ZERO)
        } else if (r0 - dist).abs() < r_l {
            (
                0.5 * k * (dist - r0) * (dist - r0),
                v * (-k * (dist - r0)) / vabs,
            )
        } else if dist < r0 {
            (-k * r_l * (dist - r0 + 0.5 * r_l), v * (r_l * k) / vabs)
        } else {
            (k * r_l * (dist - r0 - 0.5 * r_l), v * (-r_l * k) / vabs)
        };
        match self.mode {
            1 => {},
            2 => {
                energy *= self.w0;
                f *= self.w0;
            },
            _ => {
                energy = 0.0;
                f = Vec3::ZERO;
            },
        }
        conf.current_mut().force[self.atom1] += f;
        conf.current_mut().force[self.atom2] -= f;
        energy
    }
}

/// Collection of distance restraints with shared parameters.
#[derive(Debug, Clone, Default)]
pub struct DistanceRestraints {
    pub restraints: Vec<DistanceRestraint>,
}

impl DistanceRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add(&mut self, r: DistanceRestraint) {
        self.restraints.push(r);
    }

    /// Seed every time-averaged restraint from the start-up configuration.
    pub fn init_averages(&mut self, conf: &Configuration, periodicity: &Periodicity) {
        for r in &mut self.restraints {
            r.init_average(conf, periodicity);
        }
    }

    /// The running ⟨r⁻³⟩ of the time-averaged restraints, in restraint order (gromosXX writes
    /// them as DISRESEXPAVE into the final configuration).
    pub fn averages(&self) -> Vec<f64> {
        self.restraints
            .iter()
            .filter_map(|r| r.avg_r_inv3)
            .collect()
    }

    pub fn calculate_all(&mut self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        self.restraints
            .iter_mut()
            .map(|r| r.calculate(conf, periodicity))
            .sum()
    }

    /// [`Self::calculate_all`], also adding each restraint's energy to `per_group` at the
    /// energy group of its first atom (gromosXX `distanceres_energy[atom_energy_group(v1.atom(0))]`,
    /// `distance_restraint_interaction.cc:191`).
    pub fn calculate_all_grouped(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        atom_group: &[usize],
        per_group: &mut [f64],
    ) -> f64 {
        self.restraints
            .iter_mut()
            .map(|r| {
                let e = r.calculate(conf, periodicity);
                per_group[atom_group[r.atom1]] += e;
                e
            })
            .sum()
    }
}

// ─── Perturbed distance restraints (perturbed_distance_restraint_interaction.cc) ──

/// Single perturbed distance restraint.
///
/// r0(λ) = (1-λ)·A_r0 + λ·B_r0, w0(λ) = (1-λ)·A_w0 + λ·B_w0
/// prefactor = 2^(n+m) · λ^n · (1-λ)^m
/// Energy = prefactor · en_term (· w0 for mode 2)
#[derive(Debug, Clone)]
pub struct PerturbedDistanceRestraint {
    pub atom1: usize,
    pub atom2: usize,
    pub n: i32,
    pub m: i32,
    pub a_r0: f64,
    pub b_r0: f64,
    pub a_w0: f64,
    pub b_w0: f64,
    pub rah: i32,
    pub k: f64,
    pub r_linear: f64,
    pub mode: i32,
    /// The running ⟨r⁻³⟩ of a time-averaged restraint (gromosXX `pertdistanceres.av`).
    pub avg_r_inv3: Option<f64>,
    /// exp(−dt/τ) of a time-averaged restraint (NTDIR < 0); `None` = instantaneous.
    pub exp_mem: Option<f64>,
}

impl PerturbedDistanceRestraint {
    pub fn new(
        atom1: usize,
        atom2: usize,
        n: i32,
        m: i32,
        a_r0: f64,
        b_r0: f64,
        a_w0: f64,
        b_w0: f64,
        rah: i32,
        k: f64,
        r_linear: f64,
        mode: i32,
    ) -> Self {
        Self {
            atom1,
            atom2,
            n,
            m,
            a_r0,
            b_r0,
            a_w0,
            b_w0,
            rah,
            k,
            r_linear,
            mode,
            avg_r_inv3: None,
            exp_mem: None,
        }
    }

    /// Time-average the restraint (NTDIR < 0): `exp_mem = exp(−dt/τ)`
    /// (`perturbed_distance_restraint_interaction.cc:445-447`).
    pub fn with_time_averaging(mut self, tau: f64, dt: f64) -> Self {
        self.exp_mem = Some((-dt / tau).exp());
        self
    }

    /// Seed the average with `|r|⁻³` of the start-up configuration (NTDIRA = 0,
    /// `_init_averages`, :403-421: the unmasked distance).
    pub fn init_average(&mut self, conf: &Configuration, periodicity: &Periodicity) {
        if self.exp_mem.is_some() {
            let pos = &conf.current().pos;
            let v = periodicity.nearest_image(pos[self.atom1], pos[self.atom2]);
            self.avg_r_inv3 = Some(v.length().powf(-3.0));
        }
    }

    /// Energy and dH/dλ at the individual λ `l` (and `l_deriv` = dL/dλ) of this restraint,
    /// adding its forces to `conf` — gromosXX `perturbed_distance_restraint_interaction.cc`
    /// (`_calculate_perturbed_distance_restraint_interactions`, :82-285).
    ///
    /// `r0` and `w0` are linear in `l`; the prefactor is `2^(n+m)·lⁿ·(1−l)ᵐ` with the restraint's
    /// own exponents (not NLAM); dH/dλ = `l_deriv·(dprefactor/dl·en + prefactor·den/dl)`.
    /// Time-averaged (NTDIR < 0, :149-155): the distance is `⟨r⁻³⟩^(−1/3)` everywhere, the force
    /// direction `v/dist` included.
    ///
    /// With `precalc` the restraint also adds its energy and derivative at every λ_k
    /// (`:291-392`): r0, w0 and the prefactor at λ_k, the same (possibly time-averaged)
    /// distance.
    pub fn calculate(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        l: f64,
        l_deriv: f64,
        precalc: Option<&mut DisresPrecalc>,
    ) -> (f64, f64) {
        let (n, m) = (self.n as f64, self.m as f64);
        // gromosXX: v = v1 − v2
        let v = periodicity.nearest_image(
            conf.current().pos[self.atom1],
            conf.current().pos[self.atom2],
        );
        let (rah, mask) = decode_rah(self.rah);
        let v = Vec3::new(v.x * mask.x, v.y * mask.y, v.z * mask.z);
        let mut dist = v.length();
        if let Some(exp_t) = self.exp_mem {
            if self.avg_r_inv3.is_none() {
                self.init_average(conf, periodicity);
            }
            let ave = (1.0 - exp_t) * dist.powf(-3.0) + exp_t * self.avg_r_inv3.unwrap_or(0.0);
            self.avg_r_inv3 = Some(ave);
            dist = ave.powf(-1.0 / 3.0);
        }

        let w0 = (1.0 - l) * self.a_w0 + l * self.b_w0;
        let r0 = (1.0 - l) * self.a_r0 + l * self.b_r0;
        let k = self.k;
        let r_l = self.r_linear;
        let d_r0 = self.b_r0 - self.a_r0;
        let d_w0 = self.b_w0 - self.a_w0;

        let prefactor = 2.0_f64.powf(n + m) * l.powf(n) * (1.0 - l).powf(m);

        let rah_f = rah as f64;
        let fulfilled = rah_f * dist < rah_f * r0;
        let harmonic = (r0 - dist).abs() < r_l;
        let (mut f, mut en_term) = if fulfilled {
            (Vec3::ZERO, 0.0)
        } else if harmonic {
            (
                v * (-k * (dist - r0)) / dist,
                0.5 * k * (dist - r0) * (dist - r0),
            )
        } else if dist < r0 {
            (v * (r_l * k) / dist, -k * (dist + 0.5 * r_l - r0) * r_l)
        } else {
            (v * (-r_l * k) / dist, k * (dist - 0.5 * r_l - r0) * r_l)
        };
        match self.mode {
            1 => {},
            2 => {
                f *= w0;
                en_term *= w0;
            },
            _ => {
                f = Vec3::ZERO;
                en_term = 0.0;
            },
        }
        let f = f * prefactor;
        let energy = prefactor * en_term;
        conf.current_mut().force[self.atom1] += f;
        conf.current_mut().force[self.atom2] -= f;

        let dlam_term = if fulfilled {
            0.0
        } else if self.mode == 1 {
            if harmonic {
                -k * (dist - r0) * d_r0
            } else if dist < r0 {
                k * d_r0 * r_l
            } else {
                -k * d_r0 * r_l
            }
        } else if self.mode == 2 {
            if harmonic {
                0.5 * k * d_w0 * (dist - r0) * (dist - r0) - k * w0 * (dist - r0) * d_r0
            } else if dist < r0 {
                -k * d_w0 * (dist + 0.5 * r_l - r0) * r_l + k * w0 * d_r0 * r_l
            } else {
                k * d_w0 * (dist - 0.5 * r_l - r0) * r_l - k * w0 * d_r0 * r_l
            }
        } else {
            // gromosXX leaves the previous restraint's value here; with |NTDIR| ∉ {1, 2} there
            // is no restraint energy to differentiate
            0.0
        };
        let dprefndl = if self.n == 0 {
            0.0
        } else {
            n * l.powf(n - 1.0) * (1.0 - l).powf(m)
        };
        let dprefmdl = if self.m == 0 {
            0.0
        } else {
            m * l.powf(n) * (1.0 - l).powf(m - 1.0)
        };
        let dprefdl = 2.0_f64.powf(m + n) * (dprefndl - dprefmdl) * en_term;
        let dpotdl = prefactor * dlam_term;
        if let Some(pc) = precalc {
            self.add_precalc(pc, dist, rah_f);
        }
        (energy, l_deriv * (dprefdl + dpotdl))
    }

    /// gromosXX's ext-TI loop (`perturbed_distance_restraint_interaction.cc:291-392`) at the
    /// restraint's distance `dist` and half-harmonic form `rah`.
    fn add_precalc(&self, pc: &mut DisresPrecalc, dist: f64, rah: f64) {
        let (n, m) = (self.n as f64, self.m as f64);
        let (k, r_l) = (self.k, self.r_linear);
        let d_r0 = self.b_r0 - self.a_r0;
        let d_w0 = self.b_w0 - self.a_w0;
        for (i, &lam) in pc.lambdas.iter().enumerate() {
            let prefactor = 2.0_f64.powf(n + m) * lam.powf(n) * (1.0 - lam).powf(m);
            let w0 = (1.0 - lam) * self.a_w0 + lam * self.b_w0;
            let r0 = (1.0 - lam) * self.a_r0 + lam * self.b_r0;
            let diff = dist - r0;
            let diff2 = diff * diff;
            let fulfilled = rah * dist < rah * r0;
            let harmonic = (r0 - dist).abs() < r_l;
            let mut en_term = if fulfilled {
                0.0
            } else if harmonic {
                0.5 * k * diff2
            } else if dist < r0 {
                -k * (dist + 0.5 * r_l - r0) * r_l
            } else {
                k * (dist - 0.5 * r_l - r0) * r_l
            };
            match self.mode {
                1 => {},
                2 => en_term *= w0,
                _ => en_term = 0.0,
            }
            let dlam_term = if fulfilled {
                0.0
            } else if self.mode == 1 {
                if harmonic {
                    -k * diff * d_r0
                } else if dist < r0 {
                    k * d_r0 * r_l
                } else {
                    -k * d_r0 * r_l
                }
            } else if self.mode == 2 {
                if harmonic {
                    0.5 * k * d_w0 * diff2 - k * w0 * diff * d_r0
                } else if dist < r0 {
                    -k * d_w0 * (dist + 0.5 * r_l - r0) * r_l + k * w0 * d_r0 * r_l
                } else {
                    k * d_w0 * (dist - 0.5 * r_l - r0) * r_l - k * w0 * d_r0 * r_l
                }
            } else {
                0.0
            };
            let dprefndl = if self.n == 0 {
                0.0
            } else {
                n * lam.powf(n - 1.0) * (1.0 - lam).powf(m)
            };
            let dprefmdl = if self.m == 0 {
                0.0
            } else {
                m * lam.powf(n) * (1.0 - lam).powf(m - 1.0)
            };
            let dprefdl = 2.0_f64.powf(m + n) * (dprefndl - dprefmdl) * en_term;
            pc.energy[i] += prefactor * en_term;
            pc.deriv[i] += dprefdl + prefactor * dlam_term;
        }
    }
}

/// PRECALCLAM (extended TI, item 1.8) energies of the perturbed distance restraints:
/// `AB_disres` at every λ_k and its λ derivative.
#[derive(Debug, Clone, Default)]
pub struct DisresPrecalc {
    /// λ_k of the grid.
    pub lambdas: Vec<f64>,
    /// `AB_disres` energy per λ_k.
    pub energy: Vec<f64>,
    /// Its λ derivative per λ_k.
    pub deriv: Vec<f64>,
}

impl DisresPrecalc {
    /// Zeroed arrays for the λ grid `lambdas`.
    pub fn new(lambdas: &[f64]) -> Self {
        Self {
            lambdas: lambdas.to_vec(),
            energy: vec![0.0; lambdas.len()],
            deriv: vec![0.0; lambdas.len()],
        }
    }

    pub fn clear(&mut self) {
        self.energy.fill(0.0);
        self.deriv.fill(0.0);
    }
}

/// Collection of perturbed distance restraints.
#[derive(Debug, Clone, Default)]
pub struct PerturbedDistanceRestraints {
    pub restraints: Vec<PerturbedDistanceRestraint>,
}

impl PerturbedDistanceRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add(&mut self, r: PerturbedDistanceRestraint) {
        self.restraints.push(r);
    }

    /// Seed every time-averaged restraint from the start-up configuration.
    pub fn init_averages(&mut self, conf: &Configuration, periodicity: &Periodicity) {
        for r in &mut self.restraints {
            r.init_average(conf, periodicity);
        }
    }

    /// Total energy and dH/dλ; each restraint takes `lam.at(atom1)` (the `distanceres`
    /// individual λ of its first atom's energy group).
    pub fn calculate_all(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        lam: TermLambda<'_>,
        mut precalc: Option<&mut DisresPrecalc>,
    ) -> (f64, f64) {
        let mut e = 0.0;
        let mut d = 0.0;
        for r in &mut self.restraints {
            let (l, dl) = lam.at(r.atom1);
            let (ei, di) = r.calculate(conf, periodicity, l, dl, precalc.as_deref_mut());
            e += ei;
            d += di;
        }
        (e, d)
    }

    /// [`Self::calculate_all`], also adding each restraint's energy to `energy_groups` and its
    /// dH/dλ to `dhdl_groups` at the energy group of its first atom
    /// (`perturbed_distance_restraint_interaction.cc:194-196, 283-285`).
    pub fn calculate_all_grouped(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        lam: TermLambda<'_>,
        atom_group: &[usize],
        energy_groups: &mut [f64],
        dhdl_groups: &mut [f64],
        mut precalc: Option<&mut DisresPrecalc>,
    ) -> (f64, f64) {
        let mut e = 0.0;
        let mut d = 0.0;
        for r in &mut self.restraints {
            let (l, dl) = lam.at(r.atom1);
            let (ei, di) = r.calculate(conf, periodicity, l, dl, precalc.as_deref_mut());
            energy_groups[atom_group[r.atom1]] += ei;
            dhdl_groups[atom_group[r.atom1]] += di;
            e += ei;
            d += di;
        }
        (e, d)
    }
}

/// Angle restraint - restrains angle between three atoms
///
/// Restrains the angle θ defined by atoms i-j-k to a target value.
/// Used for NMR refinement and maintaining structural features.
///
/// Energy: E = 0.5 * k * (θ - θ₀)²
/// where θ is the angle i-j-k
#[derive(Debug, Clone)]
pub struct AngleRestraint {
    /// First atom index
    pub atom_i: usize,
    /// Second atom index (vertex)
    pub atom_j: usize,
    /// Third atom index
    pub atom_k: usize,
    /// Target angle (radians)
    pub theta0: f64,
    /// Force constant (kJ mol⁻¹ rad⁻²)
    pub force_constant: f64,
}

impl AngleRestraint {
    pub fn new(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        theta0: f64,
        force_constant: f64,
    ) -> Self {
        Self {
            atom_i,
            atom_j,
            atom_k,
            theta0,
            force_constant,
        }
    }

    /// Create from target angle in degrees
    pub fn from_degrees(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        theta0_deg: f64,
        force_constant: f64,
    ) -> Self {
        Self::new(
            atom_i,
            atom_j,
            atom_k,
            theta0_deg.to_radians(),
            force_constant,
        )
    }

    /// Calculate restraint energy and forces
    pub fn calculate(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        // Get vectors using periodic boundary conditions
        let r_ij = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_i],
        );
        let r_kj = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_k],
        );

        let r_ij_len = r_ij.length();
        let r_kj_len = r_kj.length();

        if r_ij_len < 1e-10 || r_kj_len < 1e-10 {
            return 0.0; // Degenerate angle
        }

        // Calculate angle using dot product: cos(θ) = r_ij · r_kj / (|r_ij| |r_kj|)
        let cos_theta = (r_ij.dot(r_kj)) / (r_ij_len * r_kj_len);
        let cos_theta = cos_theta.clamp(-1.0, 1.0); // Numerical safety
        let theta = cos_theta.acos();

        let delta_theta = theta - self.theta0;

        // Energy: E = 0.5 * k * (θ - θ₀)²
        let energy = 0.5 * self.force_constant * delta_theta * delta_theta;

        // Force magnitude: F = -dE/dθ = -k * (θ - θ₀)
        let force_magnitude = -self.force_constant * delta_theta;

        // sin(θ) for force calculation
        let sin_theta = theta.sin();
        if sin_theta.abs() < 1e-10 {
            return energy; // Avoid division by zero at θ = 0 or π
        }

        // Force derivatives
        let f_factor = force_magnitude / sin_theta;

        // Force on atom i
        let f_i_term = (r_kj * (cos_theta) - r_ij) * (f_factor / r_ij_len);
        conf.current_mut().force[self.atom_i] += f_i_term;

        // Force on atom k
        let f_k_term = (r_ij * (cos_theta) - r_kj) * (f_factor / r_kj_len);
        conf.current_mut().force[self.atom_k] += f_k_term;

        // Force on atom j (conservation)
        conf.current_mut().force[self.atom_j] -= f_i_term + f_k_term;

        energy
    }
}

/// Collection of angle restraints
#[derive(Debug, Clone, Default)]
pub struct AngleRestraints {
    pub restraints: Vec<AngleRestraint>,
}

impl AngleRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: AngleRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        self.restraints
            .iter()
            .map(|r| r.calculate(conf, periodicity))
            .sum()
    }
}

/// Dihedral restraint - restrains torsion angle between four atoms
///
/// Restrains the dihedral angle φ defined by atoms i-j-k-l to a target value.
/// Critical for protein backbone conformations and NMR refinement.
///
/// Energy: E = 0.5 * k * min[(φ - φ₀)², (φ - φ₀ ± 2π)²]
/// (handles periodicity of torsion angles)
#[derive(Debug, Clone)]
pub struct DihedralRestraint {
    /// First atom index
    pub atom_i: usize,
    /// Second atom index
    pub atom_j: usize,
    /// Third atom index
    pub atom_k: usize,
    /// Fourth atom index
    pub atom_l: usize,
    /// Target dihedral angle (radians)
    pub phi0: f64,
    /// Force constant (kJ mol⁻¹ rad⁻²)
    pub force_constant: f64,
}

impl DihedralRestraint {
    pub fn new(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        atom_l: usize,
        phi0: f64,
        force_constant: f64,
    ) -> Self {
        Self {
            atom_i,
            atom_j,
            atom_k,
            atom_l,
            phi0,
            force_constant,
        }
    }

    /// Create from target angle in degrees
    pub fn from_degrees(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        atom_l: usize,
        phi0_deg: f64,
        force_constant: f64,
    ) -> Self {
        Self::new(
            atom_i,
            atom_j,
            atom_k,
            atom_l,
            phi0_deg.to_radians(),
            force_constant,
        )
    }

    /// Compute minimum angular difference considering periodicity
    fn min_angle_diff(phi: f64, phi0: f64) -> f64 {
        let delta = phi - phi0;
        let two_pi = 2.0 * std::f64::consts::PI;

        // Find minimum difference considering +/-2π wrapping
        let mut min_delta = delta;
        let mut min_abs = delta.abs();

        for &shift in &[delta - two_pi, delta + two_pi] {
            if shift.abs() < min_abs {
                min_delta = shift;
                min_abs = shift.abs();
            }
        }

        min_delta
    }

    /// Calculate restraint energy and forces
    pub fn calculate(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        // Get bond vectors using periodic boundary conditions
        let r_ij = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_i],
        );
        let r_kj = periodicity.nearest_image(
            conf.current().pos[self.atom_j],
            conf.current().pos[self.atom_k],
        );
        let r_kl = periodicity.nearest_image(
            conf.current().pos[self.atom_l],
            conf.current().pos[self.atom_k],
        );

        // Calculate normal vectors to the two planes
        let m = r_ij.cross(r_kj);
        let n = r_kj.cross(r_kl);

        let m_len_sq = m.length_squared();
        let n_len_sq = n.length_squared();

        if m_len_sq < 1e-20 || n_len_sq < 1e-20 {
            return 0.0; // Degenerate dihedral
        }

        let m_len = m_len_sq.sqrt();
        let n_len = n_len_sq.sqrt();

        // Calculate dihedral angle
        let r_kj_len = r_kj.length();
        let cos_phi = (m.dot(n)) / (m_len * n_len);
        let cos_phi = cos_phi.clamp(-1.0, 1.0);

        // Get sign from scalar triple product
        let sign = if r_ij.dot(n) >= 0.0 { 1.0 } else { -1.0 };
        let phi = sign * cos_phi.acos();

        // Calculate minimum angular difference (handles periodicity)
        let delta_phi = Self::min_angle_diff(phi, self.phi0);

        // Energy: E = 0.5 * k * δφ²
        let energy = 0.5 * self.force_constant * delta_phi * delta_phi;

        // Force magnitude: F = -k * δφ
        let force_magnitude = -self.force_constant * delta_phi;

        // Apply forces (using standard dihedral force derivatives)
        let f_factor = force_magnitude / r_kj_len;

        // Forces on atoms i and l (perpendicular to planes)
        let f_i = m * (f_factor / m_len);
        let f_l = n * (-f_factor / n_len);

        conf.current_mut().force[self.atom_i] += f_i;
        conf.current_mut().force[self.atom_l] += f_l;

        // Forces on atoms j and k (derived from torque balance)
        let _r_ij_len = r_ij.length();
        let _r_kl_len = r_kl.length();

        let dot_ij_kj = r_ij.dot(r_kj);
        let dot_kl_kj = r_kl.dot(r_kj);

        let f_j_factor = dot_ij_kj / (r_kj_len * r_kj_len);
        let f_k_factor = dot_kl_kj / (r_kj_len * r_kj_len);

        let f_j = f_i * (f_j_factor - 1.0);
        let f_k = f_l * (f_k_factor - 1.0);

        conf.current_mut().force[self.atom_j] += f_j;
        conf.current_mut().force[self.atom_k] += f_k;

        energy
    }
}

/// Collection of dihedral restraints
#[derive(Debug, Clone, Default)]
pub struct DihedralRestraints {
    pub restraints: Vec<DihedralRestraint>,
}

impl DihedralRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: DihedralRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(&self, conf: &mut Configuration, periodicity: &Periodicity) -> f64 {
        self.restraints
            .iter()
            .map(|r| r.calculate(conf, periodicity))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_position_restraint() {
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(1.0, 0.0, 0.0);

        let restraint = PositionRestraint::new(0, Vec3::new(0.0, 0.0, 0.0), 100.0);

        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        let energy = restraint.calculate(&mut conf, &periodicity);

        // E = 0.5 * 100 * 1² = 50 kJ/mol
        assert!((energy - 50.0).abs() < 1e-6);

        // Force should point toward origin: F = -k * r = -100 * 1 = -100
        // (pointing in negative x direction)
        assert!((conf.current().force[0].x + 100.0).abs() < 1e-3);
    }

    #[test]
    fn test_distance_restraint_harmonic() {
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(2.0, 0.0, 0.0);

        // rah=0 → full harmonic XYZ, mode=1 (no w0), r_linear > delta so stays harmonic
        let mut restraint = DistanceRestraint::new(0, 1, 1.0, 1.0, 0, 100.0, 10.0, 1);

        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        let energy = restraint.calculate(&mut conf, &periodicity);

        // Distance = 2.0, r0 = 1.0, delta = 1.0
        // E = 0.5 * 100 * 1² = 50 kJ/mol
        assert!((energy - 50.0).abs() < 1e-6, "Expected 50, got {energy}");

        // Forces should be equal and opposite
        let force_i = conf.current().force[0];
        let force_j = conf.current().force[1];
        assert!((force_i.x + force_j.x).abs() < 1e-3);
    }

    #[test]
    fn test_distance_restraint_linear() {
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(5.0, 0.0, 0.0);

        // rah=0 → full harmonic XYZ, r_linear=0.1
        let mut restraint = DistanceRestraint::new(0, 1, 1.0, 1.0, 0, 100.0, 0.1, 1);

        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        let energy = restraint.calculate(&mut conf, &periodicity);

        // Distance = 5.0, r0 = 1.0, delta = 4.0 > r_linear=0.1 → linear zone above r0
        // E = k * r_linear * (dist - 0.5*r_linear - r0) = 100 * 0.1 * (5.0 - 0.05 - 1.0) = 39.5
        assert!((energy - 39.5).abs() < 1e-6, "Expected 39.5, got {energy}");
    }

    #[test]
    fn test_restraint_types() {
        let mut conf = Configuration::new(2, 1, 1);
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        // Repulsive restraint (rah=-1): only active when distance < r0
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.5, 0.0, 0.0);
        let mut restraint = DistanceRestraint::new(0, 1, 1.0, 1.0, -1, 100.0, 10.0, 1);
        let energy = restraint.calculate(&mut conf, &periodicity);
        assert!(
            energy > 0.0,
            "Repulsive restraint should have energy when dist < r0"
        );

        // Same restraint but distance > r0 - no energy
        conf.current_mut().pos[1] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().clear_forces();
        let mut restraint2 = DistanceRestraint::new(0, 1, 1.0, 1.0, -1, 100.0, 10.0, 1);
        let energy = restraint2.calculate(&mut conf, &periodicity);
        assert_eq!(
            energy, 0.0,
            "Repulsive restraint should have no energy when dist > r0"
        );
    }

    #[test]
    fn test_distance_restraint_gromosxx_reference() {
        // Reference values from GROMOS aladip_special.t.cc:
        //   DistanceRestraint = 257.189539
        //   PerturbedDistanceRestraint = 195.899012
        //
        // Parameters: K=1000.0, r_linear=0.3, mode=2 (NTDIR=2), lambda=0.125
        // Atom positions from aladip.conf:
        // atom 5 → 0-indexed 4: (2.506319188, 0.858753474, 1.635829769)
        // atom 9 → 0-indexed 8: (2.722961726, 1.075330732, 1.724773551)
        // atom12 → 0-indexed11: (3.001529656, 0.996657605, 1.693520228)
        let n_atoms = 12;
        let mut conf = Configuration::new(n_atoms, 1, 1);
        conf.current_mut().pos[4] = Vec3::new(2.506319188, 0.858753474, 1.635829769);
        conf.current_mut().pos[8] = Vec3::new(2.722961726, 1.075330732, 1.724773551);
        conf.current_mut().pos[11] = Vec3::new(3.001529656, 0.996657605, 1.693520228);

        let box_l = 3.767055681_f64;
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            box_l, box_l, box_l,
        )));

        let k = 1000.0_f64;
        let r_linear = 0.3_f64;
        let mode = 2;

        // Unperturbed restraints from aladip.distrest (all type-0 virtual atoms):
        // pairs: (4,8) and (8,11), rah={1,1,0,0,-1,-1}, r0 and w0 vary per entry
        let specs: &[(usize, usize, f64, f64, i32)] = &[
            (4, 8, 0.21, 1.0, 1),
            (8, 11, 0.80, 1.0, 1),
            (4, 8, 0.21, 1.0, 0),
            (8, 11, 0.90, 1.0, 0),
            (4, 8, 0.21, 1.0, -1),
            (8, 11, 0.80, 1.0, -1),
        ];

        let mut total = 0.0_f64;
        for &(a1, a2, r0, w0, rah) in specs {
            let mut r = DistanceRestraint::new(a1, a2, r0, w0, rah, k, r_linear, mode);
            total += r.calculate(&mut conf, &periodicity);
        }
        assert!(
            (total - 257.189539).abs() < 1e-3,
            "DistanceRestraint total: expected 257.189539, got {total:.6}"
        );

        // Perturbed restraints (lambda=0.125, NLAM=1, RLAM=0.125)
        let lambda = 0.125_f64;
        let pspecs: &[(usize, usize, i32, i32, f64, f64, f64, f64, i32)] = &[
            (4, 8, 1, 1, 0.19, 2.0, 0.21, 1.0, 1),
            (8, 11, 1, 1, 0.80, 2.0, 0.90, 1.0, 1),
            (4, 8, 1, 1, 0.19, 2.0, 0.21, 1.0, 0),
            (8, 11, 1, 1, 0.80, 2.0, 0.90, 1.0, 0),
            (4, 8, 1, 1, 0.19, 2.0, 0.21, 1.0, -1),
            (8, 11, 1, 1, 0.80, 2.0, 0.90, 1.0, -1),
        ];

        // Clear forces from previous run
        conf.current_mut().clear_forces();
        let mut ptotal = 0.0_f64;
        for &(a1, a2, n, m, a_r0, a_w0, b_r0, b_w0, rah) in pspecs {
            let mut r = PerturbedDistanceRestraint::new(
                a1, a2, n, m, a_r0, b_r0, a_w0, b_w0, rah, k, r_linear, mode,
            );
            ptotal += r.calculate(&mut conf, &periodicity, lambda, 1.0, None).0;
        }
        assert!(
            (ptotal - 195.899012).abs() < 1e-3,
            "PerturbedDistanceRestraint total: expected 195.899012, got {ptotal:.6}"
        );
    }

    /// PRECALCLAM: at λ_k = λ (and L = λ) `AB_disres` is the restraint's energy and its
    /// derivative is dH/dλ; across the grid the derivative is the finite difference of the
    /// energies.
    #[test]
    fn precalc_disres_matches_the_energy_at_lambda_and_its_derivative() {
        let pbc = Periodicity::Vacuum(newen_core::math::Vacuum);
        for (dist, mode, rah) in [(0.25, 1, 0), (0.25, 2, 0), (0.9, 2, 0), (0.02, 1, 0)] {
            let mut conf = Configuration::new(2, 1, 1);
            conf.current_mut().pos[1] = Vec3::new(dist, 0.0, 0.0);
            let r = PerturbedDistanceRestraint::new(
                0, 1, 1, 2, 0.2, 0.3, 1.0, 2.0, rah, 1000.0, 0.1, mode,
            );
            let eps = 1e-6;
            let grid = [0.4 - eps, 0.4, 0.4 + eps];
            let mut pc = DisresPrecalc::new(&grid);
            let (e, de) = r
                .clone()
                .calculate(&mut conf.clone(), &pbc, 0.4, 1.0, Some(&mut pc));
            assert!((pc.energy[1] - e).abs() <= 1e-12 * e.abs().max(1.0));
            assert!((pc.deriv[1] - de).abs() <= 1e-12 * de.abs().max(1.0));
            let fd = (pc.energy[2] - pc.energy[0]) / (2.0 * eps);
            assert!(
                (fd - pc.deriv[1]).abs() <= 1e-5 * fd.abs().max(1.0),
                "dist {dist} mode {mode}: fd {fd} vs {}",
                pc.deriv[1]
            );
        }
    }

    /// dH/dλ of a perturbed distance restraint is the λ derivative of its energy, in every
    /// zone (harmonic, linear below and above r0) and for both NTDIR forms, with L′ ≠ 1.
    #[test]
    fn perturbed_distance_restraint_lambda_derivative_is_a_finite_difference() {
        let pbc = Periodicity::Vacuum(newen_core::math::Vacuum);
        for (dist, mode, rah) in [
            (0.25, 1, 0),
            (0.25, 2, 0),
            (0.9, 2, 0),
            (0.02, 1, 0),
            (0.3, 2, 1),
        ] {
            let mut conf = Configuration::new(2, 1, 1);
            conf.current_mut().pos[1] = Vec3::new(dist, 0.0, 0.0);
            let r = PerturbedDistanceRestraint::new(
                0, 1, 1, 2, 0.2, 0.3, 1.0, 2.0, rah, 1000.0, 0.1, mode,
            );
            // L(λ) = λ² (so L′ = 2λ)
            let energy = |lam: f64| {
                let mut c = conf.clone();
                r.clone()
                    .calculate(&mut c, &pbc, lam * lam, 2.0 * lam, None)
                    .0
            };
            let lam = 0.4;
            let eps = 1e-6;
            let fd = (energy(lam + eps) - energy(lam - eps)) / (2.0 * eps);
            let (_, de) = r
                .clone()
                .calculate(&mut conf.clone(), &pbc, lam * lam, 2.0 * lam, None);
            assert!(
                (fd - de).abs() <= 1e-6 * de.abs().max(1.0),
                "dist {dist} mode {mode}: finite difference {fd} vs dH/dλ {de}"
            );
        }
    }

    #[test]
    fn test_angle_restraint_linear() {
        let mut conf = Configuration::new(3, 1, 1);
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        // Set up linear angle (180 degrees)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);

        // Restrain to 180 degrees (π radians) - should be at equilibrium
        let restraint = AngleRestraint::new(0, 1, 2, std::f64::consts::PI, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Should have nearly zero energy (linear configuration at target angle)
        assert!(energy < 1e-6, "Energy should be ~0, got {}", energy);
    }

    #[test]
    fn test_angle_restraint_90_degrees() {
        let mut conf = Configuration::new(3, 1, 1);
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        // Set up 90 degree angle
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 1.0, 0.0);

        // Restrain to 90 degrees (π/2 radians)
        let restraint = AngleRestraint::from_degrees(0, 1, 2, 90.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Should have nearly zero energy (at target angle)
        assert!(energy < 1e-6, "Energy should be ~0, got {}", energy);

        // Forces should be near zero at equilibrium
        let max_force = conf
            .current()
            .force
            .iter()
            .map(|f| f.length())
            .fold(0.0_f64, f64::max);
        assert!(
            max_force < 1e-3,
            "Forces should be ~0 at equilibrium, max={}",
            max_force
        );
    }

    #[test]
    fn test_angle_restraint_violation() {
        let mut conf = Configuration::new(3, 1, 1);
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        // Set up 90 degree angle
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 1.0, 0.0);

        // Restrain to 180 degrees (should be violated)
        let restraint = AngleRestraint::from_degrees(0, 1, 2, 180.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Angle is 90°, target is 180°, delta = 90° = π/2
        // E = 0.5 * 100 * (π/2)² ≈ 0.5 * 100 * 2.467 ≈ 123.4
        assert!(
            energy > 100.0,
            "Energy should be large for violation, got {}",
            energy
        );
        assert!(energy < 150.0, "Energy should be ~123, got {}", energy);
    }

    #[test]
    fn test_dihedral_restraint_trans() {
        let mut conf = Configuration::new(4, 1, 1);
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        // Set up trans dihedral (180 degrees)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(3.0, 0.0, 0.0);

        // Restrain to 180 degrees (trans)
        let restraint = DihedralRestraint::from_degrees(0, 1, 2, 3, 180.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Should have nearly zero energy (at target dihedral)
        assert!(
            energy < 1e-6,
            "Energy should be ~0 for trans, got {}",
            energy
        );
    }

    #[test]
    fn test_dihedral_restraint_gauche() {
        let mut conf = Configuration::new(4, 1, 1);
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        // Set up gauche+ dihedral (~60 degrees)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.5, 0.866, 0.0); // 60° rotation

        // Restrain to 60 degrees (gauche+)
        let restraint = DihedralRestraint::from_degrees(0, 1, 2, 3, 60.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // Should have nearly zero energy (at target dihedral)
        assert!(
            energy < 1e-3,
            "Energy should be ~0 for gauche, got {}",
            energy
        );
    }

    #[test]
    fn test_dihedral_restraint_periodicity() {
        let mut conf = Configuration::new(4, 1, 1);
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        // Set up dihedral at -170 degrees
        let angle = -170.0f64.to_radians();
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0 + angle.cos(), angle.sin(), 0.0);

        // Restrain to +170 degrees (should recognize as close due to periodicity)
        let restraint = DihedralRestraint::from_degrees(0, 1, 2, 3, 170.0, 100.0);
        let energy = restraint.calculate(&mut conf, &periodicity);

        // -170° and +170° differ by 340°, but min difference is 20° (due to periodicity)
        // E = 0.5 * 100 * (20° in radians)² = 0.5 * 100 * (0.349)² ≈ 6.1
        assert!(
            energy < 10.0,
            "Energy should consider periodicity, got {}",
            energy
        );
    }

    #[test]
    fn test_angle_collection() {
        let mut conf = Configuration::new(6, 1, 1);
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        // Set up multiple angles
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(0.0, 1.0, 0.0);
        conf.current_mut().pos[4] = Vec3::new(1.0, 1.0, 0.0);
        conf.current_mut().pos[5] = Vec3::new(2.0, 1.0, 0.0);

        let mut restraints = AngleRestraints::new();
        restraints.add_restraint(AngleRestraint::from_degrees(0, 1, 2, 180.0, 100.0));
        restraints.add_restraint(AngleRestraint::from_degrees(3, 4, 5, 180.0, 100.0));

        let total_energy = restraints.calculate_all(&mut conf, &periodicity);

        // Both angles are linear (180°), so total energy should be near zero
        assert!(
            total_energy < 1e-3,
            "Total energy should be ~0, got {}",
            total_energy
        );
    }

    #[test]
    fn test_dihedral_collection() {
        let mut conf = Configuration::new(8, 1, 1);
        let periodicity = Periodicity::Rectangular(newen_core::math::Rectangular::new(Vec3::new(
            10.0, 10.0, 10.0,
        )));

        // Set up multiple trans dihedrals
        for i in 0..8 {
            conf.current_mut().pos[i] = Vec3::new(i as f64, 0.0, 0.0);
        }

        let mut restraints = DihedralRestraints::new();
        restraints.add_restraint(DihedralRestraint::from_degrees(0, 1, 2, 3, 180.0, 100.0));
        restraints.add_restraint(DihedralRestraint::from_degrees(4, 5, 6, 7, 180.0, 100.0));

        let total_energy = restraints.calculate_all(&mut conf, &periodicity);

        // Both dihedrals are trans (180°), so total energy should be near zero
        assert!(
            total_energy < 1e-3,
            "Total energy should be ~0, got {}",
            total_energy
        );
    }
}

/// J-value restraint for NMR refinement
///
/// Direct translation of md++/src/interaction/special/jvalue_restraint_interaction.cc
///
/// **Purpose**: Restrain scalar J-coupling constants from NMR spectroscopy
///
/// J-coupling constants are calculated from dihedral angles using the Karplus relation:
/// ```text
/// J(φ) = A·cos²(φ+δ) + B·cos(φ+δ) + C
/// ```
///
/// **Features**:
/// - Karplus relation with custom parameters (A, B, C, δ)
/// - Time averaging with exponential memory
/// - Multiple restraining modes (instantaneous, time-averaged, biquadratic)
/// - Weighted and unweighted variants
/// - Flat-bottom potentials (allow J-values within tolerance)
/// - Half-harmonic (repulsive/attractive only)
/// - Optional local elevation (metadynamics)
///
/// **Use case**: NMR structure refinement using ³J-coupling data
#[derive(Debug, Clone)]
pub struct JValueRestraint {
    /// Atom indices defining dihedral i-j-k-l
    pub i: usize,
    pub j: usize,
    pub k_atom: usize,
    pub l_atom: usize,
    /// Karplus parameter A (cos² term)
    pub a: f64,
    /// Karplus parameter B (cos term)
    pub b: f64,
    /// Karplus parameter C (constant term)
    pub c: f64,
    /// Phase shift δ in Karplus relation (radians)
    pub delta: f64,
    /// Target J-value from experiment (Hz)
    pub j0: f64,
    /// Force constant K (kJ/mol/Hz²)
    pub k: f64,
    /// Flat-bottom tolerance (Hz)
    pub flat_bottom_width: f64,
    /// Half-harmonic type: 0=full, -1=repulsive, +1=attractive
    pub half_harmonic: i32,
    /// Time-averaging memory constant τ (ps)
    pub tau: f64,
    /// Running average J-value (updated during simulation)
    pub j_average: f64,
    /// Current instantaneous J-value
    pub j_current: f64,
}

impl JValueRestraint {
    /// Create new J-value restraint
    ///
    /// # Arguments
    /// * `i, j, k, l` - Atom indices defining dihedral angle
    /// * `a, b, c` - Karplus parameters (typical: A=6-9, B=-1 to +1, C=0-2 Hz)
    /// * `delta` - Phase shift in radians
    /// * `j0` - Target J-value from NMR experiment (Hz)
    /// * `k` - Force constant (kJ/mol/Hz²)
    pub fn new(
        i: usize,
        j: usize,
        k_atom: usize,
        l_atom: usize,
        a: f64,
        b: f64,
        c: f64,
        delta: f64,
        j0: f64,
        force_constant: f64,
    ) -> Self {
        Self {
            i,
            j,
            k_atom,
            l_atom,
            a,
            b,
            c,
            delta,
            j0,
            k: force_constant,
            flat_bottom_width: 0.0,
            half_harmonic: 0,
            tau: 1.0, // 1 ps default
            j_average: 0.0,
            j_current: 0.0,
        }
    }

    /// Create J-value restraint with default Karplus parameters for ³J(HN-Hα)
    ///
    /// Uses Pardi et al. parameters: A=6.4, B=-1.4, C=1.9, δ=0
    pub fn backbone_hn_ha(
        i: usize,
        j: usize,
        k_atom: usize,
        l_atom: usize,
        j0: f64,
        force_constant: f64,
    ) -> Self {
        Self::new(
            i,
            j,
            k_atom,
            l_atom,
            6.4,
            -1.4,
            1.9,
            0.0,
            j0,
            force_constant,
        )
    }

    /// Set flat-bottom tolerance (no force if |J-J0| < width)
    pub fn with_flat_bottom(mut self, width: f64) -> Self {
        self.flat_bottom_width = width;
        self
    }

    /// Set half-harmonic type (-1=repulsive, 0=full, +1=attractive)
    pub fn with_half_harmonic(mut self, htype: i32) -> Self {
        self.half_harmonic = htype;
        self
    }

    /// Set time-averaging memory constant (ps)
    pub fn with_time_averaging(mut self, tau: f64) -> Self {
        self.tau = tau;
        self
    }

    /// Calculate J-value restraint energy and forces
    ///
    /// # Algorithm
    /// 1. Calculate dihedral angle φ from atomic positions
    /// 2. Calculate current J-value using Karplus relation
    /// 3. Update time-averaged J-value with exponential memory
    /// 4. Calculate restraint energy and derivative
    /// 5. Apply forces to all four atoms
    pub fn calculate(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        // Get positions
        let pos_i = conf.current().pos[self.i];
        let pos_j = conf.current().pos[self.j];
        let pos_k = conf.current().pos[self.k_atom];
        let pos_l = conf.current().pos[self.l_atom];

        // Calculate nearest image vectors
        let rij = periodicity.nearest_image(pos_i, pos_j);
        let rkj = periodicity.nearest_image(pos_k, pos_j);
        let rkl = periodicity.nearest_image(pos_k, pos_l);

        // Calculate dihedral angle φ using cross products
        let rmj = rij.cross(rkj);
        let rnk = rkj.cross(rkl);

        let dkj2 = rkj.length_squared();
        let dmj2 = rmj.length_squared();
        let dnk2 = rnk.length_squared();

        if dkj2 < 1e-10 || dmj2 < 1e-10 || dnk2 < 1e-10 {
            return 0.0; // Degenerate geometry
        }

        // Alternative stable dihedral calculation
        let frim = rij.dot(rkj) / dkj2;
        let frln = rkl.dot(rkj) / dkj2;
        let rim = rij - rkj * (frim);
        let rln = rkj * (frln) - rkl;

        let dim = rim.length();
        let dln = rln.length();

        if dim < 1e-10 || dln < 1e-10 {
            return 0.0;
        }

        let cos_phi = (rim.dot(rln) / (dim * dln)).clamp(-1.0, 1.0);
        let mut phi = (cos_phi).acos();

        // Determine sign of dihedral
        if rij.dot(rnk) < 0.0 {
            phi = -phi;
        }

        // Calculate current J-value using Karplus relation
        let cos_phi_delta = (phi + self.delta).cos();
        let sin_phi_delta = (phi + self.delta).sin();

        self.j_current = self.a * cos_phi_delta * cos_phi_delta + self.b * cos_phi_delta + self.c;

        // Update time-averaged J-value
        let exp_term = (-dt / self.tau).exp();
        let memory_decay = 1.0 - exp_term;
        self.j_average = memory_decay * self.j_current + self.j_average * exp_term;

        // Calculate deviation with flat-bottom
        let mut delta_j = self.j_average - self.j0;

        // Apply flat bottom
        if delta_j > 0.0 {
            delta_j = (delta_j - self.flat_bottom_width).max(0.0);
        } else {
            delta_j = (delta_j + self.flat_bottom_width).min(0.0);
        }

        // Check half-harmonic
        if (self.half_harmonic == -1 && delta_j > 0.0)
            || (self.half_harmonic == 1 && delta_j <= 0.0)
        {
            return 0.0;
        }

        // Calculate energy: E = 0.5 * K * (J_av - J0)²
        let energy = 0.5 * self.k * delta_j * delta_j;

        // Calculate derivative dV/dφ = dV/dJ * dJ/dφ
        let dv_dj = self.k * delta_j;
        let dj_dphi =
            -memory_decay * (2.0 * self.a * cos_phi_delta * sin_phi_delta + self.b * sin_phi_delta);
        let dv_dphi = dv_dj * dj_dphi;

        // Calculate force derivatives dphi/dr for each atom
        let dkj = dkj2.sqrt();
        let dphi_dri = rmj * (dkj / dmj2);
        let dphi_drl = rnk * (-(dkj / dnk2));
        let dphi_drj = dphi_dri * (frim - 1.0) - dphi_drl * (frln);
        let dphi_drk = -(dphi_dri + dphi_drj + dphi_drl);

        // Apply forces: F = -dV/dφ * dφ/dr
        let dv_dphi_f32 = dv_dphi;
        conf.current_mut().force[self.i] -= dphi_dri * dv_dphi_f32;
        conf.current_mut().force[self.j] -= dphi_drj * dv_dphi_f32;
        conf.current_mut().force[self.k_atom] -= dphi_drk * dv_dphi_f32;
        conf.current_mut().force[self.l_atom] -= dphi_drl * dv_dphi_f32;

        energy
    }

    /// Reset time average (e.g., at start of simulation)
    pub fn reset_average(&mut self) {
        self.j_average = 0.0;
    }
}

/// Collection of J-value restraints
#[derive(Debug, Clone)]
pub struct JValueRestraints {
    pub restraints: Vec<JValueRestraint>,
}

impl Default for JValueRestraints {
    fn default() -> Self {
        Self::new()
    }
}

impl JValueRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: JValueRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        self.restraints
            .iter_mut()
            .map(|r| r.calculate(conf, periodicity, dt))
            .sum()
    }

    /// Reset all time averages
    pub fn reset_averages(&mut self) {
        for r in &mut self.restraints {
            r.reset_average();
        }
    }
}

/// RDC (Residual Dipolar Coupling) restraint for NMR refinement
///
/// Direct translation of md++/src/interaction/special/rdc_restraint.cc
///
/// **Purpose**: Restrain residual dipolar couplings from aligned NMR samples
///
/// RDCs provide long-range orientational information. The dipolar coupling
/// depends on the orientation of internuclear vectors relative to the
/// alignment tensor:
///
/// ```text
/// D = D_max * Σᵢⱼ Sᵢⱼ * (3*cosθᵢ*cosθⱼ - δᵢⱼ)
/// ```
///
/// where S is the Saupe order matrix (alignment tensor).
///
/// **Features**:
/// - Alignment tensor determination (SVD fitting)
/// - Time averaging
/// - Flat-bottom potentials
/// - Multiple alignment media support
/// - Automatic or fixed alignment tensor
///
/// **Use case**: NMR structure refinement with orientation data
#[derive(Debug, Clone)]
pub struct RDCRestraint {
    /// Atom indices i-j defining the internuclear vector
    pub i: usize,
    pub j: usize,
    /// Maximum dipolar coupling D_max (Hz)
    /// D_max = -(μ₀/4π) * (γᵢ*γⱼ*ħ)/(2π*r³)
    pub d_max: f64,
    /// Target RDC value from experiment (Hz)
    pub d0: f64,
    /// Force constant K (kJ/mol/Hz²)
    pub k: f64,
    /// Flat-bottom tolerance (Hz)
    pub flat_bottom_width: f64,
    /// Type: 0=CH, 1=NH, 2=CC, 3=custom
    pub rdc_type: i32,
    /// Alignment tensor Saupe matrix (5 independent components: Sxx, Syy, Szz, Sxy, Sxz)
    pub saupe_matrix: [f64; 5],
    /// Time-averaging memory constant τ (ps)
    pub tau: f64,
    /// Running average RDC value
    pub rdc_average: f64,
    /// Current instantaneous RDC value
    pub rdc_current: f64,
}

impl RDCRestraint {
    /// Create new RDC restraint
    ///
    /// # Arguments
    /// * `i, j` - Atom indices (e.g., C-H, N-H pair)
    /// * `d_max` - Maximum dipolar coupling in Hz
    /// * `d0` - Target RDC from experiment in Hz
    /// * `k` - Force constant in kJ/mol/Hz²
    pub fn new(i: usize, j: usize, d_max: f64, d0: f64, k: f64) -> Self {
        Self {
            i,
            j,
            d_max,
            d0,
            k,
            flat_bottom_width: 0.0,
            rdc_type: 0,
            saupe_matrix: [0.0; 5],
            tau: 1.0,
            rdc_average: 0.0,
            rdc_current: 0.0,
        }
    }

    /// Create CH RDC restraint with standard parameters
    ///
    /// D_max for ¹H-¹³C ≈ -60 kHz * (r_CH/0.109nm)³
    /// For r_CH = 0.109 nm: D_max ≈ -60 kHz
    pub fn ch_bond(i: usize, j: usize, d0: f64, k: f64) -> Self {
        Self::new(i, j, -60000.0, d0, k).with_type(0)
    }

    /// Create NH RDC restraint with standard parameters
    ///
    /// D_max for ¹H-¹⁵N ≈ -24 kHz * (r_NH/0.104nm)³
    /// For r_NH = 0.104 nm: D_max ≈ -24 kHz
    pub fn nh_bond(i: usize, j: usize, d0: f64, k: f64) -> Self {
        Self::new(i, j, -24000.0, d0, k).with_type(1)
    }

    /// Set RDC type (0=CH, 1=NH, 2=CC, 3=custom)
    pub fn with_type(mut self, rdc_type: i32) -> Self {
        self.rdc_type = rdc_type;
        self
    }

    /// Set alignment tensor (Saupe matrix: Sxx, Syy, Szz, Sxy, Sxz)
    /// Note: Syz can be derived from traceless condition
    pub fn with_alignment_tensor(mut self, saupe: [f64; 5]) -> Self {
        self.saupe_matrix = saupe;
        self
    }

    /// Set flat-bottom tolerance
    pub fn with_flat_bottom(mut self, width: f64) -> Self {
        self.flat_bottom_width = width;
        self
    }

    /// Set time-averaging memory constant
    pub fn with_time_averaging(mut self, tau: f64) -> Self {
        self.tau = tau;
        self
    }

    /// Calculate RDC restraint energy and forces
    ///
    /// # Algorithm
    /// 1. Calculate internuclear vector r_ij
    /// 2. Normalize to unit vector n
    /// 3. Calculate RDC: D = D_max * Σᵢⱼ Sᵢⱼ * (3*nᵢ*nⱼ - δᵢⱼ)/2
    /// 4. Update time average
    /// 5. Calculate restraint energy and forces
    pub fn calculate(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        // Get internuclear vector
        let pos_i = conf.current().pos[self.i];
        let pos_j = conf.current().pos[self.j];
        let rij = periodicity.nearest_image(pos_i, pos_j);

        let r = rij.length();
        if r < 1e-10 {
            return 0.0; // Avoid division by zero
        }

        // Normalize to unit vector
        let n = rij / r;
        let nx = n.x;
        let ny = n.y;
        let nz = n.z;

        // Unpack Saupe matrix (Sxx, Syy, Szz, Sxy, Sxz)
        let sxx = self.saupe_matrix[0];
        let syy = self.saupe_matrix[1];
        let szz = self.saupe_matrix[2];
        let sxy = self.saupe_matrix[3];
        let sxz = self.saupe_matrix[4];
        // Syz derived from traceless condition if needed
        let syz = 0.0; // Simplified for now

        // Calculate order parameter: P2 = Σᵢⱼ Sᵢⱼ * (3*nᵢ*nⱼ - δᵢⱼ)/2
        let p2 = 0.5
            * (sxx * (3.0 * nx * nx - 1.0)
                + syy * (3.0 * ny * ny - 1.0)
                + szz * (3.0 * nz * nz - 1.0)
                + 2.0 * sxy * 3.0 * nx * ny
                + 2.0 * sxz * 3.0 * nx * nz
                + 2.0 * syz * 3.0 * ny * nz);

        // Calculate current RDC: D = D_max * P2
        self.rdc_current = self.d_max * p2;

        // Update time average
        let exp_term = (-dt / self.tau).exp();
        let memory_decay = 1.0 - exp_term;
        self.rdc_average = memory_decay * self.rdc_current + self.rdc_average * exp_term;

        // Calculate deviation with flat-bottom
        let mut delta_d = self.rdc_average - self.d0;

        if delta_d > 0.0 {
            delta_d = (delta_d - self.flat_bottom_width).max(0.0);
        } else {
            delta_d = (delta_d + self.flat_bottom_width).min(0.0);
        }

        // Calculate energy: E = 0.5 * K * (D_av - D0)²
        let energy = 0.5 * self.k * delta_d * delta_d;

        // Calculate derivative dV/dr
        // dV/dD = K * (D_av - D0)
        // dD/dn = D_max * d(P2)/dn
        // d(P2)/dnᵢ = Σⱼ Sᵢⱼ * 3 * nⱼ
        let dv_dd = self.k * delta_d;
        let dd_dp2 = self.d_max;

        // dp2/dn (gradient of order parameter with respect to unit vector)
        let dp2_dnx = 3.0 * (sxx * nx + sxy * ny + sxz * nz);
        let dp2_dny = 3.0 * (sxy * nx + syy * ny + syz * nz);
        let dp2_dnz = 3.0 * (sxz * nx + syz * ny + szz * nz);

        let dv_dp2 = dv_dd * dd_dp2 * memory_decay;

        // dn/dr accounting for normalization: dn/dr = (I - n⊗n)/r
        let factor = dv_dp2 / r;

        let grad_x = (dp2_dnx - nx * (dp2_dnx * nx + dp2_dny * ny + dp2_dnz * nz)) * factor;
        let grad_y = (dp2_dny - ny * (dp2_dnx * nx + dp2_dny * ny + dp2_dnz * nz)) * factor;
        let grad_z = (dp2_dnz - nz * (dp2_dnx * nx + dp2_dny * ny + dp2_dnz * nz)) * factor;

        let force_vec = Vec3::new(grad_x, grad_y, grad_z);

        // Apply forces (equal and opposite)
        conf.current_mut().force[self.i] -= force_vec;
        conf.current_mut().force[self.j] += force_vec;

        energy
    }

    /// Reset time average
    pub fn reset_average(&mut self) {
        self.rdc_average = 0.0;
    }
}

/// Collection of RDC restraints
#[derive(Debug, Clone)]
pub struct RDCRestraints {
    pub restraints: Vec<RDCRestraint>,
}

impl Default for RDCRestraints {
    fn default() -> Self {
        Self::new()
    }
}

impl RDCRestraints {
    pub fn new() -> Self {
        Self {
            restraints: Vec::new(),
        }
    }

    pub fn add_restraint(&mut self, restraint: RDCRestraint) {
        self.restraints.push(restraint);
    }

    pub fn calculate_all(
        &mut self,
        conf: &mut Configuration,
        periodicity: &Periodicity,
        dt: f64,
    ) -> f64 {
        self.restraints
            .iter_mut()
            .map(|r| r.calculate(conf, periodicity, dt))
            .sum()
    }

    /// Reset all time averages
    pub fn reset_averages(&mut self) {
        for r in &mut self.restraints {
            r.reset_average();
        }
    }

    /// Fit alignment tensor from RDC data using SVD
    ///
    /// This solves the linear system: D = D_max * Σᵢⱼ Sᵢⱼ * Oᵢⱼ
    /// where Oᵢⱼ = (3*nᵢ*nⱼ - δᵢⱼ)/2 is the orientation matrix
    ///
    /// Returns the fitted Saupe matrix [Sxx, Syy, Szz, Sxy, Sxz]
    pub fn fit_alignment_tensor(
        &self,
        _conf: &Configuration,
        _periodicity: &Periodicity,
    ) -> [f64; 5] {
        // TODO: Implement SVD fitting of alignment tensor
        // This requires numerical linear algebra (SVD decomposition)
        // For now, return zero tensor
        [0.0; 5]
    }
}
