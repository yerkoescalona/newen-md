//! QM/MM Hybrid Simulations
//!
//! Direct translation of md++/src/interaction/qmmm/*.cc
//!
//! **Purpose**: Combine quantum mechanics (QM) and molecular mechanics (MM)
//! for reactive molecular dynamics simulations
//!
//! ## Overview
//!
//! QM/MM allows treating a small reactive region quantum mechanically while
//! the environment is treated classically. This enables:
//! - **Chemical reactions**: Bond breaking/formation
//! - **Electronic excitations**: Spectroscopy, photochemistry
//! - **Charge transfer**: Redox reactions
//! - **Metalloprotein catalysis**: Transition metal chemistry
//!
//! ## Architecture
//!
//! ```text
//! ┌──────────────────────────────────────────┐
//! │            Full System (MD)               │
//! │  ┌──────────────┐                        │
//! │  │  QM Region   │  ←─ Quantum Calculation│
//! │  │  (50-200 at.)│     (DFT, SCC-DFTB,   │
//! │  └──────────────┘      Semi-empirical)  │
//! │         ↕                                 │
//! │  QM/MM Interface                          │
//! │  - Electrostatic embedding               │
//! │  - Link atoms (capping)                  │
//! │  - Boundary treatment                    │
//! │         ↕                                 │
//! │  ┌──────────────────────┐                │
//! │  │   MM Region          │  ←─ Classical │
//! │  │  (thousands of atoms)│      Force    │
//! │  └──────────────────────┘      Field    │
//! └──────────────────────────────────────────┘
//! ```
//!
//! ## Partitioning Schemes
//!
//! 1. **Subtractive**: E_total = E_MM(full) + E_QM(QM) - E_MM(QM)
//! 2. **Additive**: E_total = E_QM(QM) + E_MM(MM) + E_QM-MM(interface)
//!
//! ## QM Engines Supported (md++ has 9, newen framework for all)
//!
//! 1. **xTB** (GFN1-xTB, GFN2-xTB) - Fast, accurate tight-binding DFT
//! 2. **DFTB+** - Self-Consistent Charge DFTB
//! 3. **Gaussian** - Full DFT/HF/post-HF
//! 4. **ORCA** - Modern DFT package
//! 5. **MOPAC** - Semi-empirical (AM1, PM3, PM6, PM7)
//! 6. **MNDO** - Classic semi-empirical
//! 7. **Turbomole** - Efficient DFT
//! 8. **Neural Network** - Machine learning potentials
//! 9. **Ghost** - Mock engine for testing
//!
//! ## References
//!
//! - Senn & Thiel (2009). "QM/MM methods for biomolecular systems."
//!   Angew. Chem. Int. Ed. 48:1198-1229
//! - Lin & Truhlar (2007). "QM/MM: what have we learned, where are we, and where do we go from here?"
//!   Theor. Chem. Acc. 117:185-199
//! - Groenhof (2013). "Introduction to QM/MM simulations." Methods Mol. Biol. 924:43-66

use newen_core::configuration::Configuration;
use newen_core::math::Vec3;
use newen_core::topology::Topology;
use std::io::Write;
use std::process::{Command, Stdio};

/// QM engine type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QMEngine {
    /// No QM (pure MM)
    None = 0,
    /// xTB tight-binding DFT (GFN1-xTB, GFN2-xTB)
    XTB = 1,
    /// DFTB+ self-consistent charge DFTB
    DFTBPlus = 2,
    /// Gaussian (G09, G16)
    Gaussian = 3,
    /// ORCA DFT package
    ORCA = 4,
    /// MOPAC semi-empirical
    MOPAC = 5,
    /// MNDO semi-empirical
    MNDO = 6,
    /// Turbomole DFT
    Turbomole = 7,
    /// Neural network potential (SchNet, ANI, etc.)
    NeuralNetwork = 8,
    /// Ghost engine (testing/debugging)
    Ghost = 9,
}

/// QM calculation method
#[derive(Debug, Clone)]
pub enum QMMethod {
    /// GFN1-xTB (fast, accurate for geometry and vibrations)
    GFN1XTB,
    /// GFN2-xTB (improved energies, good for reactions)
    GFN2XTB,
    /// SCC-DFTB (Self-Consistent Charge DFTB)
    SCCDFTB,
    /// AM1 semi-empirical
    AM1,
    /// PM3 semi-empirical
    PM3,
    /// PM6 semi-empirical (improved for inorganic)
    PM6,
    /// PM7 semi-empirical (latest MOPAC)
    PM7,
    /// B3LYP/6-31G* DFT
    B3LYP,
    /// PBE0/def2-SVP DFT
    PBE0,
    /// Custom method string
    Custom(String),
}

/// QM/MM coupling scheme
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CouplingScheme {
    /// Mechanical embedding: QM charges not polarized by MM
    Mechanical = 1,
    /// Electrostatic embedding: MM point charges included in QM Hamiltonian
    Electrostatic = 2,
    /// Polarized embedding: MM dipoles respond to QM density
    Polarized = 3,
}

/// Link atom (hydrogen cap) for QM/MM boundary
///
/// When a bond crosses the QM/MM boundary, it's capped with a link atom (usually H)
#[derive(Debug, Clone)]
pub struct LinkAtom {
    /// QM atom index (the atom in QM region)
    pub qm_atom: usize,
    /// MM atom index (the atom in MM region)
    pub mm_atom: usize,
    /// Link atom position (dynamically calculated)
    pub position: Vec3,
    /// Bond ratio g (link atom position: r_link = r_QM + g*(r_MM - r_QM))
    pub bond_ratio: f64,
    /// Element for link atom (usually 1 for H)
    pub element: u8,
}

impl LinkAtom {
    /// Create hydrogen link atom
    pub fn hydrogen(qm_atom: usize, mm_atom: usize) -> Self {
        Self {
            qm_atom,
            mm_atom,
            position: Vec3::ZERO,
            bond_ratio: 0.723, // C-H / C-C bond ratio
            element: 1,        // Hydrogen
        }
    }

    /// Calculate link atom position
    pub fn update_position(&mut self, conf: &Configuration) {
        let r_qm = conf.current().pos[self.qm_atom];
        let r_mm = conf.current().pos[self.mm_atom];
        self.position = r_qm + (r_mm - r_qm) * (self.bond_ratio);
    }
}

/// QM region definition
#[derive(Debug, Clone)]
pub struct QMRegion {
    /// Atoms in QM region
    pub atoms: Vec<usize>,
    /// Total charge of QM region (integer)
    pub charge: i32,
    /// Spin multiplicity (1=singlet, 2=doublet, 3=triplet, ...)
    pub multiplicity: i32,
    /// Link atoms at QM/MM boundary
    pub link_atoms: Vec<LinkAtom>,
}

impl QMRegion {
    pub fn new(atoms: Vec<usize>, charge: i32, multiplicity: i32) -> Self {
        Self {
            atoms,
            charge,
            multiplicity,
            link_atoms: Vec::new(),
        }
    }

    pub fn add_link_atom(&mut self, link: LinkAtom) {
        self.link_atoms.push(link);
    }

    /// Update all link atom positions
    pub fn update_link_atoms(&mut self, conf: &Configuration) {
        for link in &mut self.link_atoms {
            link.update_position(conf);
        }
    }
}

/// QM calculation result
#[derive(Debug, Clone)]
pub struct QMResult {
    /// QM energy (kJ/mol)
    pub energy: f64,
    /// Forces on QM atoms (kJ/mol/nm)
    pub forces: Vec<Vec3>,
    /// Mulliken charges on QM atoms (e)
    pub charges: Vec<f64>,
    /// Success flag
    pub success: bool,
    /// Error message if failed
    pub error_message: Option<String>,
}

/// QM engine interface trait
///
/// All QM engines must implement this trait
pub trait QMEngineInterface {
    /// Initialize QM engine
    fn initialize(&mut self, qm_region: &QMRegion) -> Result<(), String>;

    /// Run QM calculation
    fn calculate(
        &mut self,
        qm_region: &QMRegion,
        conf: &Configuration,
        mm_charges: &[(Vec3, f64)], // MM point charges for electrostatic embedding
    ) -> Result<QMResult, String>;

    /// Finalize QM engine (cleanup)
    fn finalize(&mut self);

    /// Get engine name
    fn name(&self) -> &str;
}

/// xTB QM engine implementation
pub struct XTBEngine {
    /// xTB executable path
    pub xtb_binary: String,
    /// QM method (GFN1-xTB or GFN2-xTB)
    pub method: QMMethod,
    /// Working directory for QM calculations
    pub work_dir: String,
    /// Verbosity level
    pub verbosity: i32,
}

impl XTBEngine {
    pub fn new(method: QMMethod) -> Self {
        Self {
            xtb_binary: "xtb".to_string(),
            method,
            work_dir: "/tmp/qmmm_xtb".to_string(),
            verbosity: 0,
        }
    }
}

impl QMEngineInterface for XTBEngine {
    fn initialize(&mut self, _qm_region: &QMRegion) -> Result<(), String> {
        // Create working directory
        std::fs::create_dir_all(&self.work_dir)
            .map_err(|e| format!("Failed to create work dir: {}", e))?;
        Ok(())
    }

    fn calculate(
        &mut self,
        qm_region: &QMRegion,
        conf: &Configuration,
        mm_charges: &[(Vec3, f64)],
    ) -> Result<QMResult, String> {
        // Write XYZ file for QM region
        let xyz_file = format!("{}/qm_region.xyz", self.work_dir);
        self.write_xyz_file(&xyz_file, qm_region, conf)?;

        // Write point charges file for electrostatic embedding
        if !mm_charges.is_empty() {
            let pc_file = format!("{}/point_charges.pc", self.work_dir);
            self.write_point_charges(&pc_file, mm_charges)?;
        }

        // Run xTB calculation
        let method_flag = match self.method {
            QMMethod::GFN1XTB => "--gfn 1",
            QMMethod::GFN2XTB => "--gfn 2",
            _ => "--gfn 2", // Default
        };

        let output = Command::new(&self.xtb_binary)
            .current_dir(&self.work_dir)
            .arg("qm_region.xyz")
            .arg(method_flag)
            .arg("--grad")
            .arg(format!("--chrg {}", qm_region.charge))
            .arg(format!("--uhf {}", qm_region.multiplicity - 1))
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()
            .map_err(|e| format!("Failed to run xTB: {}", e))?;

        if !output.status.success() {
            return Err(format!(
                "xTB calculation failed: {}",
                String::from_utf8_lossy(&output.stderr)
            ));
        }

        // Parse results
        self.parse_xtb_output(qm_region)
    }

    fn finalize(&mut self) {
        // Cleanup temporary files (optional)
    }

    fn name(&self) -> &str {
        "xTB"
    }
}

impl XTBEngine {
    fn write_xyz_file(
        &self,
        filename: &str,
        qm_region: &QMRegion,
        conf: &Configuration,
    ) -> Result<(), String> {
        let mut file = std::fs::File::create(filename)
            .map_err(|e| format!("Failed to create XYZ file: {}", e))?;

        writeln!(
            file,
            "{}",
            qm_region.atoms.len() + qm_region.link_atoms.len()
        )
        .map_err(|e| e.to_string())?;
        writeln!(file, "QM region from newen").map_err(|e| e.to_string())?;

        // Write QM atoms
        for &atom_idx in &qm_region.atoms {
            let pos = conf.current().pos[atom_idx];
            // TODO: Get element symbol from topology
            let element = "C"; // Placeholder
            writeln!(
                file,
                "{} {:12.6} {:12.6} {:12.6}",
                element,
                pos.x * 10.0,
                pos.y * 10.0,
                pos.z * 10.0
            ) // nm to Å
            .map_err(|e| e.to_string())?;
        }

        // Write link atoms
        for link in &qm_region.link_atoms {
            let pos = link.position;
            writeln!(
                file,
                "H {:12.6} {:12.6} {:12.6}",
                pos.x * 10.0,
                pos.y * 10.0,
                pos.z * 10.0
            )
            .map_err(|e| e.to_string())?;
        }

        Ok(())
    }

    fn write_point_charges(
        &self,
        filename: &str,
        mm_charges: &[(Vec3, f64)],
    ) -> Result<(), String> {
        let mut file = std::fs::File::create(filename)
            .map_err(|e| format!("Failed to create point charges file: {}", e))?;

        writeln!(file, "{}", mm_charges.len()).map_err(|e| e.to_string())?;

        for (pos, charge) in mm_charges {
            writeln!(
                file,
                "{:12.6} {:12.6} {:12.6} {:12.6}",
                pos.x * 10.0,
                pos.y * 10.0,
                pos.z * 10.0,
                charge
            )
            .map_err(|e| e.to_string())?;
        }

        Ok(())
    }

    fn parse_xtb_output(&self, qm_region: &QMRegion) -> Result<QMResult, String> {
        // Parse energy from energy file
        let energy_file = format!("{}/energy", self.work_dir);
        let energy_str = std::fs::read_to_string(&energy_file)
            .map_err(|e| format!("Failed to read energy file: {}", e))?;
        let energy: f64 = energy_str
            .trim()
            .parse()
            .map_err(|e| format!("Failed to parse energy: {}", e))?;
        let energy_kj = energy * 2625.5; // Hartree to kJ/mol

        // Parse forces from gradient file
        let grad_file = format!("{}/gradient", self.work_dir);
        let grad_content = std::fs::read_to_string(&grad_file)
            .map_err(|e| format!("Failed to read gradient file: {}", e))?;

        let mut forces = Vec::new();
        for line in grad_content.lines().skip(2) {
            // Skip header
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                let fx: f64 = parts[0].parse().unwrap_or(0.0);
                let fy: f64 = parts[1].parse().unwrap_or(0.0);
                let fz: f64 = parts[2].parse().unwrap_or(0.0);
                // Convert Hartree/Bohr to kJ/mol/nm
                forces.push(Vec3::new(fx * 49614.75, fy * 49614.75, fz * 49614.75));
            }
        }

        // Parse Mulliken charges (optional)
        let charges = vec![0.0; qm_region.atoms.len()];

        Ok(QMResult {
            energy: energy_kj,
            forces,
            charges,
            success: true,
            error_message: None,
        })
    }
}

/// DFTB+ engine (skeleton)
pub struct DFTBPlusEngine {
    pub dftb_binary: String,
    pub slako_dir: String,
}

impl QMEngineInterface for DFTBPlusEngine {
    fn initialize(&mut self, _qm_region: &QMRegion) -> Result<(), String> {
        Ok(())
    }

    fn calculate(
        &mut self,
        _qm_region: &QMRegion,
        _conf: &Configuration,
        _mm_charges: &[(Vec3, f64)],
    ) -> Result<QMResult, String> {
        Err("DFTB+ engine not yet implemented".to_string())
    }

    fn finalize(&mut self) {}

    fn name(&self) -> &str {
        "DFTB+"
    }
}

/// QM/MM calculator
///
/// Main driver for QM/MM simulations
pub struct QMMMCalculator {
    /// QM engine
    pub engine: Box<dyn QMEngineInterface>,
    /// QM region definition
    pub qm_region: QMRegion,
    /// Coupling scheme
    pub coupling: CouplingScheme,
    /// MM point charge cutoff (nm)
    pub mm_cutoff: f64,
}

impl QMMMCalculator {
    pub fn new(
        engine: Box<dyn QMEngineInterface>,
        qm_region: QMRegion,
        coupling: CouplingScheme,
    ) -> Self {
        Self {
            engine,
            qm_region,
            coupling,
            mm_cutoff: 2.0, // 2 nm cutoff for MM charges
        }
    }

    /// Calculate QM/MM energy and forces
    ///
    /// # Algorithm
    /// 1. Update link atom positions
    /// 2. Collect MM point charges within cutoff
    /// 3. Run QM calculation
    /// 4. Redistribute QM forces to MM atoms (chain rule through link atoms)
    /// 5. Calculate MM forces in MM region
    /// 6. Return total energy
    pub fn calculate(&mut self, topo: &Topology, conf: &mut Configuration) -> Result<f64, String> {
        // 1. Update link atoms
        self.qm_region.update_link_atoms(conf);

        // 2. Collect MM point charges (for electrostatic embedding)
        let mm_charges = if self.coupling == CouplingScheme::Electrostatic {
            self.collect_mm_charges(topo, conf)
        } else {
            Vec::new()
        };

        // 3. Run QM calculation
        let qm_result = self.engine.calculate(&self.qm_region, conf, &mm_charges)?;

        if !qm_result.success {
            return Err(qm_result
                .error_message
                .unwrap_or_else(|| "QM calculation failed".to_string()));
        }

        // 4. Apply QM forces to QM atoms
        for (i, &atom_idx) in self.qm_region.atoms.iter().enumerate() {
            if i < qm_result.forces.len() {
                conf.current_mut().force[atom_idx] += qm_result.forces[i];
            }
        }

        // 5. Redistribute link atom forces (chain rule)
        self.redistribute_link_forces(&qm_result, conf);

        // 6. Calculate MM forces in MM region
        // (This is done separately by the MM force field)

        Ok(qm_result.energy)
    }

    /// Collect MM point charges within cutoff of QM region
    fn collect_mm_charges(&self, topo: &Topology, conf: &Configuration) -> Vec<(Vec3, f64)> {
        let mut mm_charges = Vec::new();

        // QM center of mass
        let mut qm_com = Vec3::ZERO;
        for &atom_idx in &self.qm_region.atoms {
            qm_com += conf.current().pos[atom_idx];
        }
        qm_com /= self.qm_region.atoms.len() as f64;

        // Collect MM charges
        for i in 0..topo.num_atoms() {
            // Skip QM atoms
            if self.qm_region.atoms.contains(&i) {
                continue;
            }

            let pos = conf.current().pos[i];
            let dist = (pos - qm_com).length();

            if dist < self.mm_cutoff {
                let charge = topo.charge[i];
                mm_charges.push((pos, charge));
            }
        }

        mm_charges
    }

    /// Redistribute link atom forces to QM and MM atoms
    ///
    /// Chain rule: F_QM += (1-g)*F_link, F_MM += g*F_link
    fn redistribute_link_forces(&self, qm_result: &QMResult, conf: &mut Configuration) {
        let n_qm_atoms = self.qm_region.atoms.len();

        for (i, link) in self.qm_region.link_atoms.iter().enumerate() {
            let link_force_idx = n_qm_atoms + i;
            if link_force_idx >= qm_result.forces.len() {
                continue;
            }

            let f_link = qm_result.forces[link_force_idx];
            let g = link.bond_ratio;

            conf.current_mut().force[link.qm_atom] += f_link * (1.0 - g);
            conf.current_mut().force[link.mm_atom] += f_link * g;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qm_region() {
        let mut qm_region = QMRegion::new(vec![0, 1, 2], 0, 1);
        qm_region.add_link_atom(LinkAtom::hydrogen(2, 3));

        assert_eq!(qm_region.atoms.len(), 3);
        assert_eq!(qm_region.link_atoms.len(), 1);
        assert_eq!(qm_region.charge, 0);
        assert_eq!(qm_region.multiplicity, 1);
    }

    #[test]
    fn test_link_atom() {
        let link = LinkAtom::hydrogen(0, 1);
        assert_eq!(link.element, 1);
        assert!((link.bond_ratio - 0.723).abs() < 1e-6);
    }

    #[test]
    fn test_xtb_engine() {
        let engine = XTBEngine::new(QMMethod::GFN2XTB);
        assert_eq!(engine.name(), "xTB");
    }
}
