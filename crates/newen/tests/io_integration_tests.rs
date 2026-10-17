//! Integration tests for GROMOS I/O modules
//!
//! Tests the complete I/O workflow across all sub-crates:
//! - IMD file reading and parsing
//! - Trajectory writing (TRC)
//! - Energy writing (TRE)
//! - Force writing (TRF)
//! - Coordinate file parsing
//! - Topology file parsing

use newen::configuration::{Box as SimBox, Configuration};
use newen::io::{
    coordinate::read_coordinate_file,
    energy::{EnergyFrame, EnergyWriter},
    force::ForceWriter,
    imd::{read_imd_file, ImdParameters},
    topology::read_topology_file,
    trajectory::TrajectoryWriter,
};
use newen::math::Vec3;
use std::fs;
use std::io::Write;
use std::path::Path;

// ─── helpers ────────────────────────────────────────────────────────────────

fn tmp(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(name)
}

fn write_tmp(content: &str, name: &str) -> std::path::PathBuf {
    let path = tmp(name);
    let mut f = fs::File::create(&path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
    path
}

// ─── IMD ────────────────────────────────────────────────────────────────────

const IMD_FULL: &str = r"TITLE
  Integration test simulation
END

SYSTEM
#  NPM    NSM
     1      0
END

STEP
#  NSTLIM      T        DT
    1000    0.0     0.00200
END

BOUNDCOND
#   NTB  NDFMIN
      1       0
END

MULTIBATH
#  ALGORITHM
          1
#  NUM
      1
#  TEMP0         TAU
  300.00      0.1000
#  DOFSET
      1
#  LAST  COM-BATH  IR-BATH
      2     1         1
END

CONSTRAINT
#  NTC  NTCP  NTCP0(1)  NTCS  NTCS0(1)
    2     0         0      1         0
END

PAIRLIST
#  ALGORITHM  NSNB  RCUTP   RCUTL     SIZE  TYPE
          0     5   0.80    1.40     0.40     0
END

NONBONDED
#  NLRELE  APPAK    RCRF   EPSRF  NSLFEXCL
       1    0.0    1.40     0.0         1
END

INITIALISE
#  NTIVEL  NTISHK  NTINHT  NTINHB  NTISHI  NTIRTC  NTICOM  NTISTI
        1       0       0       0       0       0       0       0
#      IG   TEMPI
    12345  300.00
END

WRITETRAJ
#  NTWX  NTWSE  NTWV  NTWF  NTWE
    100      0     0     0    10
END

PRINTOUT
#  NTPR
     10
END
";

#[test]
fn test_imd_write_and_read() {
    let path = write_tmp(IMD_FULL, "integ_test.imd");
    let params = read_imd_file(&path).expect("Failed to read IMD file");

    assert_eq!(params.title, "Integration test simulation");
    assert_eq!(params.npm, 1);
    assert_eq!(params.nsm, 0);
    assert_eq!(params.nstlim, 1000);
    assert_eq!(params.dt, 0.002);
    assert_eq!(params.ntc, 2);
    assert_eq!(params.nlrele, 1);
    assert_eq!(params.rcutl, 1.4);
    assert_eq!(params.ntwx, 100);
    assert_eq!(params.ntwe, 10);
    assert_eq!(params.tempi, 300.0);

    fs::remove_file(path).ok();
}

#[test]
fn test_imd_parameter_defaults() {
    let params = ImdParameters::default();

    assert_eq!(params.nstlim, 1000);
    assert_eq!(params.dt, 0.002);
    assert_eq!(params.ntc, 1);
    assert_eq!(params.nlrele, 1);
    assert_eq!(params.epsrf, 0.0);
    assert_eq!(params.tempi, 300.0);
}

// ─── Trajectory ─────────────────────────────────────────────────────────────

fn make_conf(n_atoms: usize) -> Configuration {
    let mut config = Configuration::new(n_atoms, 1, 1);
    config.current_mut().pos = (0..n_atoms)
        .map(|i| Vec3::new(i as f64 * 0.1, 0.0, 0.0))
        .collect();
    config.current_mut().vel = (0..n_atoms)
        .map(|i| Vec3::new(i as f64 * 0.01, 0.0, 0.0))
        .collect();
    config.current_mut().box_config = SimBox::rectangular(3.0, 3.0, 3.0);
    config
}

#[test]
fn test_trajectory_writer_integration() {
    let path = tmp("integ_traj.trc");
    let config = make_conf(3);

    let mut writer = TrajectoryWriter::new(&path, "Integration test trajectory", true, false)
        .expect("Failed to create trajectory writer");

    for step in 0..5_usize {
        writer
            .write_frame(step, step as f64 * 0.002, &config)
            .expect("Failed to write frame");
    }
    writer.flush().expect("Failed to flush");

    assert!(path.exists());
    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("TITLE"));
    assert!(content.contains("TIMESTEP"));
    assert!(content.contains("POSITIONRED"));
    assert!(content.contains("VELOCITYRED"));
    assert!(content.contains("GENBOX"));
    assert_eq!(content.matches("TIMESTEP").count(), 5);
    assert_eq!(writer.frame_count(), 5);

    fs::remove_file(path).ok();
}

#[test]
fn test_trajectory_no_velocities() {
    let path = tmp("integ_traj_novel.trc");
    let config = make_conf(3);

    let mut writer = TrajectoryWriter::new(&path, "no-vel test", false, false).unwrap();
    writer.write_frame(0, 0.0, &config).unwrap();
    writer.flush().unwrap();

    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("POSITIONRED"));
    assert!(!content.contains("VELOCITYRED"));

    fs::remove_file(path).ok();
}

// ─── Energy ─────────────────────────────────────────────────────────────────

#[test]
fn test_energy_writer_integration() {
    let path = tmp("integ_energy.tre");

    let mut writer = EnergyWriter::new(&path, "Integration test energy")
        .expect("Failed to create energy writer");

    for step in 0..10_usize {
        let time = step as f64 * 0.002;
        let mut frame = EnergyFrame::new(time, 100.0 + step as f64, -200.0 - step as f64, 300.0);
        frame.bond = -50.0;
        frame.angle = -30.0;
        frame.lj = -80.0;
        frame.coul_real = -40.0;
        frame.update_potential();
        writer
            .write_frame(&frame)
            .expect("Failed to write energy frame");
    }

    writer.finalize().expect("Failed to finalize");

    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("TITLE"));
    // gromosXX's native layout (0.0.34): TIMESTEP + ENERGY03 + VOLUMEPRESSURE03 per frame
    assert!(content.contains("ENERGY03"));
    assert!(content.contains("VOLUMEPRESSURE03"));
    // the self-description in the header names the TIMESTEP block too: count block openers, not
    // substring hits
    assert_eq!(content.lines().filter(|l| *l == "TIMESTEP").count(), 10);
    assert_eq!(writer.frame_count(), 10);

    fs::remove_file(path).ok();
}

#[test]
fn test_energy_frame_updates() {
    let mut frame = EnergyFrame {
        kinetic: 100.0,
        bond: -50.0,
        angle: -30.0,
        lj: -80.0,
        coul_real: -40.0,
        ..Default::default()
    };

    frame.update_potential();
    frame.update_total();

    assert_eq!(frame.potential, -200.0);
    assert_eq!(frame.total, -100.0);
}

// ─── Force ──────────────────────────────────────────────────────────────────

#[test]
fn test_force_writer_integration() {
    let path = tmp("integ_force.trf");
    let forces = vec![
        Vec3::new(1.0, 2.0, 3.0),
        Vec3::new(4.0, 5.0, 6.0),
        Vec3::new(7.0, 8.0, 9.0),
    ];

    let mut writer = ForceWriter::new(&path, "Integration test forces", false)
        .expect("Failed to create force writer");

    for step in 0..3_usize {
        writer
            .write_frame(step, step as f64 * 0.002, &forces, None)
            .expect("Failed to write force frame");
    }
    writer.flush().expect("Failed to flush");

    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("TITLE"));
    assert!(content.contains("TIMESTEP"));
    assert!(content.contains("FREEFORCERED"));
    assert_eq!(content.matches("TIMESTEP").count(), 3);

    fs::remove_file(path).ok();
}

#[test]
fn test_force_writer_with_constraints() {
    let path = tmp("integ_force_cons.trf");
    let n = 3;
    let forces: Vec<Vec3> = (0..n).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
    let constraints: Vec<Vec3> = (0..n).map(|i| Vec3::new(0.0, i as f64, 0.0)).collect();

    let mut writer = ForceWriter::new(&path, "constraint force test", true).unwrap();
    writer
        .write_frame(0, 0.0, &forces, Some(&constraints))
        .unwrap();
    writer.flush().unwrap();

    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("FREEFORCERED"));
    assert!(content.contains("CONSFORCERED"));

    fs::remove_file(path).ok();
}

// ─── Coordinate ─────────────────────────────────────────────────────────────

const CONF_4ATOMS: &str = "\
TITLE
  4-atom integration test
END
POSITION
    1 RES    C         1   1.000000000  2.000000000  3.000000000
    1 RES    H         2   1.100000000  2.100000000  3.100000000
    1 RES    H         3   0.900000000  1.900000000  2.900000000
    1 RES    H         4   1.000000000  2.000000000  3.200000000
END
BOX
  10.0  10.0  10.0
END
";

#[test]
fn test_coordinate_read_integration() {
    let path = write_tmp(CONF_4ATOMS, "integ_conf.conf");
    let conf = read_coordinate_file(&path, 1, 1).expect("Failed to read coordinate file");

    assert_eq!(conf.current().pos.len(), 4);
    assert!((conf.current().pos[0].x - 1.0_f64).abs() < 1e-5);
    assert!((conf.current().pos[0].y - 2.0_f64).abs() < 1e-5);
    assert!((conf.current().pos[0].z - 3.0_f64).abs() < 1e-5);

    let dims = conf.current().box_config.dimensions();
    assert!((dims.x - 10.0_f64).abs() < 1e-5);

    fs::remove_file(path).ok();
}

#[test]
fn test_coordinate_read_missing_file() {
    let result = read_coordinate_file("/no/such/file.conf", 1, 1);
    assert!(result.is_err());
}

// ─── Topology ───────────────────────────────────────────────────────────────

const TOPO_4ATOMS: &str = "\
TITLE
  4-atom integration test topology
END
SOLUTEATOM
# ATNM MRES PANM IAC   MASS      CG   CGC INE IEXCL
4
    1    1   C     6  12.000000  0.000000  0  2  2  3
                                              0
    2    1   H     1   1.000000  0.000000  0  1  1
                                              0
    3    1   H     1   1.000000  0.000000  0  1  1
                                              0
    4    1   H     1   1.000000  0.000000  0  0
                                              0
END
BONDSTRETCHTYPE
# CB        CHB       B0
1
  1000.000  500.000   0.109
END
BOND
# IB  JB  ICB
3
  1   2   1
  1   3   1
  1   4   1
END
BONDANGLEBENDTYPE
# CT        CHT       T0
1
   50.000  100.000  109.500
END
BONDANGLE
# IT  JT  KT  ICT
2
  2   1   3   1
  2   1   4   1
END
";

#[test]
fn test_topology_read_integration() {
    let path = write_tmp(TOPO_4ATOMS, "integ_topo.topo");
    let topo = read_topology_file(&path).expect("Failed to read topology file");

    assert_eq!(topo.n_atoms, 4);
    assert_eq!(topo.bonds.len(), 3);
    assert_eq!(topo.angles.len(), 2);
    assert_eq!(topo.iac[0], 5); // file stores IAC 1-indexed; parser subtracts 1 at boundary
    assert!((topo.masses[0] - 12.0).abs() < 1e-6);

    // Bond 0: atoms 1-2 in file → 0-1 (0-indexed)
    assert_eq!(topo.bonds[0], (0, 1, 0));

    fs::remove_file(path).ok();
}

#[test]
fn test_topology_zero_indexed_atoms() {
    let path = write_tmp(TOPO_4ATOMS, "integ_topo_idx.topo");
    let topo = read_topology_file(&path).unwrap();

    for (i, j, _) in &topo.bonds {
        assert!(*i < topo.n_atoms, "bond atom i={i} out of range");
        assert!(*j < topo.n_atoms, "bond atom j={j} out of range");
    }
    for (i, j, k, _) in &topo.angles {
        assert!(*i < topo.n_atoms);
        assert!(*j < topo.n_atoms);
        assert!(*k < topo.n_atoms);
    }

    fs::remove_file(path).ok();
}

// ─── Complete workflow ───────────────────────────────────────────────────────

#[test]
fn test_complete_io_workflow() {
    let imd_content = r"TITLE
  Complete workflow test
END
SYSTEM
#  NPM    NSM
     1      0
END
STEP
#  NSTLIM      T        DT
    100      0.0     0.002
END
BOUNDCOND
#  NTB  NDFMIN
      1       0
END
WRITETRAJ
#  NTWX  NTWSE  NTWV  NTWF  NTWE
      5      0     0     0     5
END
PRINTOUT
#  NTPR
     5
END
";

    let imd_path = write_tmp(imd_content, "workflow_test.imd");
    let trc_path = tmp("workflow_test.trc");
    let tre_path = tmp("workflow_test.tre");

    // 1. Read IMD
    let params = read_imd_file(&imd_path).expect("Failed to read IMD");
    assert_eq!(params.nstlim, 100);
    assert_eq!(params.ntwx, 5);

    // 2. Create configuration
    let mut config = Configuration::new(10, 1, 1);
    for i in 0..10_usize {
        config.current_mut().pos[i] = Vec3::new(i as f64 * 0.1, 0.0, 0.0);
    }
    config.current_mut().box_config = SimBox::rectangular(5.0, 5.0, 5.0);

    // 3. Open writers
    let mut traj = TrajectoryWriter::new(&trc_path, &params.title, false, false).unwrap();
    let mut energy = EnergyWriter::new(&tre_path, &params.title).unwrap();

    // 4. Simulate steps
    for step in 0..params.nstlim {
        let time = step as f64 * params.dt;
        if step % params.ntwx == 0 {
            traj.write_frame(step, time, &config).unwrap();
        }
        if step % params.ntwe == 0 {
            let frame = EnergyFrame::new(time, 100.0, -200.0, 300.0);
            energy.write_frame(&frame).unwrap();
        }
    }

    traj.flush().unwrap();
    energy.finalize().unwrap();

    // 5. Verify
    assert!(Path::new(&trc_path).exists());
    assert!(Path::new(&tre_path).exists());
    assert_eq!(traj.frame_count(), 20); // 100 / 5
    assert_eq!(energy.frame_count(), 20); // 100 / 5

    fs::remove_file(imd_path).ok();
    fs::remove_file(trc_path).ok();
    fs::remove_file(tre_path).ok();
}
