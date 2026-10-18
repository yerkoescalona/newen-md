//! GROMOS96 (.g96) file format reader/writer
//!
//! GROMOS96 format is a block-based format with keywords like TITLE, POSITION, VELOCITY, BOX

use crate::coordinate::G96Atom;
use newen_core::math::Vec3;
use newen_core::topology::Topology;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// The 24-character label a GROMOS96 position line carries before the coordinates
/// (`out_configuration.cc:982`): residue number, residue name, atom name, atom number.
pub fn position_prefix(resnum: usize, resname: &str, atomname: &str, atomnum: usize) -> String {
    format!("{resnum:>5} {resname:<5} {atomname:<6}{atomnum:>6}")
}

/// Labels for atom `i` as gromosXX writes them: a solute atom carries its own residue
/// number, a solvent atom carries the number of the solvent molecule it belongs to
/// (`out_configuration.cc:1015`) — not its index inside that molecule.
fn atom_label(topo: &Topology, i: usize) -> (usize, &str, &str) {
    let Some(a) = topo.moltype_atom(i) else {
        return (1, "UNK", "AT");
    };
    let n_solute = topo.num_solute_atoms();
    let resnum = if i < n_solute {
        a.residue_nr
    } else {
        let per_molecule = topo.atoms_per_solvent().max(1);
        (i - n_solute) / per_molecule + 1
    };
    (resnum, a.residue_name.as_str(), a.name.as_str())
}

/// Write a GROMOS96 format file
pub struct G96Writer {
    pub writer: BufWriter<File>,
}

impl G96Writer {
    /// Create a new G96 writer
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let file = File::create(path).map_err(|e| format!("Cannot create .g96 file: {}", e))?;
        Ok(G96Writer {
            writer: BufWriter::new(file),
        })
    }

    /// Write title block
    pub fn write_title(&mut self, title: &str) -> Result<(), String> {
        writeln!(self.writer, "TITLE").map_err(|e| format!("Write error: {}", e))?;
        writeln!(self.writer, "{}", title).map_err(|e| format!("Write error: {}", e))?;
        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write the TIMESTEP block (step counter and time, ps) — gromosXX `_print_timestep`
    /// (out_configuration.cc:879-888): width 15, a space, width 14 with 9 decimals.
    pub fn write_timestep(&mut self, step: usize, time: f64) -> Result<(), String> {
        writeln!(self.writer, "TIMESTEP\n{step:15} {time:14.9}\nEND")
            .map_err(|e| format!("Write error: {}", e))
    }

    /// Write position block
    ///
    /// GROMOS96 position format:
    /// ```text
    /// POSITION
    /// %5d %-5s %-5s%7d%15.9f%15.9f%15.9f
    /// (resnum, resname, atomname, atomnum, x, y, z)
    /// END
    /// ```
    pub fn write_positions(
        &mut self,
        positions: &[Vec3],
        topology: Option<&Topology>,
    ) -> Result<(), String> {
        writeln!(self.writer, "POSITION").map_err(|e| format!("Write error: {}", e))?;
        writeln!(self.writer, "# first 24 chars ignored")
            .map_err(|e| format!("Write error: {}", e))?;

        for (i, pos) in positions.iter().enumerate() {
            let (resnum, resname, atomname) =
                topology.map_or((1, "UNK", "AT"), |topo| atom_label(topo, i));

            writeln!(
                self.writer,
                "{}{:15.9}{:15.9}{:15.9}",
                position_prefix(resnum, resname, atomname, i + 1),
                pos.x,
                pos.y,
                pos.z
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write velocity block (optional)
    pub fn write_velocities(&mut self, velocities: &[Vec3]) -> Result<(), String> {
        writeln!(self.writer, "VELOCITY").map_err(|e| format!("Write error: {}", e))?;

        for vel in velocities {
            let vx = vel.x;
            let vy = vel.y;
            let vz = vel.z;

            writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", vx, vy, vz)
                .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write box dimensions (for periodic boundary conditions)
    pub fn write_box(&mut self, box_dims: Vec3) -> Result<(), String> {
        writeln!(self.writer, "BOX").map_err(|e| format!("Write error: {}", e))?;

        let x = box_dims.x;
        let y = box_dims.y;
        let z = box_dims.z;

        writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", x, y, z)
            .map_err(|e| format!("Write error: {}", e))?;

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write the PERTDATA block (λ) — gromosXX `_print_pertdata` (out_configuration.cc:2271-2277),
    /// written into the final configuration of a perturbed run.
    pub fn write_pertdata(&mut self, lambda: f64) -> Result<(), String> {
        writeln!(self.writer, "PERTDATA\n{lambda:15.7}\nEND")
            .map_err(|e| format!("Write error: {}", e))
    }

    /// Write the DISRESEXPAVE block (the ⟨r⁻³⟩ of time-averaged distance restraints, five per
    /// line) — gromosXX `_print_distance_restraint_averages` (out_configuration.cc:2555-2578).
    pub fn write_disresexpave(&mut self, averages: &[f64]) -> Result<(), String> {
        let err = |e: std::io::Error| format!("Write error: {}", e);
        writeln!(self.writer, "DISRESEXPAVE").map_err(err)?;
        for chunk in averages.chunks(5) {
            for a in chunk {
                write!(self.writer, "{a:15.7}").map_err(err)?;
            }
            writeln!(self.writer).map_err(err)?;
        }
        writeln!(self.writer, "END").map_err(err)
    }

    /// Write the AEDSSEARCH block (the A-EDS search state a continuation run starts from) —
    /// gromosXX `_print_aedssearch` (out_configuration.cc:3985-4004): fixed format, 9 decimals,
    /// width 15; flags and counts as integers.
    pub fn write_aedssearch(
        &mut self,
        s: &newen_core::configuration::AedsSearch,
    ) -> Result<(), String> {
        let err = |e: std::io::Error| format!("Write error: {}", e);
        let w = &mut self.writer;
        writeln!(w, "AEDSSEARCH").map_err(err)?;
        writeln!(w, "{:15.9}\n{:15.9}\n{:15.9}", s.emax, s.emin, s.searchemax).map_err(err)?;
        writeln!(
            w,
            "{:15}\n{:15}\n{:15}",
            s.emaxcounts,
            s.oldstate,
            u8::from(s.fullemin)
        )
        .map_err(err)?;
        for i in 0..s.num_states() {
            writeln!(
                w,
                "{:15.9} {:15.9} {:15.9} {:15} {:15} {:15.9} {:15.9} {:15.9} {:15.9} {:15}",
                s.eir[i],
                s.lnexpde[i],
                s.statefren[i],
                u8::from(s.visited[i]),
                s.visitcounts[i],
                s.avgenergy[i],
                s.eiravgenergy[i],
                s.bigs[i],
                s.stdevenergy[i],
                s.framecounts[i]
            )
            .map_err(err)?;
        }
        writeln!(w, "END").map_err(err)
    }

    /// Write position block with full atom labels (residue number, name, atom name).
    ///
    /// Atoms are renumbered sequentially starting from 1.
    pub fn write_labeled_positions(&mut self, atoms: &[G96Atom]) -> Result<(), String> {
        writeln!(self.writer, "POSITION").map_err(|e| format!("Write error: {}", e))?;

        for (i, atom) in atoms.iter().enumerate() {
            writeln!(
                self.writer,
                "{}{:15.9}{:15.9}{:15.9}",
                position_prefix(atom.res_num, &atom.res_name, &atom.atom_name, i + 1),
                atom.pos.x,
                atom.pos.y,
                atom.pos.z
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write POSRESSPEC block — specifies which atoms are positionally restrained.
    ///
    /// Same format as POSITION but with block name POSRESSPEC.
    /// Only the atom sequence numbers matter to the MD engine (first 17 chars
    /// and coordinates are ignored), but we write full labels for readability.
    pub fn write_posresspec(&mut self, atoms: &[G96Atom]) -> Result<(), String> {
        writeln!(self.writer, "POSRESSPEC").map_err(|e| format!("Write error: {}", e))?;

        for (i, atom) in atoms.iter().enumerate() {
            writeln!(
                self.writer,
                "{}{:15.9}{:15.9}{:15.9}",
                position_prefix(atom.res_num, &atom.res_name, &atom.atom_name, i + 1),
                atom.pos.x,
                atom.pos.y,
                atom.pos.z
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write REFPOSITION block — reference positions for position restraints.
    ///
    /// Same format as POSITION but with block name REFPOSITION.
    /// Must contain all atoms; the MD engine uses coordinates of atoms
    /// listed in the POSRESSPEC file as restraint targets.
    pub fn write_refposition(&mut self, atoms: &[G96Atom]) -> Result<(), String> {
        writeln!(self.writer, "REFPOSITION").map_err(|e| format!("Write error: {}", e))?;

        for (i, atom) in atoms.iter().enumerate() {
            writeln!(
                self.writer,
                "{}{:15.9}{:15.9}{:15.9}",
                position_prefix(atom.res_num, &atom.res_name, &atom.atom_name, i + 1),
                atom.pos.x,
                atom.pos.y,
                atom.pos.z
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Flush and close the writer
    pub fn close(mut self) -> Result<(), String> {
        self.writer
            .flush()
            .map_err(|e| format!("Flush error: {}", e))
    }
}

/// Write a configuration to a GROMOS96 file
pub fn write_g96<P: AsRef<Path>>(
    path: P,
    title: &str,
    positions: &[Vec3],
    velocities: Option<&[Vec3]>,
    box_dims: Option<Vec3>,
    topology: Option<&Topology>,
) -> Result<(), String> {
    let mut writer = G96Writer::new(path)?;

    writer.write_title(title)?;
    writer.write_positions(positions, topology)?;

    if let Some(vels) = velocities {
        writer.write_velocities(vels)?;
    }

    if let Some(box_size) = box_dims {
        writer.write_box(box_size)?;
    }

    writer.close()
}

/// Generate a .por (position restraint specification) file from coordinate data.
///
/// The .por file contains a POSRESSPEC block with only the solute atoms
/// (first `n_solute_atoms` atoms). The MD engine uses only the atom indices
/// from this block to determine which atoms are restrained.
pub fn write_por<P: AsRef<Path>>(
    path: P,
    atoms: &[G96Atom],
    n_solute_atoms: usize,
) -> Result<(), String> {
    let solute = &atoms[..n_solute_atoms];
    let mut writer = G96Writer::new(path)?;
    writer.write_title("Position restraint specification: solute atoms")?;
    writer.write_posresspec(solute)?;
    writer.close()
}

/// Generate a .rpr (reference positions) file from coordinate data.
///
/// The .rpr file contains a REFPOSITION block with all atoms.
/// The MD engine reads reference coordinates for restrained atoms from this file.
pub fn write_rpr<P: AsRef<Path>>(path: P, atoms: &[G96Atom]) -> Result<(), String> {
    let mut writer = G96Writer::new(path)?;
    writer.write_title("Reference positions for position restraints")?;
    writer.write_refposition(atoms)?;
    writer.close()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_g96_write() {
        let positions = vec![Vec3::new(1.0, 2.0, 3.0), Vec3::new(4.0, 5.0, 6.0)];

        let result = write_g96(
            "/tmp/test.g96",
            "Test structure",
            &positions,
            None,
            Some(Vec3::new(10.0, 10.0, 10.0)),
            None,
        );

        assert!(result.is_ok());
    }

    #[test]
    fn pertdata_is_written_as_gromosxx_does() {
        let path = std::env::temp_dir().join(format!("newen_pertdata_{}.cnf", std::process::id()));
        let mut w = G96Writer::new(&path).unwrap();
        w.write_pertdata(0.25).unwrap();
        w.close().unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(text, "PERTDATA\n      0.2500000\nEND\n");
        std::fs::remove_file(path).ok();
    }

    /// The block gromosXX wrote for a FORM 6 search (aladip, three states), byte for byte, and
    /// read back.
    #[test]
    fn aedssearch_is_written_as_gromosxx_does_and_reads_back() {
        let block = "AEDSSEARCH
    5.000000000
   -5.000000000
    0.000000000
              7
              0
              0
    0.494146008    -0.695482463     1.734757904               1             157     0.160923377    -0.491641558   109.985544149     0.839663944             250
    1.594612468    -1.315188147     3.280504045               0              18    -1.104058078    -2.490014747    13.835621966     0.902141960              58
   -2.088758476     0.330672389    -0.824803745               1             125    -0.976024228     1.080319529    14.158800520     0.337911046             170
END
";
        let lines: Vec<String> = block
            .lines()
            .skip(1)
            .take_while(|l| *l != "END")
            .map(|l| l.trim().to_string())
            .collect();
        let s = crate::coordinate::parse_aedssearch(&lines, 3).unwrap();
        assert_eq!((s.emaxcounts, s.oldstate, s.fullemin), (7, 0, false));
        assert_eq!(s.visited, vec![true, false, true]);
        assert_eq!(s.framecounts, vec![250, 58, 170]);
        assert!(!s.initialise);
        let path = std::env::temp_dir().join(format!("newen_aedss_{}.cnf", std::process::id()));
        let mut w = G96Writer::new(&path).unwrap();
        w.write_aedssearch(&s).unwrap();
        w.close().unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), block);
        std::fs::remove_file(path).ok();
        // a state line short, or a NUMSTATES that does not match the block
        assert!(crate::coordinate::parse_aedssearch(&lines, 2).is_err());
        let mut bad = lines.clone();
        bad[7] = "1.0 2.0".to_string();
        let e = crate::coordinate::parse_aedssearch(&bad, 3).unwrap_err();
        assert!(
            e.to_string().contains("bad line in AEDSSEARCH block"),
            "{e}"
        );
    }

    #[test]
    fn disresexpave_has_five_values_per_line() {
        let path =
            std::env::temp_dir().join(format!("newen_disresexpave_{}.cnf", std::process::id()));
        let mut w = G96Writer::new(&path).unwrap();
        w.write_disresexpave(&[32.304151, 39.9546641, 1.0, 2.0, 3.0, 4.0])
            .unwrap();
        w.close().unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(
            text,
            "DISRESEXPAVE\n     32.3041510     39.9546641      1.0000000      2.0000000      \
             3.0000000\n      4.0000000\nEND\n"
        );
        std::fs::remove_file(path).ok();
    }
}
