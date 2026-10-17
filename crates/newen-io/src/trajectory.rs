//! GROMOS trajectory file writer (.trc/.trj)
//!
//! Writes coordinates (and optionally velocities and forces) at regular intervals.
//!
//! # Format
//! ```text
//! TITLE
//!   System trajectory
//! END
//! TIMESTEP
//!        1     0.000000
//! END
//! POSITIONRED
//!     1 RES    ATOM      1    x        y        z
//!     2 RES    ATOM      2    x        y        z
//! END
//! VELOCITYRED (optional)
//!     1 RES    ATOM      1    vx       vy       vz
//! END
//! FORCERED (optional, added in newen)
//!     1 RES    ATOM      1    fx       fy       fz
//! #         10  (comment every 10 atoms)
//! END
//! LATTICESHIFTS (optional)
//!     1    0    0    0
//! END
//! GENBOX
//!     lx   ly   lz
//! END
//! ```
//!
//! Note: For dedicated force trajectory output with FREEFORCERED/CONSFORCERED blocks,
//! use the `ForceWriter` from `io::force` module instead.

use crate::IoError;
use newen_core::configuration::Configuration;
use newen_core::math::Vec3;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// GROMOS trajectory writer
pub struct TrajectoryWriter {
    writer: BufWriter<File>,
    write_velocities: bool,
    write_forces: bool,
    step_count: usize,
}

impl TrajectoryWriter {
    /// Create a new trajectory writer
    pub fn new<P: AsRef<Path>>(
        path: P,
        title: &str,
        write_velocities: bool,
        write_forces: bool,
    ) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "TITLE")?;
        writeln!(writer, "  {}", title)?;
        writeln!(writer, "END")?;

        Ok(Self {
            writer,
            write_velocities,
            write_forces,
            step_count: 0,
        })
    }

    /// Write a trajectory frame from a Configuration's *current* state.
    pub fn write_frame(
        &mut self,
        step: usize,
        time: f64,
        config: &Configuration,
    ) -> Result<(), IoError> {
        self.write_state(step, time, config.current())
    }

    /// Write the frame a step's *energies* describe — the `old` buffer after the step's
    /// double-buffer swap. gromosXX writes the configuration and the energies of a step at the
    /// same point, so its trajectory starts at the input configuration and a frame's structure
    /// belongs with the energies of the same frame. Writing `current` instead put the two one
    /// step out of step with each other.
    pub fn write_frame_of_step(
        &mut self,
        step: usize,
        time: f64,
        config: &Configuration,
    ) -> Result<(), IoError> {
        self.write_state(step, time, config.old())
    }

    fn write_state(
        &mut self,
        step: usize,
        time: f64,
        state: &newen_core::configuration::State,
    ) -> Result<(), IoError> {
        let dims = state.box_config.dimensions();
        let box_opt = if dims.x > 0.0 { Some(dims) } else { None };
        self.write_trc_frame(step, time, &state.pos, box_opt)?;

        if self.write_velocities {
            writeln!(self.writer, "VELOCITYRED")?;
            for vel in &state.vel {
                writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", vel.x, vel.y, vel.z)?;
            }
            writeln!(self.writer, "END")?;
        }

        if self.write_forces {
            writeln!(self.writer, "FREEFORCERED")?;
            for force in &state.force {
                writeln!(
                    self.writer,
                    "{:18.9}{:18.9}{:18.9}",
                    force.x, force.y, force.z
                )?;
            }
            writeln!(self.writer, "END")?;
        }

        Ok(())
    }

    /// Write a velocity frame: `TIMESTEP` + `VELOCITYRED` (+ `GENBOX`), the layout gromosXX
    /// writes into the file named by `@trv`.
    pub fn write_velocity_frame(
        &mut self,
        step: usize,
        time: f64,
        velocities: &[Vec3],
        box_dims: Option<Vec3>,
    ) -> Result<(), IoError> {
        writeln!(self.writer, "TIMESTEP")?;
        writeln!(self.writer, "{:>15}{:20.9}", step, time)?;
        writeln!(self.writer, "END")?;
        writeln!(self.writer, "VELOCITYRED")?;
        for v in velocities {
            writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", v.x, v.y, v.z)?;
        }
        writeln!(self.writer, "END")?;
        if let Some(b) = box_dims {
            writeln!(self.writer, "GENBOX")?;
            writeln!(self.writer, "{:>8}", 1)?;
            writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", b.x, b.y, b.z)?;
            writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", 90.0, 90.0, 90.0)?;
            writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0)?;
            writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0)?;
            writeln!(self.writer, "END")?;
        }
        self.step_count += 1;
        Ok(())
    }

    /// Write a frame from raw positions (no Configuration needed).
    ///
    /// Writes TIMESTEP + POSITIONRED in standard 3-column GROMOS format,
    /// plus GENBOX if `box_dims` is Some.
    pub fn write_trc_frame(
        &mut self,
        step: usize,
        time: f64,
        positions: &[Vec3],
        box_dims: Option<Vec3>,
    ) -> Result<(), IoError> {
        writeln!(self.writer, "TIMESTEP")?;
        writeln!(self.writer, "{:>15}{:20.9}", step, time)?;
        writeln!(self.writer, "END")?;

        writeln!(self.writer, "POSITIONRED")?;
        for pos in positions {
            writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", pos.x, pos.y, pos.z)?;
        }
        writeln!(self.writer, "END")?;

        if let Some(b) = box_dims {
            if b.x > 0.0 {
                // gromosXX GENBOX: NTB (1 rectangular), lengths, angles, Euler angles, origin
                writeln!(self.writer, "GENBOX")?;
                writeln!(self.writer, "{:>8}", 1)?;
                writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", b.x, b.y, b.z)?;
                writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", 90.0, 90.0, 90.0)?;
                writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0)?;
                writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0)?;
                writeln!(self.writer, "END")?;
            }
        }

        self.step_count += 1;
        Ok(())
    }

    /// Flush buffered data to disk
    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush()?;
        Ok(())
    }

    /// Get number of frames written
    pub fn frame_count(&self) -> usize {
        self.step_count
    }
}

/// Write a single trajectory frame to a file (convenience function)
pub fn write_trajectory_frame<P: AsRef<Path>>(
    path: P,
    step: usize,
    time: f64,
    config: &Configuration,
    title: &str,
) -> Result<(), IoError> {
    let mut writer = TrajectoryWriter::new(path, title, false, false)?;
    writer.write_frame(step, time, config)?;
    writer.flush()?;
    Ok(())
}

/// A single frame from a trajectory file
#[derive(Debug, Clone)]
pub struct TrajectoryFrame {
    /// Simulation step number
    pub step: usize,
    /// Simulation time (ps)
    pub time: f64,
    /// Atomic positions (nm)
    pub positions: Vec<Vec3>,
    /// Atomic velocities (nm/ps), if present
    pub velocities: Option<Vec<Vec3>>,
    /// Atomic forces (kJ/(mol·nm)), if present
    pub forces: Option<Vec<Vec3>>,
    /// Lattice shifts (integer vectors), if present
    pub lattice_shifts: Option<Vec<(i32, i32, i32)>>,
    /// Box dimensions (nm)
    pub box_dims: Vec3,
}

/// GROMOS trajectory reader (.trc/.trj files)
///
/// Reads trajectory files frame by frame, parsing POSITIONRED, VELOCITYRED, FORCERED,
/// LATTICESHIFTS, and GENBOX blocks.
pub struct TrajectoryReader {
    reader: crate::gz::TextReader,
    title: String,
    frames_read: usize,
    buffer: String,
}

impl TrajectoryReader {
    /// Open a trajectory file for reading
    ///
    /// # Arguments
    /// * `path` - Path to the .trc or .trj file
    ///
    /// # Returns
    /// * `Ok(TrajectoryReader)` - Successfully opened the file
    /// * `Err(IoError)` - Failed to open or parse the header
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, IoError> {
        let mut reader = crate::gz::TextReader::open(path)?;
        let mut buffer = String::new();

        // Read TITLE block
        let title = Self::read_title_block(&mut reader, &mut buffer)?;

        Ok(Self {
            reader,
            title,
            frames_read: 0,
            buffer,
        })
    }

    /// Get the title from the trajectory file
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Get the number of frames read so far
    pub fn frames_read(&self) -> usize {
        self.frames_read
    }

    /// Read the next frame from the trajectory
    ///
    /// # Returns
    /// * `Ok(Some(TrajectoryFrame))` - Successfully read a frame
    /// * `Ok(None)` - End of file reached
    /// * `Err(IoError)` - Error reading the frame
    pub fn read_frame(&mut self) -> Result<Option<TrajectoryFrame>, IoError> {
        // Try to read TIMESTEP block; a configuration file may start its frame with the
        // position block directly (no TIMESTEP), as gromos++'s InG96 accepts.
        let mut header_seen = false;
        let (step, time) = match Self::read_timestep_block(
            &mut self.reader,
            &mut self.buffer,
            &mut header_seen,
        )? {
            Some(st) => st,
            None => return Ok(None), // End of file
        };

        // Read POSITION/POSITIONRED block (required)
        let positions = Self::read_position_block(&mut self.reader, &mut self.buffer, header_seen)?;

        // Try to read optional blocks
        let velocities = Self::try_read_velocity_block(&mut self.reader, &mut self.buffer)?;
        let forces = Self::try_read_force_block(&mut self.reader, &mut self.buffer)?;
        let lattice_shifts = Self::try_read_lattice_shifts(&mut self.reader, &mut self.buffer)?;

        // Read GENBOX block (optional — vacuum systems may omit it)
        let box_dims = Self::try_read_genbox_block(&mut self.reader, &mut self.buffer)?;

        self.frames_read += 1;

        Ok(Some(TrajectoryFrame {
            step,
            time,
            positions,
            velocities,
            forces,
            lattice_shifts,
            box_dims,
        }))
    }

    /// Read all frames from the trajectory
    ///
    /// # Returns
    /// * `Ok(Vec<TrajectoryFrame>)` - All frames successfully read
    /// * `Err(IoError)` - Error reading frames
    pub fn read_all_frames(&mut self) -> Result<Vec<TrajectoryFrame>, IoError> {
        let mut frames = Vec::new();
        while let Some(frame) = self.read_frame()? {
            frames.push(frame);
        }
        Ok(frames)
    }

    // Helper functions for reading blocks

    fn read_title_block(
        reader: &mut crate::gz::TextReader,
        buffer: &mut String,
    ) -> Result<String, IoError> {
        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("TITLE") {
            return Err(IoError::FormatError(
                "Expected TITLE block at start of trajectory".to_string(),
            ));
        }

        // The title runs until END and may span several lines (gromosXX writes a blank line
        // and "final configuration" into its TITLE block).
        let mut title_lines = Vec::new();
        loop {
            buffer.clear();
            if reader.read_line(buffer)? == 0 {
                return Err(IoError::FormatError("Expected END after TITLE".to_string()));
            }
            if buffer.trim() == "END" {
                break;
            }
            title_lines.push(buffer.trim().to_string());
        }
        let title = title_lines.join("\n").trim().to_string();

        Ok(title)
    }

    fn read_timestep_block(
        reader: &mut crate::gz::TextReader,
        buffer: &mut String,
        position_header_seen: &mut bool,
    ) -> Result<Option<(usize, f64)>, IoError> {
        buffer.clear();
        let bytes_read = reader.read_line(buffer)?;
        if bytes_read == 0 {
            return Ok(None); // EOF
        }

        let line = buffer.trim();
        if line.is_empty() {
            return Self::read_timestep_block(reader, buffer, position_header_seen);
            // Skip empty lines
        }

        if line.starts_with("POSITION") {
            // a frame without TIMESTEP: step 0, time 0 (gromos++ InG96)
            *position_header_seen = true;
            return Ok(Some((0, 0.0)));
        }

        if !line.starts_with("TIMESTEP") {
            // some other block (a configuration file's STOCHINT, PERTDATA, … after GENBOX): skip it
            loop {
                buffer.clear();
                if reader.read_line(buffer)? == 0 {
                    return Ok(None);
                }
                if buffer.trim() == "END" {
                    break;
                }
            }
            return Self::read_timestep_block(reader, buffer, position_header_seen);
        }

        // Read step and time
        buffer.clear();
        reader.read_line(buffer)?;
        let parts: Vec<&str> = buffer.split_whitespace().collect();
        if parts.len() < 2 {
            return Err(IoError::FormatError(
                "TIMESTEP data should have step and time".to_string(),
            ));
        }

        let step = parts[0]
            .parse::<usize>()
            .map_err(|e| IoError::FormatError(format!("Invalid step number: {}", e)))?;
        let time = parts[1]
            .parse::<f64>()
            .map_err(|e| IoError::FormatError(format!("Invalid time: {}", e)))?;

        // Read END
        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("END") {
            return Err(IoError::FormatError(
                "Expected END after TIMESTEP".to_string(),
            ));
        }

        Ok(Some((step, time)))
    }

    fn read_position_block(
        reader: &mut crate::gz::TextReader,
        buffer: &mut String,
        header_seen: bool,
    ) -> Result<Vec<Vec3>, IoError> {
        if !header_seen {
            buffer.clear();
            reader.read_line(buffer)?;
            if !buffer.trim().starts_with("POSITION") {
                return Err(IoError::FormatError(
                    "Expected POSITION or POSITIONRED block".to_string(),
                ));
            }
        }

        let mut positions = Vec::new();
        loop {
            buffer.clear();
            reader.read_line(buffer)?;
            let line = buffer.trim();

            if line.starts_with("END") {
                break;
            }

            if line.starts_with('#') || line.is_empty() {
                continue; // Skip comments and empty lines
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            // Standard GROMOS POSITIONRED: 3 columns (x y z)
            // Legacy 7-column format: atom_num RES ATOM serial x y z
            let (xi, yi, zi) = if parts.len() == 3 {
                (0, 1, 2)
            } else {
                (4, 5, 6)
            };
            if parts.len() > zi {
                let x = parts[xi]
                    .parse::<f64>()
                    .map_err(|e| IoError::FormatError(format!("Invalid x: {e}")))?;
                let y = parts[yi]
                    .parse::<f64>()
                    .map_err(|e| IoError::FormatError(format!("Invalid y: {e}")))?;
                let z = parts[zi]
                    .parse::<f64>()
                    .map_err(|e| IoError::FormatError(format!("Invalid z: {e}")))?;
                positions.push(Vec3::new(x, y, z));
            }
        }

        Ok(positions)
    }

    fn try_read_velocity_block(
        reader: &mut crate::gz::TextReader,
        buffer: &mut String,
    ) -> Result<Option<Vec<Vec3>>, IoError> {
        if reader.read_line(buffer)? == 0 {
            return Ok(None);
        }
        if !buffer.trim().starts_with("VELOCITY") {
            // not our block — put the line back for the next reader
            reader.unread_line(&buffer.clone());
            return Ok(None);
        }

        let mut velocities = Vec::new();
        loop {
            buffer.clear();
            reader.read_line(buffer)?;
            let line = buffer.trim();

            if line.starts_with("END") {
                break;
            }

            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            // reduced (3 columns) or labelled (7 columns) block
            let (xi, yi, zi) = if parts.len() == 3 {
                (0, 1, 2)
            } else {
                (4, 5, 6)
            };
            if parts.len() > zi {
                let vx = parts[xi]
                    .parse::<f64>()
                    .map_err(|e| IoError::FormatError(format!("Invalid vx velocity: {}", e)))?;
                let vy = parts[yi]
                    .parse::<f64>()
                    .map_err(|e| IoError::FormatError(format!("Invalid vy velocity: {}", e)))?;
                let vz = parts[zi]
                    .parse::<f64>()
                    .map_err(|e| IoError::FormatError(format!("Invalid vz velocity: {}", e)))?;
                velocities.push(Vec3::new(vx, vy, vz));
            }
        }

        Ok(Some(velocities))
    }

    fn try_read_force_block(
        reader: &mut crate::gz::TextReader,
        buffer: &mut String,
    ) -> Result<Option<Vec<Vec3>>, IoError> {
        if reader.read_line(buffer)? == 0 {
            return Ok(None);
        }
        if !buffer.trim().starts_with("FORCE") {
            // not our block — put the line back for the next reader
            reader.unread_line(&buffer.clone());
            return Ok(None);
        }

        let mut forces = Vec::new();
        loop {
            buffer.clear();
            reader.read_line(buffer)?;
            let line = buffer.trim();

            if line.starts_with("END") {
                break;
            }

            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            // reduced (3 columns) or labelled (7 columns) block
            let (xi, yi, zi) = if parts.len() == 3 {
                (0, 1, 2)
            } else {
                (4, 5, 6)
            };
            if parts.len() > zi {
                let fx = parts[xi]
                    .parse::<f64>()
                    .map_err(|e| IoError::FormatError(format!("Invalid fx force: {}", e)))?;
                let fy = parts[yi]
                    .parse::<f64>()
                    .map_err(|e| IoError::FormatError(format!("Invalid fy force: {}", e)))?;
                let fz = parts[zi]
                    .parse::<f64>()
                    .map_err(|e| IoError::FormatError(format!("Invalid fz force: {}", e)))?;
                forces.push(Vec3::new(fx, fy, fz));
            }
        }

        Ok(Some(forces))
    }

    fn try_read_lattice_shifts(
        reader: &mut crate::gz::TextReader,
        buffer: &mut String,
    ) -> Result<Option<Vec<(i32, i32, i32)>>, IoError> {
        if reader.read_line(buffer)? == 0 {
            return Ok(None);
        }
        if !buffer.trim().starts_with("LATTICESHIFTS") {
            // not our block — put the line back for the next reader
            reader.unread_line(&buffer.clone());
            return Ok(None);
        }

        let mut shifts = Vec::new();
        loop {
            buffer.clear();
            reader.read_line(buffer)?;
            let line = buffer.trim();

            if line.starts_with("END") {
                break;
            }

            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                let sx = parts[1]
                    .parse::<i32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid shift x: {}", e)))?;
                let sy = parts[2]
                    .parse::<i32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid shift y: {}", e)))?;
                let sz = parts[3]
                    .parse::<i32>()
                    .map_err(|e| IoError::FormatError(format!("Invalid shift z: {}", e)))?;
                shifts.push((sx, sy, sz));
            }
        }

        Ok(Some(shifts))
    }

    // TODO: wire up when GENBOX block reading is needed for box reshaping
    #[allow(dead_code)]
    fn read_genbox_block(
        reader: &mut crate::gz::TextReader,
        buffer: &mut String,
    ) -> Result<Vec3, IoError> {
        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("GENBOX") {
            return Err(IoError::FormatError("Expected GENBOX block".to_string()));
        }

        buffer.clear();
        reader.read_line(buffer)?;
        let parts: Vec<&str> = buffer.split_whitespace().collect();
        if parts.len() < 3 {
            return Err(IoError::FormatError(
                "GENBOX should have 3 dimensions".to_string(),
            ));
        }

        let lx = parts[0]
            .parse::<f64>()
            .map_err(|e| IoError::FormatError(format!("Invalid box x: {}", e)))?;
        let ly = parts[1]
            .parse::<f64>()
            .map_err(|e| IoError::FormatError(format!("Invalid box y: {}", e)))?;
        let lz = parts[2]
            .parse::<f64>()
            .map_err(|e| IoError::FormatError(format!("Invalid box z: {}", e)))?;

        buffer.clear();
        reader.read_line(buffer)?;
        if !buffer.trim().starts_with("END") {
            return Err(IoError::FormatError(
                "Expected END after GENBOX".to_string(),
            ));
        }

        Ok(Vec3::new(lx, ly, lz))
    }

    /// Try to read GENBOX; return Vec3::ZERO if missing (vacuum trajectory).
    fn try_read_genbox_block(
        reader: &mut crate::gz::TextReader,
        buffer: &mut String,
    ) -> Result<Vec3, IoError> {
        if reader.read_line(buffer)? == 0 {
            return Ok(Vec3::ZERO);
        }
        if !buffer.trim().starts_with("GENBOX") {
            // not GENBOX — put the line back and report a zero box (vacuum)
            reader.unread_line(&buffer.clone());
            return Ok(Vec3::ZERO);
        }
        // gromosXX writes NTB, then the box lengths, angles, Euler angles and origin; older
        // newen files carried the lengths alone. Read up to END either way.
        let mut lengths: Option<Vec3> = None;
        let mut first = true;
        loop {
            buffer.clear();
            if reader.read_line(buffer)? == 0 {
                break;
            }
            let t = buffer.trim();
            if t == "END" {
                break;
            }
            if t.is_empty() || t.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = t.split_whitespace().collect();
            if first {
                first = false;
                if parts.len() == 1 {
                    continue; // NTB line
                }
            }
            if lengths.is_none() && parts.len() >= 3 {
                let v: Vec<f64> = parts
                    .iter()
                    .take(3)
                    .map(|x| x.parse::<f64>().unwrap_or(0.0))
                    .collect();
                lengths = Some(Vec3::new(v[0], v[1], v[2]));
            }
        }
        Ok(lengths.unwrap_or(Vec3::ZERO))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use newen_core::math::Vec3;

    fn round_trip(positions: Vec<Vec3>, box_dims: Option<Vec3>) -> TrajectoryFrame {
        let tmp =
            std::env::temp_dir().join(format!("newen_trc_test_{}_{}.trc", std::process::id(), {
                static C: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
                C.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
            }));
        {
            let mut w = TrajectoryWriter::new(&tmp, "test", false, false).unwrap();
            w.write_trc_frame(42, 0.084, &positions, box_dims).unwrap();
            w.flush().unwrap();
        }
        let mut r = TrajectoryReader::new(&tmp).unwrap();
        let frame = r.read_frame().unwrap().expect("expected one frame");
        std::fs::remove_file(&tmp).ok();
        frame
    }

    #[test]
    fn write_trc_frame_produces_3column_positionred() {
        let tmp = std::env::temp_dir().join(format!("newen_fmt_{}_{}.trc", std::process::id(), {
            static C: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
            C.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        }));
        {
            let mut w = TrajectoryWriter::new(&tmp, "fmt test", false, false).unwrap();
            w.write_trc_frame(1, 0.002, &[Vec3::new(1.0, 2.0, 3.0)], None)
                .unwrap();
            w.flush().unwrap();
        }
        let content = std::fs::read_to_string(&tmp).unwrap();
        std::fs::remove_file(&tmp).ok();

        // POSITIONRED line must have exactly 3 floats (standard GROMOS format)
        let pos_line = content
            .lines()
            .skip_while(|l| !l.starts_with("POSITIONRED"))
            .nth(1)
            .unwrap()
            .trim()
            .to_string();
        let cols: Vec<&str> = pos_line.split_whitespace().collect();
        assert_eq!(
            cols.len(),
            3,
            "POSITIONRED must have 3 columns, got: {pos_line}"
        );
        assert!((cols[0].parse::<f64>().unwrap() - 1.0).abs() < 1e-9);
        assert!((cols[1].parse::<f64>().unwrap() - 2.0).abs() < 1e-9);
        assert!((cols[2].parse::<f64>().unwrap() - 3.0).abs() < 1e-9);
    }

    #[test]
    fn round_trip_vacuum() {
        let positions = vec![
            Vec3::new(0.1, 0.2, 0.3),
            Vec3::new(1.5, -0.5, 2.1),
            Vec3::new(0.0, 0.0, 0.0),
        ];
        let frame = round_trip(positions.clone(), None);
        assert_eq!(frame.step, 42);
        assert!((frame.time - 0.084).abs() < 1e-9);
        assert_eq!(frame.positions.len(), 3);
        for (a, b) in positions.iter().zip(&frame.positions) {
            assert!((a.x - b.x).abs() < 1e-9, "x mismatch: {a:?} vs {b:?}");
            assert!((a.y - b.y).abs() < 1e-9);
            assert!((a.z - b.z).abs() < 1e-9);
        }
    }

    #[test]
    fn round_trip_with_box() {
        let positions = vec![Vec3::new(0.5, 0.5, 0.5), Vec3::new(1.0, 1.0, 1.0)];
        let box_dims = Vec3::new(2.0, 2.0, 2.0);
        let frame = round_trip(positions.clone(), Some(box_dims));
        assert!((frame.box_dims.x - 2.0).abs() < 1e-9);
        assert!((frame.box_dims.y - 2.0).abs() < 1e-9);
        assert!((frame.box_dims.z - 2.0).abs() < 1e-9);
        for (a, b) in positions.iter().zip(&frame.positions) {
            assert!((a.x - b.x).abs() < 1e-9);
        }
    }

    #[test]
    fn reads_3column_gromos_format() {
        // Standard GROMOS .trc with 3-column POSITIONRED (as written by gromosXX)
        let tmp = std::env::temp_dir().join(format!("newen_3col_{}_{}.trc", std::process::id(), {
            static C: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
            C.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        }));
        std::fs::write(
            &tmp,
            "TITLE\n  test\nEND\nTIMESTEP\n              5    0.010000000\nEND\n\
             POSITIONRED\n    1.234567890    2.345678901    3.456789012\nEND\n\
             GENBOX\n    4.000000000    4.000000000    4.000000000\nEND\n",
        )
        .unwrap();
        let mut r = TrajectoryReader::new(&tmp).unwrap();
        let frame = r.read_frame().unwrap().unwrap();
        std::fs::remove_file(&tmp).ok();
        assert_eq!(frame.step, 5);
        assert!((frame.time - 0.01).abs() < 1e-9);
        assert_eq!(frame.positions.len(), 1);
        assert!((frame.positions[0].x - 1.234567890).abs() < 1e-9);
        assert!((frame.box_dims.x - 4.0).abs() < 1e-9);
    }

    #[test]
    fn reads_7column_legacy_format() {
        // Legacy 7-column format (written by older newen)
        let tmp = std::env::temp_dir().join(format!("newen_7col_{}_{}.trc", std::process::id(), {
            static C: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
            C.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        }));
        std::fs::write(
            &tmp,
            "TITLE\n  legacy\nEND\nTIMESTEP\n              1    0.002000000\nEND\n\
             POSITIONRED\n     1    RES   ATOM      1    1.100000    2.200000    3.300000\nEND\n\
             GENBOX\n    3.000000000    3.000000000    3.000000000\nEND\n",
        )
        .unwrap();
        let mut r = TrajectoryReader::new(&tmp).unwrap();
        let frame = r.read_frame().unwrap().unwrap();
        std::fs::remove_file(&tmp).ok();
        assert!((frame.positions[0].x - 1.1).abs() < 1e-6);
        assert!((frame.positions[0].y - 2.2).abs() < 1e-6);
        assert!((frame.positions[0].z - 3.3).abs() < 1e-6);
    }

    #[test]
    fn multi_frame_round_trip() {
        let tmp =
            std::env::temp_dir().join(format!("newen_multi_{}_{}.trc", std::process::id(), {
                static C: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
                C.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
            }));
        let frames_in = vec![
            (0usize, 0.0_f64, vec![Vec3::new(0.1, 0.2, 0.3)]),
            (1, 0.002, vec![Vec3::new(0.4, 0.5, 0.6)]),
            (2, 0.004, vec![Vec3::new(0.7, 0.8, 0.9)]),
        ];
        {
            let mut w = TrajectoryWriter::new(&tmp, "multi", false, false).unwrap();
            for (step, time, ref pos) in &frames_in {
                w.write_trc_frame(*step, *time, pos, None).unwrap();
            }
            w.flush().unwrap();
        }
        let mut r = TrajectoryReader::new(&tmp).unwrap();
        let frames_out = r.read_all_frames().unwrap();
        std::fs::remove_file(&tmp).ok();
        assert_eq!(frames_out.len(), 3);
        for (i, ((_step, time, pos), frame)) in frames_in.iter().zip(&frames_out).enumerate() {
            assert!((frame.time - time).abs() < 1e-9, "frame {i} time mismatch");
            assert!(
                (frame.positions[0].x - pos[0].x).abs() < 1e-9,
                "frame {i} x mismatch"
            );
        }
    }
}

// (merged into the tests module above)
