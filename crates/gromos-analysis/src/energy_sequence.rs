//! A list of `.tre` or `.trg` files read as one stream of frames — gromos++'s `@en_files` /
//! `@fr_files` loop (`ene_ana.cc`, `ext_ti_ana.cc`): when a file runs out, the next one is opened,
//! its preamble read and its layout bound against the library before reading on.

use gromos_io::energy_traj::{read_preamble, EnergyTraj};
use gromos_io::gz::TextReader;

/// The files of one kind (`"ENERTRJ"` or `"FRENERTRJ"`), in the order given.
pub struct EnergySequence {
    files: Vec<String>,
    next: usize,
    reader: Option<TextReader>,
    kind: &'static str,
    lib_name: String,
}

impl EnergySequence {
    /// `kind` is the block a frame starts with; `lib_name` only labels messages.
    pub fn new(files: Vec<String>, kind: &'static str, lib_name: &str) -> Self {
        EnergySequence {
            files,
            next: 0,
            reader: None,
            kind,
            lib_name: lib_name.to_string(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.files.is_empty()
    }

    /// Read the next frame into `etrj`, opening the next file when the current one runs out;
    /// `Ok(false)` once every file is exhausted. Reader warnings go to stderr.
    pub fn read_frame(&mut self, etrj: &mut EnergyTraj) -> Result<bool, String> {
        loop {
            if self.reader.is_none() {
                let Some(path) = self.files.get(self.next) else {
                    return Ok(false);
                };
                self.next += 1;
                let mut reader = TextReader::open(path)
                    .map_err(|e| format!("{}: cannot open {path}: {e}", self.kind))?;
                let preamble = read_preamble(&mut reader).map_err(|e| format!("{path}: {e}"))?;
                etrj.bind(self.kind, &preamble, &self.lib_name, path)?;
                for w in etrj.drain_warnings() {
                    eprintln!("{w}");
                }
                self.reader = Some(reader);
            }
            let reader = self.reader.as_mut().unwrap();
            let file = &self.files[self.next - 1];
            let more = etrj
                .read_frame(reader, self.kind)
                .map_err(|e| format!("{file}: {e}"))?;
            for w in etrj.drain_warnings() {
                eprintln!("{w}");
            }
            if more {
                return Ok(true);
            }
            self.reader = None;
        }
    }
}
