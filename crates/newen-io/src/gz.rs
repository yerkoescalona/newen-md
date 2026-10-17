//! Transparent gzip input, the way every GROMOS reader gets it.
//!
//! gromosXX opens each of its inputs through `io::igzstream` and gromos++ through
//! `gio::Ginstream`; both sit on zlib's `gzopen`, which reads an uncompressed file just as
//! happily as a compressed one. So the decision is made on the file's content — the two-byte
//! gzip magic `1f 8b` — not on its name: a `.trc` that happens to be compressed reads, and a
//! `.gz` name on a plain file does not lie to us.
//!
//! This matters for the analysis programs: `mk_script` appends `gzip <trajectory>` to every run
//! script it writes (`mk_script.cc:3853`), so a production run leaves `.trc.gz` / `.tre.gz`
//! behind and every analysis input in the GROMOS tutorials is compressed.

use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// The gzip magic number at the start of a member header (RFC 1952 §2.3.1).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Open `path` for line-oriented reading, decompressing it if it is gzipped.
///
/// `MultiGzDecoder` is used rather than `GzDecoder` because concatenated members — what
/// `cat a.trc.gz b.trc.gz > all.trc.gz` produces — are a single valid gzip stream that must
/// read as the concatenation of both, which is also what zlib's `gzread` does.
pub fn open_text<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    if n == 2 && magic == GZIP_MAGIC {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// A line reader with one line of pushback, over a possibly-gzipped file.
///
/// GROMOS block readers peek the next block's name and put the line back when the block is not
/// theirs (a `.trc` frame's `VELOCITY`/`FORCE`/`GENBOX` blocks are all optional). A plain file
/// can be rewound with `seek`; a gzip stream cannot, so the line is remembered instead.
pub struct TextReader {
    inner: Box<dyn BufRead>,
    pending: Option<String>,
}

impl TextReader {
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(TextReader {
            inner: open_text(path)?,
            pending: None,
        })
    }

    /// Read the next line into `buf` (cleared first); `Ok(0)` at end of file.
    pub fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        buf.clear();
        if let Some(line) = self.pending.take() {
            buf.push_str(&line);
            return Ok(buf.len());
        }
        self.inner.read_line(buf)
    }

    /// Put `line` back; the next `read_line` returns it again.
    pub fn unread_line(&mut self, line: &str) {
        self.pending = Some(line.to_string());
    }
}

/// Whether `path` holds a gzip stream (by its magic number, not its name).
pub fn is_gzipped<P: AsRef<Path>>(path: P) -> bool {
    File::open(path)
        .and_then(|mut f| {
            let mut magic = [0u8; 2];
            let n = f.read(&mut magic)?;
            Ok(n == 2 && magic == GZIP_MAGIC)
        })
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    fn tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("newen_gz_{}_{name}", std::process::id()))
    }

    #[test]
    fn reads_plain_and_gzipped_the_same() {
        let text = "TITLE\nhello\nEND\n";
        let plain = tmp("plain.txt");
        std::fs::write(&plain, text).unwrap();

        let gz = tmp("compressed.txt.gz");
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(text.as_bytes()).unwrap();
        std::fs::write(&gz, enc.finish().unwrap()).unwrap();

        for p in [&plain, &gz] {
            let lines: Vec<String> = open_text(p).unwrap().lines().map(|l| l.unwrap()).collect();
            assert_eq!(lines, ["TITLE", "hello", "END"], "{p:?}");
        }
        assert!(!is_gzipped(&plain));
        assert!(is_gzipped(&gz));
        std::fs::remove_file(plain).ok();
        std::fs::remove_file(gz).ok();
    }

    /// `cat a.gz b.gz` is one stream of two members; zlib reads both, so we must too.
    #[test]
    fn reads_concatenated_members() {
        let gz = tmp("multi.txt.gz");
        let mut bytes = Vec::new();
        for part in ["one\n", "two\n"] {
            let mut enc = GzEncoder::new(Vec::new(), Compression::default());
            enc.write_all(part.as_bytes()).unwrap();
            bytes.extend(enc.finish().unwrap());
        }
        std::fs::write(&gz, bytes).unwrap();
        let lines: Vec<String> = open_text(&gz)
            .unwrap()
            .lines()
            .map(|l| l.unwrap())
            .collect();
        assert_eq!(lines, ["one", "two"]);
        std::fs::remove_file(gz).ok();
    }
}
