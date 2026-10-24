//! Read a GROMOS `.imd` and write it back through `write_imd` — the tool
//! `scripts/roundtrip_imd_gromosXX.py` uses to prove gromosXX accepts what newen writes
//! (3.9 assumption A7).
//!
//! Usage: `imd_roundtrip <input.imd> <output.imd> [n_atoms]`

use std::process::exit;

use newen_io::imd::read_imd_file;
use newen_io::imd_write::write_imd_file;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: imd_roundtrip <input.imd> <output.imd> [n_atoms]");
        exit(2);
    }
    let n_atoms = args.get(3).and_then(|s| s.parse::<usize>().ok());
    let params = match read_imd_file(&args[1]) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("error reading {}: {e}", args[1]);
            exit(1);
        },
    };
    if let Err(e) = write_imd_file(&args[2], &params, n_atoms) {
        eprintln!("error writing {}: {e}", args[2]);
        exit(1);
    }
}
