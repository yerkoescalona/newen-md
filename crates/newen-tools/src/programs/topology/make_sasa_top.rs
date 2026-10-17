//! make_sasa_top — append a SASAPARAMETERS block (SASA implicit-solvent parameters per solute
//! atom) to a topology (gromos++ `make_sasa_top`).
//!
//! Usage: make_sasa_top @topo <top> @sasaspec <SASASPEC file> [@noH]
//!
//! The specification block lists `radius probability sigma n iac…` per line; every solute atom
//! (hydrogens — mass 1.008 — excluded with `@noH`) gets the parameters of its atom type. The
//! topology is written unchanged, followed by the new block.

use crate::programs::io_err;
use newen_io::args::{Arguments, Program, Status};
use newen_io::topology::{build_topology, read_topology_file, write_parsed_topology};
use std::collections::HashMap;
use std::io::Write;

pub const USAGE: &str = "# make_sasa_top
\t@topo       <molecular topology file>
\t@sasaspec   <sasa specification file>
\t[@noH       <do not include hydrogen atoms (default: include)";

const H_MASS: f64 = 1.008;

pub const KNOWN: &[&str] = &["topo", "sasaspec", "noH"];

pub const PROGRAM: Program = Program {
    name: "make_sasa_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_path = args.value("topo")?;
    let text = std::fs::read_to_string(topo_path).map_err(|e| format!("@topo {topo_path}: {e}"))?;
    let title: Vec<&str> = text
        .lines()
        .skip_while(|l| l.trim() != "TITLE")
        .skip(1)
        .take_while(|l| l.trim() != "END")
        .collect();
    let title = format!("{}\n", title.join("\n"));
    let parsed = read_topology_file(topo_path).map_err(|e| e.to_string())?;
    let no_h = args.count("noH") != -1;
    if args.count("sasaspec") != 1 {
        return Err("No sasa specification file".into());
    }
    let spec_path = args.value("sasaspec")?;
    let spec_text =
        std::fs::read_to_string(spec_path).map_err(|e| format!("@sasaspec {spec_path}: {e}"))?;
    let lines: Vec<&str> = spec_text
        .lines()
        .map(|l| l.split('#').next().unwrap_or("").trim())
        .filter(|l| !l.is_empty())
        .collect();
    let start = lines
        .iter()
        .position(|l| *l == "SASASPEC")
        .ok_or("sasa specification file does not contain a SASASPEC block!")?;
    let mut params: HashMap<usize, (f64, f64, f64)> = HashMap::new();
    for line in &lines[start + 1..] {
        if *line == "END" {
            break;
        }
        let tok: Vec<&str> = line.split_whitespace().collect();
        let num = |i: usize| -> Result<f64, String> {
            tok.get(i)
                .and_then(|t| t.parse::<f64>().ok())
                .ok_or_else(|| "bad line in SASASPEC block!".to_string())
        };
        let (radius, probability, sigma) = (num(0)?, num(1)?, num(2)?);
        let n = num(3)? as usize;
        for k in 0..n {
            let iac = tok
                .get(4 + k)
                .and_then(|t| t.parse::<usize>().ok())
                .ok_or_else(|| {
                    format!("bad line in SASASPEC block: could not read {n} IACs from line.")
                })?;
            params.insert(iac - 1, (radius, probability, sigma));
        }
    }
    let topo = build_topology(parsed.clone());
    let n_solute = topo.num_solute_atoms();
    let molecules: Vec<std::ops::Range<usize>> = topo
        .molecules
        .iter()
        .filter(|r| r.start < n_solute)
        .cloned()
        .collect();
    let include = |i: usize| !no_h || parsed.masses[i] != H_MASS;
    let count = molecules
        .iter()
        .flat_map(|r| r.clone())
        .filter(|&i| include(i))
        .count();
    let mut text = write_parsed_topology(&parsed, &title);
    text.push_str(&format!(
        "SASAPARAMETERS\n#NRSASAA\n{count}\n#ISASA    RADI      PI     SIGMAI\n"
    ));
    for r in &molecules {
        for (k, i) in r.clone().enumerate() {
            if !include(i) {
                continue;
            }
            let iac = parsed.iac[i];
            let (radius, probability, sigma) = params
                .get(&iac)
                .ok_or_else(|| format!("No SASA parameters for atom type: {}", iac + 1))?;
            // gromos++ leaves cout in fixed notation after the topology: three decimals
            text.push_str(&format!(
                "{:>6}{:>8.3}{:>8.3}{:>11.3}\n",
                k + 1,
                radius,
                probability,
                sigma
            ));
        }
    }
    text.push_str("END\n");
    write!(out, "{text}").map_err(io_err)?;
    Ok(())
}
