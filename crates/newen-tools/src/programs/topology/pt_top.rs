//! pt_top — apply a perturbation topology to a topology (gromos++ `pt_top`).
//!
//! Usage: pt_top @topo <top> @pttopo <ptp> @type TOPO [@npt <n>] [@firstatom <atom spec>]
//!
//! State B (the perturbation's second state) replaces atom IAC, mass and charge, atom-pair
//! exclusions / 1-4 pairs and the bonded types listed in the perturbation topology; the result
//! is written as a topology (`@type PERTTOPO` — writing the inverse perturbation — is not
//! implemented). `@npt` other than 1 (multiple perturbations) is not supported. The title
//! gains one line per perturbed molecule with its net charge in state A and B (1.8).
//! Typed verb: [`pt_top`].

use crate::programs::io_err;
use newen_core::topology::PairInteraction;
use newen_core::NetCharge;
use newen_io::args::{Arguments, Program, Status};
use newen_io::ptp::{read_pttopo, PerturbedTopology};
use newen_io::topology::{read_topology_file, write_parsed_topology, ParsedTopology};
use std::io::Write;

pub const USAGE: &str = "# pt_top
\t@topo      <molecular topology file>
\t@pttopo    <perturbation topology>
\t@type      <output format: TOPO, PERTTOPO>
\t@npt       <sequence number of the perturbation to apply, default 1 (state B)>
\t@firstatom <first atom to which the perturbation will be applied>";

pub const KNOWN: &[&str] = &["topo", "pttopo", "firstatom", "npt", "type"];

pub const PROGRAM: Program = Program {
    name: "pt_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let ptp_file = args.value("pttopo")?;
    let topo = read_topology_file(topo_file).map_err(|e| format!("@topo: {e}"))?;
    let pt = read_pttopo(ptp_file).map_err(|e| format!("@pttopo: {e}"))?;
    let npt: usize = args.get("npt", 1)?;
    if npt != 1 {
        return Err("@npt: only the first perturbation (state B) is supported".into());
    }
    let start: usize = match args.values("firstatom").first() {
        None => 0,
        Some(s) => {
            let sel = newen_core::selection::AtomSelection::from_string(
                s,
                &newen_io::topology::build_topology(
                    read_topology_file(topo_file).map_err(|e| e.to_string())?,
                ),
            )
            .map_err(|e| format!("@firstatom: {e}"))?;
            *sel.indices().first().ok_or("@firstatom selects nothing")?
        },
    };
    let kind = args
        .value("type")?
        .chars()
        .next()
        .map(|c| c.to_ascii_lowercase())
        .unwrap_or('t');
    if kind == 'p' {
        return Err("@type PERTTOPO is not implemented; TOPO writes the state-B topology".into());
    }
    if kind != 't' {
        return Err("Unkown output format.".into());
    }
    let (topo, charges) = pt_top(topo, &pt, start)?;
    let mut title = format!(
        "Topology based on\n{topo_file} and \n{ptp_file} (perturbation {npt}{})",
        if start > 0 {
            format!("; shifted to start at atomnumber {}", start + 1)
        } else {
            String::new()
        }
    );
    for (m, q) in charges {
        title += &format!("\n{}", net_charge_line(m, &q));
    }
    write!(out, "{}", write_parsed_topology(&topo, &title)).map_err(io_err)?;
    Ok(())
}

/// State B of `topo`: the perturbation `pt` applied with its atom 1 at atom `start + 1`
/// (IAC, mass and charge of every `PERTATOMPARAM` atom, the `PERTATOMPAIR` exclusions and 1-4
/// pairs, the perturbed bond, angle and dihedral types), and the net charge in A and B of each
/// molecule holding a perturbed atom (molecule index 0-based).
pub fn pt_top(
    mut topo: ParsedTopology,
    pt: &PerturbedTopology,
    start: usize,
) -> Result<(ParsedTopology, Vec<(usize, NetCharge)>), String> {
    let ps = &pt.perturbed_solute;
    let mut charge_a = topo.charges.clone();
    let mut perturbed = vec![false; topo.n_atoms];
    for a in &ps.atoms {
        let i = a.seq + start;
        if i >= topo.n_atoms {
            return Err(format!(
                "perturbed atom {} beyond the topology ({} atoms)",
                i + 1,
                topo.n_atoms
            ));
        }
        charge_a[i] = a.a_charge;
        perturbed[i] = true;
        topo.iac[i] = a.b_iac;
        topo.masses[i] = a.b_mass;
        topo.charges[i] = a.b_charge;
    }
    for p in &ps.atom_pairs {
        let (i, j) = (p.i + start, p.j + start);
        let excl = &mut topo.exclusions[i];
        let one4 = &mut topo.one_four_pairs[i];
        excl.retain(|&x| x != j);
        one4.retain(|&x| x != j);
        match p.b_type {
            PairInteraction::Excluded => excl.push(j),
            PairInteraction::OneFour => one4.push(j),
            PairInteraction::Normal => {},
        }
        excl.sort();
        one4.sort();
    }
    for b in &ps.bonds {
        for t in topo.bonds.iter_mut() {
            if t.0 == b.i + start && t.1 == b.j + start {
                t.2 = b.b_type;
            }
        }
    }
    for a in &ps.angles {
        for t in topo.angles.iter_mut() {
            if t.0 == a.i + start && t.1 == a.j + start && t.2 == a.k + start {
                t.3 = a.b_type;
            }
        }
    }
    for d in &ps.proper_dihedrals {
        for t in topo.proper_dihedrals.iter_mut() {
            if (t.0, t.1, t.2, t.3) == (d.i + start, d.j + start, d.k + start, d.l + start) {
                t.4 = d.b_type;
            }
        }
    }
    let mut molecules = Vec::new();
    let mut first = 0;
    for &last in &topo.solute_molecules {
        molecules.push(first..last);
        first = last;
    }
    let charges =
        NetCharge::of_perturbed_molecules(&molecules, &charge_a, &topo.charges, &perturbed);
    Ok((topo, charges))
}

/// `net charge of molecule <m>: A <q> B <q> (ΔQ = <q>)`, `m` 0-based, printed 1-based.
pub(crate) fn net_charge_line(m: usize, q: &NetCharge) -> String {
    format!("net charge of molecule {}: {q}", m + 1)
}
