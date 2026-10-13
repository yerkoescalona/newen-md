//! Write GROMOS `.imd` input files — the inverse of [`crate::imd::read_imd_file`].
//!
//! Every block the parser models is regenerated from `ImdParameters` fields in the positional
//! GROMOS layout (comment line with the field names, then the values), in the block order the
//! reference inputs use. Blocks the parser does **not** model but found in the file
//! (`raw_blocks` entries such as GAMD, QMMM, …) are re-emitted verbatim after the modelled
//! ones, so a read → write round trip loses nothing.
//!
//! Two rules learnt from feeding the output to the real gromosXX (3.9 A7):
//! - values are always separated by whitespace, never packed into fixed-width columns — a
//!   number wider than its column (`1e-10` prints as `0.0000000001`) would fuse with its
//!   neighbour and gromosXX would read `20.0000000001` for `NA2CLC`;
//! - optional blocks whose content is "off" are **not** written unless the file that was read
//!   carried them: gromosXX validates a DISTANCERES block's `DIR0`/`TAUDIR` even with `NTDIR=0`,
//!   and an absent block is what a recipe built in memory means anyway.
//!
//! Numbers are printed with Rust's shortest round-trip formatting (`{}`), so
//! `parse(write(p)) == p` holds exactly for every finite value.
//!
//! `n_atoms` is needed only for blocks whose GROMOS layout carries atom counts that a recipe
//! built in memory does not know (`MULTIBATH` DOFSET `LAST`, `FORCE` `NRE`); a file that was
//! read keeps its own values and ignores it.

use std::fs;
use std::path::Path;

use crate::imd::{ImdParameters, PressureParameters, TempBathParameters};
use crate::IoError;

/// Names of the blocks the parser/writer model; anything else in `raw_blocks` is passthrough.
pub const MODELLED_BLOCKS: &[&str] = &[
    "TITLE",
    "SYSTEM",
    "INITIALISE",
    "STEP",
    "BOUNDCOND",
    "MULTIBATH",
    "PRESSURESCALE",
    "COMTRANSROT",
    "PRINTOUT",
    "WRITETRAJ",
    "CONSTRAINT",
    "COVALENTFORM",
    "FORCE",
    "PAIRLIST",
    "NONBONDED",
    "POSITIONRES",
    "ROTTRANS",
    "DISTANCERES",
    "PERTURBATION",
    "LAMBDAS",
    "PRECALCLAM",
    "ENERGYMIN",
    "EDS",
];

/// One data line: indented, values separated by two spaces.
fn row(o: &mut String, values: &[String]) {
    o.push_str("    ");
    o.push_str(&values.join("  "));
    o.push('\n');
}

macro_rules! vals {
    ($($v:expr),* $(,)?) => { &[$($v.to_string()),*][..] };
}

fn block(out: &mut String, name: &str, body: impl FnOnce(&mut String)) {
    out.push_str(name);
    out.push('\n');
    body(out);
    out.push_str("END\n");
}

/// Render `params` as an `.imd` file. See the module docs for what `n_atoms` is for.
pub fn write_imd(params: &ImdParameters, n_atoms: Option<usize>) -> String {
    let mut out = String::new();
    let p = params;
    let present = |name: &str| p.raw_blocks.contains_key(name);

    block(&mut out, "TITLE", |o| {
        if p.title.is_empty() {
            o.push_str("GROMOS simulation\n");
        }
        for line in p.title.lines() {
            o.push_str(line);
            o.push('\n');
        }
    });

    block(&mut out, "SYSTEM", |o| {
        o.push_str("#      NPM      NSM\n");
        row(o, vals![p.npm, p.nsm]);
    });

    block(&mut out, "INITIALISE", |o| {
        o.push_str("#  NTIVEL  NTISHK  NTINHT  NTINHB\n");
        row(o, vals![p.ntivel, p.ntishk, p.ntinht, p.ntinhb]);
        o.push_str("#  NTISHI  NTIRTC  NTICOM\n");
        row(o, vals![p.ntishi, p.ntirtc, p.nticom]);
        o.push_str("#  NTISTI\n");
        row(o, vals![p.ntisti]);
        o.push_str("#      IG   TEMPI\n");
        row(o, vals![p.ig, p.tempi]);
    });

    block(&mut out, "STEP", |o| {
        o.push_str("#   NSTLIM         T        DT\n");
        row(o, vals![p.nstlim, p.t0, p.dt]);
    });

    block(&mut out, "BOUNDCOND", |o| {
        o.push_str("#      NTB    NDFMIN\n");
        row(o, vals![p.ntb, p.ndfmin]);
    });

    // MULTIBATH: absent means "no temperature coupling" in gromosXX — only written when a
    // bath exists. A bath without DOF sets gets one covering every atom.
    if let Some(bath) = p.temp_bath.first() {
        block(&mut out, "MULTIBATH", |o| write_multibath(o, bath, n_atoms));
    }

    if let Some(pp) = &p.pressure_parameters {
        block(&mut out, "PRESSURESCALE", |o| {
            write_pressurescale(o, p.couple_pressure, pp)
        });
    }

    // Absent COMTRANSROT == NSCM 0 (gromosXX parameter.h); written when read or non-default.
    if present("COMTRANSROT") || p.nscm != 0 {
        block(&mut out, "COMTRANSROT", |o| {
            o.push_str("#     NSCM\n");
            row(o, vals![p.nscm]);
        });
    }

    block(&mut out, "PRINTOUT", |o| {
        o.push_str("#     NTPR      NTPP\n");
        row(o, vals![p.ntpr, p.ntpp]);
    });

    block(&mut out, "WRITETRAJ", |o| {
        o.push_str("#     NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB\n");
        row(
            o,
            vals![p.ntwx, p.ntwse, p.ntwv, p.ntwf, p.ntwe, p.ntwg, p.ntwb],
        );
    });

    // Only written when it says something other than the GROMOS defaults, as gromosXX's own
    // inputs do: an absent COVALENTFORM means quartic bonds, cosine-harmonic angles, arbitrary
    // dihedral phase shifts.
    if p.ntbbh != 0 || p.ntbah != 0 || p.ntbdn != 0 {
        block(&mut out, "COVALENTFORM", |o| {
            o.push_str("#   NTBBH    NTBAH    NTBDN\n");
            row(o, vals![p.ntbbh, p.ntbah, p.ntbdn]);
        });
    }

    block(&mut out, "CONSTRAINT", |o| {
        o.push_str("#       NTC\n");
        row(o, vals![p.ntc]);
        o.push_str("#       NTCP\n");
        row(o, vals![constraint_algorithm(p.ntcp)]);
        // SHAKE takes a tolerance, LINCS an expansion order; "off" and SETTLE take no parameter,
        // and gromosXX reads none for them — writing one would shift every value after it.
        if p.ntcp == 2 {
            o.push_str("#       NTCP0(1)\n");
            row(o, vals![p.lincs_order_solute]);
        } else if p.ntcp == 1 {
            o.push_str("#       NTCP0(1)\n");
            row(o, vals![p.shake_tol]);
        }
        o.push_str("#       NTCS\n");
        row(o, vals![constraint_algorithm(p.ntcs)]);
        // SETTLE and "off" take no NTCS0 parameter — the reference inputs omit the line.
        if p.ntcs == 2 {
            o.push_str("#       NTCS0(1)\n");
            row(o, vals![p.lincs_order_solvent]);
        } else if p.ntcs == 1 {
            o.push_str("#       NTCS0(1)\n");
            row(o, vals![p.ntcs0]);
        }
    });

    block(&mut out, "FORCE", |o| {
        o.push_str("#      NTF array\n");
        o.push_str("# bonds    angles    imp.     dihe     charge nonbonded\n");
        row(
            o,
            vals![p.ntf[0], p.ntf[1], p.ntf[2], p.ntf[3], p.ntf[4], p.ntf[5]],
        );
        o.push_str("# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)\n");
        if p.nre.is_empty() {
            // One energy group spanning the whole system, if we know its size.
            match n_atoms {
                Some(n) => row(o, vals![1, n]),
                None => row(o, vals![p.negr.max(1)]),
            }
        } else {
            let mut line = vec![p.negr.to_string()];
            line.extend(p.nre.iter().map(|n| n.to_string()));
            row(o, &line);
        }
    });

    block(&mut out, "PAIRLIST", |o| {
        o.push_str("#       ALGORITHM       NSNB    RCUTP   RCUTL   SIZE    TYPE\n");
        let algorithm = match p.algorithm {
            0 => "standard".to_string(),
            1 => "grid".to_string(),
            2 => "grid_cell".to_string(),
            other => other.to_string(),
        };
        let size = if p.size == 0.0 {
            "auto".to_string()
        } else {
            p.size.to_string()
        };
        let type_ = match p.type_ {
            0 => "chargegroup".to_string(),
            1 => "atomic".to_string(),
            other => other.to_string(),
        };
        row(o, vals![algorithm, p.nsnb, p.rcutp, p.rcutl, size, type_]);
    });

    block(&mut out, "NONBONDED", |o| {
        let x = &p.nonbonded_extra;
        o.push_str("#   NLRELE\n");
        row(o, vals![p.nlrele]);
        o.push_str("#    APPAK      RCRF     EPSRF  NSLFEXCL\n");
        row(o, vals![p.appak, p.rcrf, p.epsrf, p.nslfexcl]);
        o.push_str("#   NSHAPE    ASHAPE    NA2CLC     TOLA2    EPSLS\n");
        row(o, vals![x.nshape, x.ashape, x.na2clc, x.tola2, x.epsls]);
        o.push_str("#      NKX       NKY       NKZ       NK2\n");
        row(o, vals![p.grid_x, p.grid_y, p.grid_z, x.nk2]);
        o.push_str("#      NGX       NGY       NGZ    NASORD    NFDORD    NALIAS    NSPORD\n");
        row(
            o,
            vals![x.ng[0], x.ng[1], x.ng[2], x.nasord, x.nfdord, x.nalias, x.nspord],
        );
        o.push_str("#   NQEVAL    FACCUR    NRDGRD    NWRGDR\n");
        row(o, vals![x.nqeval, x.faccur, x.nrdgrd, x.nwrgdr]);
        o.push_str("#    NLRLJ    SLVDNS\n");
        row(o, vals![x.nlrlj, x.slvdns]);
    });

    if present("POSITIONRES") || p.ntpor != 0 || p.cpor != 0.0 {
        block(&mut out, "POSITIONRES", |o| {
            o.push_str("#   NTPOR  NTPORB  NTPORS    CPOR\n");
            row(o, vals![p.ntpor, p.ntporb, p.ntpors, p.cpor]);
        });
    }

    if present("ROTTRANS") || p.rtc != 0 {
        block(&mut out, "ROTTRANS", |o| {
            o.push_str(
                "#     RTC  RTCLAST
",
            );
            row(o, vals![p.rtc, p.rtclast]);
        });
    }

    // gromosXX validates DIR0 > 0 and TAUDIR > 0 even when NTDIR = 0, so an "off" block must
    // not be invented: written only when read or when it says something.
    let distanceres_set = p.ntdir != 0
        || p.ntdira != 0
        || p.cdir != 0.0
        || p.dir0 != 0.0
        || p.taudir != 0.0
        || p.forcescale != 0
        || p.vdir != 0
        || p.ntwdir != 0;
    if present("DISTANCERES") || distanceres_set {
        block(&mut out, "DISTANCERES", |o| {
            o.push_str("#   NTDIR  NTDIRA    CDIR    DIR0  TAUDIR  FORCESCALE  VDIR  NTWDIR\n");
            row(
                o,
                vals![
                    p.ntdir,
                    p.ntdira,
                    p.cdir,
                    p.dir0,
                    p.taudir,
                    p.forcescale,
                    p.vdir,
                    p.ntwdir
                ],
            );
        });
    }

    let perturbation_set = p.ntg != 0
        || p.nrdgl != 0
        || p.rlam != 0.0
        || p.dlamt != 0.0
        || p.alphlj != 0.0
        || p.alphc != 0.0
        || p.nlam != 1
        || p.nscale != 0;
    if present("PERTURBATION") || perturbation_set {
        block(&mut out, "PERTURBATION", |o| {
            o.push_str("#    NTG   NRDGL   RLAM   DLAMT\n");
            row(o, vals![p.ntg, p.nrdgl, p.rlam, p.dlamt]);
            o.push_str("#  ALPHLJ   ALPHC   NLAM   NSCALE\n");
            row(o, vals![p.alphlj, p.alphc, p.nlam, p.nscale]);
        });
    }

    if let Some(l) = &p.lambdas {
        block(&mut out, "LAMBDAS", |o| {
            o.push_str("#    NTIL\n");
            row(o, vals![if l.enabled { "on" } else { "off" }]);
            o.push_str("#    NTLI  NILG1  NILG2  ALI  BLI  CLI  DLI  ELI\n");
            for r in &l.rows {
                let [a, b, c, d, e] = r.coefficients;
                row(o, vals![r.kind.gromos_name(), r.g1, r.g2, a, b, c, d, e]);
            }
        });
    }

    if let Some(b) = &p.precalclam {
        block(&mut out, "PRECALCLAM", |o| {
            o.push_str("#    NRLAM   MINLAM   MAXLAM\n");
            row(o, vals![b.nr_lambdas, b.min_lambda, b.max_lambda]);
        });
    }

    if present("ENERGYMIN") || p.ntem != 0 {
        block(&mut out, "ENERGYMIN", |o| {
            o.push_str("#    NTEM    NCYC    DELE    DX0    DXM   NMIN   FLIM\n");
            row(
                o,
                vals![p.ntem, p.ncyc, p.dele, p.dx0, p.dxm, p.nmin, p.flim],
            );
        });
    }

    if let Some(b) = &p.eds {
        block(&mut out, "EDS", |o| {
            o.push_str("#    EDS\n");
            row(o, vals![i32::from(b.enabled)]);
            o.push_str("#    FORM  NUMSTATES\n");
            row(o, vals![b.form, b.num_states]);
            if !b.enabled {
                return;
            }
            if b.form == 3 {
                o.push_str("#    i  j  S\n");
                for pr in &b.pairs {
                    row(o, vals![pr.i, pr.j, pr.s]);
                }
            } else {
                o.push_str("#    S\n");
                row(o, &b.s.iter().map(|v| v.to_string()).collect::<Vec<_>>());
            }
            o.push_str("#    EIR\n");
            row(o, &b.eir.iter().map(|v| v.to_string()).collect::<Vec<_>>());
        });
    }

    // Passthrough: blocks the parser does not model, verbatim, in name order.
    let mut extra: Vec<(&String, &Vec<String>)> = p
        .raw_blocks
        .iter()
        .filter(|(name, _)| !MODELLED_BLOCKS.contains(&name.as_str()))
        .collect();
    extra.sort_by(|a, b| a.0.cmp(b.0));
    for (name, lines) in extra {
        block(&mut out, name, |o| {
            for line in lines {
                o.push_str(line);
                o.push('\n');
            }
        });
    }

    out
}

/// Write `params` to `path`.
pub fn write_imd_file<P: AsRef<Path>>(
    path: P,
    params: &ImdParameters,
    n_atoms: Option<usize>,
) -> Result<(), IoError> {
    fs::write(path.as_ref(), write_imd(params, n_atoms))
        .map_err(|e| IoError::WriteError(format!("{}: {e}", path.as_ref().display())))
}

fn constraint_algorithm(code: i32) -> String {
    match code {
        1 => "shake".to_string(),
        2 => "lincs".to_string(),
        3 => "settle".to_string(),
        other => other.to_string(),
    }
}

fn write_multibath(o: &mut String, bath: &TempBathParameters, n_atoms: Option<usize>) {
    match bath.algorithm {
        0 => row(o, vals!["weak-coupling"]),
        1 => row(o, vals!["nose-hoover"]),
        n => row(o, vals!["nose-hoover-chains", n.max(2)]),
    }
    o.push_str("#   NBATHS\n");
    let n_baths = bath.num_bath_groups.max(bath.temp0.len()).max(1);
    row(o, vals![n_baths]);
    o.push_str("#   TEMP0  TAU\n");
    let mut line = Vec::new();
    for i in 0..n_baths {
        line.push(bath.temp0.get(i).copied().unwrap_or(300.0).to_string());
        line.push(bath.tau.get(i).copied().unwrap_or(-1.0).to_string());
    }
    row(o, &line);
    o.push_str("#   DOFSET\n");
    let dof_sets: Vec<[usize; 3]> = if bath.dof_sets.is_empty() {
        // One DOF set covering every atom, coupled to bath 1 (COM and IR alike).
        vec![[n_atoms.unwrap_or(0), 1, 1]]
    } else {
        bath.dof_sets.clone()
    };
    row(o, vals![dof_sets.len()]);
    o.push_str("#   LAST   COM-BATH  IR-BATH\n");
    for set in &dof_sets {
        row(o, vals![set[0], set[1], set[2]]);
    }
}

fn write_pressurescale(o: &mut String, couple_pressure: bool, pp: &PressureParameters) {
    o.push_str("#   COUPLE   SCALE    COMP    TAUP  VIRIAL\n");
    let couple = if !couple_pressure {
        "off"
    } else {
        match pp.couple {
            0 => "off",
            1 => "calc",
            _ => "scale",
        }
    };
    let scale = match pp.algorithm {
        0 => "off",
        1 => "iso",
        2 => "aniso",
        3 => "full",
        4 => "semianiso",
        _ => "iso",
    };
    let virial = match pp.virial {
        0 => "none",
        1 => "atomic",
        2 => "molecular",
        _ => "none",
    };
    row(
        o,
        vals![couple, scale, pp.compressibility[0][0], pp.tau_p, virial],
    );
    o.push_str("#   SEMI\n");
    row(o, vals![pp.semi[0], pp.semi[1], pp.semi[2]]);
    o.push_str("#   reference pressure\n");
    for r in &pp.pressure0 {
        row(o, vals![r[0], r[1], r[2]]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::imd::{parse_imd_str, read_imd_file};

    const REFERENCE_IN: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../gromos-md/tests/gromosXX_references/aladip_solvated/aladip_solvated.in"
    );

    #[test]
    fn a_reference_input_survives_write_then_read() {
        let p = read_imd_file(REFERENCE_IN).expect("reference input");
        let text = write_imd(&p, Some(72));
        let q =
            parse_imd_str(&text).unwrap_or_else(|e| panic!("re-read what we wrote: {e}\n{text}"));
        assert_eq!(q.nstlim, p.nstlim);
        assert_eq!(q.dt, p.dt);
        assert_eq!(q.ntc, p.ntc);
        assert_eq!(q.ntf, p.ntf);
        assert_eq!(q.rcutp, p.rcutp);
        assert_eq!(q.rcutl, p.rcutl);
        assert_eq!(q.epsrf, p.epsrf);
        assert_eq!(q.nslfexcl, p.nslfexcl);
        assert_eq!(q.rlam, p.rlam);
        assert_eq!(q.temp_bath.len(), p.temp_bath.len());
        if let (Some(a), Some(b)) = (q.temp_bath.first(), p.temp_bath.first()) {
            assert_eq!(a.temp0, b.temp0);
            assert_eq!(a.tau, b.tau);
            assert_eq!(a.dof_sets, b.dof_sets);
        }
    }

    #[test]
    fn defaults_round_trip_and_carry_the_gromosxx_nslfexcl() {
        let p = ImdParameters::default();
        let q = parse_imd_str(&write_imd(&p, None)).expect("defaults re-read");
        assert_eq!(q.nslfexcl, 1);
        assert_eq!(q.nstlim, p.nstlim);
        assert_eq!(q.ntc, p.ntc);
    }
}
