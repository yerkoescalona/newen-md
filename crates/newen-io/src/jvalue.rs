//! J-value restraint files: the `JVALRESSPEC` specification block (gromos++ `jval`, `jepot`,
//! gromosXX `@jval`), the `JVALUERESEPS` local-elevation block of a final configuration and the
//! restraint trajectory (`.trs`: `TIMESTEP` + `JVALUERESEPS` per frame, gromos++ `RestrTraj`).
//! `#` starts a comment and blank lines are ignored, as in gromos++'s `Ginstream`.

use crate::IoError;

/// One `JVALRESSPEC` line: `i j k l weight J0 delta A B C` (atoms 1-based in the file).
#[derive(Debug, Clone)]
pub struct JvalRestraint {
    /// global atom indices, 0-based
    pub atoms: [usize; 4],
    pub weight: f64,
    pub j0: f64,
    pub delta: f64,
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

/// The `JVALUERESEPS` entry of one restraint in a restraint trajectory frame.
#[derive(Debug, Clone)]
pub struct JvalueEps {
    /// atom numbers as written (1-based)
    pub atoms: [usize; 4],
    pub epsilon: Vec<f64>,
}

/// One frame of a restraint trajectory.
#[derive(Debug, Clone, Default)]
pub struct RestraintFrame {
    pub step: i64,
    pub time: f64,
    pub jvalue_eps: Vec<JvalueEps>,
}

/// `(name, content lines)` of every block, comments stripped and blank lines dropped.
fn blocks(text: &str) -> Vec<(String, Vec<String>)> {
    let mut out: Vec<(String, Vec<String>)> = Vec::new();
    let mut current: Option<(String, Vec<String>)> = None;
    for raw in text.lines() {
        let line = raw.split('#').next().unwrap_or("").trim_end();
        if line.trim().is_empty() {
            continue;
        }
        match current.as_mut() {
            None => current = Some((line.trim().to_string(), Vec::new())),
            Some(b) => {
                if line.trim() == "END" {
                    out.push(current.take().unwrap());
                } else {
                    b.1.push(line.to_string());
                }
            },
        }
    }
    out
}

fn numbers(line: &str, what: &str) -> Result<Vec<f64>, IoError> {
    line.split_whitespace()
        .map(|t| {
            t.parse::<f64>().map_err(|_| {
                IoError::FormatError(format!("{what}: bad number '{t}' in line\n{line}"))
            })
        })
        .collect()
}

/// Read the `JVALRESSPEC` block of a J-value specification file.
pub fn read_jvalrespec(path: &str) -> Result<Vec<JvalRestraint>, IoError> {
    let text =
        std::fs::read_to_string(path).map_err(|_| IoError::FileNotFound(path.to_string()))?;
    let block = blocks(&text)
        .into_iter()
        .find(|b| b.0 == "JVALRESSPEC")
        .ok_or_else(|| {
            IoError::FormatError("jval file does not contain an JVALRESSPEC block!".into())
        })?;
    block
        .1
        .iter()
        .map(|line| {
            let r = numbers(line, "jval file")?;
            if r.len() < 10 {
                return Err(IoError::FormatError(format!(
                    "Bad line in jval-file\n{line}"
                )));
            }
            Ok(JvalRestraint {
                atoms: [
                    r[0] as usize - 1,
                    r[1] as usize - 1,
                    r[2] as usize - 1,
                    r[3] as usize - 1,
                ],
                weight: r[4],
                j0: r[5],
                delta: r[6],
                a: r[7],
                b: r[8],
                c: r[9],
            })
        })
        .collect()
}

/// The epsilon lines of the first `JVALUERESEPS` block of a configuration file (one line per
/// restraint, as gromos++ `jepot @fin` reads it).
pub fn read_jvalue_eps_block(path: &str) -> Result<Vec<Vec<f64>>, IoError> {
    let text =
        std::fs::read_to_string(path).map_err(|_| IoError::FileNotFound(path.to_string()))?;
    let block = blocks(&text)
        .into_iter()
        .find(|b| b.0 == "JVALUERESEPS")
        .ok_or_else(|| IoError::FormatError(format!("{path}: no JVALUERESEPS block")))?;
    block.1.iter().map(|l| numbers(l, "JVALUERESEPS")).collect()
}

/// Read a restraint trajectory: every `TIMESTEP` starts a frame, the `JVALUERESEPS` blocks that
/// follow belong to it (an atom line `i j k l` and an epsilon line per restraint).
pub fn read_restraint_trajectory(path: &str) -> Result<Vec<RestraintFrame>, IoError> {
    let text =
        std::fs::read_to_string(path).map_err(|_| IoError::FileNotFound(path.to_string()))?;
    let mut frames: Vec<RestraintFrame> = Vec::new();
    for (name, lines) in blocks(&text) {
        match name.as_str() {
            "TITLE" | "XRAYRVALUE" => {},
            "TIMESTEP" => {
                let first = lines
                    .first()
                    .ok_or_else(|| IoError::FormatError(format!("{path}: empty TIMESTEP block")))?;
                let v = numbers(first, "TIMESTEP")?;
                if v.len() < 2 {
                    return Err(IoError::FormatError(format!(
                        "{path}: bad line in TIMESTEP block. Got\n{first}"
                    )));
                }
                frames.push(RestraintFrame {
                    step: v[0] as i64,
                    time: v[1],
                    jvalue_eps: Vec::new(),
                });
            },
            "JVALUERESEPS" => {
                let frame = frames.last_mut().ok_or_else(|| {
                    IoError::FormatError(format!(
                        "{path}: trajectory does not contain a TIMESTEP block before JVALUERESEPS"
                    ))
                })?;
                let mut it = lines.iter();
                while let Some(atom_line) = it.next() {
                    let a = numbers(atom_line, "JVALUERESEPS")?;
                    if a.len() < 4 {
                        return Err(IoError::FormatError(format!(
                            "Bad line in JVALUERESEPS block:\n{atom_line}\nTrying to read atom numbers"
                        )));
                    }
                    let eps_line = it.next().ok_or_else(|| {
                        IoError::FormatError("Not enough lines in JVALUERESEPS block\n".into())
                    })?;
                    frame.jvalue_eps.push(JvalueEps {
                        atoms: [a[0] as usize, a[1] as usize, a[2] as usize, a[3] as usize],
                        epsilon: numbers(eps_line, "JVALUERESEPS")?,
                    });
                }
            },
            other => {
                return Err(IoError::FormatError(format!(
                    "Block {other} is unknown in a restraint trajectory file"
                )))
            },
        }
    }
    Ok(frames)
}
