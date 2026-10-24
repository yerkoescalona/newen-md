//! gromos++ `utils::Time`: `@time <t0> <dt>` overrides the times in the trajectory (t0 + n·dt);
//! without it the frame's own time is used. Printed fixed with nine decimals in a 15-wide field.

use crate::args::Arguments;

#[derive(Debug, Clone)]
pub struct Time {
    t0: f64,
    dt: f64,
    read: bool,
    steps: usize,
    current: f64,
}

impl Time {
    pub fn from_args(args: &Arguments) -> Result<Self, String> {
        let v = args.values("time");
        let read = v.is_empty();
        let t0 = v
            .first()
            .map(|s| {
                s.parse::<f64>()
                    .map_err(|_| "@time: starting time is not numeric".to_string())
            })
            .transpose()?
            .unwrap_or(0.0);
        let dt = v
            .get(1)
            .map(|s| {
                s.parse::<f64>()
                    .map_err(|_| "@time: dt is not numeric".to_string())
            })
            .transpose()?
            .unwrap_or(1.0);
        Ok(Time {
            t0,
            dt,
            read,
            steps: 0,
            current: 0.0,
        })
    }

    /// Advance with a frame: its own time, or t0 + n·dt when `@time` was given.
    pub fn next(&mut self, frame_time: f64) -> f64 {
        self.current = if self.read {
            frame_time
        } else {
            self.t0 + self.steps as f64 * self.dt
        };
        self.steps += 1;
        self.current
    }

    pub fn current(&self) -> f64 {
        self.current
    }

    /// gromos++ `Time::print`: `setw(15) << setprecision(9)`, fixed.
    pub fn fmt(t: f64) -> String {
        format!("{:15.9}", t)
    }
}

/// gromos++ `@timespec ALL|EVERY|SPEC` with `@timepts`: which frames of a series a program
/// evaluates (`computeJval`, `computeJepot`, …). Frames are counted from one.
#[derive(Debug, Clone)]
pub struct TimeSpec {
    spec: String,
    pts: Vec<usize>,
    n: usize,
    /// frames selected so far
    pub written: usize,
    /// every `SPEC` frame has been seen — the caller may stop reading
    pub done: bool,
}

impl TimeSpec {
    pub fn from_args(args: &Arguments) -> Result<Self, String> {
        let spec = args
            .values("timespec")
            .first()
            .cloned()
            .unwrap_or_else(|| "ALL".into());
        if !["ALL", "EVERY", "SPEC"].contains(&spec.as_str()) {
            return Err(format!("timespec format {spec} unknown.\n"));
        }
        let mut pts = Vec::new();
        if spec != "ALL" {
            for p in args.values("timepts") {
                pts.push(
                    p.parse::<usize>()
                        .map_err(|_| format!("@timepts: '{p}' is not a number"))?,
                );
            }
            if pts.is_empty() {
                return Err("if you give EVERY or SPEC you have to use @timepts as well".into());
            }
            if spec == "EVERY" && pts.len() != 1 {
                return Err(
                    "if you give EVERY you have to give exactly one number with @timepts".into(),
                );
            }
        }
        Ok(TimeSpec {
            spec,
            pts,
            n: 0,
            written: 0,
            done: false,
        })
    }

    /// Register the next frame; true when it is selected.
    pub fn select(&mut self) -> bool {
        self.n += 1;
        let hit = match self.spec.as_str() {
            "ALL" => true,
            "EVERY" => self.n.is_multiple_of(self.pts[0]),
            _ => self.pts.contains(&self.n),
        };
        if hit {
            self.written += 1;
            if self.spec == "SPEC" && self.written == self.pts.len() {
                self.done = true;
            }
        }
        hit
    }
}
