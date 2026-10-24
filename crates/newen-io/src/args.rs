//! gromos++ `Arguments`: `@key value value …` on the command line (or in an `@f` file, expanded
//! by `newen_io::gromos_args`), every key multi-valued, absent / flag / valued distinguished.

use std::collections::BTreeMap;

#[derive(Debug, Clone, Default)]
pub struct Arguments {
    program: String,
    values: BTreeMap<String, Vec<String>>,
    order: Vec<String>,
}

impl Arguments {
    /// Parse the process arguments; `known` are the accepted keys (an unknown key is an error,
    /// as in gromos++). `--help`/`-h` returns `Err(usage)`.
    pub fn parse(known: &[&str], usage: &str) -> Result<Self, String> {
        Self::from_vec(crate::gromos_args(), known, usage)
    }

    pub fn from_vec(raw: Vec<String>, known: &[&str], usage: &str) -> Result<Self, String> {
        let program = raw.first().cloned().unwrap_or_default();
        let mut a = Arguments {
            program,
            ..Default::default()
        };
        let mut key: Option<String> = None;
        for tok in raw.iter().skip(1) {
            if tok == "--help" || tok == "-h" {
                return Err(usage.to_string());
            }
            if let Some(k) = tok.strip_prefix("--").or_else(|| tok.strip_prefix('@')) {
                if !known.contains(&k) {
                    return Err(format!("unknown argument @{k}\n{usage}"));
                }
                a.values.entry(k.to_string()).or_default();
                if !a.order.contains(&k.to_string()) {
                    a.order.push(k.to_string());
                }
                key = Some(k.to_string());
            } else if let Some(k) = &key {
                a.values.get_mut(k).unwrap().push(tok.clone());
            } else {
                return Err(format!("stray argument '{tok}'\n{usage}"));
            }
        }
        Ok(a)
    }

    pub fn program(&self) -> &str {
        &self.program
    }

    /// gromos++ `count`: −1 when the key is absent, else the number of values (0 for a flag).
    pub fn count(&self, key: &str) -> i64 {
        self.values.get(key).map_or(-1, |v| v.len() as i64)
    }

    pub fn has(&self, key: &str) -> bool {
        self.values.contains_key(key)
    }

    /// All values of `key` (empty when absent).
    pub fn values(&self, key: &str) -> &[String] {
        self.values.get(key).map_or(&[], |v| v.as_slice())
    }

    /// The first value of `key`, or an error naming it (gromos++ `args[key]` / `check`).
    pub fn value(&self, key: &str) -> Result<&str, String> {
        self.values(key)
            .first()
            .map(|s| s.as_str())
            .ok_or_else(|| format!("@{key} is required"))
    }

    /// The first value parsed as `T`, or `default` when absent.
    pub fn get<T: std::str::FromStr>(&self, key: &str, default: T) -> Result<T, String> {
        match self.values(key).first() {
            None => Ok(default),
            Some(v) => v
                .parse::<T>()
                .map_err(|_| format!("@{key}: cannot parse '{v}'")),
        }
    }

    /// The first value parsed as `T`, required.
    pub fn require<T: std::str::FromStr>(&self, key: &str) -> Result<T, String> {
        let v = self.value(key)?;
        v.parse::<T>()
            .map_err(|_| format!("@{key}: cannot parse '{v}'"))
    }
}

/// Print `err` the gromos++ way and exit 1 (usage goes out unchanged).
pub fn fail(err: &str) -> ! {
    eprintln!("{err}");
    std::process::exit(1)
}

/// How faithful a program is to the gromos++ program of the same name.
///
/// The audit tables in the crates' `CONTEXT.md` as data: a caller (Python included) can see what
/// they are calling. Porting a stand-in changes its status, not its API.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// A reference test compares our output with gromos++'s.
    Verified,
    /// Ported from gromos++, no reference test yet.
    Ported,
    /// Same name, simplified or different semantics — a utility, not a port.
    StandIn,
    /// Recognised and parsed, does not do the work yet.
    Stub,
    /// No gromos++ program of this name.
    NoCounterpart,
}

impl Status {
    pub fn as_str(&self) -> &'static str {
        match self {
            Status::Verified => "verified",
            Status::Ported => "ported",
            Status::StandIn => "stand_in",
            Status::Stub => "stub",
            Status::NoCounterpart => "no_counterpart",
        }
    }
}

/// One gromos++-style program, callable in process.
///
/// The body lives in its crate's library (`<crate>::programs::<name>`), so the binary and the
/// Python binding run the same code: `run` takes an already-[`crate::normalize`]d argv and writes
/// what the program prints to `out`. It never exits and never touches the process stdout.
pub struct Program {
    pub name: &'static str,
    /// The `@key`s the program accepts, for `Arguments::from_vec` and for the registry.
    pub known: &'static [&'static str],
    pub usage: &'static str,
    pub status: Status,
    pub run: fn(&[String], &mut dyn std::io::Write) -> Result<(), String>,
}

impl std::fmt::Debug for Program {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.debug_struct("Program")
            .field("name", &self.name)
            .field("known", &self.known)
            .field("status", &self.status)
            .finish_non_exhaustive()
    }
}

/// Run `program` as a command-line binary: the whole `fn main` of every program binary.
///
/// The one place allowed to read the process arguments, lock stdout and exit — everything a
/// `Program`'s `run` must not do.
pub fn run_main(program: &Program) -> ! {
    use std::io::Write;

    let argv = match crate::normalize(std::env::args().collect()) {
        Ok(a) => a,
        Err(e) => fail(&e),
    };

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    let result = (program.run)(&argv, &mut out);
    let flushed = out.flush();
    drop(out);

    match result {
        Ok(()) => match flushed {
            Ok(()) => std::process::exit(0),
            Err(e) => fail(&format!("{}: cannot write output: {e}", program.name)),
        },
        Err(e) => fail(&e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keys_are_multi_valued_and_flags_count_zero() {
        let a = Arguments::from_vec(
            [
                "p", "--topo", "t.top", "--traj", "a.trc", "b.trc", "--nots", "@time", "0", "0.5",
            ]
            .iter()
            .map(|s| s.to_string())
            .collect(),
            &["topo", "traj", "nots", "time"],
            "usage",
        )
        .unwrap();
        assert_eq!(a.count("traj"), 2);
        assert_eq!(a.count("nots"), 0);
        assert_eq!(a.count("absent"), -1);
        assert_eq!(a.value("topo").unwrap(), "t.top");
        assert_eq!(a.get::<f64>("time", 1.0).unwrap(), 0.0);
        assert_eq!(a.values("time"), &["0", "0.5"]);
    }

    #[test]
    fn unknown_keys_are_rejected() {
        assert!(Arguments::from_vec(vec!["p".into(), "--nope".into()], &["topo"], "u").is_err());
    }
}
