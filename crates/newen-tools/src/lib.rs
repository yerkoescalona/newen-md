//! newen system construction tools, as library functions.
//!
//! Every program has its body here — `programs::<name>::run(argv, out)` — not in its binary's
//! `main`, so the command line and an in-process caller (the Python binding) run the same code.
//! A [`Program`] record per program carries the name, the accepted `@key`s, the usage text and
//! how faithful the program is to gromos++'s; [`PROGRAMS`] is the registry.
//!
//! ```no_run
//! let mut out = Vec::new();
//! let argv: Vec<String> = ["duplicate", "--topo", "a.top", "--pos", "a.cnf", "--pbc", "r"]
//!     .iter()
//!     .map(|s| s.to_string())
//!     .collect();
//! (newen_tools::find("duplicate").unwrap().run)(&argv, &mut out).unwrap();
//! ```

pub mod programs;

pub use newen_io::args::{Program, Status};
pub use programs::{find, PROGRAMS};
