//! An energy trajectory read according to a library file's declaration of it — gromos++'s
//! `utils::EnergyTraj` plus the `gmath::Expression` evaluator it leans on.
//!
//! `ene_ana`, `ext_ti_ana` and friends do not know the layout of a `.tre`/`.trg` at compile
//! time. A library file declares it:
//!
//! ```text
//! ENERTRJ
//!   block ENERGY03
//!     subblock ENER 50 1
//!     size NUM_BATHS
//!     subblock KINENER NUM_BATHS 3
//!     subblock NONBONDED matrix_NUM_ENERGY_GROUPS 6
//! END
//! VARIABLES
//!   totpot = ENER[3]
//!   solute_intra = BONDED[1][1] + BONDED[1][2]
//! END
//! ```
//!
//! `block` names a GROMOS block in the file; the `subblock`s that follow consume its values in
//! order. `size` reads one integer from the stream and binds it as a dimension for later
//! subblocks — and simultaneously binds `matrix_<NAME>` to N(N+1)/2, which is how the
//! energy-group pair table is declared. `VARIABLES` then names arithmetic over those tables.
//! That indirection is what lets one library file span GROMOS versions, and why the trajectory
//! carries an `ENEVERSION` block to check against the library's.
//!
//! **Quirk reproduced deliberately:** gromos++ splits an operator chain at the *rightmost*
//! `+`/`-`, but at the *leftmost* `*`/`/` (`Expression.cc:310`), so `*` and `/` come out
//! right-associative: `a / b / c` evaluates as `a / (b / c)`, not `(a / b) / c`. The shipped
//! library's only multi-operator entry is `densit = MASS[1] * 1.66056 / VOLUME[1]`, where the
//! two groupings differ only in the last bits, but reproducing the rule keeps every digit of
//! our output equal to gromos++'s.

use crate::gz::TextReader;
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fmt::Write as _;

/// Version of the library profile (`docs/src/reference/energy-library.md`): the canonical
/// form, the `# energy-schema` / `# energy-layout` lines and the reader's obligations.
pub const PROFILE_VERSION: u32 = 1;

/// A subblock dimension: a literal count, or a `size` bound while reading the frame.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Dim {
    Fixed(usize),
    Named(String),
}

impl Dim {
    fn parse(s: &str, declared: &HashMap<String, usize>) -> Result<Dim, String> {
        if let Ok(n) = s.parse::<usize>() {
            return Ok(Dim::Fixed(n));
        }
        if declared.contains_key(s)
            || declared.contains_key(s.strip_prefix("matrix_").unwrap_or(""))
        {
            Ok(Dim::Named(s.to_string()))
        } else {
            Err(format!("variable name {s} in block size not defined"))
        }
    }
}

impl std::fmt::Display for Dim {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Dim::Fixed(n) => write!(f, "{n}"),
            Dim::Named(n) => f.write_str(n),
        }
    }
}

/// One line of an `ENERTRJ`/`FRENERTRJ` declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Entry {
    /// A GROMOS block in the trajectory file; the entries after it consume its values.
    Block(String),
    /// Read one integer and bind it (and `matrix_<name>`) as a dimension.
    Size(String),
    /// A named `rows × cols` table of values.
    Sub { name: String, rows: Dim, cols: Dim },
}

impl std::fmt::Display for Entry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Entry::Block(b) => write!(f, "block {b}"),
            Entry::Size(s) => write!(f, "size {s}"),
            Entry::Sub { name, rows, cols } => write!(f, "subblock {name} {rows} {cols}"),
        }
    }
}

/// The declaration of one file type — the body of an `ENERTRJ` or `FRENERTRJ` block — in
/// declaration order. Two schemas that place every value in the same slot have the same
/// [`fingerprint`](Schema::fingerprint), whatever they call the tables.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Schema {
    file_type: String,
    entries: Vec<Entry>,
}

impl Schema {
    /// Parse declaration lines (`block`, `size`, `subblock`); comments and blank lines are
    /// skipped, so the body of a library block can be passed as read.
    pub fn parse<'a>(
        file_type: &str,
        lines: impl IntoIterator<Item = &'a str>,
    ) -> Result<Schema, String> {
        let mut entries = Vec::new();
        let mut declared: HashMap<String, usize> = HashMap::new();
        for raw in lines {
            let line = raw.split('#').next().unwrap_or("").trim();
            let t: Vec<&str> = line.split_whitespace().collect();
            if t.is_empty() {
                continue;
            }
            match t[0].to_ascii_lowercase().as_str() {
                "block" => {
                    if t.len() != 2 {
                        return Err(format!("'block' takes one name: {line}"));
                    }
                    entries.push(Entry::Block(t[1].to_string()));
                },
                "size" => {
                    if t.len() != 2 {
                        return Err(format!("'size' takes one name: {line}"));
                    }
                    declared.insert(t[1].to_string(), 0);
                    entries.push(Entry::Size(t[1].to_string()));
                },
                "subblock" => {
                    if t.len() != 4 {
                        return Err(format!("'subblock' takes a name and two sizes: {line}"));
                    }
                    entries.push(Entry::Sub {
                        name: t[1].to_string(),
                        rows: Dim::parse(t[2], &declared)?,
                        cols: Dim::parse(t[3], &declared)?,
                    });
                },
                other => return Err(format!("unknown library keyword '{other}': {line}")),
            }
        }
        Ok(Schema {
            file_type: file_type.to_string(),
            entries,
        })
    }

    /// `ENERTRJ` or `FRENERTRJ`.
    pub fn file_type(&self) -> &str {
        &self.file_type
    }

    /// The declarations with their names, one per line, single-spaced: what `# energy-layout`
    /// lines carry and what a library block lists.
    pub fn declaration_lines(&self) -> Vec<String> {
        self.entries.iter().map(|e| e.to_string()).collect()
    }

    /// The declarations with names dropped and sizes numbered by first occurrence (profile R2):
    /// the text the fingerprint hashes.
    pub fn canonical_form(&self) -> String {
        fn dim(d: &Dim, numbered: &[&str]) -> String {
            match d {
                Dim::Fixed(n) => n.to_string(),
                Dim::Named(n) => {
                    let (base, matrix) = match n.strip_prefix("matrix_") {
                        Some(b) => (b, true),
                        None => (n.as_str(), false),
                    };
                    let k = numbered
                        .iter()
                        .position(|s| *s == base)
                        .map_or(0, |i| i + 1);
                    if matrix {
                        format!("matrix(${k})")
                    } else {
                        format!("${k}")
                    }
                },
            }
        }
        let mut numbered: Vec<&str> = Vec::new();
        let mut out = String::new();
        for e in &self.entries {
            match e {
                Entry::Block(b) => {
                    let _ = writeln!(out, "block {b}");
                },
                Entry::Size(s) => {
                    if !numbered.contains(&s.as_str()) {
                        numbered.push(s);
                    }
                    let _ = writeln!(out, "size {}", dim(&Dim::Named(s.clone()), &numbered));
                },
                Entry::Sub { rows, cols, .. } => {
                    let _ = writeln!(
                        out,
                        "subblock {} {}",
                        dim(rows, &numbered),
                        dim(cols, &numbered)
                    );
                },
            }
        }
        out
    }

    /// `sha256:` + hex SHA-256 of the canonical form.
    pub fn fingerprint(&self) -> String {
        format!(
            "sha256:{:x}",
            Sha256::digest(self.canonical_form().as_bytes())
        )
    }

    /// The block as a library file lists it, in md++'s indentation.
    pub fn library_block(&self) -> String {
        let mut out = format!("{}\n", self.file_type);
        for e in &self.entries {
            let indent = if matches!(e, Entry::Block(_)) {
                "  "
            } else {
                "    "
            };
            let _ = writeln!(out, "{indent}{e}");
        }
        out.push_str("END\n");
        out
    }

    /// Where two schemas differ, one line per declaration slot, e.g.
    /// `  ENERGY03 ENER:   library 26 x 2,   file 52 x 1`. Empty when the layouts agree
    /// (names aside).
    pub fn diff(&self, other: &Schema, labels: (&str, &str)) -> Vec<String> {
        fn by_block(s: &Schema) -> Vec<(String, Vec<&Entry>)> {
            let mut out: Vec<(String, Vec<&Entry>)> = Vec::new();
            for e in &s.entries {
                match e {
                    Entry::Block(b) => out.push((b.clone(), Vec::new())),
                    _ => match out.last_mut() {
                        Some((_, v)) => v.push(e),
                        None => out.push((String::new(), vec![e])),
                    },
                }
            }
            out
        }
        fn shape(e: Option<&&Entry>) -> String {
            match e {
                None => "absent".into(),
                Some(Entry::Sub { rows, cols, .. }) => format!("{rows} x {cols}"),
                Some(Entry::Size(_)) => "size".into(),
                Some(Entry::Block(_)) => unreachable!(),
            }
        }
        fn same_shape(a: &Entry, b: &Entry) -> bool {
            match (a, b) {
                (Entry::Size(_), Entry::Size(_)) => true,
                (
                    Entry::Sub {
                        rows: r, cols: c, ..
                    },
                    Entry::Sub {
                        rows: r2, cols: c2, ..
                    },
                ) => r == r2 && c == c2,
                _ => false,
            }
        }
        let (mine, theirs) = (by_block(self), by_block(other));
        let mut blocks: Vec<&str> = mine.iter().map(|(b, _)| b.as_str()).collect();
        for (b, _) in &theirs {
            if !blocks.contains(&b.as_str()) {
                blocks.push(b);
            }
        }
        let mut out = Vec::new();
        for block in blocks {
            let a = mine
                .iter()
                .find(|(b, _)| b == block)
                .map(|(_, v)| v.as_slice());
            let b = theirs
                .iter()
                .find(|(b, _)| b == block)
                .map(|(_, v)| v.as_slice());
            let (Some(a), Some(b)) = (a, b) else {
                let (has, lacks) = if a.is_some() {
                    labels
                } else {
                    (labels.1, labels.0)
                };
                out.push(format!(
                    "  {block}:   {has} declares the block,   {lacks} does not"
                ));
                continue;
            };
            for i in 0..a.len().max(b.len()) {
                let (x, y) = (a.get(i), b.get(i));
                if let (Some(x), Some(y)) = (x, y) {
                    if same_shape(x, y) {
                        continue;
                    }
                }
                let name = match x.or(y) {
                    Some(Entry::Sub { name, .. }) | Some(Entry::Size(name)) => name.as_str(),
                    _ => "?",
                };
                out.push(format!(
                    "  {block} {name}:   {} {},   {} {}",
                    labels.0,
                    shape(x),
                    labels.1,
                    shape(y)
                ));
            }
        }
        out
    }
}

/// The upstream layouts gromos-rs knows, by `ENEVERSION` and file type. A new upstream layout
/// is a new entry under its own version; an existing entry is never edited (profile R9).
const LAYOUTS: &[(&str, &str, &str)] = &[
    (
        "2023-04-15",
        "ENERTRJ",
        "block TIMESTEP
         subblock TIME 2 1
         block ENERGY03
         subblock ENER 52 1
         size NUM_BATHS
         subblock KINENER NUM_BATHS 3
         size NUM_ENERGY_GROUPS
         subblock BONDED NUM_ENERGY_GROUPS 5
         subblock NONBONDED matrix_NUM_ENERGY_GROUPS 6
         subblock SPECIAL NUM_ENERGY_GROUPS 13
         size NUM_EDS_STATES
         subblock EDS NUM_EDS_STATES 6
         size NUM_GAMD_GROUPS
         subblock GAMD NUM_GAMD_GROUPS 5
         size NUM_LAMBDAS
         subblock PRECALCLAM NUM_LAMBDAS 12
         subblock ABDIH 1 2
         block VOLUMEPRESSURE03
         subblock MASS 1 1
         size NUM_BATHS
         subblock TEMPERATURE NUM_BATHS 4
         subblock VOLUME 10 1
         subblock PRESSURE 30 1",
    ),
    (
        "2023-04-15",
        "FRENERTRJ",
        "block TIMESTEP
         subblock TIME 2 1
         block FREEENERDERIVS03
         subblock RLAM 1 1
         subblock FREEENER 52 1
         size NUM_BATHS
         subblock FREEKINENER NUM_BATHS 3
         size NUM_ENERGY_GROUPS
         subblock FREEBONDED NUM_ENERGY_GROUPS 5
         subblock FREENONBONDED matrix_NUM_ENERGY_GROUPS 6
         subblock FREESPECIAL NUM_ENERGY_GROUPS 13
         size NUM_EDS_STATES
         subblock FREEEDS NUM_EDS_STATES 6
         size NUM_GAMD_GROUPS
         subblock FREEGAMD NUM_GAMD_GROUPS 5
         size NUM_LAMBDAS
         subblock FREEPRECALCLAM NUM_LAMBDAS 12
         subblock FREEABDIH 1 2",
    ),
];

/// The `ENEVERSION` strings with a built-in layout, newest first.
pub const KNOWN_VERSIONS: &[&str] = &["2023-04-15"];

/// The built-in layout of `file_type` (`ENERTRJ`/`FRENERTRJ`) for an upstream `ENEVERSION`.
pub fn builtin_schema(version: &str, file_type: &str) -> Option<Schema> {
    LAYOUTS
        .iter()
        .find(|(v, t, _)| *v == version && *t == file_type)
        .map(|(_, t, decl)| Schema::parse(t, decl.lines()).expect("built-in layout parses"))
}

/// The library `gromos-io` ships: schema sections generated from the built-in layouts,
/// `VARIABLES` taken verbatim from md++'s `ene_ana.md++.lib`. `ene_ana` reads it when no
/// `@library` is given; the checked-in copy is tested equal to [`official_library_text`].
pub const OFFICIAL_LIBRARY: &str = include_str!("../data/ene_ana.md++.lib");

/// The `ENEVERSION` the official library describes.
pub const OFFICIAL_LIBRARY_VERSION: &str = "2023-04-15";

/// A complete library file: `TITLE`, `ENEVERSION`, one `# energy-schema` line per schema, the
/// schema blocks, then `variables_block` (the lines between `VARIABLES` and `END`, verbatim).
pub fn library_text(
    title: &str,
    version: &str,
    schemas: &[Schema],
    variables_block: &str,
) -> String {
    let mut out = format!("TITLE\n{title}\nEND\nENEVERSION\n  {version}\nEND\n");
    for s in schemas {
        let _ = writeln!(
            out,
            "# energy-schema {PROFILE_VERSION} {} {}",
            s.file_type(),
            s.fingerprint()
        );
    }
    for s in schemas {
        out.push_str(&s.library_block());
    }
    out.push_str("VARIABLES\n");
    out.push_str(variables_block);
    if !variables_block.ends_with('\n') {
        out.push('\n');
    }
    out.push_str("END\n");
    out
}

/// What [`OFFICIAL_LIBRARY`] must equal, byte for byte.
pub fn official_library_text() -> String {
    let schemas: Vec<Schema> = ["ENERTRJ", "FRENERTRJ"]
        .iter()
        .map(|t| builtin_schema(OFFICIAL_LIBRARY_VERSION, t).expect("official layout"))
        .collect();
    library_text(
        "  ene_ana library for the md++ 2023-04-15 energy and free-energy trajectory layout\n  \
         schema blocks generated by gromos-io from its layout definition; VARIABLES from md++",
        OFFICIAL_LIBRARY_VERSION,
        &schemas,
        &official_variables(),
    )
}

/// The `VARIABLES` section of the official library, comments included, as md++ wrote it.
pub fn official_variables() -> String {
    raw_block(OFFICIAL_LIBRARY, "VARIABLES").unwrap_or_default()
}

/// The lines of block `name` in a GROMOS text, comments included, as written.
fn raw_block(text: &str, name: &str) -> Option<String> {
    let mut lines = text.lines().skip_while(|l| l.trim() != name).skip(1);
    let mut out = String::new();
    let mut ended = false;
    for l in lines.by_ref() {
        if l.trim() == "END" {
            ended = true;
            break;
        }
        out.push_str(l);
        out.push('\n');
    }
    ended.then_some(out)
}

/// `# energy-schema N TYPE sha256:…` lines of a library (profile R4), as `(N, type, fingerprint)`.
fn schema_lines(text: &str) -> Result<Vec<(u32, String, String)>, String> {
    let mut out = Vec::new();
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("# energy-schema") else {
            continue;
        };
        let t: Vec<&str> = rest.split_whitespace().collect();
        if t.len() != 3 {
            return Err(format!("malformed self-description line: {line}"));
        }
        let n: u32 = t[0]
            .parse()
            .map_err(|_| format!("malformed self-description line: {line}"))?;
        if n > PROFILE_VERSION {
            return Err(format!(
                "energy-schema {n} is newer than this reader (profile {PROFILE_VERSION}); \
                 update gromos-rs"
            ));
        }
        out.push((n, t[1].to_string(), t[2].to_string()));
    }
    Ok(out)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Op {
    Add,
    Sub,
    Mul,
    Div,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Fun {
    Sin,
    Cos,
    Exp,
    Log,
}

#[derive(Debug, Clone)]
enum Expr {
    Const(f64),
    /// `SUBBLOCK[i]` / `SUBBLOCK[i][j]`, or a bare name resolved at evaluation time.
    Name {
        name: String,
        i: usize,
        j: usize,
    },
    Fun(Fun, Box<Expr>),
    Bin(Op, Box<Expr>, Box<Expr>),
}

/// The `# energy-schema` / `# energy-layout` lines a writer puts after `ENEVERSION`
/// (profile R3): the file's own verifiable statement of its layout.
/// The energy-group pair a row of a `matrix_NUM_ENERGY_GROUPS` subblock holds, 0-based.
/// gromosXX writes them `for i in 0..n { for j in i..n }` — upper triangle, row by row,
/// diagonal included (`out_configuration.cc:3148`). `None` past the last pair.
pub fn group_pair(row: usize, n_groups: usize) -> Option<(usize, usize)> {
    let mut r = row;
    for i in 0..n_groups {
        let width = n_groups - i;
        if r < width {
            return Some((i, i + r));
        }
        r -= width;
    }
    None
}

pub fn self_description_lines(schema: &Schema) -> Vec<String> {
    let mut out = vec![format!(
        "# energy-schema {PROFILE_VERSION} {} {}",
        schema.file_type(),
        schema.fingerprint()
    )];
    out.extend(
        schema
            .declaration_lines()
            .into_iter()
            .map(|l| format!("# energy-layout {l}")),
    );
    out
}

/// What a writer records about the run behind a trajectory (profile R8).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Provenance {
    /// The program and version, e.g. `gromos-rs md 0.0.47`.
    pub writer: String,
    /// `(sha256:<hex>, path as given)` of the topology.
    pub topology: Option<(String, String)>,
    /// The energy groups as inclusive 1-based atom ranges, in order.
    pub energy_groups: Vec<(usize, usize)>,
}

impl Provenance {
    /// The `# energy-provenance` lines.
    pub fn lines(&self) -> Vec<String> {
        let mut out = Vec::new();
        if !self.writer.is_empty() {
            out.push(format!("# energy-provenance writer {}", self.writer));
        }
        if let Some((digest, path)) = &self.topology {
            out.push(format!("# energy-provenance topology {digest}   {path}"));
        }
        if !self.energy_groups.is_empty() {
            let ranges: Vec<String> = self
                .energy_groups
                .iter()
                .map(|(a, b)| format!("{a}-{b}"))
                .collect();
            out.push(format!(
                "# energy-provenance energy-groups {}",
                ranges.join(" ")
            ));
        }
        out
    }
}

/// `sha256:<hex>` of a file's bytes, for [`Provenance::topology`].
pub fn sha256_file(path: impl AsRef<std::path::Path>) -> std::io::Result<String> {
    let mut file = std::fs::File::open(path)?;
    let mut hasher = Sha256::new();
    std::io::copy(&mut file, &mut hasher)?;
    Ok(format!("sha256:{:x}", hasher.finalize()))
}

/// What a trajectory says about itself before its first frame: `TITLE`, `ENEVERSION`, and
/// the `# energy-…` lines a gromos-rs writer adds (profile R3, R8).
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Preamble {
    pub title: Vec<String>,
    pub version: Option<String>,
    /// The layout the file declares for itself, fingerprint verified.
    pub schema: Option<Schema>,
    /// `# energy-provenance <key> <rest>` lines, in file order.
    pub provenance: Vec<(String, String)>,
}

impl Preamble {
    /// The `writer` provenance entry, e.g. `gromos-rs md 0.0.47`.
    pub fn writer(&self) -> Option<&str> {
        self.provenance
            .iter()
            .find(|(k, _)| k == "writer")
            .map(|(_, v)| v.as_str())
    }

    /// `ENEVERSION 2023-04-15, written by gromos-rs md 0.0.47` — for messages.
    fn describe(&self) -> String {
        let mut s = format!("ENEVERSION {}", self.version.as_deref().unwrap_or("NONE"));
        if let Some(w) = self.writer() {
            let _ = write!(s, ", written by {w}");
        }
        s
    }
}

/// Read a trajectory's preamble, leaving the reader at its first data block. gromos++
/// consumes `TITLE` and `ENEVERSION` in `Ginstream::open`; the profile lines sit between
/// `ENEVERSION`'s `END` and the first `TIMESTEP` and are comments to every other reader.
pub fn read_preamble(reader: &mut TextReader) -> Result<Preamble, String> {
    let mut pre = Preamble::default();
    let mut schema_line: Option<(String, String)> = None;
    let mut layout: Vec<String> = Vec::new();
    let mut buf = String::new();
    loop {
        buf.clear();
        if reader.read_line(&mut buf).map_err(|e| e.to_string())? == 0 {
            break;
        }
        if let Some(rest) = buf.strip_prefix("# energy-") {
            let (key, value) = rest
                .trim_end()
                .split_once(char::is_whitespace)
                .unwrap_or((rest.trim_end(), ""));
            let value = value.trim();
            match key {
                "schema" => {
                    let t: Vec<&str> = value.split_whitespace().collect();
                    if t.len() != 3 {
                        return Err(format!(
                            "trajectory self-description is corrupt: {}",
                            buf.trim_end()
                        ));
                    }
                    let n: u32 = t[0].parse().map_err(|_| {
                        format!("trajectory self-description is corrupt: {}", buf.trim_end())
                    })?;
                    if n > PROFILE_VERSION {
                        return Err(format!(
                            "energy-schema {n} is newer than this reader (profile \
                             {PROFILE_VERSION}); update gromos-rs"
                        ));
                    }
                    if schema_line.is_some() {
                        return Err(
                            "trajectory self-description is corrupt: two energy-schema lines"
                                .into(),
                        );
                    }
                    schema_line = Some((t[1].to_string(), t[2].to_string()));
                },
                "layout" => layout.push(value.to_string()),
                "provenance" => {
                    let (k, v) = value.split_once(char::is_whitespace).unwrap_or((value, ""));
                    pre.provenance.push((k.to_string(), v.trim().to_string()));
                },
                _ => {},
            }
            continue;
        }
        let name = buf.split('#').next().unwrap_or("").trim();
        if name.is_empty() {
            continue;
        }
        if name != "TITLE" && name != "ENEVERSION" {
            reader.unread_line(&buf);
            break;
        }
        let name = name.to_string();
        let mut body: Vec<String> = Vec::new();
        loop {
            buf.clear();
            if reader.read_line(&mut buf).map_err(|e| e.to_string())? == 0 {
                return Err(format!("no END in {name} block"));
            }
            let l = buf.split('#').next().unwrap_or("").trim();
            if l == "END" {
                break;
            }
            if !l.is_empty() {
                body.push(l.to_string());
            }
        }
        if name == "ENEVERSION" {
            pre.version = Some(body.concat().split_whitespace().collect());
        } else {
            pre.title = body;
        }
    }
    match (schema_line, layout.is_empty()) {
        (None, true) => {},
        (None, false) => {
            return Err(
                "trajectory self-description is corrupt: energy-layout lines without an \
                 energy-schema line"
                    .into(),
            )
        },
        (Some((file_type, claimed)), _) => {
            let schema = Schema::parse(&file_type, layout.iter().map(String::as_str))
                .map_err(|e| format!("trajectory self-description is corrupt: {e}"))?;
            let actual = schema.fingerprint();
            if actual != claimed {
                return Err(format!(
                    "trajectory self-description is corrupt: energy-schema says {claimed}, \
                     the energy-layout lines hash to {actual}"
                ));
            }
            pre.schema = Some(schema);
        },
    }
    Ok(pre)
}

/// How a trajectory's layout was established before reading it (profile R5).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Tier {
    /// The file carries its own verified layout; the library agreed with it.
    SelfDescribed,
    /// The file's `ENEVERSION` has a built-in layout; the library agreed with it.
    KnownVersion(String),
    /// Nothing to check the library against; its declaration was taken on trust.
    Unverified,
}

/// An energy trajectory plus the library that describes it.
#[derive(Debug, Default)]
pub struct EnergyTraj {
    /// File type (`ENERTRJ`, `FRENERTRJ`) → its declaration.
    layouts: HashMap<String, Schema>,
    /// Subblock name → `[rows][cols]` of the frame currently held.
    data: HashMap<String, Vec<Vec<f64>>>,
    /// `size` bindings of the frame currently held, including the `matrix_` forms.
    sizes: HashMap<String, usize>,
    /// The sizes of the first frame read, per file type: a later frame must repeat them.
    first_sizes: HashMap<(String, String), usize>,
    constants: HashMap<String, f64>,
    variables: HashMap<String, Expr>,
    version: Option<String>,
    /// File types bound but not yet read: their first frame gets the invariant checks.
    unchecked: Vec<(String, Tier)>,
    warnings: Vec<String>,
}

impl EnergyTraj {
    /// Read a library file: the `ENERTRJ`/`FRENERTRJ` declarations, `VARIABLES` and
    /// `ENEVERSION`. A `# energy-schema` line the library carries must match its own
    /// declaration (profile R4).
    pub fn from_library(text: &str) -> Result<Self, String> {
        let mut traj = EnergyTraj::default();
        // gromos++ `set_standards`
        traj.add_constant("BOLTZ", gromos_core::units::kB);

        for (name, body) in blocks(text) {
            match name.as_str() {
                "ENERTRJ" | "FRENERTRJ" => {
                    let schema = Schema::parse(&name, body.iter().map(String::as_str))?;
                    traj.set_schema(schema);
                },
                "ENEVERSION" => {
                    let v: String = body.concat().split_whitespace().collect();
                    if traj.version.is_some() {
                        return Err("library has two ENEVERSION blocks".into());
                    }
                    traj.version = Some(v);
                },
                "VARIABLES" => {
                    for (var, expr) in variable_definitions(&body) {
                        let parsed = parse_expression(&expr)
                            .map_err(|e| format!("library variable '{var} = {expr}': {e}"))?;
                        traj.variables.insert(var, parsed);
                    }
                },
                _ => {},
            }
        }
        // a library may carry only VARIABLES; the layout then comes from the trajectory or
        // the built-in one for its version, at `bind`
        for (_, file_type, claimed) in schema_lines(text)? {
            let Some(schema) = traj.layouts.get(&file_type) else {
                return Err(format!(
                    "library schema section edited after generation: \
                     # energy-schema names {file_type}, which the library does not declare"
                ));
            };
            let actual = schema.fingerprint();
            if actual != claimed {
                return Err(format!(
                    "library schema section edited after generation: \
                     {file_type} declares {actual}, # energy-schema says {claimed}"
                ));
            }
        }
        Ok(traj)
    }

    /// Establish the layout of a trajectory about to be read as `file_type` and check the
    /// library against it (profile R5). Returns the tier that applied; the error is the R5
    /// diff. `lib_name`/`traj_name` only label the messages.
    pub fn bind(
        &mut self,
        file_type: &str,
        preamble: &Preamble,
        lib_name: &str,
        traj_name: &str,
    ) -> Result<Tier, String> {
        // gromos++ `check_version`, verbatim: a mismatch is only ever a warning there
        let (lv, tv) = (self.version.as_deref(), preamble.version.as_deref());
        if lv != tv {
            let what = if lv.is_some() && tv.is_some() {
                "mismatch"
            } else {
                "missing"
            };
            self.warnings.push(format!(
                "WARNING: Version number {what}!\n         Library {lib_name} version: {}\n         \
                 Energy Trajectory {traj_name} version: {}",
                lv.unwrap_or("NONE"),
                tv.unwrap_or("NONE")
            ));
        }
        let footer = |traj: &EnergyTraj| {
            format!(
                "  library:    {lib_name}  (ENEVERSION {})\n  trajectory: {traj_name}  ({})",
                traj.version.as_deref().unwrap_or("NONE"),
                preamble.describe()
            )
        };
        if let Some(theirs) = &preamble.schema {
            if theirs.file_type() != file_type {
                return Err(format!(
                    "{traj_name} describes itself as {}, but is being read as {file_type}",
                    theirs.file_type()
                ));
            }
            match self.layouts.get(file_type) {
                Some(mine) if mine.fingerprint() != theirs.fingerprint() => {
                    return Err(format!(
                        "energy library does not describe this trajectory\n{}\n{}",
                        mine.diff(theirs, ("library", "file")).join("\n"),
                        footer(self)
                    ));
                },
                Some(_) => {},
                None => self.set_schema(theirs.clone()),
            }
            self.schedule_check(file_type, Tier::SelfDescribed);
            return Ok(Tier::SelfDescribed);
        }
        if let Some(builtin) = tv.and_then(|v| builtin_schema(v, file_type)) {
            let version = tv.unwrap_or_default().to_string();
            match self.layouts.get(file_type) {
                Some(mine) if mine.fingerprint() != builtin.fingerprint() => {
                    return Err(format!(
                        "energy library does not match the {version} layout\n{}\n{}",
                        mine.diff(&builtin, ("library", "built-in")).join("\n"),
                        footer(self)
                    ));
                },
                Some(_) => {},
                None => self.set_schema(builtin),
            }
            self.schedule_check(file_type, Tier::KnownVersion(version.clone()));
            return Ok(Tier::KnownVersion(version));
        }
        if !self.layouts.contains_key(file_type) {
            return Err(format!(
                "{traj_name}: no way to establish the layout — the library declares no \
                 {file_type}, the trajectory carries no self-description and its \
                 {} has no built-in layout",
                preamble.describe()
            ));
        }
        self.warnings.push(format!(
            "WARNING: layout unverified: {traj_name} ({}) is read as {lib_name} declares it",
            preamble.describe()
        ));
        self.schedule_check(file_type, Tier::Unverified);
        Ok(Tier::Unverified)
    }

    fn schedule_check(&mut self, file_type: &str, tier: Tier) {
        self.unchecked.retain(|(t, _)| t != file_type);
        self.unchecked.push((file_type.to_string(), tier));
    }

    /// Warnings raised so far (version mismatches, unverified layouts, failed invariants),
    /// oldest first; the caller prints them.
    pub fn drain_warnings(&mut self) -> Vec<String> {
        std::mem::take(&mut self.warnings)
    }

    /// The declaration of `file_type`, if the library carried one.
    pub fn schema(&self, file_type: &str) -> Option<&Schema> {
        self.layouts.get(file_type)
    }

    /// Install (or replace) the declaration of a file type.
    pub fn set_schema(&mut self, schema: Schema) {
        for e in &schema.entries {
            match e {
                Entry::Block(_) => {},
                // gromos++ registers the plain and the matrix form together
                Entry::Size(n) => {
                    self.sizes.entry(n.clone()).or_insert(0);
                    self.sizes.entry(format!("matrix_{n}")).or_insert(0);
                },
                Entry::Sub { name, .. } => {
                    self.data.entry(name.clone()).or_default();
                },
            }
        }
        self.layouts.insert(schema.file_type.clone(), schema);
    }

    /// A named constant available to every expression (`ene_ana` adds `MASS` and `NUMMOL`
    /// from the topology).
    pub fn add_constant(&mut self, name: &str, value: f64) {
        self.constants.insert(name.to_string(), value);
    }

    pub fn version(&self) -> Option<&str> {
        self.version.as_deref()
    }

    /// The subblock currently held, as `[row][column]` — the shape the library declares, so
    /// `NONBONDED` is one row per energy-group *pair* in gromosXX's order ([`group_pair`]).
    /// `None` if the library does not declare it.
    pub fn table(&self, subblock: &str) -> Option<&[Vec<f64>]> {
        self.data.get(subblock).map(|t| t.as_slice())
    }

    /// The subblock names the library declares, sorted.
    pub fn subblock_names(&self) -> Vec<&str> {
        let mut v: Vec<&str> = self.data.keys().map(|s| s.as_str()).collect();
        v.sort_unstable();
        v
    }

    /// The value of a `size` in the frame currently held (`NUM_ENERGY_GROUPS`, `NUM_BATHS`,
    /// and the `matrix_` forms).
    pub fn size(&self, name: &str) -> Option<usize> {
        self.sizes.get(name).copied()
    }

    /// The names a library declared, for `@library print`.
    pub fn variable_names(&self) -> Vec<&str> {
        let mut v: Vec<&str> = self.variables.keys().map(|s| s.as_str()).collect();
        v.sort_unstable();
        v
    }

    fn dim(&self, d: &Dim) -> usize {
        match d {
            Dim::Fixed(n) => *n,
            Dim::Named(n) => self.sizes.get(n).copied().unwrap_or(0),
        }
    }

    /// Read the next frame of `file_type` (`ENERTRJ`/`FRENERTRJ`); `Ok(false)` at end of file.
    pub fn read_frame(&mut self, reader: &mut TextReader, file_type: &str) -> Result<bool, String> {
        let layout = self
            .layouts
            .get(file_type)
            .ok_or_else(|| format!("don't know how to read from file type {file_type}"))?
            .entries
            .clone();

        let mut values: std::vec::IntoIter<f64> = Vec::new().into_iter();
        let mut block = String::new();
        let mut block_len = 0usize;
        // gromosXX writes `# totals`, `# baths`, … between the tables of a block; a declared
        // table that straddles one is mis-sized, whatever the library says
        let mut markers: Vec<(usize, String)> = Vec::new();
        // gromos++ guards against a library that does not describe the file it is reading
        // (`EnergyTraj::read_block`); without it a misaligned `size` silently yields empty
        // subblocks and the mismatch only surfaces as a missing value much later.
        let mut check_leftover: Option<String> = None;
        for entry in &layout {
            match entry {
                Entry::Block(name) => {
                    let Some(read) = read_gromos_block(reader)? else {
                        return Ok(false);
                    };
                    if read.name != *name {
                        return Err(format!("Block {name} expected. Got {}", read.name));
                    }
                    let mut nums = Vec::new();
                    for tok in read.body.iter().flat_map(|l| l.split_whitespace()) {
                        nums.push(
                            tok.parse::<f64>()
                                .map_err(|_| format!("{name}: '{tok}' is not a number"))?,
                        );
                    }
                    if let Some(leftover) = check_leftover.take() {
                        return Err(leftover);
                    }
                    block_len = nums.len();
                    values = nums.into_iter();
                    markers = read.markers;
                    block = name.clone();
                },
                Entry::Size(name) => {
                    let raw = values
                        .next()
                        .ok_or_else(|| format!("Tried to read an integer for {name}"))?;
                    // a misaligned library lands an energy in a size slot; refuse rather than
                    // allocate a table of that many rows
                    let x = size_value(raw)
                        .ok_or_else(|| format!("size {name}: expected an integer, got {raw}"))?;
                    if x > values.len() {
                        return Err(format!(
                            "size {name}: {x} rows, but only {} values left in the block",
                            values.len()
                        ));
                    }
                    // per file type: md++ writes `numaccelgroups` 1 in the .tre and 0 in the
                    // .trg of the same run, and gromos++ never compares the two
                    let key = (file_type.to_string(), name.clone());
                    match self.first_sizes.get(&key) {
                        None => {
                            self.first_sizes.insert(key, x);
                        },
                        Some(&first) if first != x => {
                            return Err(format!(
                                "size {name}: {x} in this frame, {first} in the first"
                            ));
                        },
                        Some(_) => {},
                    }
                    self.sizes.insert(name.clone(), x);
                    self.sizes.insert(format!("matrix_{name}"), x * (x + 1) / 2);
                },
                Entry::Sub { name, rows, cols } => {
                    let (si, sj) = (self.dim(rows), self.dim(cols));
                    if si.saturating_mul(sj) > values.len() {
                        return Err(format!(
                            "subblock {name}: needs {si} x {sj} values, block has {}",
                            values.len()
                        ));
                    }
                    let (start, end) =
                        (block_len - values.len(), block_len - values.len() + si * sj);
                    if let Some((at, marker)) = markers.iter().find(|(p, _)| start < *p && *p < end)
                    {
                        return Err(format!(
                            "\"{marker}\" falls inside subblock {name} ({si} x {sj} values from \
                             value {} of {block}; the marker is at value {})",
                            start + 1,
                            at + 1
                        ));
                    }
                    let mut table = vec![vec![0.0; sj]; si];
                    for row in table.iter_mut() {
                        for cell in row.iter_mut() {
                            *cell = values
                                .next()
                                .ok_or_else(|| format!("Not enough values in {name}"))?;
                        }
                    }
                    self.data.insert(name.clone(), table);
                },
            }
            if !block.is_empty() {
                let rest = values.len();
                check_leftover = (rest > 0).then(|| {
                    format!(
                        "Block definition does not agree with trajectory data for {block}: {rest} leftover values"
                    )
                });
            }
        }
        if let Some(leftover) = check_leftover {
            return Err(leftover);
        }
        if let Some(i) = self.unchecked.iter().position(|(t, _)| t == file_type) {
            let (_, tier) = self.unchecked.remove(i);
            self.check_invariants(file_type, &tier);
        }
        Ok(true)
    }

    /// The identities every md++ frame satisfies, checked on the first frame after a bind. A
    /// reshaped table keeps the value count and passes every other check, but breaks these;
    /// with a verified layout a failure means the writer composed the totals inconsistently.
    fn check_invariants(&mut self, file_type: &str, tier: &Tier) {
        let get = |name: &str, i: usize| {
            self.data
                .get(name)
                .and_then(|t| t.get(i))
                .and_then(|r| r.first())
                .copied()
        };
        let mut failed = Vec::new();
        if let (Some(total), Some(kin), Some(pot), Some(special)) = (
            get("ENER", 0),
            get("ENER", 1),
            get("ENER", 2),
            get("ENER", 20),
        ) {
            let sum = kin + pot + special;
            if (total - sum).abs() > 1e-6 * total.abs().max(1.0) {
                failed.push(format!(
                    "ENER[1] is not ENER[2] + ENER[3] + ENER[21]: {total} vs {sum}"
                ));
            }
            if kin < 0.0 {
                failed.push(format!("ENER[2] (kinetic energy) is negative: {kin}"));
            }
        }
        if let Some(vol) = get("VOLUME", 0) {
            if vol <= 0.0 {
                failed.push(format!("VOLUME[1] is not positive: {vol}"));
            }
        }
        let because = match tier {
            Tier::Unverified => "the library's layout is unverified and may mislabel the values",
            _ => "the layout is verified, so the writer composed these totals inconsistently",
        };
        for f in failed {
            self.warnings
                .push(format!("WARNING: {file_type} first frame: {f}; {because}"));
        }
    }

    /// The value of a property: `SUBBLOCK[i]`, `SUBBLOCK[i][j]`, a library variable, or a
    /// constant. Indices are 1-based, as in the library.
    pub fn value(&self, prop: &str) -> Result<f64, String> {
        let expr = parse_expression(prop)
            .map_err(|e| format!("Trying to access an unknown variable {prop}: {e}"))?;
        self.eval(&expr, &mut Vec::new())
    }

    fn eval(&self, e: &Expr, stack: &mut Vec<String>) -> Result<f64, String> {
        match e {
            Expr::Const(v) => Ok(*v),
            Expr::Name { name, i, j } => {
                if let Some(table) = self.data.get(name) {
                    return table
                        .get(*i)
                        .and_then(|r| r.get(*j))
                        .copied()
                        .ok_or_else(|| {
                            format!("{name}[{}][{}] is outside the block read", i + 1, j + 1)
                        });
                }
                if let Some(v) = self.constants.get(name) {
                    return Ok(*v);
                }
                if let Some(v) = self.sizes.get(name) {
                    return Ok(*v as f64);
                }
                if let Some(sub) = self.variables.get(name) {
                    if stack.iter().any(|s| s == name) {
                        return Err(format!("variable {name} is defined in terms of itself"));
                    }
                    stack.push(name.clone());
                    let v = self.eval(sub, stack);
                    stack.pop();
                    return v;
                }
                Err(format!("Trying to access an unknown variable {name}"))
            },
            Expr::Fun(f, a) => {
                let v = self.eval(a, stack)?;
                Ok(match f {
                    Fun::Sin => v.sin(),
                    Fun::Cos => v.cos(),
                    Fun::Exp => v.exp(),
                    Fun::Log => v.ln(),
                })
            },
            Expr::Bin(op, a, b) => {
                let (a, b) = (self.eval(a, stack)?, self.eval(b, stack)?);
                Ok(match op {
                    Op::Add => a + b,
                    Op::Sub => a - b,
                    Op::Mul => a * b,
                    Op::Div => a / b,
                })
            },
        }
    }
}

/// A `size` slot as an integer: finite, non-negative and whole.
fn size_value(raw: f64) -> Option<usize> {
    (raw.is_finite() && raw >= 0.0 && raw.fract() == 0.0 && raw < u32::MAX as f64)
        .then_some(raw as usize)
}

/// Split a GROMOS text file into `(block name, body lines)`; comments and blank lines dropped.
fn blocks(text: &str) -> Vec<(String, Vec<String>)> {
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut body: Vec<String> = Vec::new();
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if line == "END" {
            if let Some(n) = name.take() {
                out.push((n, std::mem::take(&mut body)));
            }
            continue;
        }
        match &name {
            None => name = Some(line.to_string()),
            Some(_) => body.push(line.to_string()),
        }
    }
    out
}

/// A GROMOS block as read: its data lines, and the one-word comment lines gromosXX puts
/// between tables (`# baths`), each with the number of values that precede it.
struct RawBlock {
    name: String,
    body: Vec<String>,
    markers: Vec<(usize, String)>,
}

/// Read the next GROMOS block from a stream: `Ok(None)` at end of file.
fn read_gromos_block(reader: &mut TextReader) -> Result<Option<RawBlock>, String> {
    let mut buf = String::new();
    let name = loop {
        if reader.read_line(&mut buf).map_err(|e| e.to_string())? == 0 {
            return Ok(None);
        }
        let line = buf.split('#').next().unwrap_or("").trim();
        if !line.is_empty() {
            break line.to_string();
        }
    };
    let mut body: Vec<String> = Vec::new();
    let mut markers = Vec::new();
    let mut count = 0usize;
    loop {
        if reader.read_line(&mut buf).map_err(|e| e.to_string())? == 0 {
            return Err(format!("no END in {name} block"));
        }
        let trimmed = buf.trim();
        if let Some(comment) = trimmed.strip_prefix('#') {
            let word = comment.trim();
            if !word.is_empty() && !word.contains(char::is_whitespace) {
                markers.push((count, format!("# {word}")));
            }
            continue;
        }
        let line = trimmed.split('#').next().unwrap_or("").trim();
        if line == "END" {
            break;
        }
        if !line.is_empty() {
            count += line.split_whitespace().count();
            body.push(line.to_string());
        }
    }
    Ok(Some(RawBlock {
        name,
        body,
        markers,
    }))
}

/// `name = expression` pairs from a `VARIABLES` body. gromos++ reads the block as one token
/// stream and cuts each definition at the token before the next `=`.
fn variable_definitions(body: &[String]) -> Vec<(String, String)> {
    let tokens: Vec<&str> = body.iter().flat_map(|l| l.split_whitespace()).collect();
    let eq: Vec<usize> = (0..tokens.len()).filter(|&i| tokens[i] == "=").collect();
    let mut out = Vec::new();
    for (k, &i) in eq.iter().enumerate() {
        if i == 0 {
            continue;
        }
        // up to the name of the next definition, or the end
        let end = eq.get(k + 1).map_or(tokens.len(), |&next| next - 1);
        let rhs = tokens[i + 1..end].join(" ");
        if !rhs.is_empty() {
            out.push((tokens[i - 1].to_string(), rhs));
        }
    }
    out
}

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Op(Op),
    Fun(Fun),
    Open,
    Close,
    Atom(String),
}

fn tokenize(s: &str) -> Vec<Tok> {
    s.split_whitespace()
        .map(|t| match t {
            "+" => Tok::Op(Op::Add),
            "-" => Tok::Op(Op::Sub),
            "*" => Tok::Op(Op::Mul),
            "/" => Tok::Op(Op::Div),
            "(" => Tok::Open,
            ")" => Tok::Close,
            "sin" => Tok::Fun(Fun::Sin),
            "cos" => Tok::Fun(Fun::Cos),
            "exp" => Tok::Fun(Fun::Exp),
            "log" => Tok::Fun(Fun::Log),
            other => Tok::Atom(other.to_string()),
        })
        .collect()
}

fn parse_expression(s: &str) -> Result<Expr, String> {
    let toks = tokenize(s);
    if toks.is_empty() {
        return Err("empty expression".into());
    }
    parse_tokens(&toks)
}

/// gromos++ `Expression::calc`: split at the rightmost top-level `+`/`-`, else at the
/// *leftmost* top-level `*`/`/` — see the module note on associativity.
fn parse_tokens(toks: &[Tok]) -> Result<Expr, String> {
    if toks.is_empty() {
        return Err("empty expression".into());
    }
    let mut depth = 0i32;
    let mut add_sub: Option<usize> = None;
    let mut mul_div: Option<usize> = None;
    for (i, t) in toks.iter().enumerate() {
        match t {
            Tok::Open => depth += 1,
            Tok::Close => depth -= 1,
            Tok::Op(op) if depth == 0 && i > 0 => match op {
                Op::Add | Op::Sub => add_sub = Some(i), // keep the rightmost
                Op::Mul | Op::Div => {
                    if mul_div.is_none() {
                        mul_div = Some(i) // keep the leftmost
                    }
                },
            },
            _ => {},
        }
    }
    if depth != 0 {
        return Err("unbalanced brackets".into());
    }
    if let Some(k) = add_sub.or(mul_div) {
        let Tok::Op(op) = toks[k] else { unreachable!() };
        return Ok(Expr::Bin(
            op,
            Box::new(parse_tokens(&toks[..k])?),
            Box::new(parse_tokens(&toks[k + 1..])?),
        ));
    }
    parse_atom(toks)
}

fn parse_atom(toks: &[Tok]) -> Result<Expr, String> {
    match toks.first() {
        None => Err("empty expression".into()),
        Some(Tok::Op(_)) => Err("an expression cannot start with an operator".into()),
        Some(Tok::Fun(f)) => Ok(Expr::Fun(*f, Box::new(parse_atom(&toks[1..])?))),
        Some(Tok::Open) => {
            if toks.last() != Some(&Tok::Close) {
                return Err("unbalanced brackets".into());
            }
            parse_tokens(&toks[1..toks.len() - 1])
        },
        Some(Tok::Close) => Err("unbalanced brackets".into()),
        Some(Tok::Atom(a)) => {
            if toks.len() > 1 {
                return Err(format!(
                    "two values in a row ({a}) — an operator is missing"
                ));
            }
            if let Ok(v) = a.parse::<f64>() {
                return Ok(Expr::Const(v));
            }
            parse_reference(a)
        },
    }
}

/// `NAME`, `NAME[i]` or `NAME[i][j]`, with 1-based indices as the library writes them.
fn parse_reference(s: &str) -> Result<Expr, String> {
    let Some(open) = s.find('[') else {
        return Ok(Expr::Name {
            name: s.to_string(),
            i: 0,
            j: 0,
        });
    };
    let name = s[..open].to_string();
    let mut idx = Vec::new();
    for part in s[open..].split('[').skip(1) {
        let close = part
            .find(']')
            .ok_or_else(|| format!("missing ']' in {s}"))?;
        let n: usize = part[..close]
            .trim()
            .parse()
            .map_err(|_| format!("index in {s} is not an integer"))?;
        if n == 0 {
            return Err(format!("Indexes start at 1, invalid variable {s}"));
        }
        idx.push(n - 1);
    }
    match idx.len() {
        1 => Ok(Expr::Name {
            name,
            i: idx[0],
            j: 0,
        }),
        2 => Ok(Expr::Name {
            name,
            i: idx[0],
            j: idx[1],
        }),
        _ => Err(format!("{s} has more than two indices")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LIB: &str = "\
TITLE
a library
END
ENEVERSION
  2023-04-15
END
ENERTRJ
  block TIMESTEP
    subblock TIME 2 1
  block ENERGY03
    subblock ENER 4 1
    size NUM_BATHS
    subblock KINENER NUM_BATHS 3
    size NUM_ENERGY_GROUPS
    subblock NONBONDED matrix_NUM_ENERGY_GROUPS 2
END
VARIABLES
time = TIME[2]
totene = ENER[1]
totpot = ENER[3]
sum = ENER[1] + ENER[3]
scaled = ENER[1] * 2 / 4
nested = totene + totpot
END
";

    /// A unique scratch path: the tests run in parallel and each deletes its own file.
    fn scratch(tag: &str) -> std::path::PathBuf {
        static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        std::env::temp_dir().join(format!(
            "gromos_etrj_{tag}_{}_{}.tre",
            std::process::id(),
            N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
        ))
    }

    fn frame() -> (EnergyTraj, TextReader) {
        // two baths (3 values each), two energy groups → 3 pairs of 2 values
        let text = "\
TIMESTEP
   1  0.002
END
ENERGY03
  10.0 20.0 30.0 40.0
  2
  1 2 3   4 5 6
  2
  7 8   9 10   11 12
END
";
        let path = scratch("frame");
        std::fs::write(&path, text).unwrap();
        let reader = TextReader::open(&path).unwrap();
        std::fs::remove_file(&path).ok();
        (EnergyTraj::from_library(LIB).unwrap(), reader)
    }

    #[test]
    fn library_declares_the_layout_and_the_names() {
        let traj = EnergyTraj::from_library(LIB).unwrap();
        assert_eq!(traj.version(), Some("2023-04-15"));
        assert_eq!(
            traj.variable_names(),
            ["nested", "scaled", "sum", "time", "totene", "totpot"]
        );
    }

    #[test]
    fn a_frame_fills_the_declared_subblocks() {
        let (mut traj, mut reader) = frame();
        assert!(traj.read_frame(&mut reader, "ENERTRJ").unwrap());

        assert_eq!(traj.value("TIME[2]").unwrap(), 0.002);
        assert_eq!(traj.value("totene").unwrap(), 10.0);
        assert_eq!(traj.value("totpot").unwrap(), 30.0);
        // the size read from the stream sets the following subblock's row count
        assert_eq!(traj.value("KINENER[2][3]").unwrap(), 6.0);
        // matrix_N expands 2 energy groups to 3 pairs
        assert_eq!(traj.value("NONBONDED[3][2]").unwrap(), 12.0);
        assert_eq!(traj.value("sum").unwrap(), 40.0);
        assert_eq!(traj.value("nested").unwrap(), 40.0);

        // one frame only
        assert!(!traj.read_frame(&mut reader, "ENERTRJ").unwrap());
    }

    #[test]
    fn indices_are_one_based_and_bounds_checked() {
        let (mut traj, mut reader) = frame();
        traj.read_frame(&mut reader, "ENERTRJ").unwrap();
        assert!(traj.value("ENER[0]").unwrap_err().contains("start at 1"));
        assert!(traj.value("ENER[5]").unwrap_err().contains("outside"));
        assert!(traj.value("NOSUCH").unwrap_err().contains("unknown"));
    }

    #[test]
    fn a_block_out_of_order_is_an_error() {
        let text = "ENERGY03\n 1 2 3 4\n 1\n 1 2 3\n 1\n 1 2\nEND\n";
        let path = scratch("order");
        std::fs::write(&path, text).unwrap();
        let mut reader = TextReader::open(&path).unwrap();
        std::fs::remove_file(&path).ok();
        let mut traj = EnergyTraj::from_library(LIB).unwrap();
        let err = traj.read_frame(&mut reader, "ENERTRJ").unwrap_err();
        assert!(err.contains("TIMESTEP expected"), "{err}");
    }

    /// gromos++ splits at the rightmost `+`/`-` (left-associative) but the leftmost `*`/`/`
    /// (right-associative). Both are reproduced.
    #[test]
    fn operator_associativity_matches_gromospp() {
        let ev = |s: &str| {
            EnergyTraj::default()
                .eval(&parse_expression(s).unwrap(), &mut Vec::new())
                .unwrap()
        };
        assert_eq!(ev("10 - 4 - 3"), 3.0); // (10-4)-3, not 10-(4-3)
        assert_eq!(ev("2 / 4 / 8"), 4.0); // 2/(4/8), the gromos++ grouping
        assert_eq!(ev("1 + 2 * 3"), 7.0); // + binds loosest
        assert_eq!(ev("( 1 + 2 ) * 3"), 9.0);
        assert_eq!(ev("exp ( 0 ) + 4"), 5.0);
    }

    // ---- profile: preamble, tiers, markers, invariants ----

    fn reader_of(text: &str) -> TextReader {
        let path = scratch("pre");
        std::fs::write(&path, text).unwrap();
        let reader = TextReader::open(&path).unwrap();
        std::fs::remove_file(&path).ok();
        reader
    }

    const XX_HEADER: &str = "TITLE\nMD++\nEND\nENEVERSION\n\t2023-04-15\nEND\n";

    /// The R3/R8 lines a gromos-rs writer adds, for the md++ ENERTRJ layout.
    fn self_description() -> String {
        let mut s = String::new();
        for l in self_description_lines(&md_schema("ENERTRJ")) {
            let _ = writeln!(s, "{l}");
        }
        let prov = Provenance {
            writer: "gromos-rs md 0.0.47".into(),
            ..Default::default()
        };
        for l in prov.lines() {
            let _ = writeln!(s, "{l}");
        }
        s
    }

    /// One md++-shaped frame: one bath, one energy group, no EDS/GaMD/λ tables.
    fn md_frame(total: f64, kin: f64, pot: f64) -> String {
        let mut s = String::from("TIMESTEP\n 1 0.002\nEND\nENERGY03\n# totals\n");
        let mut ener = vec![0.0; 52];
        ener[0] = total;
        ener[1] = kin;
        ener[2] = pot;
        for v in ener {
            let _ = writeln!(s, " {v}");
        }
        s.push_str("# baths\n1\n 1 2 3\n# bonded\n1\n 1 2 3 4 5\n# nonbonded\n 1 2 3 4 5 6\n");
        s.push_str("# special\n 1 2 3 4 5 6 7 8 9 10 11 12 13\n# eds\n# numstates\n0\n# gamd\n# numaccelgroups\n0\n");
        s.push_str("# precalclam\n# nr_lambdas\n0\n# ABdih\n 0 0\nEND\n");
        s.push_str("VOLUMEPRESSURE03\n# mass\n 100\n# temperature\n1\n 300 300 300 1\n# volume\n 27 3 0 0 0 3 0 0 0 3\n# pressure\n");
        for _ in 0..30 {
            s.push_str(" 0\n");
        }
        s.push_str("END\n");
        s
    }

    #[test]
    fn preamble_of_a_gromosxx_file_has_version_only() {
        let mut r = reader_of(&format!("{XX_HEADER}{}", md_frame(3.0, 1.0, 2.0)));
        let pre = read_preamble(&mut r).unwrap();
        assert_eq!(pre.version.as_deref(), Some("2023-04-15"));
        assert_eq!(pre.title, ["MD++"]);
        assert!(pre.schema.is_none() && pre.provenance.is_empty());
        // the reader is left at the first frame
        let mut traj = EnergyTraj::from_library(OFFICIAL_LIBRARY).unwrap();
        assert!(traj.read_frame(&mut r, "ENERTRJ").unwrap());
        assert_eq!(traj.value("totene").unwrap(), 3.0);
    }

    #[test]
    fn preamble_self_description_is_parsed_and_verified() {
        let mut r = reader_of(&format!(
            "{XX_HEADER}{}{}",
            self_description(),
            md_frame(3.0, 1.0, 2.0)
        ));
        let pre = read_preamble(&mut r).unwrap();
        assert_eq!(pre.schema.as_ref(), Some(&md_schema("ENERTRJ")));
        assert_eq!(pre.writer(), Some("gromos-rs md 0.0.47"));

        let forged = self_description().replace("subblock ENER 52 1", "subblock ENER 26 2");
        let err = read_preamble(&mut reader_of(&format!("{XX_HEADER}{forged}"))).unwrap_err();
        assert!(
            err.starts_with("trajectory self-description is corrupt"),
            "{err}"
        );

        let headless: String = self_description()
            .lines()
            .skip(1)
            .map(|l| format!("{l}\n"))
            .collect();
        let err = read_preamble(&mut reader_of(&format!("{XX_HEADER}{headless}"))).unwrap_err();
        assert!(err.contains("without an energy-schema line"), "{err}");

        let newer = self_description().replace("energy-schema 1", "energy-schema 7");
        let err = read_preamble(&mut reader_of(&format!("{XX_HEADER}{newer}"))).unwrap_err();
        assert!(err.starts_with("energy-schema 7 is newer"), "{err}");
    }

    /// The official library without its `# energy-schema` lines: free to edit.
    fn unsigned_library() -> String {
        OFFICIAL_LIBRARY
            .lines()
            .filter(|l| !l.starts_with("# energy-schema"))
            .map(|l| format!("{l}\n"))
            .collect()
    }

    /// The LiveCoMS tutorial's copy: an older layout under the same `ENEVERSION`.
    fn stale_library() -> String {
        unsigned_library()
            .replace("subblock ENER 52 1", "subblock ENER 50 1")
            .replace(
                "subblock SPECIAL NUM_ENERGY_GROUPS 13",
                "subblock SPECIAL NUM_ENERGY_GROUPS 12",
            )
    }

    #[test]
    fn tier_a_the_file_is_the_authority() {
        let mut r = reader_of(&format!(
            "{XX_HEADER}{}{}",
            self_description(),
            md_frame(3.0, 1.0, 2.0)
        ));
        let pre = read_preamble(&mut r).unwrap();

        let mut lib = EnergyTraj::from_library(OFFICIAL_LIBRARY).unwrap();
        assert_eq!(
            lib.bind("ENERTRJ", &pre, "official", "md.tre").unwrap(),
            Tier::SelfDescribed
        );
        assert!(lib.drain_warnings().is_empty());

        let mut stale = EnergyTraj::from_library(&stale_library()).unwrap();
        let err = stale
            .bind("ENERTRJ", &pre, "tutorial/ene_ana.md++.lib", "md.tre")
            .unwrap_err();
        assert_eq!(
            err,
            "energy library does not describe this trajectory\n\
             \x20 ENERGY03 ENER:   library 50 x 1,   file 52 x 1\n\
             \x20 ENERGY03 SPECIAL:   library NUM_ENERGY_GROUPS x 12,   file NUM_ENERGY_GROUPS x 13\n\
             \x20 library:    tutorial/ene_ana.md++.lib  (ENEVERSION 2023-04-15)\n\
             \x20 trajectory: md.tre  (ENEVERSION 2023-04-15, written by gromos-rs md 0.0.47)"
        );

        // a VARIABLES-only library takes the file's layout
        let mut vars = EnergyTraj::from_library("VARIABLES\n kin = ENER[2]\nEND\n").unwrap();
        assert_eq!(
            vars.bind("ENERTRJ", &pre, "vars.lib", "md.tre").unwrap(),
            Tier::SelfDescribed
        );
        assert!(vars.read_frame(&mut r, "ENERTRJ").unwrap());
        assert_eq!(vars.value("kin").unwrap(), 1.0);
        // the version warning is gromos++'s, verbatim
        assert_eq!(
            vars.drain_warnings(),
            ["WARNING: Version number missing!\n         Library vars.lib version: NONE\n         Energy Trajectory md.tre version: 2023-04-15"]
        );
    }

    #[test]
    fn tier_b_a_known_version_is_checked_against_the_builtin_layout() {
        let mut r = reader_of(&format!("{XX_HEADER}{}", md_frame(3.0, 1.0, 2.0)));
        let pre = read_preamble(&mut r).unwrap();

        let mut lib = EnergyTraj::from_library(OFFICIAL_LIBRARY).unwrap();
        assert_eq!(
            lib.bind("ENERTRJ", &pre, "official", "md.tre").unwrap(),
            Tier::KnownVersion("2023-04-15".into())
        );
        assert!(lib.read_frame(&mut r, "ENERTRJ").unwrap());
        assert!(lib.drain_warnings().is_empty());

        let mut stale = EnergyTraj::from_library(&stale_library()).unwrap();
        let err = stale
            .bind("ENERTRJ", &pre, "old.lib", "md.tre")
            .unwrap_err();
        assert!(err.starts_with("energy library does not match the 2023-04-15 layout\n  ENERGY03 ENER:   library 50 x 1,   built-in 52 x 1\n"), "{err}");

        let mut vars = EnergyTraj::from_library(
            "ENEVERSION\n2023-04-15\nEND\nVARIABLES\n kin = ENER[2]\nEND\n",
        )
        .unwrap();
        assert!(vars.bind("ENERTRJ", &pre, "vars.lib", "md.tre").is_ok());
        assert_eq!(vars.schema("ENERTRJ"), Some(&md_schema("ENERTRJ")));
    }

    #[test]
    fn tier_c_takes_the_library_on_trust_and_says_so() {
        let old = format!(
            "TITLE\nold\nEND\nENEVERSION\n 2011-04-05\nEND\n{}",
            md_frame(3.0, 1.0, 2.0)
        );
        let mut r = reader_of(&old);
        let pre = read_preamble(&mut r).unwrap();
        let mut lib = EnergyTraj::from_library(OFFICIAL_LIBRARY).unwrap();
        assert_eq!(
            lib.bind("ENERTRJ", &pre, "official", "old.tre").unwrap(),
            Tier::Unverified
        );
        assert!(lib.read_frame(&mut r, "ENERTRJ").unwrap());
        let w = lib.drain_warnings();
        assert!(
            w[0].starts_with("WARNING: Version number mismatch!"),
            "{}",
            w[0]
        );
        assert_eq!(w[1], "WARNING: layout unverified: old.tre (ENEVERSION 2011-04-05) is read as official declares it");
        assert_eq!(w.len(), 2, "{w:?}");

        let mut vars = EnergyTraj::from_library("VARIABLES\n kin = ENER[2]\nEND\n").unwrap();
        let err = vars
            .bind("ENERTRJ", &pre, "vars.lib", "old.tre")
            .unwrap_err();
        assert!(err.contains("no way to establish the layout"), "{err}");
    }

    #[test]
    fn a_marker_inside_a_table_is_a_mis_sized_table() {
        // the stale library reads 50 totals, then takes `1` (the bath count) and the first
        // bath value as totals 51 and 52 — but `# baths` sits between them
        let mut r = reader_of(&format!("TITLE\nx\nEND\n{}", md_frame(3.0, 1.0, 2.0)));
        let pre = read_preamble(&mut r).unwrap();
        let mut stale = EnergyTraj::from_library(
            &unsigned_library().replace("subblock ENER 52 1", "subblock ENER 53 1"),
        )
        .unwrap();
        stale.bind("ENERTRJ", &pre, "stale.lib", "x.tre").unwrap();
        let err = stale.read_frame(&mut r, "ENERTRJ").unwrap_err();
        assert_eq!(err, "\"# baths\" falls inside subblock ENER (53 x 1 values from value 1 of ENERGY03; the marker is at value 53)");
    }

    #[test]
    fn first_frame_invariants_catch_a_reshape() {
        // ENER 26 x 2 reads the same 52 values; ENER[2] is now the old ENER[3]
        let mut r = reader_of(&format!("TITLE\nx\nEND\n{}", md_frame(3.0, 1.0, 2.0)));
        let pre = read_preamble(&mut r).unwrap();
        let mut lib = EnergyTraj::from_library(
            &unsigned_library().replace("subblock ENER 52 1", "subblock ENER 26 2"),
        )
        .unwrap();
        lib.bind("ENERTRJ", &pre, "reshaped.lib", "x.tre").unwrap();
        assert!(lib.read_frame(&mut r, "ENERTRJ").unwrap());
        let w = lib.drain_warnings();
        assert!(
            w.iter().any(
                |w| w.contains("ENER[1] is not ENER[2] + ENER[3] + ENER[21]: 3 vs 2")
                    && w.contains("unverified")
            ),
            "{w:?}"
        );

        // a verified layout whose totals do not add up blames the writer
        let mut r = reader_of(&format!("{XX_HEADER}{}", md_frame(3.0, 1.0, 2.5)));
        let pre = read_preamble(&mut r).unwrap();
        let mut lib = EnergyTraj::from_library(OFFICIAL_LIBRARY).unwrap();
        lib.bind("ENERTRJ", &pre, "official", "x.tre").unwrap();
        assert!(lib.read_frame(&mut r, "ENERTRJ").unwrap());
        let w = lib.drain_warnings();
        assert_eq!(w.len(), 1);
        assert!(
            w[0].contains("3 vs 3.5") && w[0].contains("writer"),
            "{}",
            w[0]
        );
    }

    #[test]
    fn a_malformed_expression_is_rejected() {
        assert!(parse_expression("+ 1").is_err());
        assert!(parse_expression("1 2").is_err());
        assert!(parse_expression("( 1 + 2").is_err());
    }

    fn read_text(lib: &str, text: &str) -> Result<EnergyTraj, String> {
        let path = scratch("read");
        std::fs::write(&path, text).unwrap();
        let mut reader = TextReader::open(&path).unwrap();
        std::fs::remove_file(&path).ok();
        let mut traj = EnergyTraj::from_library(lib)?;
        traj.read_frame(&mut reader, "ENERTRJ")?;
        Ok(traj)
    }

    #[test]
    fn a_size_slot_must_hold_an_integer() {
        // the ENER table is one row short, so a bath energy lands in the NUM_BATHS slot
        let short = LIB.replace("subblock ENER 4 1", "subblock ENER 3 1");
        let text = "TIMESTEP\n 1 0.002\nEND\nENERGY03\n 10 20 30 40.5\n 2\n 1 2 3 4 5 6\n 2\n 7 8 9 10 11 12\nEND\n";
        let err = read_text(&short, text).unwrap_err();
        assert_eq!(err, "size NUM_BATHS: expected an integer, got 40.5");

        // an integer, but far more rows than the block holds: refused before allocating
        let text = "TIMESTEP\n 1 0.002\nEND\nENERGY03\n 10 20 30 40\n 2000000000\n 1 2 3\nEND\n";
        let err = read_text(LIB, text).unwrap_err();
        assert!(err.starts_with("size NUM_BATHS: 2000000000 rows"), "{err}");

        // a table wider than what is left in the block
        let wide = LIB.replace(
            "subblock KINENER NUM_BATHS 3",
            "subblock KINENER NUM_BATHS 30",
        );
        let text = "TIMESTEP\n 1 0.002\nEND\nENERGY03\n 10 20 30 40\n 2\n 1 2 3 4 5 6\n 2\n 7 8 9 10 11 12\nEND\n";
        let err = read_text(&wide, text).unwrap_err();
        assert_eq!(err, "subblock KINENER: needs 2 x 30 values, block has 13");
    }

    #[test]
    fn a_size_must_not_change_between_frames() {
        let one = "TIMESTEP\n 1 0.002\nEND\nENERGY03\n 10 20 30 40\n 1\n 1 2 3\n 1\n 7 8\nEND\n";
        let two =
            "TIMESTEP\n 2 0.004\nEND\nENERGY03\n 10 20 30 40\n 2\n 1 2 3 4 5 6\n 1\n 7 8\nEND\n";
        let path = scratch("frames");
        std::fs::write(&path, format!("{one}{two}")).unwrap();
        let mut reader = TextReader::open(&path).unwrap();
        std::fs::remove_file(&path).ok();
        let mut traj = EnergyTraj::from_library(LIB).unwrap();
        assert!(traj.read_frame(&mut reader, "ENERTRJ").unwrap());
        let err = traj.read_frame(&mut reader, "ENERTRJ").unwrap_err();
        assert_eq!(err, "size NUM_BATHS: 2 in this frame, 1 in the first");
    }

    // ---- profile: canonical form, fingerprint, official library ----

    fn md_schema(file_type: &str) -> Schema {
        builtin_schema("2023-04-15", file_type).unwrap()
    }

    #[test]
    fn canonical_form_drops_names_and_numbers_sizes() {
        let canon = md_schema("ENERTRJ").canonical_form();
        assert!(
            canon.starts_with(
                "block TIMESTEP\nsubblock 2 1\nblock ENERGY03\nsubblock 52 1\nsize $1\n\
                 subblock $1 3\nsize $2\nsubblock $2 5\nsubblock matrix($2) 6\nsubblock $2 13\n"
            ),
            "{canon}"
        );
        // the second `size NUM_BATHS`, in VOLUMEPRESSURE03, keeps its number
        assert!(canon.contains("block VOLUMEPRESSURE03\nsubblock 1 1\nsize $1\nsubblock $1 4\n"));
        assert!(
            !canon.contains("NUM_BATHS") && !canon.contains("ENER "),
            "{canon}"
        );
    }

    #[test]
    fn fingerprint_is_blind_to_names_and_sensitive_to_shape() {
        let base = md_schema("ENERTRJ");
        let edit = |from: &str, to: &str| -> Schema {
            let decl = base.declaration_lines().join("\n");
            assert!(decl.contains(from), "{from} not in layout");
            Schema::parse("ENERTRJ", decl.replace(from, to).lines()).unwrap()
        };
        let renamed = edit("subblock ENER 52 1", "subblock TOTALS 52 1");
        let renamed_size = Schema::parse(
            "ENERTRJ",
            base.declaration_lines()
                .join("\n")
                .replace("NUM_BATHS", "NBATH")
                .lines(),
        )
        .unwrap();
        assert_eq!(renamed.fingerprint(), base.fingerprint());
        assert_eq!(renamed_size.fingerprint(), base.fingerprint());
        assert!(renamed.diff(&base, ("a", "b")).is_empty());

        // the same 52 values in a different shape, one value more, a size gone
        let reshaped = edit("subblock ENER 52 1", "subblock ENER 26 2");
        let shifted = edit("subblock ENER 52 1", "subblock ENER 53 1");
        let stale = edit(
            "subblock SPECIAL NUM_ENERGY_GROUPS 13",
            "subblock SPECIAL NUM_ENERGY_GROUPS 12",
        );
        let no_size = Schema::parse(
            "ENERTRJ",
            base.declaration_lines()
                .join("\n")
                .replace("size NUM_LAMBDAS\nsubblock PRECALCLAM NUM_LAMBDAS 12\n", "")
                .lines(),
        )
        .unwrap();
        let prints: Vec<String> = [&base, &reshaped, &shifted, &stale, &no_size]
            .iter()
            .map(|s| s.fingerprint())
            .collect();
        for (i, a) in prints.iter().enumerate() {
            assert!(a.starts_with("sha256:") && a.len() == 7 + 64);
            for b in &prints[i + 1..] {
                assert_ne!(a, b);
            }
        }
        assert_eq!(
            reshaped.diff(&base, ("library", "file")),
            ["  ENERGY03 ENER:   library 26 x 2,   file 52 x 1"]
        );
        assert_eq!(
            no_size.diff(&base, ("library", "file")),
            [
                "  ENERGY03 ABDIH:   library 1 x 2,   file size",
                "  ENERGY03 PRECALCLAM:   library absent,   file NUM_LAMBDAS x 12",
                "  ENERGY03 ABDIH:   library absent,   file 1 x 2",
            ]
        );
    }

    /// The fingerprint of the layout in code equals the fingerprint of md++'s own library —
    /// the file the tutorials tell users to copy. Skipped when the reference checkout is absent.
    #[test]
    fn builtin_layout_matches_mdpp_library() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../.local/gromosXX/md++/data/ene_ana.md++.lib"
        );
        let Ok(text) = std::fs::read_to_string(path) else {
            eprintln!("skipped: {path} not present");
            return;
        };
        let theirs = EnergyTraj::from_library(&text).unwrap();
        for t in ["ENERTRJ", "FRENERTRJ"] {
            let ours = md_schema(t);
            let theirs = theirs.schema(t).unwrap();
            assert!(
                ours.diff(theirs, ("gromos-rs", "md++")).is_empty(),
                "{t}: {:#?}",
                ours.diff(theirs, ("gromos-rs", "md++"))
            );
            assert_eq!(ours.fingerprint(), theirs.fingerprint(), "{t}");
        }
    }

    /// The checked-in library is what the code generates. `UPDATE_OFFICIAL_LIBRARY=1` rewrites
    /// it (after a layout change, or to take a new md++ `VARIABLES` section).
    #[test]
    fn official_library_is_generated_from_the_layout() {
        let expected = official_library_text();
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/data/ene_ana.md++.lib");
        if std::env::var_os("UPDATE_OFFICIAL_LIBRARY").is_some() {
            std::fs::write(path, &expected).unwrap();
        }
        assert_eq!(
            OFFICIAL_LIBRARY, expected,
            "run with UPDATE_OFFICIAL_LIBRARY=1 to regenerate {path}"
        );

        let lib = EnergyTraj::from_library(OFFICIAL_LIBRARY).unwrap();
        assert_eq!(lib.version(), Some("2023-04-15"));
        for t in ["ENERTRJ", "FRENERTRJ"] {
            assert_eq!(lib.schema(t).unwrap(), &md_schema(t));
        }
        assert!(lib.variable_names().contains(&"totene"));
        assert!(lib.variable_names().contains(&"densit"));
    }

    #[test]
    fn an_edited_schema_section_is_refused() {
        let edited = OFFICIAL_LIBRARY.replace("subblock ENER 52 1", "subblock ENER 50 1");
        let err = EnergyTraj::from_library(&edited).unwrap_err();
        assert!(
            err.starts_with("library schema section edited after generation"),
            "{err}"
        );

        let newer =
            OFFICIAL_LIBRARY.replace("# energy-schema 1 ENERTRJ", "# energy-schema 2 ENERTRJ");
        let err = EnergyTraj::from_library(&newer).unwrap_err();
        assert!(
            err.starts_with("energy-schema 2 is newer than this reader"),
            "{err}"
        );

        // without the self-description an edit is a different, valid, library
        let unsigned: String = edited
            .lines()
            .filter(|l| !l.starts_with("# energy-schema"))
            .map(|l| format!("{l}\n"))
            .collect();
        assert!(EnergyTraj::from_library(&unsigned).is_ok());
    }
}
