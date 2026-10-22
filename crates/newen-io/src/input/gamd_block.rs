//! The gromosXX GAMD input: the `GAMD` block of the `.imd` file (gromosXX `read_GAMD`,
//! `in_parameter.cc:4148-4320`), the `GAMDATOMS`/`GAMDGROUPS` blocks of the `@gamd` file
//! (`In_GAMD::read`, `in_gamd.cc:76-238`) and the `GAMDSTAT` block of a configuration
//! (`_read_gamdstat`, `in_configuration.cc:3754-3785`; `_print_gamdstat`,
//! `out_configuration.cc:4013-4034`).
//!
//! The `GAMD` block is one token stream, read in gromosXX's order:
//!
//! ```text
//! GAMD
//! # GAMD
//!   1
//! # SEARCH  FORM  THRESH  NTIGAMDS
//!     0       1      1        0
//! # AGROUPS
//!   2
//! # IGROUPS
//!   2
//! # DIHSTD  TOTSTD
//!   24.79   24.79
//! # ED   (IGROUPS values, likewise ET, KD, KT)
//!   10.0   100.0
//! # ET
//!   -10.0   20.0
//! # KD
//!   10.0   100.0
//! # KT
//!   -10.0   20.0
//! # EQSTEPS
//!  5000
//! # GWINDOW
//!  0
//! END
//! ```
//!
//! SEARCH 0 = fixed parameters, 1 = cMD search, 2 = GaMD search; FORM 1 = dual boost,
//! 2 = potential-energy boost, 3 = dihedral boost; THRESH 1 = lower bound, 2 = upper bound.
use newen_integrators::gamd::{BoostForm, GamdSettings, SearchMode, ThresholdType};
use std::collections::BTreeMap;

/// The `GAMD` block as gromosXX reads it. The per-group vectors hold IGROUPS values (the input
/// order); [`Self::settings`] prepends gromosXX's non-accelerated group 0.
#[derive(Debug, Clone, PartialEq)]
pub struct GamdBlock {
    /// GAMD: 1 = Gaussian accelerated MD on.
    pub enabled: bool,
    /// SEARCH (0, 1, 2).
    pub search: i32,
    /// FORM (1, 2, 3).
    pub form: i32,
    /// THRESH: gromosXX accepts 1-4; 3 and 4 select nothing and leave the lower bound.
    pub thresh: i32,
    /// NTIGAMDS: 1 = start the search statistics from scratch, 0 = from `GAMDSTAT`.
    pub ntigamds: bool,
    /// AGROUPS: atom groups of the `@gamd` file.
    pub agroups: usize,
    /// IGROUPS: accelerated interaction groups.
    pub igroups: usize,
    pub dihstd: f64,
    pub totstd: f64,
    pub ed: Vec<f64>,
    pub et: Vec<f64>,
    pub kd: Vec<f64>,
    pub kt: Vec<f64>,
    /// EQSTEPS.
    pub eqsteps: i64,
    /// GWINDOW.
    pub gwindow: i64,
}

/// gromosXX's `Block` token reader for one block: its error messages
/// (`blockinput.h:354-416`, `blockinput.cc:270-288`).
struct BlockTokens<'a> {
    block: &'a str,
    values: Vec<String>,
    pos: usize,
}

impl<'a> BlockTokens<'a> {
    fn new(block: &'a str, lines: &[String]) -> Self {
        let values = lines
            .iter()
            .map(|l| l.split('#').next().unwrap_or_default())
            .flat_map(|l| l.split_whitespace().map(str::to_string))
            .collect();
        Self {
            block,
            values,
            pos: 0,
        }
    }

    fn next(&mut self, name: &str) -> Result<String, String> {
        let v = self.values.get(self.pos).cloned().ok_or_else(|| {
            format!(
                "{} block reached END before {name} could be read!",
                self.block
            )
        })?;
        self.pos += 1;
        Ok(v)
    }

    fn int(&mut self, name: &str) -> Result<i64, String> {
        let v = self.next(name)?;
        v.parse()
            .map_err(|_| format!("{} block: wrong value type for {name}: {v}", self.block))
    }

    fn real(&mut self, name: &str) -> Result<f64, String> {
        let v = self.next(name)?;
        v.parse()
            .map_err(|_| format!("{} block: wrong value type for {name}: {v}", self.block))
    }

    fn out_of_range(&self, name: &str, v: impl std::fmt::Display, should: &str) -> String {
        format!(
            "{} block: {name} parameter is {v}, should be: {should}",
            self.block
        )
    }

    /// `allowed` ("0,1,2"): the value must be one of them.
    fn one_of(&mut self, name: &str, allowed: &[i64]) -> Result<i64, String> {
        let v = self.int(name)?;
        if allowed.contains(&v) {
            Ok(v)
        } else {
            let list: Vec<String> = allowed.iter().map(i64::to_string).collect();
            Err(self.out_of_range(name, v, &format!(" one of: {}", list.join(","))))
        }
    }

    fn leftover(&mut self) -> Result<(), String> {
        match self.values.get(self.pos) {
            Some(extra) => Err(format!(
                "Left-over parameters in {} block: {extra}\n",
                self.block
            )),
            None => Ok(()),
        }
    }
}

impl GamdBlock {
    /// Parse the block's data lines (between `GAMD` and `END`, comments allowed), every
    /// parameter in gromosXX's order with its range, whatever GAMD says.
    pub fn from_lines(lines: &[String]) -> Result<Self, String> {
        let mut t = BlockTokens::new("GAMD", lines);
        let enabled = t.one_of("GAMD", &[0, 1])? == 1;
        let search = t.one_of("SEARCH", &[0, 1, 2])? as i32;
        let form = t.one_of("FORM", &[1, 2, 3])? as i32;
        let thresh = t.one_of("THRESH", &[1, 2, 3, 4])? as i32;
        let ntigamds = t.one_of("NTIGAMD", &[0, 1])? == 1;
        let agroups = t.int("AGROUPS")?;
        if agroups <= 0 {
            return Err(t.out_of_range("AGROUPS", agroups, " >0"));
        }
        // gromosXX reads IGROUPS under the name AGROUPS (`in_parameter.cc:4223`)
        let igroups = t.int("AGROUPS")?;
        if igroups <= 0 {
            return Err(t.out_of_range("AGROUPS", igroups, " >0"));
        }
        let mut positive = |name: &str| -> Result<f64, String> {
            let v = t.real(name)?;
            if v > 0.0 {
                Ok(v)
            } else {
                Err(t.out_of_range(name, v, " >0"))
            }
        };
        let dihstd = positive("DIHSTD")?;
        let totstd = positive("TOTSTD")?;
        let n = igroups as usize;
        let mut per_group = |name: &str| -> Result<Vec<f64>, String> {
            (0..n).map(|i| t.real(&format!("{name}[{i}]"))).collect()
        };
        let ed = per_group("ED")?;
        let et = per_group("ET")?;
        let kd = per_group("KD")?;
        let kt = per_group("KT")?;
        let mut non_negative = |name: &str| -> Result<i64, String> {
            let v = t.int(name)?;
            if v >= 0 {
                Ok(v)
            } else {
                Err(t.out_of_range(name, v, " >=0"))
            }
        };
        let eqsteps = non_negative("EQSTEPS")?;
        let gwindow = non_negative("GWINDOW")?;
        t.leftover()?;
        Ok(Self {
            enabled,
            search,
            form,
            thresh,
            ntigamds,
            agroups: agroups as usize,
            igroups: n,
            dihstd,
            totstd,
            ed,
            et,
            kd,
            kt,
            eqsteps,
            gwindow,
        })
    }

    /// gromosXX's interaction-group count, the non-accelerated group 0 included
    /// (`param.gamd.igroups = IGROUPS + 1`): the rows of the `# gamd` energy section.
    pub fn n_interaction_groups(&self) -> usize {
        self.igroups + 1
    }

    pub fn search_mode(&self) -> SearchMode {
        match self.search {
            1 => SearchMode::CmdSearch,
            2 => SearchMode::GamdSearch,
            _ => SearchMode::NoSearch,
        }
    }

    pub fn boost_form(&self) -> BoostForm {
        match self.form {
            2 => BoostForm::TotalBoost,
            3 => BoostForm::DihedralBoost,
            _ => BoostForm::DualBoost,
        }
    }

    /// THRESH 2 is the upper bound; 1, 3 and 4 the lower bound (gromosXX's default, which its
    /// switch leaves in place for 3 and 4).
    pub fn threshold(&self) -> ThresholdType {
        if self.thresh == 2 {
            ThresholdType::UpperBound
        } else {
            ThresholdType::LowerBound
        }
    }

    /// The algorithm's settings, with group 0 (all zero) in front of the input's groups.
    pub fn settings(&self) -> GamdSettings {
        let with_zero = |v: &[f64]| std::iter::once(0.0).chain(v.iter().copied()).collect();
        GamdSettings {
            search: self.search_mode(),
            form: self.boost_form(),
            thresh: self.threshold(),
            ntisearch: self.ntigamds,
            dihstd: self.dihstd,
            totstd: self.totstd,
            equilibration: self.eqsteps,
            window: self.gwindow,
            ed: with_zero(&self.ed),
            et: with_zero(&self.et),
            kd: with_zero(&self.kd),
            kt: with_zero(&self.kt),
        }
    }

    /// The block in gromosXX's layout (its example block), parsed back unchanged.
    pub fn render(&self) -> String {
        let list = |v: &[f64]| {
            v.iter()
                .map(|x| format!("{x:?}"))
                .collect::<Vec<_>>()
                .join("   ")
        };
        format!(
            "GAMD\n# GAMD\n  {}\n# SEARCH  FORM  THRESH  NTIGAMDS\n    {}       {}      {}        {}\n\
             # AGROUPS  IGROUPS\n  {}    {}\n# DIHSTD  TOTSTD\n  {:?}   {:?}\n# ED\n  {}\n# ET\n  {}\n\
             # KD\n  {}\n# KT\n  {}\n# EQSTEPS\n  {}\n# GWINDOW\n  {}\nEND\n",
            u8::from(self.enabled),
            self.search,
            self.form,
            self.thresh,
            u8::from(self.ntigamds),
            self.agroups,
            self.igroups,
            self.dihstd,
            self.totstd,
            list(&self.ed),
            list(&self.et),
            list(&self.kd),
            list(&self.kt),
            self.eqsteps,
            self.gwindow,
        )
    }
}

impl std::fmt::Display for GamdBlock {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.render())
    }
}

/// The `@gamd` file: which atom group each atom is in (`GAMDATOMS`, unlisted atoms in group 0)
/// and which interaction group each pair of atom groups is booked to (`GAMDGROUPS`, unlisted
/// pairs in group 0).
#[derive(Debug, Clone, PartialEq)]
pub struct GamdGroups {
    /// Atom group per atom.
    pub atom_group: Vec<u32>,
    /// Interaction group per (atom group, atom group), both orders.
    pub pairs: BTreeMap<(u32, u32), u32>,
}

/// The data lines (between the name and `END`, comments dropped) of every block of a GROMOS
/// file.
fn blocks(text: &str) -> BTreeMap<String, Vec<String>> {
    let mut out = BTreeMap::new();
    let mut current: Option<(String, Vec<String>)> = None;
    for raw in text.lines() {
        let line = raw.split('#').next().unwrap_or_default().trim();
        if line.is_empty() {
            continue;
        }
        match current.as_mut() {
            None => current = Some((line.to_string(), Vec::new())),
            Some((name, data)) => {
                if line == "END" {
                    out.entry(std::mem::take(name))
                        .or_insert_with(|| std::mem::take(data));
                    current = None;
                } else {
                    data.push(line.to_string());
                }
            },
        }
    }
    out
}

impl GamdGroups {
    /// Read `GAMDATOMS` and `GAMDGROUPS` for a topology of `num_atoms` atoms against the
    /// block's AGROUPS and IGROUPS, with gromosXX's errors (`in_gamd.cc:83-235`).
    pub fn parse(
        text: &str,
        num_atoms: usize,
        agroups: usize,
        igroups: usize,
    ) -> Result<Self, String> {
        let blocks = blocks(text);
        let lines = |name: &str| {
            blocks
                .get(name)
                .filter(|l| !l.is_empty())
                .ok_or_else(|| format!("no {name} block in GAMD definition file"))
        };
        let inconsistent = || "Inconsistent number of groups between imd file and gamd file";

        let atoms = lines("GAMDATOMS")?;
        let numgroups: usize = atoms[0]
            .split_whitespace()
            .next()
            .and_then(|v| v.parse().ok())
            .ok_or("bad line in GAMDATOMS block; expected atom GROUPS\n")?;
        if numgroups != agroups {
            return Err(inconsistent().to_string());
        }
        let mut atom_group = vec![0u32; num_atoms];
        let mut assigned = vec![false; num_atoms];
        for line in &atoms[1..] {
            let v: Vec<u64> = line
                .split_whitespace()
                .take(3)
                .map_while(|x| x.parse().ok())
                .collect();
            let [iatom, fatom, agroup] = v[..] else {
                return Err(
                    "bad line in GAMDATOMS block; expected format:\nGAMDATOMS\n#  INATOM   \
                            FINATOM   AGROUP\nEND\n"
                        .to_string(),
                );
            };
            // gromosXX checks `iatom > num_atoms` after the 0-based shift (`in_gamd.cc:138`),
            // which lets FINATOM = num_atoms + 1 write past the end; newen refuses it
            if iatom == 0 || fatom == 0 || iatom > num_atoms as u64 || fatom > num_atoms as u64 {
                return Err(format!(
                    "GAMD atoms ({iatom}-{fatom}) atom indices out of range."
                ));
            }
            if agroup > numgroups as u64 {
                return Err("GAMD atoms acceleration group out of range.".to_string());
            }
            for i in (iatom - 1)..fatom {
                let i = i as usize;
                if assigned[i] {
                    return Err("Duplicated atoms in the GAMDATOMS.".to_string());
                }
                assigned[i] = true;
                atom_group[i] = agroup as u32;
            }
        }

        let groups = lines("GAMDGROUPS")?;
        let anumgroups: usize = groups[0]
            .split_whitespace()
            .next()
            .and_then(|v| v.parse().ok())
            .ok_or("bad line in GAMDATOMS block; expected acceleration GROUPS\n")?;
        if anumgroups != igroups {
            return Err(inconsistent().to_string());
        }
        let mut pairs = BTreeMap::new();
        for line in &groups[1..] {
            let v: Vec<u32> = line
                .split_whitespace()
                .take(3)
                .map_while(|x| x.parse().ok())
                .collect();
            let [g1, g2, ig] = v[..] else {
                return Err(
                    "bad line in GAMDGROUPS block; expected format:\nGAMDGROUPS\n#  \
                            GROUP_1   GROUP_2   IGROUP\nEND\n"
                        .to_string(),
                );
            };
            if ig as usize > anumgroups {
                return Err("Acceleration groups in GAMDGROUPS are out of range".to_string());
            }
            // std::map::insert keeps the first entry of a repeated pair
            pairs.entry((g1, g2)).or_insert(ig);
            pairs.entry((g2, g1)).or_insert(ig);
        }
        Ok(Self { atom_group, pairs })
    }

    /// The interaction group of a pair of atom groups (gromosXX `gamd_interaction_group`).
    pub fn interaction_group(&self, a: u32, b: u32) -> u32 {
        self.pairs.get(&(a, b)).copied().unwrap_or(0)
    }

    /// The one interaction group every interaction of the run is booked to — every pair of the
    /// atom groups that occur maps to it — or `None` when they split into several.
    pub fn single_interaction_group(&self) -> Option<u32> {
        let present: std::collections::BTreeSet<u32> = self.atom_group.iter().copied().collect();
        let mut groups = present
            .iter()
            .flat_map(|&a| present.iter().map(move |&b| self.interaction_group(a, b)));
        let first = groups.next()?;
        groups.all(|g| g == first).then_some(first)
    }

    /// The one interaction group of the reaction-field excluded pairs of the first `num_solute`
    /// atoms: gromosXX looks it up with the key (0, group of j), its `gamdi` being shadowed
    /// (`nonbonded_innerloop.cc:1694, 1711-1723`); `None` when the solute's groups give several.
    pub fn excluded_pair_group(&self, num_solute: usize) -> Option<u32> {
        let present: std::collections::BTreeSet<u32> =
            self.atom_group.iter().take(num_solute).copied().collect();
        let mut groups = present.iter().map(|&b| self.interaction_group(0, b));
        let Some(first) = groups.next() else {
            return Some(0);
        };
        groups.all(|g| g == first).then_some(first)
    }
}

/// The `GAMDSTAT` block of a configuration: `stepsdone`, then one row per interaction group
/// (group 0 included) of M2D M2T VmaxD VmaxT VmeanD VmeanT VminD VminT sigmaVD sigmaVT.
#[derive(Debug, Clone, PartialEq)]
pub struct GamdStat {
    pub stepsdone: i64,
    pub rows: Vec<[f64; 10]>,
}

impl GamdStat {
    /// The `GAMDSTAT` block of configuration text, for `n_groups` interaction groups (group 0
    /// included); `Ok(None)` without one.
    pub fn from_configuration(text: &str, n_groups: usize) -> Result<Option<Self>, String> {
        let Some(lines) = blocks(text).remove("GAMDSTAT") else {
            return Ok(None);
        };
        let bad = || "bad line in GAMDSTAT block".to_string();
        let stepsdone = lines
            .first()
            .and_then(|l| l.split_whitespace().next())
            .and_then(|v| v.parse().ok())
            .ok_or_else(bad)?;
        let mut rows = Vec::with_capacity(n_groups);
        for g in 0..n_groups {
            let v: Vec<f64> = lines
                .get(g + 1)
                .ok_or_else(bad)?
                .split_whitespace()
                .take(10)
                .map(|x| x.parse().map_err(|_| bad()))
                .collect::<Result<_, _>>()?;
            rows.push(<[f64; 10]>::try_from(v).map_err(|_| bad())?);
        }
        Ok(Some(Self { stepsdone, rows }))
    }

    /// The block as gromosXX prints it (`%15.9f`, `out_configuration.cc:4017-4033`).
    pub fn render(&self) -> String {
        let mut s = format!("GAMDSTAT\n{:>15}\n", self.stepsdone);
        for row in &self.rows {
            let cells: Vec<String> = row.iter().map(|v| format!("{v:>15.9}")).collect();
            s.push_str(&cells.join(" "));
            s.push('\n');
        }
        s.push_str("END\n");
        s
    }

    /// Append the block to the configuration file at `path` (it is the last block gromosXX
    /// writes, `out_configuration.cc:619-621`).
    pub fn append_to(&self, path: &str) -> Result<(), String> {
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new()
            .append(true)
            .open(path)
            .map_err(|e| format!("could not open {path}: {e}"))?;
        f.write_all(self.render().as_bytes())
            .map_err(|e| format!("could not write GAMDSTAT to {path}: {e}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn data_lines(block: &str) -> Vec<String> {
        let mut lines = block.lines().map(str::to_string).collect::<Vec<_>>();
        assert_eq!(lines.first().map(|l| l.trim()), Some("GAMD"));
        lines.remove(0);
        assert_eq!(
            lines.pop().map(|l| l.trim().to_string()),
            Some("END".into())
        );
        lines
    }

    /// gromosXX's example block (`in_parameter.cc:4155-4205`).
    const EXAMPLE: &str = "GAMD
# GAMD      0,1
#              0: no gaussian accelerated molecular dynamics [default]
#              1: aussian accelerated molecular dynamics
# SEARCH       0-2
# FORM         1-3
# THRESH       1-2: type of energy threshold to use
# NTIGAMDS   0,1
#
# GAMD
  1
# SEARCH  FORM  THRESH  NTIGAMDS
    0       1      1        0
# AGROUPS
  2
# IGROUPS
  2
# DIHSTD  TOTSTD
  24.79   24.79
# ED
  10.0   100.0
# ET
  -10.0   20.0
# KD
  10.0   100.0
# KT
  -10.0   20.0
# EQSTEPS
 5000
# GWINDOW
 0
END";

    #[test]
    fn parses_the_gromosxx_example_block() {
        let b = GamdBlock::from_lines(&data_lines(EXAMPLE)).unwrap();
        assert!(b.enabled);
        assert_eq!((b.search, b.form, b.thresh), (0, 1, 1));
        assert!(!b.ntigamds);
        assert_eq!((b.agroups, b.igroups), (2, 2));
        assert_eq!((b.dihstd, b.totstd), (24.79, 24.79));
        assert_eq!(b.ed, vec![10.0, 100.0]);
        assert_eq!(b.et, vec![-10.0, 20.0]);
        assert_eq!(b.kd, vec![10.0, 100.0]);
        assert_eq!(b.kt, vec![-10.0, 20.0]);
        assert_eq!((b.eqsteps, b.gwindow), (5000, 0));
        assert_eq!(b.boost_form(), BoostForm::DualBoost);
        assert_eq!(b.search_mode(), SearchMode::NoSearch);
        assert_eq!(b.threshold(), ThresholdType::LowerBound);
        assert_eq!(b.n_interaction_groups(), 3);
        let s = b.settings();
        assert_eq!(s.ed, vec![0.0, 10.0, 100.0]);
        assert_eq!(s.kt, vec![0.0, -10.0, 20.0]);
        assert_eq!(s.equilibration, 5000);
    }

    /// The reference input `aladip_gamd.in`: FORM 3 is the dihedral boost.
    #[test]
    fn parses_the_aladip_gamd_reference_input() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../newen-md/tests/gromosXX_references/aladip_gamd/aladip_gamd.in"
        );
        let imd = crate::imd::read_imd_file(path).unwrap();
        let b = GamdBlock::from_lines(&imd.raw_blocks["GAMD"]).unwrap();
        assert!(b.enabled);
        assert_eq!(b.search_mode(), SearchMode::NoSearch);
        assert_eq!(b.boost_form(), BoostForm::DihedralBoost);
        assert_eq!((b.agroups, b.igroups), (1, 1));
        assert_eq!(
            (b.ed[0], b.et[0], b.kd[0], b.kt[0]),
            (30.0, -100.0, 0.005, 0.0005)
        );
        assert_eq!((b.eqsteps, b.gwindow), (0, 0));
    }

    #[test]
    fn form_and_search_codes_follow_gromosxx() {
        let block = |line: &str| {
            let text = EXAMPLE.replace("    0       1      1        0", line);
            GamdBlock::from_lines(&data_lines(&text)).unwrap()
        };
        assert_eq!(block("0 2 1 0").boost_form(), BoostForm::TotalBoost);
        assert_eq!(block("0 3 1 0").boost_form(), BoostForm::DihedralBoost);
        assert_eq!(block("1 3 2 1").search_mode(), SearchMode::CmdSearch);
        assert_eq!(block("2 3 2 1").search_mode(), SearchMode::GamdSearch);
        assert_eq!(block("2 3 2 1").threshold(), ThresholdType::UpperBound);
        // accepted by gromosXX's range check, mapped by none of its cases
        assert_eq!(block("0 1 4 0").threshold(), ThresholdType::LowerBound);
    }

    #[test]
    fn errors_are_gromosxx_messages() {
        let err = |from: &str, to: &str| {
            let text = EXAMPLE.replace(from, to);
            GamdBlock::from_lines(&data_lines(&text)).unwrap_err()
        };
        assert_eq!(
            err(
                "    0       1      1        0",
                "    0       4      1        0"
            ),
            "GAMD block: FORM parameter is 4, should be:  one of: 1,2,3"
        );
        assert_eq!(
            err("# IGROUPS\n  2", "# IGROUPS\n  0"),
            "GAMD block: AGROUPS parameter is 0, should be:  >0"
        );
        assert_eq!(
            err("24.79   24.79", "24.79   0.0"),
            "GAMD block: TOTSTD parameter is 0, should be:  >0"
        );
        assert_eq!(
            err(" 5000\n", " x\n"),
            "GAMD block: wrong value type for EQSTEPS: x"
        );
        assert_eq!(
            err("# KT\n  -10.0   20.0", "# KT\n  -10.0   x"),
            "GAMD block: wrong value type for KT[1]: x"
        );
        // a missing value takes the next one: KT[1] = 5000, EQSTEPS = 0, GWINDOW missing
        assert_eq!(
            err("# KT\n  -10.0   20.0", "# KT\n  -10.0"),
            "GAMD block reached END before GWINDOW could be read!"
        );
        assert_eq!(
            err(" 0\nEND", " 0 7\nEND"),
            "Left-over parameters in GAMD block: 7\n"
        );
    }

    #[test]
    fn render_round_trips() {
        let b = GamdBlock::from_lines(&data_lines(EXAMPLE)).unwrap();
        let text = b.to_string();
        assert_eq!(
            GamdBlock::from_lines(&data_lines(text.trim_end())).unwrap(),
            b
        );
        // and through the IMD writer, which passes the block through verbatim
        let imd = crate::imd::parse_imd_str(&format!("SYSTEM\n  1 0\nEND\n{text}")).unwrap();
        let again = crate::imd::parse_imd_str(&crate::imd_write::write_imd(&imd, None)).unwrap();
        assert_eq!(GamdBlock::from_lines(&again.raw_blocks["GAMD"]).unwrap(), b);
    }

    const GROUPS: &str =
        "TITLE\ntwo groups\nEND\nGAMDATOMS\n2\n# INATOM FINATOM AGROUP\n  1 3 1\n  4 5 \
                          2\nEND\nGAMDGROUPS\n2\n  1 1 1\n  2 2 2\n  1 1 2\nEND\n";

    #[test]
    fn reads_the_gamd_groups_file() {
        let g = GamdGroups::parse(GROUPS, 6, 2, 2).unwrap();
        assert_eq!(g.atom_group, vec![1, 1, 1, 2, 2, 0]);
        // the first entry of a repeated pair wins
        assert_eq!(g.interaction_group(1, 1), 1);
        assert_eq!(g.interaction_group(2, 2), 2);
        assert_eq!(g.interaction_group(1, 2), 0);
        assert_eq!(g.single_interaction_group(), None);
        let one = "GAMDATOMS\n1\n 1 6 1\nEND\nGAMDGROUPS\n1\n 1 1 1\nEND\n";
        let g = GamdGroups::parse(one, 6, 1, 1).unwrap();
        assert_eq!(g.single_interaction_group(), Some(1));
        assert_eq!(g.excluded_pair_group(4), Some(0));
        let with_zero = "GAMDATOMS\n1\n 1 6 1\nEND\nGAMDGROUPS\n1\n 1 1 1\n 0 1 1\nEND\n";
        let g = GamdGroups::parse(with_zero, 6, 1, 1).unwrap();
        assert_eq!(g.excluded_pair_group(4), Some(1));
        assert_eq!(
            GamdGroups::parse(GROUPS, 6, 2, 3).unwrap_err(),
            "Inconsistent number of groups between imd file and gamd file"
        );
        assert_eq!(
            GamdGroups::parse(GROUPS, 4, 2, 2).unwrap_err(),
            "GAMD atoms (4-5) atom indices out of range."
        );
        assert_eq!(
            GamdGroups::parse("GAMDATOMS\n1\n1 2 1\nEND\n", 6, 1, 1).unwrap_err(),
            "no GAMDGROUPS block in GAMD definition file"
        );
    }

    #[test]
    fn gamdstat_round_trips_in_gromosxx_format() {
        let s = GamdStat {
            stepsdone: 100,
            rows: vec![
                [0.0; 10],
                [
                    1026.631365723,
                    0.0,
                    16.956357753,
                    0.0,
                    8.608420333,
                    0.0,
                    3.919509044,
                    0.0,
                    3.204108871,
                    0.0,
                ],
            ],
        };
        let text = s.render();
        assert!(text.starts_with("GAMDSTAT\n            100\n    0.000000000     0.000000000 "));
        assert!(text.contains("\n 1026.631365723     0.000000000    16.956357753 "));
        let conf = format!("TITLE\nx\nEND\nPOSITION\n# c\nEND\n{text}");
        assert_eq!(GamdStat::from_configuration(&conf, 2).unwrap(), Some(s));
        assert_eq!(
            GamdStat::from_configuration("TITLE\nx\nEND\n", 2).unwrap(),
            None
        );
    }
}
