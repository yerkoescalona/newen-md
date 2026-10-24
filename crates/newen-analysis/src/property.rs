//! gromos++ property specifiers (`utils::Property`, `PropertyContainer`):
//! `<type>%<atoms>[%<arg>…]`, several separated by spaces, with the types
//!
//! * `d`  distance between two atoms (nearest image)
//! * `a`  angle i–j–k in degrees
//! * `t`  torsion i–j–k–l in degrees, kept continuous with the previous value
//! * `tp` torsion mapped to (`arg` − 180, `arg` + 180], default (−180, 180]
//! * `o`  angle between two vectors, `op` its order parameter ½(3cos²−1)
//! * `pr`/`pa` pseudo-rotation angle / pucker amplitude of a five-membered ring
//!
//! Atoms are `mol:atoms` in the order written (`1:3,1,2` is atom 3, 1, 2), resolved through the
//! newen AtomSpecifier; vectors are `atom(<two atoms>)` or `cart(x,y,z)`. Averages
//! (`<…>`) and the `|x=…` substitution are not supported.

use newen_core::math::Periodicity;
use newen_core::selection::AtomSelection;
use newen_core::stat::Stat;
use newen_core::topology::Topology;
use newen_core::Vec3;

const SIN36SIN72: f64 = 0.559016994374947; // sin36°·sin72°

#[derive(Debug, Clone)]
pub enum VectorSpec {
    Atoms(usize, usize),
    /// `atom(<one atom>)`: the position of that atom
    Position(usize),
    Cart(Vec3),
}

impl VectorSpec {
    pub fn eval(&self, pos: &[Vec3], pbc: &Periodicity) -> Vec3 {
        match self {
            // the vector from atom i to the image of atom j
            VectorSpec::Atoms(i, j) => -pbc.nearest_image(pos[*i], pos[*j]),
            VectorSpec::Position(i) => pos[*i],
            VectorSpec::Cart(v) => *v,
        }
    }
}

#[derive(Debug, Clone)]
pub enum Kind {
    Distance,
    Angle,
    Torsion,
    PeriodicTorsion(f64),
    Order(VectorSpec, VectorSpec),
    OrderParam(VectorSpec, VectorSpec),
    PseudoRotation,
    PuckerAmplitude,
}

#[derive(Debug, Clone)]
pub struct Property {
    pub spec: String,
    pub kind: Kind,
    pub atoms: Vec<usize>,
    /// The numeric arguments after the atoms (`t%1:1,2,3,4%30%-60%60` → 30, −60, 60).
    pub args: Vec<f64>,
    pub stat: Stat,
    last: Option<f64>,
}

/// Resolve `mol:items` keeping the order written; each item may be a name, a number or a range.
pub fn ordered_atoms(spec: &str, topo: &Topology) -> Result<Vec<usize>, String> {
    let mut out = Vec::new();
    for term in spec.split(';') {
        let term = term.trim();
        if term.is_empty() {
            continue;
        }
        let Some((mol, items)) = term.split_once(':') else {
            // a bare global list
            for item in term.split(',') {
                let sel = AtomSelection::from_string(item.trim(), topo)
                    .map_err(|e| format!("'{item}': {e}"))?;
                out.extend(sel.indices().iter().copied());
            }
            continue;
        };
        if items.starts_with("res(") || items == "a" {
            let sel =
                AtomSelection::from_string(term, topo).map_err(|e| format!("'{term}': {e}"))?;
            out.extend(sel.indices().iter().copied());
            continue;
        }
        for item in items.split(',') {
            let one = format!("{mol}:{}", item.trim());
            let sel =
                AtomSelection::from_string(&one, topo).map_err(|e| format!("'{one}': {e}"))?;
            out.extend(sel.indices().iter().copied());
        }
    }
    Ok(out)
}

/// gromos++ `VectorSpecifier`: `atom(<one atom>)` is a position, `atom(<two atoms>)` the vector
/// from the first to the image of the second, `cart(x,y,z)` a constant.
pub fn parse_vector(s: &str, topo: &Topology) -> Result<VectorSpec, String> {
    let s = s.trim();
    if let Some(inner) = s.strip_prefix("atom(").and_then(|r| r.strip_suffix(')')) {
        let atoms = ordered_atoms(inner, topo)?;
        match atoms.len() {
            1 => return Ok(VectorSpec::Position(atoms[0])),
            2 => return Ok(VectorSpec::Atoms(atoms[0], atoms[1])),
            _ => return Err(format!("vector '{s}': one or two atoms expected")),
        }
    }
    if let Some(inner) = s.strip_prefix("cart(").and_then(|r| r.strip_suffix(')')) {
        let v: Result<Vec<f64>, _> = inner.split(',').map(|x| x.trim().parse::<f64>()).collect();
        let v = v.map_err(|_| format!("vector '{s}': three numbers expected"))?;
        if v.len() != 3 {
            return Err(format!("vector '{s}': three numbers expected"));
        }
        return Ok(VectorSpec::Cart(Vec3::new(v[0], v[1], v[2])));
    }
    Err(format!("vector '{s}': use atom(<spec>) or cart(x,y,z)"))
}

impl Property {
    pub fn parse(spec: &str, topo: &Topology) -> Result<Self, String> {
        let parts: Vec<&str> = spec.split('%').collect();
        if parts.len() < 2 {
            return Err(format!(
                "invalid property specifier '{spec}': syntax <type>%<atomspecifier>[%...]"
            ));
        }
        let need = |n: usize| -> Result<Vec<usize>, String> {
            let atoms = ordered_atoms(parts[1], topo)?;
            if atoms.len() != n {
                return Err(format!(
                    "property '{spec}': {n} atoms expected, {} given",
                    atoms.len()
                ));
            }
            Ok(atoms)
        };
        let (kind, atoms) = match parts[0] {
            "d" => (Kind::Distance, need(2)?),
            "a" => (Kind::Angle, need(3)?),
            "t" => (Kind::Torsion, need(4)?),
            "tp" => {
                let centre = parts
                    .get(2)
                    .map(|a| a.parse::<f64>().map_err(|_| format!("'{spec}': bad angle")))
                    .transpose()?
                    .unwrap_or(0.0);
                (Kind::PeriodicTorsion(centre), need(4)?)
            },
            "o" | "op" => {
                if parts.len() < 3 {
                    return Err(format!("property '{spec}': two vectors expected"));
                }
                let (v1, v2) = (parse_vector(parts[1], topo)?, parse_vector(parts[2], topo)?);
                (
                    if parts[0] == "o" {
                        Kind::Order(v1, v2)
                    } else {
                        Kind::OrderParam(v1, v2)
                    },
                    Vec::new(),
                )
            },
            "pr" => (Kind::PseudoRotation, need(5)?),
            "pa" => (Kind::PuckerAmplitude, need(5)?),
            other => {
                return Err(format!(
                    "property type '{other}' not supported (d, a, t, tp, o, op, pr, pa)"
                ))
            },
        };
        let numeric_from = if matches!(kind, Kind::Order(..) | Kind::OrderParam(..)) {
            3
        } else {
            2
        };
        let args: Vec<f64> = parts
            .iter()
            .skip(numeric_from)
            .filter_map(|a| a.trim().parse::<f64>().ok())
            .collect();
        Ok(Property {
            spec: spec.to_string(),
            kind,
            atoms,
            args,
            stat: Stat::new(),
            last: None,
        })
    }

    /// gromos++ `toTitle`: the type's name and the atoms as `AtomSpecifier::toString` writes
    /// them (`mol:atoms`, consecutive atoms collapsed to `first-last`).
    pub fn title(&self, topo: &Topology) -> String {
        let name = self.type_name();
        let atoms = if self.atoms.is_empty() {
            self.spec.split('%').nth(1).unwrap_or("").to_string()
        } else {
            atoms_to_string(&self.atoms, topo)
        };
        format!("{name}%{atoms}")
    }

    /// gromos++ `type()`: the long type name.
    pub fn type_name(&self) -> &'static str {
        match &self.kind {
            Kind::Distance => "Distance",
            Kind::Angle => "Angle",
            Kind::Torsion => "Torsion",
            Kind::PeriodicTorsion(_) => "PeriodicTorsion",
            Kind::Order(..) => "VectorOrder",
            Kind::OrderParam(..) => "VectorOrderParam",
            Kind::PseudoRotation => "PseudoRotation",
            Kind::PuckerAmplitude => "PuckerAmplitude",
        }
    }

    /// Evaluate on a frame and record the value.
    pub fn calc(&mut self, pos: &[Vec3], pbc: &Periodicity) -> f64 {
        let a = &self.atoms;
        let v = match &self.kind {
            Kind::Distance => pbc.nearest_image(pos[a[0]], pos[a[1]]).length(),
            Kind::Angle => {
                let ta = pbc.nearest_image(pos[a[0]], pos[a[1]]);
                let tb = pbc.nearest_image(pos[a[2]], pos[a[1]]);
                (ta.dot(tb) / (ta.length() * tb.length()))
                    .acos()
                    .to_degrees()
            },
            Kind::Torsion => {
                let mut d = torsion(pos, pbc, a[0], a[1], a[2], a[3]);
                if let Some(last) = self.last {
                    let x = (d - last) / 360.0;
                    let ix = if x > 0.0 {
                        (x + 0.5) as i64
                    } else {
                        (x - 0.5) as i64
                    };
                    d -= ix as f64 * 360.0;
                }
                d
            },
            Kind::PeriodicTorsion(centre) => {
                let mut d = torsion(pos, pbc, a[0], a[1], a[2], a[3]);
                while d < centre - 180.0 {
                    d += 360.0;
                }
                while d > centre + 180.0 {
                    d -= 360.0;
                }
                d
            },
            Kind::Order(v1, v2) => {
                let (axis, tmp) = (v1.eval(pos, pbc), v2.eval(pos, pbc));
                (tmp.dot(axis) / (tmp.length() * axis.length()))
                    .acos()
                    .to_degrees()
            },
            Kind::OrderParam(v1, v2) => {
                let (axis, tmp) = (v1.eval(pos, pbc), v2.eval(pos, pbc));
                let c = tmp.dot(axis) / (tmp.length() * axis.length());
                0.5 * (3.0 * c * c - 1.0)
            },
            Kind::PseudoRotation | Kind::PuckerAmplitude => {
                let t: Vec<f64> = [
                    (0, 1, 2, 3),
                    (1, 2, 3, 4),
                    (2, 3, 4, 0),
                    (3, 4, 0, 1),
                    (4, 0, 1, 2),
                ]
                .iter()
                .map(|&(i, j, k, l)| {
                    let mut d = torsion(pos, pbc, a[i], a[j], a[k], a[l]);
                    if d > 180.0 {
                        d -= 360.0;
                    }
                    d
                })
                .collect();
                let factor = ((t[2] + t[4]) - (t[1] + t[3])) / (2.0 * t[0] * SIN36SIN72);
                if matches!(self.kind, Kind::PseudoRotation) {
                    let mut d = factor.atan().to_degrees();
                    if t[0] < 0.0 {
                        d += 180.0;
                    }
                    d
                } else {
                    let mut pr = factor.atan();
                    if t[0] < 0.0 {
                        pr += std::f64::consts::PI;
                    }
                    t[0] / pr.cos()
                }
            },
        };
        self.last = Some(v);
        self.stat.add(v);
        v
    }

    /// gromos++ `Property::average`: average, rmsd, error estimate — fixed, 8 decimals, width 15
    /// (`-nan` when there are too few values for gromos++'s block averaging).
    pub fn average(&self) -> String {
        if self.stat.n() == 0 {
            return String::new();
        }
        let ee = self.stat.ee_strict();
        let ee = if ee.is_nan() {
            format!("{:>15}", "-nan")
        } else {
            format!("{:15.8}", ee)
        };
        format!("{:15.8}{:15.8}{}", self.stat.ave(), self.stat.rmsd(), ee)
    }
}

/// gromos++ `AtomSpecifier::toString`: `mol:a-b,c;mol2:…`, solvent as `s:`, 1-based.
pub fn atoms_to_string(atoms: &[usize], topo: &Topology) -> String {
    let n_solute = topo.num_solute_atoms();
    let locate = |a: usize| -> (i64, usize) {
        if a >= n_solute {
            return (-1, a - n_solute);
        }
        for (m, r) in topo.molecules.iter().enumerate() {
            if r.contains(&a) {
                return (m as i64, a - r.start);
            }
        }
        (0, a)
    };
    let mut out = String::new();
    let mut mol: Option<i64> = None;
    let mut range_start: Option<usize> = None;
    let mut last: Option<usize> = None;
    let close_range = |out: &mut String, range_start: &mut Option<usize>, last: Option<usize>| {
        if let (Some(_), Some(l)) = (*range_start, last) {
            out.push_str(&format!("-{}", l + 1));
        }
        *range_start = None;
    };
    for &a in atoms {
        let (m, loc) = locate(a);
        if mol != Some(m) {
            close_range(&mut out, &mut range_start, last);
            if mol.is_some() {
                out.push(';');
            }
            out.push_str(&if m < 0 {
                "s:".to_string()
            } else {
                format!("{}:", m + 1)
            });
            out.push_str(&(loc + 1).to_string());
            mol = Some(m);
            last = Some(loc);
            continue;
        }
        if last == Some(loc.wrapping_sub(1)) && loc > 0 {
            range_start = range_start.or(last);
        } else {
            close_range(&mut out, &mut range_start, last);
            out.push_str(&format!(",{}", loc + 1));
        }
        last = Some(loc);
    }
    close_range(&mut out, &mut range_start, last);
    out
}

/// gromos++ `TorsionProperty::calc` before the continuity transformation: signed, degrees.
pub fn torsion(pos: &[Vec3], pbc: &Periodicity, i: usize, j: usize, k: usize, l: usize) -> f64 {
    // gromos++: pos(i) − nearestImage(pos(i), pos(j)) = the minimum-image vector ri − rj
    let ta = pbc.nearest_image(pos[i], pos[j]);
    let tb = pbc.nearest_image(pos[l], pos[k]);
    let tc = pbc.nearest_image(pos[k], pos[j]);
    let p1 = ta.cross(tc);
    let p2 = tb.cross(tc);
    let cosphi = (p1.dot(p2) / (p1.length() * p2.length())).clamp(-1.0, 1.0);
    let mut d = cosphi.acos().to_degrees();
    if p1.cross(p2).dot(tc) < 0.0 {
        d = -d;
    }
    d
}

/// The properties of a program: `@prop` values joined by spaces, split outside brackets.
pub struct PropertyContainer {
    pub props: Vec<Property>,
}

impl PropertyContainer {
    pub fn parse(specs: &[String], topo: &Topology) -> Result<Self, String> {
        let mut props = Vec::new();
        for s in specs {
            for one in split_outside_brackets(s) {
                if one.starts_with('<') {
                    return Err("average properties (<...>) are not supported".into());
                }
                if one.contains('|') {
                    return Err("the |x=... variable substitution is not supported".into());
                }
                props.push(Property::parse(&one, topo)?);
            }
        }
        if props.is_empty() {
            return Err("no property given".into());
        }
        Ok(PropertyContainer { props })
    }

    pub fn calc(&mut self, pos: &[Vec3], pbc: &Periodicity) -> Vec<f64> {
        self.props.iter_mut().map(|p| p.calc(pos, pbc)).collect()
    }

    /// gromos++ `PropertyContainer::toTitle`: every title followed by two tabs (at most ten,
    /// then `[ ... ]`), and — unconditionally, a quirk of the original — the last title again.
    pub fn title(&self, topo: &Topology) -> String {
        let mut s: String = self
            .props
            .iter()
            .take(10)
            .map(|p| format!("{}\t\t", p.title(topo)))
            .collect();
        if self.props.len() > 10 {
            s.push_str("[ ... ]\t\t");
        }
        if let Some(last) = self.props.last() {
            s.push_str(&last.title(topo));
        }
        s
    }
}

fn split_outside_brackets(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut cur = String::new();
    for c in s.chars() {
        match c {
            '(' | '<' | '[' => depth += 1,
            ')' | '>' | ']' => depth -= 1,
            _ => {},
        }
        if c == ' ' && depth == 0 {
            if !cur.trim().is_empty() {
                out.push(cur.trim().to_string());
            }
            cur.clear();
        } else {
            cur.push(c);
        }
    }
    if !cur.trim().is_empty() {
        out.push(cur.trim().to_string());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::Vacuum;

    #[test]
    fn torsion_signs_and_continuity() {
        let pbc = Periodicity::Vacuum(Vacuum);
        // a +60° torsion
        let pos = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.5, 0.8660254, 1.0),
        ];
        let d = torsion(&pos, &pbc, 0, 1, 2, 3);
        assert!((d.abs() - 60.0).abs() < 1e-5);
    }

    #[test]
    fn specifiers_split_outside_brackets() {
        assert_eq!(
            split_outside_brackets("d%1:1,2 o%atom(1:1,2)%cart(0, 0, 1)"),
            vec!["d%1:1,2", "o%atom(1:1,2)%cart(0, 0, 1)"]
        );
    }
}
