//! `@pbc <r|t|c|v> [gather method]` and the gromos++ gather methods on a frame.
//!
//! The trajectory frames carry the box lengths; rectangular (`r`) and vacuum (`v`) are what the
//! analysis programs see. The default gather method is gromos++'s `cog`: molecule 1 is made whole
//! from its first atom (nearest image to the origin), every other solute molecule and every
//! solvent molecule is made whole starting from its image nearest to molecule 1's centre of
//! geometry. `nog` leaves the coordinates alone.

use newen_core::math::{Periodicity, Rectangular, Vacuum};
use newen_core::topology::Topology;
use newen_core::Vec3;

use crate::args::Arguments;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Gather {
    Cog,
    None,
}

#[derive(Debug, Clone)]
pub struct Pbc {
    pub kind: char,
    pub gather: Gather,
}

impl Pbc {
    pub fn from_args(args: &Arguments) -> Result<Self, String> {
        let v = args.values("pbc");
        let kind = v.first().and_then(|s| s.chars().next()).unwrap_or('v');
        if !"rtcv".contains(kind) {
            return Err(format!("@pbc: unknown boundary type '{kind}'"));
        }
        if kind == 't' || kind == 'c' {
            return Err(format!("@pbc {kind}: only rectangular (r) and vacuum (v) boxes are supported by the analysis programs"));
        }
        let gather = match v.get(1).map(|s| s.as_str()) {
            None | Some("cog") | Some("0") => Gather::Cog,
            Some("nog") | Some("nogather") => Gather::None,
            Some(other) => {
                return Err(format!(
                    "@pbc: gather method '{other}' not supported (cog, nog)"
                ))
            },
        };
        Ok(Pbc { kind, gather })
    }

    /// The periodicity of one frame.
    pub fn periodicity(&self, box_dims: Vec3) -> Periodicity {
        if self.kind == 'r' && box_dims.x > 0.0 {
            Periodicity::Rectangular(Rectangular::new(box_dims))
        } else {
            Periodicity::Vacuum(Vacuum)
        }
    }

    /// gromos++ `Boundary::coggather` on `positions` (solute molecules from the topology's
    /// molecule ranges, then the solvent in molecules of `atoms_per_solvent`).
    pub fn gather(&self, topo: &Topology, positions: &mut [Vec3], periodicity: &Periodicity) {
        if self.gather == Gather::None || matches!(periodicity, Periodicity::Vacuum(_)) {
            return;
        }
        let mols: Vec<std::ops::Range<usize>> = if topo.molecules.is_empty() {
            vec![std::ops::Range {
                start: 0,
                end: topo.num_solute_atoms().min(positions.len()),
            }]
        } else {
            topo.molecules.clone()
        };
        let Some(first) = mols.first().filter(|m| !m.is_empty()) else {
            return;
        };
        chain(positions, first.clone(), Vec3::ZERO, periodicity);
        let n = first.len() as f64;
        let cog = positions[first.clone()]
            .iter()
            .fold(Vec3::ZERO, |a, p| a + *p)
            / n;
        for m in mols.iter().skip(1) {
            chain(positions, m.clone(), cog, periodicity);
        }
        let n_solute = topo.num_solute_atoms().min(positions.len());
        let per = topo.atoms_per_solvent().max(1);
        let mut i = n_solute;
        while i + per <= positions.len() {
            chain(positions, i..i + per, cog, periodicity);
            i += per;
        }
    }
}

/// First atom to its image nearest to `reference`, every following atom nearest to its predecessor.
fn chain(
    positions: &mut [Vec3],
    range: std::ops::Range<usize>,
    reference: Vec3,
    periodicity: &Periodicity,
) {
    let mut prev = reference;
    for i in range {
        // gromos++ nearestImage(prev, pos): the image position = prev − (prev − pos)_min
        positions[i] = prev - periodicity.nearest_image(prev, positions[i]);
        prev = positions[i];
    }
}
