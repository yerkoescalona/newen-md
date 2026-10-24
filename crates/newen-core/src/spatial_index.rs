//! Spatial-index query service — general "neighbors within r" queries.
//!
//! The MD pairlist ([`crate::pairlist::PairlistContainer`]) is shaped for one purpose:
//! gromosXX's own twin-range, charge-group-based, exclusion-filtered nonbonded loop, built
//! once per run for a fixed cutoff pair. It cannot double as a general query service —
//! an ML model's radial graph or a QM-zone gather may ask for a different cutoff at query
//! time, and wants the literal geometric neighbor graph, not GROMOS's electrostatic
//! exclusions. `SpatialIndex` is therefore a separate, atom-level distance query, not a
//! wrapper around the pairlist's internals (see `architecture.md` "Supporting shapes" and
//! P3 (`newen-forces/.claude/CONTEXT.md`): one query-based service meant to serve the MD
//! pairlist, an ML radial graph, and the QM-zone gather from the same shared positions).

use crate::configuration::Configuration;
use crate::math::{Periodicity, Vec3};
use crate::selection::AtomSelection;

/// Query-based neighbor service: "give me atom pairs within `r` of this selection."
///
/// Pairs-only for now; graph/triplet queries (P3) are additive extensions of
/// this trait later, not a redesign.
pub trait SpatialIndex {
    /// Atom-index pairs `(i, j)` with `i < j`, where at least one of `i, j` is in `atoms`
    /// and the nearest-image distance is within `cutoff`. `shift` is the periodic-image
    /// shift applied to atom `j` so that `pos[i] - (pos[j] + shift)` equals the returned
    /// minimum-image vector (zero in vacuum / non-periodic boundaries). Carrying the shift
    /// avoids a breaking signature change once a periodic ML/QM provider needs it (GROMACS's
    /// NNPot interface passes PBC shift vectors alongside pair indices for the same reason).
    fn neighbor_pairs(&self, atoms: &AtomSelection, cutoff: f64) -> Vec<(usize, usize, Vec3)>;
}

/// Brute-force O(|selection| x N) neighbor query over a configuration's current positions.
///
/// Deliberately atom-level and exclusion-free (unlike the MD pairlist): a generic consumer
/// wants the geometric neighbor graph at whatever cutoff it asks for. Fine for the region
/// sizes providers query (an ML/QM zone, not the whole system); if this becomes a hot path
/// for large regions, swap in a cell-list-backed implementation behind the same trait.
pub struct ConfigurationSpatialIndex<'a> {
    positions: &'a [Vec3],
    periodicity: &'a Periodicity,
}

impl<'a> ConfigurationSpatialIndex<'a> {
    /// Borrow `conf`'s current positions and `periodicity`'s boundary condition for queries.
    pub fn new(conf: &'a Configuration, periodicity: &'a Periodicity) -> Self {
        Self {
            positions: &conf.current().pos,
            periodicity,
        }
    }
}

impl SpatialIndex for ConfigurationSpatialIndex<'_> {
    fn neighbor_pairs(&self, atoms: &AtomSelection, cutoff: f64) -> Vec<(usize, usize, Vec3)> {
        let cutoff2 = cutoff * cutoff;
        let n = self.positions.len();

        let mut in_selection = vec![false; n];
        for i in atoms.iter() {
            in_selection[i] = true;
        }

        let mut pairs = Vec::new();
        for &i in atoms.indices() {
            #[allow(clippy::needless_range_loop)] // one index over several arrays
            for j in 0..n {
                if j == i {
                    continue;
                }
                // Both i and j in the selection: only emit each pair once (from the lower index).
                if in_selection[j] && j < i {
                    continue;
                }
                let (lo, hi) = if i < j { (i, j) } else { (j, i) };
                let r = self
                    .periodicity
                    .nearest_image(self.positions[lo], self.positions[hi]);
                if r.length_squared() <= cutoff2 {
                    let shift = (self.positions[lo] - self.positions[hi]) - r;
                    pairs.push((lo, hi, shift));
                }
            }
        }
        pairs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{Rectangular, Vacuum};

    fn conf_with_positions(positions: Vec<Vec3>) -> Configuration {
        let n = positions.len();
        let mut conf = Configuration::new(n, 1, 1);
        conf.current_mut().pos = positions;
        conf
    }

    #[test]
    fn vacuum_finds_pairs_within_cutoff() {
        let conf = conf_with_positions(vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
        ]);
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let region = AtomSelection::from_indices(vec![0], 3).unwrap();

        let pairs = index.neighbor_pairs(&region, 1.0);

        assert_eq!(pairs.len(), 1);
        assert_eq!((pairs[0].0, pairs[0].1), (0, 1));
        assert_eq!(pairs[0].2, Vec3::ZERO); // no periodicity, no shift
    }

    #[test]
    fn rectangular_wraps_and_reports_shift() {
        // Box of 1.0 nm; atom 0 near the low edge, atom 1 near the high edge —
        // the minimum image is across the periodic boundary, not the direct distance.
        let box_size = Vec3::new(1.0, 1.0, 1.0);
        let conf = conf_with_positions(vec![Vec3::new(0.05, 0.5, 0.5), Vec3::new(0.95, 0.5, 0.5)]);
        let periodicity = Periodicity::Rectangular(Rectangular::new(box_size));
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let region = AtomSelection::all(2);

        let pairs = index.neighbor_pairs(&region, 0.2);

        assert_eq!(pairs.len(), 1);
        let (i, j, shift) = pairs[0];
        assert_eq!((i, j), (0, 1));
        // direct distance is 0.9 nm (outside cutoff); minimum-image distance is 0.1 nm.
        let direct = conf.current().pos[i] - conf.current().pos[j];
        let minimum_image = direct - shift;
        assert!((minimum_image.length() - 0.1).abs() < 1e-9);
    }
}
