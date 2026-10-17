//! Radial distribution function, as gromos++ `utils::RDF` computes it[^gromospp].
//!
//! For each frame and each centre atom $c$, the nearest-image distances to every with-atom are
//! binned over $[0, r_\mathrm{cut})$ in `grid` bins of width $\Delta r$, and bin $k$ (centre
//! $r_k$) contributes
//!
//! $$
//! g_c(r_k) = \frac{N_c(k)}{4\pi r_k^2 \, \Delta r \, \rho_c},
//! \qquad \rho_c = \frac{N_\mathrm{with} - \delta_c}{V},
//! $$
//!
//! where $\delta_c = 1$ when $c$ is itself a with-atom. The result is the mean over all frames
//! and centres. The shell volume is $4\pi r_k^2 \Delta r$, not the exact
//! $\tfrac{4}{3}\pi(r_\mathrm{out}^3 - r_\mathrm{in}^3)$; this is gromos++'s choice and is kept.
//!
//! **gromos++ quirk (kept):** with `nointra`, pairs in the same molecule are skipped, and every
//! solvent atom belongs to one shared pseudo-molecule (`AtomSpecifier::mol` = −1), so all
//! solvent–solvent pairs are skipped, not just those within one solvent molecule.
//!
//! [^gromospp]: `gromos++/src/utils/RDF.cc`, `calculateAll` and `calculateInter`; the program is
//! `gromos++/programs/rdf.cc`.

use crate::distribution::Distribution;
use newen_core::math::Periodicity;
use newen_core::Vec3;

/// A g(r) accumulated over frames.
#[derive(Debug, Clone)]
pub struct Rdf {
    cut: f64,
    grid: usize,
    empty: Distribution,
    sum: Vec<f64>,
    centre_frames: usize,
}

impl Rdf {
    /// `grid` bins over $[0, \mathrm{cut})$.
    pub fn new(cut: f64, grid: usize) -> Result<Self, String> {
        Ok(Rdf {
            cut,
            grid,
            empty: Distribution::new(0.0, cut, grid)?,
            sum: vec![0.0; grid],
            centre_frames: 0,
        })
    }

    /// Add one frame. `molecule_of[a]` is atom `a`'s molecule (−1 for every solvent atom) and is
    /// read only when `nointra`; `volume` is the box volume.
    #[allow(clippy::too_many_arguments)]
    pub fn add_frame(
        &mut self,
        positions: &[Vec3],
        centres: &[usize],
        withs: &[usize],
        molecule_of: &[i64],
        periodicity: &Periodicity,
        volume: f64,
        nointra: bool,
    ) {
        let mut is_with = vec![false; positions.len()];
        for &w in withs {
            is_with[w] = true;
        }
        let correct = 4.0 * std::f64::consts::PI * self.cut / self.grid as f64;
        for &c in centres {
            let mut dist = self.empty.clone();
            for &w in withs {
                let skip = if nointra {
                    molecule_of[w] == molecule_of[c]
                } else {
                    w == c
                };
                if !skip {
                    dist.add(
                        periodicity
                            .nearest_image(positions[w], positions[c])
                            .length(),
                    );
                }
            }
            let inwith = usize::from(is_with[c]);
            let dens = (withs.len() - inwith) as f64 / volume;
            for (k, s) in self.sum.iter_mut().enumerate() {
                let r = dist.centre(k);
                *s += dist.count(k) as f64 / (dens * correct * r * r);
            }
        }
        self.centre_frames += centres.len();
    }

    /// g(r) per bin: the sum over frames and centres divided by (frames × centres).
    pub fn g(&self) -> Vec<f64> {
        let n = self.centre_frames as f64;
        self.sum.iter().map(|s| s / n).collect()
    }

    /// The centre of bin `k`, as gromos++ prints it.
    pub fn r(&self, k: usize) -> f64 {
        (k as f64 + 0.5) * self.cut / self.grid as f64
    }

    pub fn grid(&self) -> usize {
        self.grid
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::{Rectangular, Vacuum};

    #[test]
    fn nearest_image_counts_across_the_boundary() {
        // two atoms 0.1 apart through the x boundary of a 2 nm box
        let pos = [Vec3::new(0.05, 1.0, 1.0), Vec3::new(1.95, 1.0, 1.0)];
        let pbc = Periodicity::Rectangular(Rectangular::new(Vec3::splat(2.0)));
        let mut rdf = Rdf::new(0.5, 5).unwrap();
        rdf.add_frame(&pos, &[0, 1], &[0, 1], &[-1, -1], &pbc, 8.0, false);
        let g = rdf.g();
        // bin 1 (r = 0.15): one neighbour, density 1/8, 4π·0.1·0.15²
        let want = 1.0 / ((1.0 / 8.0) * 4.0 * std::f64::consts::PI * 0.1 * 0.15 * 0.15);
        assert!((g[1] - want).abs() < 1e-12 * want);
        assert_eq!(g.iter().filter(|x| **x != 0.0).count(), 1);
    }

    #[test]
    fn nointra_skips_one_molecule_and_all_solvent() {
        let pos = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::new(0.0, 0.25, 0.0),
            Vec3::new(0.0, 0.0, 0.35),
        ];
        let pbc = Periodicity::Vacuum(Vacuum);
        let mol = [0, 0, -1, -1];
        let mut rdf = Rdf::new(0.5, 5).unwrap();
        rdf.add_frame(&pos, &[0], &[1, 2, 3], &mol, &pbc, 1.0, true);
        let g = rdf.g();
        assert_eq!((g[1], g[2] > 0.0, g[3] > 0.0), (0.0, true, true));
        let mut solvent = Rdf::new(0.5, 5).unwrap();
        solvent.add_frame(&pos, &[2], &[2, 3], &mol, &pbc, 1.0, true);
        assert!(solvent.g().iter().all(|x| *x == 0.0));
    }
}
