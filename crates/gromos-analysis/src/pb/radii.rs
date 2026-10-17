//! PB radii from the topology's Lennard-Jones parameters: gromos++ `utils::compute_atomic_radii_vdw`
//! (`AtomicRadii.cc:64-86`) and the per-atom post-processing of `dGslv_pbsolv.cc:689-716`.

use gromos_core::Topology;

/// gromos++ treats a C6 below this as no dispersion and gives the atom radius 0
/// (`AtomicRadii.cc:65`).
const SMALL_C6: f64 = 1.0E-20;

/// Which Lennard-Jones length the radius is (`@rminORsigma`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RadiusDefinition {
    /// The minimum of the probe–atom potential, rmin = (2 C12/C6)^(1/6).
    Rmin,
    /// σ = rmin / 2^(1/6).
    Sigma,
}

/// rmin(atom, probe) − `probe_radius` for every solute atom of `topo`, from the LJ pair
/// (IAC of the atom, `probe_iac`); 0 where C6 < 1e-20. `probe_iac` is 0-based.
pub fn atomic_radii_vdw(
    topo: &Topology,
    probe_iac: usize,
    probe_radius: f64,
) -> Result<Vec<f64>, String> {
    let n_types = topo.num_atom_types();
    if probe_iac >= n_types {
        return Err(format!(
            "probe IAC {} is not an atom type of the topology ({n_types} types)",
            probe_iac + 1
        ));
    }
    (0..topo.num_solute_atoms())
        .map(|a| {
            let iac = topo.iac[a];
            if iac >= n_types {
                return Err(format!(
                    "atom {} has IAC {} beyond the LJ table",
                    a + 1,
                    iac + 1
                ));
            }
            let lj = topo.lj_parameter(iac, probe_iac);
            Ok(if lj.c6 >= SMALL_C6 {
                ((2.0 * lj.c12 / lj.c6).ln() / 6.0).exp() - probe_radius
            } else {
                0.0
            })
        })
        .collect()
}

/// The radius `dGslv_pbsolv` gives one `@atoms` atom (`dGslv_pbsolv.cc:695-716`): a radius of
/// (almost) zero becomes `hydrogen_radius`, unscaled; otherwise rmin is turned into σ on request
/// and scaled by `radscal`.
pub fn pb_radius(
    r: f64,
    definition: RadiusDefinition,
    probe_radius: f64,
    hydrogen_radius: f64,
    radscal: f64,
    tiny_real: f64,
) -> f64 {
    if r.abs() < tiny_real {
        hydrogen_radius
    } else {
        match definition {
            RadiusDefinition::Sigma => {
                let tmprad = (r + probe_radius) / ((1.0 / 6.0) * 2.0f64.ln()).exp() - probe_radius;
                tmprad * radscal
            },
            RadiusDefinition::Rmin => r * radscal,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_io::topology::{build_topology, read_topology_file};

    fn aladip() -> Topology {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gromos-md/tests/gromosXX_references/shared/aladip.topo");
        build_topology(read_topology_file(path).unwrap())
    }

    /// rmin(atom, OW) − 0.14 nm for the aladip atom types, OW being IAC 4.
    #[test]
    fn aladip_radii_against_water_oxygen() {
        let topo = aladip();
        let r = atomic_radii_vdw(&topo, 3, 0.14).unwrap();
        // CB(CH3) C O N H CA(CH1) CB(CH3) C O N H CB(CH3)
        let expected = [
            0.247, 0.226, 0.198, 0.228, 0.0, 0.307, 0.247, 0.226, 0.198, 0.228, 0.0, 0.247,
        ];
        assert_eq!(r.len(), expected.len());
        for (i, (&got, &want)) in r.iter().zip(&expected).enumerate() {
            assert!((got - want).abs() < 5e-4, "atom {}: {got} vs {want}", i + 1);
        }
        // Hydrogens have no C6 with OW and become the hydrogen radius; the rest are scaled.
        let tiny = 1e-10;
        assert_eq!(
            pb_radius(r[4], RadiusDefinition::Rmin, 0.14, 0.05, 1.1, tiny),
            0.05
        );
        assert_eq!(
            pb_radius(r[2], RadiusDefinition::Rmin, 0.14, 0.05, 1.1, tiny),
            r[2] * 1.1
        );
        let sigma = pb_radius(r[2], RadiusDefinition::Sigma, 0.14, 0.05, 1.0, tiny);
        assert!((sigma - ((r[2] + 0.14) / 2f64.powf(1.0 / 6.0) - 0.14)).abs() < 1e-12);
        assert!(atomic_radii_vdw(&topo, 1000, 0.14).is_err());
    }
}
