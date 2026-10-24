//! C FFI interface for integration with C++ GROMOS code
//!
//! # Status
//! These entry points define the stable C ABI. The bodies are stubs pending
//! the `newen-forces` nonbonded innerloop being wired into a C-shim crate.
//! `newen-core` cannot depend on `newen-forces` (would be circular).
//!
//! TODO: move this module to a dedicated `newen-ffi` crate that depends on
//! both `newen-core` and `newen-forces`, then rewire the bodies.

use std::slice;

// ABI-equivalent aliases — no libc dependency needed.
type CDouble = f64;
type CFloat = f32;
type CUint = u32;

/// C-compatible representation of 3D vector
#[repr(C)]
pub struct CVec3 {
    pub x: CFloat,
    pub y: CFloat,
    pub z: CFloat,
}

/// C-compatible LJ parameters
#[repr(C)]
pub struct CLJParameters {
    pub c6: CDouble,
    pub c12: CDouble,
}

/// C FFI: LJ + CRF inner loop
///
/// # Safety
/// Caller must ensure all pointers are valid and have correct length.
///
/// # Parameters
/// * `positions`      – flat `[x0,y0,z0,…]` of length n_atoms×3
/// * `charges`        – length n_atoms
/// * `iac`            – integer atom codes (atom types), length n_atoms
/// * `n_atoms`        – number of atoms
/// * `pairlist`       – flat `[i0,j0,i1,j1,…]` of length n_pairs×2
/// * `n_pairs`        – number of pairs
/// * `lj_params`      – flat `[c6,c12,…]` of length n_types×n_types×2
/// * `n_types`        – number of atom types
/// * `box_data`       – box vectors (3 for rectangular, 9 for triclinic)
/// * `boundary_type`  – 0=vacuum, 1=rectangular, 2=triclinic
/// * `crf_cut`        – CRF cutoff
/// * `crf_2cut3i`     – CRF parameter 2/r_cut³
/// * `crf_cut3i`      – CRF parameter 1/r_cut³
/// * `forces`         – output forces, flat length n_atoms×3
/// * `energies`       – output `[e_lj, e_crf]`
/// * `virial`         – output virial tensor, flat 3×3 = 9 elements
#[no_mangle]
pub unsafe extern "C" fn rust_lj_crf_innerloop(
    positions: *const CFloat,
    charges: *const CFloat,
    iac: *const CUint,
    n_atoms: CUint,
    pairlist: *const CUint,
    n_pairs: CUint,
    lj_params: *const CDouble,
    n_types: CUint,
    box_data: *const CDouble,
    boundary_type: CUint,
    crf_cut: CDouble,
    crf_2cut3i: CDouble,
    crf_cut3i: CDouble,
    forces: *mut CFloat,
    energies: *mut CDouble,
    virial: *mut CDouble,
) {
    // Suppress unused-parameter warnings until the body is wired up.
    let _ = (
        positions,
        charges,
        iac,
        n_atoms,
        pairlist,
        n_pairs,
        lj_params,
        n_types,
        box_data,
        boundary_type,
        crf_cut,
        crf_2cut3i,
        crf_cut3i,
    );

    if forces.is_null() || energies.is_null() || virial.is_null() {
        return;
    }

    let n = n_atoms as usize;
    let forces_slice = slice::from_raw_parts_mut(forces, n * 3);
    let energies_slice = slice::from_raw_parts_mut(energies, 2);
    let virial_slice = slice::from_raw_parts_mut(virial, 9);

    forces_slice.fill(0.0);
    energies_slice.fill(0.0);
    virial_slice.fill(0.0);
}

/// C FFI: parallel version of the LJ + CRF inner loop (stub — delegates to serial for now)
#[no_mangle]
pub unsafe extern "C" fn rust_lj_crf_innerloop_parallel(
    positions: *const CFloat,
    charges: *const CFloat,
    iac: *const CUint,
    n_atoms: CUint,
    pairlist: *const CUint,
    n_pairs: CUint,
    lj_params: *const CDouble,
    n_types: CUint,
    box_data: *const CDouble,
    boundary_type: CUint,
    crf_cut: CDouble,
    crf_2cut3i: CDouble,
    crf_cut3i: CDouble,
    forces: *mut CFloat,
    energies: *mut CDouble,
    virial: *mut CDouble,
) {
    rust_lj_crf_innerloop(
        positions,
        charges,
        iac,
        n_atoms,
        pairlist,
        n_pairs,
        lj_params,
        n_types,
        box_data,
        boundary_type,
        crf_cut,
        crf_2cut3i,
        crf_cut3i,
        forces,
        energies,
        virial,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // Ignored until newen-ffi wires the real innerloop; keeps the ABI test
    // in place so it's easy to un-ignore once the body is implemented.
    #[test]
    #[ignore = "stub: rewire body before enabling (see ffi.rs module doc)"]
    fn test_ffi_simple() {
        let positions: Vec<f32> = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let charges: Vec<f32> = vec![0.5, -0.5];
        let iac: Vec<u32> = vec![0, 0];
        let pairlist: Vec<u32> = vec![0, 1];
        let lj_params: Vec<f64> = vec![0.001, 0.0001];
        let box_data: Vec<f64> = vec![10.0, 10.0, 10.0];
        let mut forces: Vec<f32> = vec![0.0; 6];
        let mut energies: Vec<f64> = vec![0.0; 2];
        let mut virial: Vec<f64> = vec![0.0; 9];

        unsafe {
            rust_lj_crf_innerloop(
                positions.as_ptr(),
                charges.as_ptr(),
                iac.as_ptr(),
                2,
                pairlist.as_ptr(),
                1,
                lj_params.as_ptr(),
                1,
                box_data.as_ptr(),
                0,
                1.4,
                0.364431,
                0.182215,
                forces.as_mut_ptr(),
                energies.as_mut_ptr(),
                virial.as_mut_ptr(),
            );
        }

        assert!((forces[0] + forces[3]).abs() < 1e-5);
        assert!(energies[0] != 0.0);
        assert!(energies[1] != 0.0);
    }
}
