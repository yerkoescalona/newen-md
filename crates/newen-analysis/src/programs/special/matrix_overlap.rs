//! matrix_overlap — overlap of two covariance matrices (gromos++ `matrix_overlap`, Hess 2002):
//! Ω = 1 − sqrt(tr((M1^½ − M2^½)²)) / sqrt(tr M1 + tr M2).
//!
//! Usage: matrix_overlap @m1 <file> @m2 <file> @dim <n>

use crate::args::Arguments;
use crate::distribution::trim_float;
use crate::programs::io_err;
use nalgebra::DMatrix;
use newen_io::args::{Program, Status};
use newen_io::table::read_columns;
use std::io::Write;

pub const USAGE: &str = "# matrix_overlap
\t@m1           <matrix1>
\t@m2             <matrix2>
\t@dim             <integer: dimension>";

pub const KNOWN: &[&str] = &["m1", "m2", "dim"];

pub const PROGRAM: Program = Program {
    name: "matrix_overlap",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

fn read_matrix(path: &str, dim: usize) -> Result<DMatrix<f64>, String> {
    let rows = read_columns(path).map_err(|e| format!("{path}: {e}"))?;
    let flat: Vec<f64> = rows.into_iter().flatten().collect();
    if flat.len() < dim * dim {
        return Err(format!(
            "{path}: {} numbers, {} needed for dim {dim}",
            flat.len(),
            dim * dim
        ));
    }
    Ok(DMatrix::from_row_slice(dim, dim, &flat[..dim * dim]))
}

/// V·diag(sqrt|λ|)·Vᵀ — the matrix square root gromos++ builds from the eigendecomposition.
fn sqrt_sym(m: &DMatrix<f64>) -> DMatrix<f64> {
    let eig = m.clone().symmetric_eigen();
    let d = DMatrix::from_diagonal(&eig.eigenvalues.map(|l| l.abs().sqrt()));
    &eig.eigenvectors * d * eig.eigenvectors.transpose()
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let dim: usize = args.require("dim")?;
    if dim == 0 {
        return Err("dim must be bigger than zero".into());
    }
    let m1 = read_matrix(args.value("m1")?, dim)?;
    let m2 = read_matrix(args.value("m2")?, dim)?;
    let tr1 = m1.trace();
    let tr2 = m2.trace();
    if tr1 < 0.0 {
        return Err(" trace of the covariance matrix 1 is negative. This will make problems with the rest of the calculations.".into());
    }
    if tr2 < 0.0 {
        return Err(" trace of the covariance matrix 2 is negative. This will make problems with the rest of the calculations.".into());
    }
    writeln!(out, "TRACE M1:  {}", trim_float(tr1)).map_err(io_err)?;
    writeln!(out, "TRACE M2:  {}", trim_float(tr2)).map_err(io_err)?;
    let diff = sqrt_sym(&m1) - sqrt_sym(&m2);
    let trace = (&diff * &diff).trace();
    let d = trace.sqrt();
    writeln!(out, "ABSOLUTE DIFFERENCE: {}", trim_float(d)).map_err(io_err)?;
    writeln!(out, "Normalizing...").map_err(io_err)?;
    writeln!(
        out,
        "NORMALIZED OVERLAP: {}",
        trim_float(1.0 - d / (tr1 + tr2).sqrt())
    )
    .map_err(io_err)?;
    Ok(())
}
