//! Python extension module for newen
//!
//! This is the main entry point for the Python `newen` package.
//! It uses the pyo3-newen crate for the low-level bindings.

use pyo3::prelude::*;

// Re-export bindings from pyo3-newen
use pyo3_newen::register_bindings;

/// newen Python Module
///
/// High-performance molecular dynamics engine written in Rust.
#[pymodule]
fn newen(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Register all bindings from pyo3-newen
    register_bindings(m)?;

    // Add module metadata
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", "Yerko Escalona")?;

    Ok(())
}
