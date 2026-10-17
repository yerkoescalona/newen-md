# Python Bindings Tutorial

A walkthrough of how the Rust → Python binding layer is structured, for
developers who want to add new classes or extend existing ones.

## The stack

```
Python:  newen.System(topo, conf)
              │
              ▼
PyO3:    #[pyclass] PySystem  ←  crates/pyo3-newen/src/system.rs
              │
              ▼
Rust:    newen_core::Topology + newen_io coordinate types
```

[PyO3](https://pyo3.rs) is the Rust library that generates the CPython extension.
[Maturin](https://github.com/PyO3/maturin) is the build tool that packages it.

## Anatomy of a PyO3 class

Here is the `System` class, simplified:

```rust
// crates/pyo3-newen/src/system.rs

use pyo3::prelude::*;
use crate::topology::PyTopology;
use crate::py_conf::PyConfiguration;

#[pyclass(name = "System")]   // exposed to Python as `newen.System`
pub struct PySystem {
    pub(crate) topology: PyTopology,
    pub(crate) configuration: PyConfiguration,
}

#[pymethods]
impl PySystem {
    #[new]                     // maps to Python __init__
    fn new(topology: &Bound<'_, PyTopology>,
           configuration: &Bound<'_, PyConfiguration>) -> PyResult<Self> {
        // validate, then construct
    }

    #[staticmethod]
    fn from_files(topo_file: &str, conf_file: &str) -> PyResult<Self> { … }

    #[getter]                  // Python property
    fn n_atoms(&self) -> usize { … }

    fn write(&self, path: &str) -> PyResult<()> { … }
}

pub fn register_system(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySystem>()?;
    Ok(())
}
```

Key PyO3 attributes:

| Attribute | Effect |
|-----------|--------|
| `#[pyclass(name = "X")]` | Exposes the struct to Python as class `X` |
| `#[pymethods]` | Marks the `impl` block containing Python-callable methods |
| `#[new]` | Maps to `__init__` / `__new__` |
| `#[staticmethod]` | Creates a Python `@staticmethod` |
| `#[getter]` | Creates a read-only Python `@property` |
| `PyResult<T>` | Return type that maps Rust errors to Python exceptions |

## Registering a new class

1. Create `crates/pyo3-newen/src/myfeature.rs`.
2. Add a `register_myfeature(m)` function in it.
3. Call it from `crates/pyo3-newen/src/lib.rs`:

```rust
mod myfeature;

#[pymodule]
fn newen(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // existing registrations …
    myfeature::register_myfeature(m)?;
    Ok(())
}
```

4. Declare `pub mod myfeature;` in `lib.rs`.
5. Re-export the new Python name in `py-newen/python/newen/__init__.py`.
6. Add the type signature to `py-newen/python/newen/newen.pyi`.

## Error handling

Rust errors become Python exceptions via `PyResult`:

```rust
// Raise ValueError
Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
    format!("Topology has {} atoms but configuration has {}", n_topo, n_conf)
))

// Raise IOError
Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
```

Common exception types: `PyValueError`, `PyIOError`, `PyRuntimeError`, `PyTypeError`.

## Returning NumPy arrays

Use the `numpy` crate (already a dependency):

```rust
use numpy::{PyArray2, IntoPyArray};
use pyo3::prelude::*;

#[getter]
fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let flat: Vec<f64> = self.pos_data.iter()
        .flat_map(|v| [v.x, v.y, v.z])
        .collect();
    let arr = ndarray::Array2::from_shape_vec((self.pos_data.len(), 3), flat)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
    Ok(arr.into_pyarray(py).into())
}
```

The returned array is a copy of the Rust data. True zero-copy would require a
structure-of-arrays (SoA) memory layout — board card V1.

## Optional arguments

Use `#[pyo3(signature = (arg1, arg2, arg3=None))]`:

```rust
#[new]
#[pyo3(signature = (arg1, arg2, arg3=None))]
fn new(arg1: &Bound<'_, PyAny>,
       arg2: &Bound<'_, PyAny>,
       arg3: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
    match arg3 {
        None    => { /* two-arg form */ }
        Some(a) => { /* three-arg form */ }
    }
}
```

## Testing without the Python interpreter

PyO3 classes cannot be instantiated in `#[test]` without special linker setup.
Extract the core logic into a plain Rust function and test that:

```rust
// In system.rs
pub(crate) fn validate_atom_count_match(n_topo: usize, n_conf: usize) -> Result<(), String> {
    if n_topo != n_conf {
        Err(format!("Topology has {} atoms but configuration has {}", n_topo, n_conf))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_count_mismatch_error() {
        let err = validate_atom_count_match(648, 3).unwrap_err();
        assert!(err.contains("must match"));
    }
}
```

Integration tests that need a real Python interpreter live in
`py-newen/tests/` and run via `make test-python`.

## Type stubs

Every public Python name needs a signature in `newen.pyi`:

```python
class System:
    def __init__(self, topology: Topology, configuration: Configuration) -> None: ...
    @staticmethod
    def from_files(topo_file: str, conf_file: str) -> System: ...
    @property
    def n_atoms(self) -> int: ...
```

Stubs enable IDE completion and `mypy` type checking without needing the compiled
extension.
