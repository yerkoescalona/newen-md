//! `SchNetInteraction` — a real SchNetPack SchNet model as a [`PotentialProvider`].
//!
//! Feature-gated (`ml`): loads a TorchScript module (see `scripts/export_toy_schnet.py`)
//! via `tch` (libtorch bindings) and runs it **in-process** — no Python subprocess, no
//! embedded interpreter, unlike gromosXX's `nn_worker`, which drives SchNetPack through CPython.
//!
//! **Version pinning is load-bearing, not incidental.** `tch`'s C++ glue (`torch-sys`) is
//! generated against one specific libtorch version; a mismatch is a build/link-time
//! failure, not a silent bug — this crate pins `tch = "0.24.0"`, which requires
//! **libtorch 2.11.0** exactly (`torch-sys`'s own build script checks this). Verified in
//! this repo's sandbox: a conda-forge MKL-linked libtorch failed to even load in Python
//! with an unrelated "cannot enable executable stack" restriction, and the official PyPI
//! CPU wheel works *only* at the exact pinned version — anything else and `torch-sys`
//! refuses to build rather than link something ABI-incompatible. `torch==2.11.0` (CPU wheel)
//! and `schnetpack` (for `scripts/export_toy_schnet.py` only) are declared as the `ml`
//! dependency-group in `py-newen/pyproject.toml`, pinned to the official PyTorch CPU index —
//! reproducible via `uv`, no hand-built throwaway venv. To build with `--features ml`:
//! ```text
//! cd py-newen && uv sync --group ml && cd ..
//! source py-newen/.venv/bin/activate
//! export LIBTORCH_USE_PYTORCH=1
//! export LD_LIBRARY_PATH="$(python3 -c 'import torch,os;print(os.path.dirname(torch.__file__))')/lib:$LD_LIBRARY_PATH"
//! cargo test -p newen-forces --features ml schnet::
//! ```
//! (`uv sync --group ml` alone drops `dev`/`notebooks` from the venv — use
//! `uv sync --all-groups --group ml` to keep everything installed together.)
//!
//! **A real SchNetPack model, not a reimplementation — this changed mid-session.** The
//! first version of this file hand-rolled a single cfconv block to prove *something*
//! TorchScript-loadable could work, because GROMOS's own real trained tutorial model
//! (SchNetPack v1 — see `.local/gromos_tutorial_livecoms/tutorial_files/t_06`) fails to
//! `torch.jit.script` outright: `schnetpack/nn/neighbors.py` has a Python-style
//! conditionally-typed return TorchScript's static type system rejects. That's a genuine
//! instance of P2 ("not every architecture exports cleanly") — but it turned out
//! to be a **v1 library problem, not an inherent SchNet problem**: verified directly that
//! **SchNetPack 2.x's `schnetpack.model.NeuralNetworkPotential` (real `SchNet`
//! representation + a `Forces` output module) scripts cleanly**, once `Forces(calc_forces=
//! true)` is included (a bare `Atomwise`-only model fails TorchScript's type inference on
//! an empty `required_derivatives` list — `Forces` is what populates it). So this provider
//! now loads a real SchNetPack 2 model, built from the actual library — still randomly
//! initialized/untrained (no chemistry-accuracy claim; that's orthogonal to the seam this
//! proves), but the genuine architecture and code, not a hand-rolled stand-in.
//!
//! **Input/output contract is `Dict[str, Tensor]`, matching `schnetpack.properties`, not
//! positional tensors** (SchNetPack 2's own convention, confirmed by inspection and by
//! running it):
//! - in: `_positions` [N,3], `_atomic_numbers` [N] long, `_idx_i`/`_idx_j` [E] long,
//!   `_offsets` [E,3], `_n_atoms` [1] long, `_idx_m` [N] long (which system each atom
//!   belongs to — SchNetPack batches multiple systems together in general; single-system
//!   inference here uses all zeros).
//! - out: `energy` [1], `forces` [N,3] — **both already computed**, unlike the old
//!   hand-rolled model. SchNetPack's `Forces` module runs the energy→forces autograd
//!   *inside* the scripted graph, so this provider no longer calls `.backward()` itself.
//!
//! Critically, `schnetpack.atomistic.PairwiseDistances` does **not** compute its own
//! neighbor list — it expects `_idx_i`/`_idx_j`/`_offsets` as input. That's exactly the
//! [`SpatialIndex`] design this codebase already committed to (P3: one spatial
//! service, host-supplied, serving the MD pairlist and an ML radial graph alike) — `_offsets`
//! is literally the same periodic-image shift `SpatialIndex::neighbor_pairs` already
//! returns, not a coincidence discovered after the fact.
//!
//! **Checked against the real GROMOS tutorial**, not just designed from docs — Tutorial 6
//! ("BuRNN", cloned into `.local/gromos_tutorial_livecoms/tutorial_files/t_06`, gitignored,
//! same pattern as `.local/gromosXX`): confirmed the real `.qmm` file's `QMUNIT` unit-
//! conversion factors and `QMZONE`'s atomic-number-indexed `QMZ` column match this design
//! (`elements` below is indexed by atomic number, `padding_idx=0`, matching the tutorial's
//! own trained model's `Embedding(100, 256, padding_idx=0)`), and ran that real trained
//! model end to end in Python for a sanity check before concluding it needed replacing for
//! the Rust seam specifically because of the scripting failure above.

use newen_core::configuration::Configuration;
use newen_core::math::{Mat3, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::SpatialIndex;
use newen_core::topology::Topology;
use tch::{CModule, Device, IValue, Kind, Tensor};

use crate::provider::{Contribution, Embedding, PotentialProvider, ProviderError, ProviderExtra};

/// A loaded SchNetPack `NeuralNetworkPotential` TorchScript model, run in-process via `tch`.
pub struct SchNetInteraction {
    model: CModule,
    cutoff: f64,
    /// Per-(global-atom-index) atomic number. Mapping `Topology` -> atomic numbers is the
    /// caller's job; this provider just forwards it as `_atomic_numbers`.
    elements: Vec<i64>,
    embedding: Embedding,
}

impl SchNetInteraction {
    /// Load a TorchScript module from `path` (see `scripts/export_toy_schnet.py`).
    /// `elements[i]` is the atomic number for global atom `i`.
    ///
    /// Defaults to [`Embedding::None`]: a SchNetPack `NeuralNetworkPotential`'s input contract
    /// (`_positions`/`_atomic_numbers`/`_idx_i`/`_idx_j`/`_offsets`/…) has no channel for
    /// environment point charges, so the model *cannot* be polarized by them — declaring
    /// anything else would be a lie about what the model computes. Use
    /// [`Self::with_embedding`] once a model that accepts embedding inputs exists.
    pub fn load(path: &str, cutoff: f64, elements: Vec<i64>) -> Result<Self, ProviderError> {
        let model = CModule::load(path).map_err(|e| {
            ProviderError::ComputationFailed(format!("failed to load TorchScript model: {e}"))
        })?;
        Ok(Self {
            model,
            cutoff,
            elements,
            embedding: Embedding::None,
        })
    }

    /// Declare a non-default embedding scheme.
    ///
    /// [`Embedding::Electrostatic`] is accepted here but rejected at `contribute()` time with a
    /// clear error rather than silently ignored — see P2.7 Step 3, which implements the
    /// point-charge input and the force path back onto MM atoms.
    pub fn with_embedding(mut self, embedding: Embedding) -> Self {
        self.embedding = embedding;
        self
    }
}

/// Find `key`'s `Tensor` value in a scripted model's `IValue::GenericDict` output.
fn dict_get(entries: &[(IValue, IValue)], key: &str) -> Result<Tensor, ProviderError> {
    entries
        .iter()
        .find_map(|(k, v)| match (k, v) {
            (IValue::String(s), IValue::Tensor(t)) if s == key => Some(t.shallow_clone()),
            _ => None,
        })
        .ok_or_else(|| ProviderError::ComputationFailed(format!("model output missing '{key}'")))
}

impl PotentialProvider for SchNetInteraction {
    fn contribute(
        &mut self,
        region: &AtomSelection,
        _topo: &Topology,
        conf: &Configuration,
        neigh: &dyn SpatialIndex,
    ) -> Result<Contribution, ProviderError> {
        let atom_indices: Vec<usize> = region.indices().to_vec();
        let n = atom_indices.len();
        if n == 0 {
            return Err(ProviderError::InvalidRegion("empty region".to_string()));
        }
        if atom_indices.iter().any(|&g| g >= self.elements.len()) {
            return Err(ProviderError::InvalidRegion(
                "elements table shorter than region's atom indices".to_string(),
            ));
        }

        // Global atom index -> local (0..n) index within this region, for the model's tensors.
        let local_of = |global: usize| atom_indices.iter().position(|&g| g == global).unwrap();

        let mut flat_pos = Vec::with_capacity(n * 3);
        let mut elements_local = Vec::with_capacity(n);
        for &g in &atom_indices {
            let p = conf.current().pos[g];
            flat_pos.extend_from_slice(&[p.x, p.y, p.z]);
            elements_local.push(self.elements[g]);
        }

        // Both edge directions, each with its own offset (sign flips with direction) — see
        // module docs for the offset-sign derivation against `SpatialIndex`'s shift convention.
        //
        // `SpatialIndex::neighbor_pairs` returns pairs where *either* endpoint is in `region`
        // (by design — P5's embedding case needs that), so cross-boundary pairs must
        // be handled per this provider's declared `Embedding` policy (P2.7 Step 2).
        // They are dropped below under `Embedding::None`, which is what this model's input
        // contract actually supports.
        match self.embedding {
            Embedding::None => {},
            Embedding::Mechanical | Embedding::Electrostatic => {
                return Err(ProviderError::ComputationFailed(format!(
                    "{:?} embedding is not implemented for SchNetInteraction — a SchNetPack \
                     NeuralNetworkPotential has no environment/point-charge input channel \
                     (see 2.7 Step 3)",
                    self.embedding
                )))
            },
        }
        let pairs = neigh.neighbor_pairs(region, self.cutoff);
        let mut idx_i = Vec::with_capacity(pairs.len() * 2);
        let mut idx_j = Vec::with_capacity(pairs.len() * 2);
        let mut offsets = Vec::with_capacity(pairs.len() * 2 * 3);
        for (i, j, shift) in &pairs {
            if !atom_indices.contains(i) || !atom_indices.contains(j) {
                continue;
            }
            let (li, lj) = (local_of(*i) as i64, local_of(*j) as i64);
            idx_i.push(li);
            idx_j.push(lj);
            offsets.extend_from_slice(&[shift.x, shift.y, shift.z]);
            idx_i.push(lj);
            idx_j.push(li);
            offsets.extend_from_slice(&[-shift.x, -shift.y, -shift.z]);
        }
        let n_edges = idx_i.len();

        let positions = Tensor::from_slice(&flat_pos)
            .to_kind(Kind::Float)
            .reshape([n as i64, 3]);
        let atomic_numbers = Tensor::from_slice(&elements_local).to_kind(Kind::Int64);
        let idx_i_t = Tensor::from_slice(&idx_i).to_kind(Kind::Int64);
        let idx_j_t = Tensor::from_slice(&idx_j).to_kind(Kind::Int64);
        let offsets_t = if n_edges == 0 {
            Tensor::zeros([0, 3], (Kind::Float, Device::Cpu))
        } else {
            Tensor::from_slice(&offsets)
                .to_kind(Kind::Float)
                .reshape([n_edges as i64, 3])
        };
        let n_atoms_t = Tensor::from_slice(&[n as i64]).to_kind(Kind::Int64);
        let idx_m_t = Tensor::zeros([n as i64], (Kind::Int64, Device::Cpu));

        let inputs = IValue::GenericDict(vec![
            (
                IValue::String("_positions".to_string()),
                IValue::Tensor(positions),
            ),
            (
                IValue::String("_atomic_numbers".to_string()),
                IValue::Tensor(atomic_numbers),
            ),
            (
                IValue::String("_idx_i".to_string()),
                IValue::Tensor(idx_i_t),
            ),
            (
                IValue::String("_idx_j".to_string()),
                IValue::Tensor(idx_j_t),
            ),
            (
                IValue::String("_offsets".to_string()),
                IValue::Tensor(offsets_t),
            ),
            (
                IValue::String("_n_atoms".to_string()),
                IValue::Tensor(n_atoms_t),
            ),
            (
                IValue::String("_idx_m".to_string()),
                IValue::Tensor(idx_m_t),
            ),
        ]);

        let output = self
            .model
            .forward_is(&[inputs])
            .map_err(|e| ProviderError::ComputationFailed(format!("SchNet forward failed: {e}")))?;

        let entries = match output {
            IValue::GenericDict(entries) => entries,
            other => {
                return Err(ProviderError::ComputationFailed(format!(
                    "expected a Dict output, got {other:?}"
                )))
            },
        };

        let energy_t = dict_get(&entries, "energy")?;
        let forces_t = dict_get(&entries, "forces")?;

        let mut forces = Vec::with_capacity(n);
        for (local_i, &global_i) in atom_indices.iter().enumerate() {
            let fx = forces_t.double_value(&[local_i as i64, 0]);
            let fy = forces_t.double_value(&[local_i as i64, 1]);
            let fz = forces_t.double_value(&[local_i as i64, 2]);
            forces.push((global_i, Vec3::new(fx, fy, fz)));
        }

        Ok(Contribution {
            energy: energy_t.double_value(&[0]), // shape [1]: one system
            forces,
            virial: Mat3::ZERO, // not computed by this proof-of-concept provider
            extra: ProviderExtra::default(),
        })
    }

    fn name(&self) -> &str {
        "SchNet (SchNetPack 2, untrained)"
    }

    fn embedding(&self) -> Embedding {
        self.embedding
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::Periodicity;
    use newen_core::spatial_index::ConfigurationSpatialIndex;
    use newen_core::topology::Topology;

    fn model_path() -> String {
        std::env::var("TOY_SCHNET_MODEL").unwrap_or_else(|_| "/tmp/toy_schnet.pt".to_string())
    }

    /// P2.7 Step 2: an unsupported embedding scheme must fail loudly, not be
    /// silently ignored. Before this, cross-boundary neighbors were dropped unconditionally
    /// with no way for a caller to learn the environment was being ignored.
    #[test]
    fn unsupported_embedding_is_rejected_not_silently_ignored() {
        let path = model_path();
        if !std::path::Path::new(&path).exists() {
            eprintln!("skipping: no model at {path} (run scripts/export_toy_schnet.py first)");
            return;
        }

        let mut interaction = SchNetInteraction::load(&path, 1.0, vec![6, 8, 1])
            .expect("model should load")
            .with_embedding(Embedding::Electrostatic);
        assert_eq!(
            PotentialProvider::embedding(&interaction),
            Embedding::Electrostatic,
            "declared policy should be readable through the trait"
        );

        let topo = Topology::new();
        let mut conf = Configuration::new(3, 1, 1);
        conf.current_mut().pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.15, 0.0, 0.0),
            Vec3::new(0.0, 0.2, 0.0),
        ];
        let region = AtomSelection::all(3);
        let periodicity = Periodicity::Vacuum(newen_core::math::Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let err = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect_err("electrostatic embedding is not implemented and must not silently pass");
        assert!(
            format!("{err}").contains("not implemented"),
            "error should say the scheme is unimplemented, got: {err}"
        );
    }

    /// The default policy must stay [`Embedding::None`] — the only honest description of a
    /// model whose input contract has no environment channel.
    #[test]
    fn default_embedding_is_none() {
        let path = model_path();
        if !std::path::Path::new(&path).exists() {
            eprintln!("skipping: no model at {path} (run scripts/export_toy_schnet.py first)");
            return;
        }
        let interaction =
            SchNetInteraction::load(&path, 1.0, vec![6, 8, 1]).expect("model should load");
        assert_eq!(PotentialProvider::embedding(&interaction), Embedding::None);
    }

    /// The seam's actual validation tier (P8): no GROMOS oracle exists for an
    /// ML potential, so instead check the model's own (internally autograd-computed)
    /// forces are self-consistent with its energy surface via central finite differences —
    /// the standard way any force provider is validated when there's no independent
    /// reference.
    #[test]
    fn model_forces_match_finite_differences() {
        let path = model_path();
        if !std::path::Path::new(&path).exists() {
            eprintln!("skipping: no model at {path} (run scripts/export_toy_schnet.py first)");
            return;
        }

        // Atomic numbers (C, O, H), matching the real embedding convention — see module docs.
        let mut interaction =
            SchNetInteraction::load(&path, 1.0, vec![6, 8, 1]).expect("model should load");

        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.15, 0.0, 0.0),
            Vec3::new(0.0, 0.2, 0.0),
        ];
        let topo = Topology::new();
        let mut conf = Configuration::new(3, 1, 1);
        conf.current_mut().pos = positions.clone();
        let region = AtomSelection::all(3);
        let periodicity = Periodicity::Vacuum(newen_core::math::Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let contribution = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect("forward pass should succeed");
        assert_eq!(contribution.forces.len(), 3);

        let mut energy_at = |positions: &[Vec3]| -> f64 {
            let mut conf = Configuration::new(3, 1, 1);
            conf.current_mut().pos = positions.to_vec();
            let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
            interaction
                .contribute(&region, &topo, &conf, &index)
                .unwrap()
                .energy
        };

        let h = 1e-4;
        for atom in 0..3 {
            for axis in 0..3 {
                let mut plus = positions.clone();
                let mut minus = positions.clone();
                let delta = match axis {
                    0 => Vec3::new(h, 0.0, 0.0),
                    1 => Vec3::new(0.0, h, 0.0),
                    _ => Vec3::new(0.0, 0.0, h),
                };
                plus[atom] += delta;
                minus[atom] -= delta;

                let finite_diff_force = -(energy_at(&plus) - energy_at(&minus)) / (2.0 * h);
                let (_, model_force) = contribution.forces[atom];
                let model_component = match axis {
                    0 => model_force.x,
                    1 => model_force.y,
                    _ => model_force.z,
                };

                // Absolute, not relative: SchNetPack 2's real interaction block (RBF expansion,
                // filter network, cutoff function) has more layers than the old hand-rolled
                // single-block model, so float32 rounding accumulates more through the forward
                // pass. Observed noise floor is ~1.4e-3 absolute across all 9 components on this
                // fixture; small-magnitude force components can look like a large *relative*
                // miss purely from that floor, so this checks absolute agreement.
                assert!(
                    (finite_diff_force - model_component).abs() < 5e-3,
                    "atom {atom} axis {axis}: finite-diff {finite_diff_force} vs model {model_component}"
                );
            }
        }
    }
}
