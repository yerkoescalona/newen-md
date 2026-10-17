//! Algorithm trait and sequence - the GROMOS Algorithm_Sequence pattern
//!
//! In GROMOS, the MD loop is driven by an ordered sequence of algorithms.
//! Each algorithm implements `apply()` which operates on the shared state
//! (topology, configuration, simulation parameters). The main loop simply
//! iterates over the sequence each step.
//!
//! This module provides the Rust equivalent using trait objects.

use crate::configuration::Configuration;
use crate::topology::Topology;

/// Simulation parameters passed to each algorithm per step.
///
/// Equivalent to GROMOS's `simulation::Simulation` - holds time step,
/// current step number, and other per-step state.
#[derive(Debug, Clone)]
pub struct SimulationState {
    /// Current time step size (ps)
    pub dt: f64,
    /// Current step number
    pub step: usize,
    /// Total number of steps
    pub n_steps: usize,
    /// Current simulation time (ps)
    pub time: f64,
}

impl SimulationState {
    /// Create a new state with the given timestep (ps) and total number of steps.
    pub fn new(dt: f64, n_steps: usize) -> Self {
        Self {
            dt,
            step: 0,
            n_steps,
            time: 0.0,
        }
    }

    /// Advance by one step, updating step counter and simulation time.
    pub fn advance(&mut self) {
        self.step += 1;
        self.time = self.step as f64 * self.dt;
    }
}

/// The core Algorithm trait - equivalent to GROMOS's `algorithm::Algorithm`.
///
/// Each algorithm in the MD sequence implements this trait.
/// The `apply()` method is called once per step in the order defined
/// by the `AlgorithmSequence`.
pub trait Algorithm {
    /// Apply this algorithm for one step.
    ///
    /// Returns Ok(()) on success, Err with description on failure.
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String>;

    /// Initialize the algorithm (called once before the main loop).
    fn init(
        &mut self,
        _topo: &Topology,
        _conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        Ok(())
    }

    /// Algorithm name for logging/debugging.
    fn name(&self) -> &str;

    /// Downcast hook: algorithms that callers may need to configure after the sequence is
    /// built (the force field's MPI hooks) return `Some(self)`. Default: not downcastable.
    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        None
    }

    /// Opt in to accumulating an internal sub-phase breakdown.
    ///
    /// Default: no-op. Algorithms that expose [`Self::sub_timings`] should honour
    /// this, since timing costs two clock reads per phase per step.
    fn enable_timing(&mut self, _enabled: bool) {}

    /// Internal sub-phase breakdown, as `(label, accumulated wall time)`.
    ///
    /// Default: empty. Reported nested under this algorithm's own total, matching
    /// the gromosXX `TIMING` block's sub-entries.
    fn sub_timings(&self) -> Vec<(&'static str, std::time::Duration)> {
        Vec::new()
    }
}

/// An ordered sequence of algorithms that defines the MD step.
///
/// Equivalent to GROMOS's `algorithm::Algorithm_Sequence`.
/// The main MD loop calls `run_step()` which iterates through all algorithms.
pub struct AlgorithmSequence {
    algorithms: Vec<Box<dyn Algorithm>>,
    /// Per-algorithm accumulated wall time, parallel to `algorithms`.
    /// Only written when `timing_enabled`; see [`Self::enable_timing`].
    timings: Vec<std::time::Duration>,
    timing_enabled: bool,
}

impl AlgorithmSequence {
    /// Create an empty algorithm sequence.
    pub fn new() -> Self {
        Self {
            algorithms: Vec::new(),
            timings: Vec::new(),
            timing_enabled: false,
        }
    }

    /// Add an algorithm to the end of the sequence.
    /// The first algorithm of concrete type `T` that opted into [`Algorithm::as_any_mut`].
    pub fn find_mut<T: Algorithm + 'static>(&mut self) -> Option<&mut T> {
        self.algorithms
            .iter_mut()
            .find_map(|a| a.as_any_mut().and_then(|any| any.downcast_mut::<T>()))
    }

    /// Append an algorithm to the end of the sequence.
    pub fn push(&mut self, alg: Box<dyn Algorithm>) {
        self.algorithms.push(alg);
        self.timings.push(std::time::Duration::ZERO);
    }

    /// Insert an algorithm at `index` (0 = first), shifting the rest back.
    ///
    /// For callers that must place an algorithm at a specific point of an already-built
    /// sequence — e.g. a provider-orchestrator term immediately after `Forcefield`, whose
    /// forces it adds to. Panics if `index > len()`, like `Vec::insert`.
    pub fn insert(&mut self, index: usize, alg: Box<dyn Algorithm>) {
        self.algorithms.insert(index, alg);
        self.timings.insert(index, std::time::Duration::ZERO);
    }

    /// Enable per-algorithm wall-clock accumulation.
    ///
    /// Off by default: each timed algorithm costs two `Instant::now()` reads per
    /// step, which is pure overhead in a production run. Enable it to produce the
    /// breakdown reported by [`Self::timings`], equivalent to the gromosXX
    /// `TIMING` block.
    pub fn enable_timing(&mut self, enabled: bool) {
        self.timing_enabled = enabled;
        for alg in &mut self.algorithms {
            alg.enable_timing(enabled);
        }
    }

    /// Accumulated wall time per algorithm, as `(name, duration)` pairs.
    ///
    /// All-zero unless [`Self::enable_timing`] was called.
    pub fn timings(&self) -> Vec<(&str, std::time::Duration)> {
        self.algorithms
            .iter()
            .map(|a| a.name())
            .zip(self.timings.iter().copied())
            .collect()
    }

    /// Per-algorithm totals together with each algorithm's own sub-phase breakdown.
    pub fn detailed_timings(
        &self,
    ) -> Vec<(
        &str,
        std::time::Duration,
        Vec<(&'static str, std::time::Duration)>,
    )> {
        self.algorithms
            .iter()
            .zip(self.timings.iter().copied())
            .map(|(a, d)| (a.name(), d, a.sub_timings()))
            .collect()
    }

    /// Initialize all algorithms in the sequence.
    pub fn init(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        for alg in &mut self.algorithms {
            alg.init(topo, conf, sim)?;
        }
        Ok(())
    }

    /// Run one complete step: apply all algorithms in order.
    pub fn run_step(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        if self.timing_enabled {
            for (alg, slot) in self.algorithms.iter_mut().zip(self.timings.iter_mut()) {
                let t = std::time::Instant::now();
                alg.apply(topo, conf, sim)?;
                *slot += t.elapsed();
            }
        } else {
            for alg in &mut self.algorithms {
                alg.apply(topo, conf, sim)?;
            }
        }
        Ok(())
    }

    /// Number of algorithms in the sequence.
    pub fn len(&self) -> usize {
        self.algorithms.len()
    }

    /// Returns `true` if no algorithms have been added.
    pub fn is_empty(&self) -> bool {
        self.algorithms.is_empty()
    }

    /// Get algorithm names for debugging.
    pub fn algorithm_names(&self) -> Vec<&str> {
        self.algorithms.iter().map(|a| a.name()).collect()
    }
}

impl Default for AlgorithmSequence {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A test algorithm that records when it's called.
    struct RecordingAlgorithm {
        label: &'static str,
        call_log: std::rc::Rc<std::cell::RefCell<Vec<&'static str>>>,
    }

    impl Algorithm for RecordingAlgorithm {
        fn apply(
            &mut self,
            _topo: &Topology,
            _conf: &mut Configuration,
            _sim: &SimulationState,
        ) -> Result<(), String> {
            self.call_log.borrow_mut().push(self.label);
            Ok(())
        }

        fn name(&self) -> &str {
            self.label
        }
    }

    fn make_trivial_system() -> (Topology, Configuration) {
        let topo = Topology::new();
        let conf = Configuration::new(0, 1, 1);
        (topo, conf)
    }

    #[test]
    fn test_sequence_runs_algorithms_in_order() {
        let log = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));

        let mut seq = AlgorithmSequence::new();
        seq.push(Box::new(RecordingAlgorithm {
            label: "Forcefield",
            call_log: log.clone(),
        }));
        seq.push(Box::new(RecordingAlgorithm {
            label: "Velocity",
            call_log: log.clone(),
        }));
        seq.push(Box::new(RecordingAlgorithm {
            label: "Position",
            call_log: log.clone(),
        }));
        seq.push(Box::new(RecordingAlgorithm {
            label: "Energy",
            call_log: log.clone(),
        }));

        let (topo, mut conf) = make_trivial_system();
        let sim = SimulationState::new(0.002, 10);

        seq.run_step(&topo, &mut conf, &sim).unwrap();

        assert_eq!(
            *log.borrow(),
            vec!["Forcefield", "Velocity", "Position", "Energy"]
        );
    }

    #[test]
    fn test_sequence_names() {
        let log = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));

        let mut seq = AlgorithmSequence::new();
        seq.push(Box::new(RecordingAlgorithm {
            label: "A",
            call_log: log.clone(),
        }));
        seq.push(Box::new(RecordingAlgorithm {
            label: "B",
            call_log: log.clone(),
        }));

        assert_eq!(seq.algorithm_names(), vec!["A", "B"]);
        assert_eq!(seq.len(), 2);
        assert!(!seq.is_empty());
    }

    #[test]
    fn test_sequence_propagates_errors() {
        struct FailAlgorithm;
        impl Algorithm for FailAlgorithm {
            fn apply(
                &mut self,
                _: &Topology,
                _: &mut Configuration,
                _: &SimulationState,
            ) -> Result<(), String> {
                Err("deliberate failure".into())
            }
            fn name(&self) -> &str {
                "Fail"
            }
        }

        let mut seq = AlgorithmSequence::new();
        seq.push(Box::new(FailAlgorithm));

        let (topo, mut conf) = make_trivial_system();
        let sim = SimulationState::new(0.002, 10);

        let result = seq.run_step(&topo, &mut conf, &sim);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "deliberate failure");
    }

    #[test]
    fn test_sequence_stops_at_first_error() {
        let log = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));

        struct FailAlgorithm;
        impl Algorithm for FailAlgorithm {
            fn apply(
                &mut self,
                _: &Topology,
                _: &mut Configuration,
                _: &SimulationState,
            ) -> Result<(), String> {
                Err("fail".into())
            }
            fn name(&self) -> &str {
                "Fail"
            }
        }

        let mut seq = AlgorithmSequence::new();
        seq.push(Box::new(RecordingAlgorithm {
            label: "Before",
            call_log: log.clone(),
        }));
        seq.push(Box::new(FailAlgorithm));
        seq.push(Box::new(RecordingAlgorithm {
            label: "After",
            call_log: log.clone(),
        }));

        let (topo, mut conf) = make_trivial_system();
        let sim = SimulationState::new(0.002, 10);

        let _ = seq.run_step(&topo, &mut conf, &sim);
        // "After" should NOT have been called
        assert_eq!(*log.borrow(), vec!["Before"]);
    }

    #[test]
    fn test_simulation_state_advance() {
        let mut sim = SimulationState::new(0.002, 100);
        assert_eq!(sim.step, 0);
        assert_eq!(sim.time, 0.0);

        sim.advance();
        assert_eq!(sim.step, 1);
        assert!((sim.time - 0.002).abs() < 1e-15);

        sim.advance();
        assert_eq!(sim.step, 2);
        assert!((sim.time - 0.004).abs() < 1e-15);
    }

    #[test]
    fn test_empty_sequence() {
        let mut seq = AlgorithmSequence::new();
        assert!(seq.is_empty());
        assert_eq!(seq.len(), 0);

        let (topo, mut conf) = make_trivial_system();
        let sim = SimulationState::new(0.002, 10);
        seq.run_step(&topo, &mut conf, &sim).unwrap();
    }

    #[test]
    fn test_multi_step_sequence() {
        let log = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));

        let mut seq = AlgorithmSequence::new();
        seq.push(Box::new(RecordingAlgorithm {
            label: "V",
            call_log: log.clone(),
        }));
        seq.push(Box::new(RecordingAlgorithm {
            label: "P",
            call_log: log.clone(),
        }));

        let (topo, mut conf) = make_trivial_system();
        let mut sim = SimulationState::new(0.002, 3);

        for _ in 0..3 {
            seq.run_step(&topo, &mut conf, &sim).unwrap();
            sim.advance();
        }

        assert_eq!(*log.borrow(), vec!["V", "P", "V", "P", "V", "P"]);
        assert_eq!(sim.step, 3);
    }
}
