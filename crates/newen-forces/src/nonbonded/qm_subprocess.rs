//! Shared subprocess mechanics for file-in/file-out QM engine providers ([`super::xtb`],
//! [`super::mopac`]). Crate-private — not a new public trait. Each engine's actual physics (input
//! format, output parsing, units) stays engine-specific in its own module; only the mechanical
//! parts (work-dir setup, stale-output removal, spawn, exit-status check) are common enough to
//! share. See `architecture.md` rule 3 ("a small, stable trait taxonomy — resist proliferation")
//! for why this is a helper module and not a `QmEngine` trait mirroring gromosXX's `QM_Worker`
//! base class.

use crate::provider::ProviderError;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Create `work_dir` if missing (idempotent — the directory is reused, overwritten, on every
/// call, no temp-directory dependency).
pub(crate) fn ensure_work_dir(work_dir: impl Into<PathBuf>) -> Result<PathBuf, ProviderError> {
    let work_dir = work_dir.into();
    std::fs::create_dir_all(&work_dir)
        .map_err(|e| ProviderError::ComputationFailed(format!("failed to create work_dir: {e}")))?;
    Ok(work_dir)
}

/// Remove a possibly-stale output file before running, so a crashed engine can't leave a prior
/// result silently readable as this call's answer. Best-effort: absence is not an error.
pub(crate) fn remove_stale(path: &Path) {
    let _ = std::fs::remove_file(path);
}

/// Run `binary` with `args` in `work_dir`; a spawn failure or nonzero exit status becomes a
/// `ProviderError` carrying stderr.
/// `run_subprocess` with a wall-clock bound: the child is killed once `timeout` elapses and the
/// call fails with a clear error instead of hanging the MD step (a QM engine that does not
/// converge, or a lost license server, must not stall a run). Its stdout/stderr go to
/// `subprocess.stdout`/`subprocess.stderr` in `work_dir` (a piped read cannot be bounded).
pub(crate) fn run_subprocess_with_timeout(
    binary: &str,
    work_dir: &Path,
    args: &[&str],
    timeout: Option<std::time::Duration>,
) -> Result<(), ProviderError> {
    let Some(timeout) = timeout else {
        return run_subprocess(binary, work_dir, args);
    };
    let stdout_path = work_dir.join("subprocess.stdout");
    let stderr_path = work_dir.join("subprocess.stderr");
    let open = |path: &Path| {
        std::fs::File::create(path).map_err(|e| {
            ProviderError::ComputationFailed(format!("cannot create {}: {e}", path.display()))
        })
    };
    let mut child = Command::new(binary)
        .current_dir(work_dir)
        .args(args)
        .stdout(std::process::Stdio::from(open(&stdout_path)?))
        .stderr(std::process::Stdio::from(open(&stderr_path)?))
        .spawn()
        .map_err(|e| ProviderError::ComputationFailed(format!("failed to run {binary}: {e}")))?;
    let started = std::time::Instant::now();
    loop {
        match child.try_wait() {
            Ok(Some(status)) => {
                if !status.success() {
                    let stderr = std::fs::read_to_string(&stderr_path).unwrap_or_default();
                    return Err(ProviderError::ComputationFailed(format!(
                        "{binary} exited with {status}: {stderr}"
                    )));
                }
                return Ok(());
            },
            Ok(None) => {
                if started.elapsed() > timeout {
                    let _ = child.kill();
                    let _ = child.wait();
                    return Err(ProviderError::ComputationFailed(format!(
                        "{binary} exceeded the {:.0} s timeout and was killed",
                        timeout.as_secs_f64()
                    )));
                }
                std::thread::sleep(std::time::Duration::from_millis(20));
            },
            Err(e) => {
                return Err(ProviderError::ComputationFailed(format!(
                    "failed to wait for {binary}: {e}"
                )))
            },
        }
    }
}

pub(crate) fn run_subprocess(
    binary: &str,
    work_dir: &Path,
    args: &[&str],
) -> Result<(), ProviderError> {
    let output = Command::new(binary)
        .current_dir(work_dir)
        .args(args)
        .output()
        .map_err(|e| ProviderError::ComputationFailed(format!("failed to run {binary}: {e}")))?;
    if !output.status.success() {
        return Err(ProviderError::ComputationFailed(format!(
            "{binary} exited with {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        )));
    }
    Ok(())
}
