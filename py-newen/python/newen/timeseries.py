"""
EnergyTimeseries — a thin wrapper around the numpy array returned by
`Simulation.run()`.

Column layout matches `newen_io::energy::EnergyFrame` (the same order the
`md` binary writes to `.tre`), so a `run()` array and a `.tre` dump of the
same trajectory line up column-for-column:

    [time, kinetic, potential, total, volume, pressure,
     bond, angle, improper, dihedral, lj, coulomb]

Dataframe and plot backends are pluggable (`to_dataframe(backend=...)`,
`plot(*components, backend=...)`). Defaults are set on the module-level
`config` (polars for dataframes, plotly for plots) — override per call or
globally via `newen.timeseries.config.dataframe_backend = "pandas"`.
"""

from dataclasses import dataclass
from typing import Union

import numpy as np

COLUMNS = (
    "time",
    "kinetic",
    "potential",
    "total",
    "volume",
    "pressure",
    "bond",
    "angle",
    "improper",
    "dihedral",
    "lj",
    "coulomb",
    "dhdl",
)

DATAFRAME_BACKENDS = ("polars", "pandas", "dict")
PLOT_BACKENDS = ("plotly", "matplotlib")


@dataclass
class Config:
    """Default backends for `EnergyTimeseries.to_dataframe()` / `.plot()`."""

    dataframe_backend: str = "polars"
    plot_backend: str = "plotly"


config = Config()


class EnergyTimeseries:
    """Wraps the `(n_frames, 13)` array from `Simulation.run()`.

    Examples
    --------
    >>> energies = sim.run(1000, ene_freq=100)
    >>> ts = EnergyTimeseries(energies)
    >>> ts.total
    array([...])
    >>> ts.block_average("total", block_size=5)
    (mean, error)
    >>> ts.plot("bond", "angle")  # plotly by default
    >>> ts.plot("total", backend="matplotlib")
    >>> ts.to_dataframe()  # polars by default
    """

    def __init__(self, array: np.ndarray):
        if array.ndim != 2 or array.shape[1] != len(COLUMNS):
            raise ValueError(
                f"expected an (n_frames, {len(COLUMNS)}) array, got shape {array.shape}"
            )
        self._array = array

    def __len__(self) -> int:
        return int(self._array.shape[0])

    def __getattr__(self, name: str) -> np.ndarray:
        if name in COLUMNS:
            return self._array[:, COLUMNS.index(name)]
        raise AttributeError(name)

    @property
    def array(self) -> np.ndarray:
        """The raw `(n_frames, 13)` numpy array."""
        return self._array

    def to_dataframe(self, backend: str | None = None) -> Union["object", dict[str, np.ndarray]]:
        """Return a dataframe of the timeseries.

        `backend` is one of `"polars"`, `"pandas"`, `"dict"` (default:
        `config.dataframe_backend`, i.e. `"polars"`).
        """
        backend = backend or config.dataframe_backend
        if backend not in DATAFRAME_BACKENDS:
            raise ValueError(
                f"unknown dataframe backend: {backend!r}, expected one of {DATAFRAME_BACKENDS}"
            )

        if backend == "polars":
            try:
                import polars as pl
            except ImportError as e:
                raise ImportError(
                    "polars is required for to_dataframe(backend='polars'); install it with "
                    "`pip install polars` or use backend='pandas'/'dict'"
                ) from e
            return pl.DataFrame(self._array, schema=list(COLUMNS))

        if backend == "pandas":
            try:
                import pandas as pd
            except ImportError as e:
                raise ImportError(
                    "pandas is required for to_dataframe(backend='pandas'); install it with "
                    "`pip install pandas` or use backend='polars'/'dict'"
                ) from e
            frame: object = pd.DataFrame(self._array, columns=list(COLUMNS))
            return frame

        return {col: self._array[:, i] for i, col in enumerate(COLUMNS)}

    def plot(self, *components: str, backend: str | None = None):
        """Plot one or more energy components against time.

        `backend` is one of `"plotly"`, `"matplotlib"` (default:
        `config.plot_backend`, i.e. `"plotly"`). Returns the figure/axes
        object so it renders in a notebook cell as the last expression.
        """
        backend = backend or config.plot_backend
        if backend not in PLOT_BACKENDS:
            raise ValueError(f"unknown plot backend: {backend!r}, expected one of {PLOT_BACKENDS}")

        if not components:
            components = ("total",)
        for component in components:
            if component not in COLUMNS:
                raise ValueError(f"unknown energy component: {component!r}")

        time = self.time
        if backend == "plotly":
            try:
                import plotly.graph_objects as go
            except ImportError as e:
                raise ImportError(
                    "plotly is required for plot(backend='plotly'); install it with "
                    "`pip install plotly` or use backend='matplotlib'"
                ) from e
            fig = go.Figure()
            for component in components:
                fig.add_trace(
                    go.Scatter(
                        x=time,
                        y=self._array[:, COLUMNS.index(component)],
                        mode="lines",
                        name=component,
                    )
                )
            fig.update_layout(xaxis_title="time (ps)", yaxis_title="energy (kJ/mol)")
            return fig

        try:
            import matplotlib.pyplot as plt
        except ImportError as e:
            raise ImportError(
                "matplotlib is required for plot(backend='matplotlib'); install it with "
                "`pip install matplotlib` or use backend='plotly'"
            ) from e
        fig, ax = plt.subplots()
        for component in components:
            ax.plot(time, self._array[:, COLUMNS.index(component)], label=component)
        ax.set_xlabel("time (ps)")
        ax.set_ylabel("energy (kJ/mol)")
        ax.legend()
        return ax

    def block_average(self, component: str, block_size: int) -> tuple[float, float]:
        """Block-average `component`, returning `(mean, standard_error)`.

        The trajectory is split into non-overlapping blocks of `block_size`
        frames (a trailing partial block is dropped); the standard error is
        the standard deviation of block means over sqrt(n_blocks).
        """
        if component not in COLUMNS:
            raise ValueError(f"unknown energy component: {component!r}")
        if block_size < 1:
            raise ValueError("block_size must be >= 1")

        values = self._array[:, COLUMNS.index(component)]
        n_blocks = len(values) // block_size
        if n_blocks < 1:
            raise ValueError(f"not enough frames ({len(values)}) for block_size={block_size}")

        blocks = values[: n_blocks * block_size].reshape(n_blocks, block_size)
        block_means = blocks.mean(axis=1)
        mean = float(block_means.mean())
        error = float(block_means.std(ddof=1) / np.sqrt(n_blocks)) if n_blocks > 1 else 0.0
        return mean, error
