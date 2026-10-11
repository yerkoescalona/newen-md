"""Teaching layer: toy systems, experiments with checks, a live viewer. No physics lives here."""

from __future__ import annotations

import importlib
from typing import Any

from . import experiments, theory
from .experiments import Check, Result
from .systems import argon_pair, diatomic, lj_gas, with_bath

__all__ = [
    "Check",
    "Recording",
    "Result",
    "SphereView",
    "argon_pair",
    "diatomic",
    "experiments",
    "lj_gas",
    "record",
    "theory",
    "view",
    "with_bath",
]

# The viewer needs `anywidget` (the `education` extra); it loads on first use.
_VIEWER = ("Recording", "SphereView", "record", "view")


def __getattr__(name: str) -> Any:
    if name in _VIEWER:
        module = importlib.import_module(f"{__name__}.viewer")
        globals().update({key: getattr(module, key) for key in _VIEWER})
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
