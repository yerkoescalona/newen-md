"""
newen.exceptions — the errors the engine raises, one class per kind.

All subclass the builtin exception the binding raised before they existed, so existing
``except ValueError`` / ``except RuntimeError`` clauses keep working:

- ``RecipeError(ValueError)``       — a recipe, term or algorithm is invalid, or an input
                                      cannot be represented as a run (unknown parameter block,
                                      typo'd field, atom-count mismatch, missing restraint file…)
- ``PlanError(ValueError)``         — an algorithm plan violates the GROMOS step-order invariants
- ``MissingFeatureError(RuntimeError)`` — a term needs a cargo feature this build lacks (``ml``)
- ``RunError(RuntimeError)``        — the engine failed to initialise or to take a step
- ``ProgramError(RuntimeError)``    — a ``newen.tools`` program failed; ``.stdout`` carries what
                                      it had written by then (``check_top`` prints its whole
                                      report before failing)
"""

from .newen import MissingFeatureError, PlanError, ProgramError, RecipeError, RunError

__all__ = [
    "RecipeError",
    "PlanError",
    "MissingFeatureError",
    "RunError",
    "ProgramError",
]
