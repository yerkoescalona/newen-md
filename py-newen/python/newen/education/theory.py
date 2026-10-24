"""The textbook results the experiments compare against, as plain functions.

Units are GROMOS units throughout: nm, ps, u, kJ/mol, K. The Boltzmann constant is
:data:`newen.education.systems.K_B`. Every function accepts NumPy arrays where a quantity can
vary.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from .systems import K_B

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

__all__ = [
    "canonical_kinetic_energy_pdf",
    "canonical_kinetic_energy_relative_width",
    "harmonic_step_limit",
    "lj_energy",
    "lj_minimum",
    "maxwell_boltzmann_mean_speed",
    "maxwell_boltzmann_speed_pdf",
    "temperature_from_kinetic",
]


def _xlogy(x: float, y: np.ndarray) -> np.ndarray:
    """:math:`x \\log y`, taken as 0 when :math:`x = 0` (so :math:`y^0 = 1` at :math:`y = 0`)."""
    if x == 0:
        return np.zeros_like(y)
    with np.errstate(divide="ignore"):
        return x * np.log(y)


def lj_energy(r: ArrayLike, c6: float, c12: float) -> np.ndarray:
    """Lennard-Jones pair energy in the GROMOS C6/C12 form.

    .. math::

        V(r) = \\frac{C_{12}}{r^{12}} - \\frac{C_6}{r^6}

    Parameters
    ----------
    r : array_like
        Pair distance (nm).
    c6, c12 : float
        Dispersion (kJ mol⁻¹ nm⁶) and repulsion (kJ mol⁻¹ nm¹²) coefficients.

    Returns
    -------
    numpy.ndarray
        Energy (kJ/mol).
    """
    r = np.asarray(r, dtype=float)
    return c12 / r**12 - c6 / r**6


def lj_minimum(c6: float, c12: float) -> tuple[float, float]:
    """Position and depth of the Lennard-Jones minimum.

    .. math::

        r_{\\min} = \\left(\\frac{2 C_{12}}{C_6}\\right)^{1/6}, \\qquad
        \\varepsilon = -V(r_{\\min}) = \\frac{C_6^2}{4 C_{12}}

    Parameters
    ----------
    c6, c12 : float
        Dispersion and repulsion coefficients.

    Returns
    -------
    r_min : float
        Distance of the minimum (nm), :math:`2^{1/6}\\sigma`.
    epsilon : float
        Well depth (kJ/mol), positive.
    """
    return (2.0 * c12 / c6) ** (1.0 / 6.0), c6**2 / (4.0 * c12)


def maxwell_boltzmann_speed_pdf(v: ArrayLike, T: float, m: float, dim: int = 3) -> np.ndarray:
    """Probability density of the speed of a classical particle at temperature ``T``.

    With :math:`a = m / k_B T` and :math:`d` dimensions,

    .. math::

        f(v) = \\frac{2}{\\Gamma(d/2)} \\left(\\frac{a}{2}\\right)^{d/2}
               v^{d-1} e^{-a v^2 / 2},

    which is :math:`a v e^{-a v^2/2}` in two dimensions and
    :math:`\\sqrt{2/\\pi}\\, a^{3/2} v^2 e^{-a v^2/2}` in three.

    Parameters
    ----------
    v : array_like
        Speed (nm/ps), non-negative.
    T : float
        Temperature (K).
    m : float
        Mass (u).
    dim : int
        Number of dimensions the particle moves in.

    Returns
    -------
    numpy.ndarray
        Density (ps/nm); integrates to 1 over :math:`v \\ge 0`.
    """
    v = np.asarray(v, dtype=float)
    a = m / (K_B * T)
    log_norm = math.log(2.0) - math.lgamma(dim / 2.0) + (dim / 2.0) * math.log(a / 2.0)
    return np.exp(log_norm + _xlogy(dim - 1, v) - a * v**2 / 2.0)


def maxwell_boltzmann_mean_speed(T: float, m: float, dim: int = 3) -> float:
    """Mean speed of the Maxwell-Boltzmann distribution.

    .. math::

        \\langle v \\rangle = \\sqrt{\\frac{2 k_B T}{m}}\\,
                              \\frac{\\Gamma\\left(\\frac{d+1}{2}\\right)}{\\Gamma(d/2)},

    which is :math:`\\sqrt{\\pi k_B T / 2m}` in two dimensions and
    :math:`\\sqrt{8 k_B T / \\pi m}` in three.

    Parameters
    ----------
    T : float
        Temperature (K).
    m : float
        Mass (u).
    dim : int
        Number of dimensions.

    Returns
    -------
    float
        Mean speed (nm/ps).
    """
    ratio = math.exp(math.lgamma((dim + 1) / 2.0) - math.lgamma(dim / 2.0))
    return math.sqrt(2.0 * K_B * T / m) * ratio


def canonical_kinetic_energy_pdf(ke: ArrayLike, T: float, dof: int) -> np.ndarray:
    """Probability density of the total kinetic energy in the canonical ensemble.

    The kinetic energy of :math:`N_{df}` quadratic degrees of freedom follows a gamma
    distribution with shape :math:`N_{df}/2` and scale :math:`k_B T`:

    .. math::

        P(K) = \\frac{K^{N_{df}/2 - 1}\\, e^{-K / k_B T}}
                    {\\Gamma(N_{df}/2)\\, (k_B T)^{N_{df}/2}}

    Parameters
    ----------
    ke : array_like
        Kinetic energy (kJ/mol), positive.
    T : float
        Temperature (K).
    dof : int
        Number of degrees of freedom :math:`N_{df}`.

    Returns
    -------
    numpy.ndarray
        Density (mol/kJ).
    """
    ke = np.asarray(ke, dtype=float)
    shape, scale = dof / 2.0, K_B * T
    return np.exp(
        _xlogy(shape - 1.0, ke) - ke / scale - math.lgamma(shape) - shape * math.log(scale)
    )


def canonical_kinetic_energy_relative_width(dof: int) -> float:
    """Relative fluctuation of the canonical kinetic energy.

    .. math::

        \\frac{\\sigma_K}{\\langle K \\rangle} = \\sqrt{\\frac{2}{N_{df}}}

    Parameters
    ----------
    dof : int
        Number of degrees of freedom.

    Returns
    -------
    float
        Standard deviation over mean.
    """
    return math.sqrt(2.0 / dof)


def harmonic_step_limit(period: float) -> float:
    """Largest time step at which leap-frog keeps a harmonic oscillator bounded.

    Leap-frog on :math:`\\ddot{x} = -\\omega^2 x` is stable for :math:`\\omega \\Delta t < 2`:

    .. math::

        \\Delta t_{\\max} = \\frac{2}{\\omega} = \\frac{T}{\\pi}

    Parameters
    ----------
    period : float
        Oscillation period :math:`T`, in any time unit.

    Returns
    -------
    float
        :math:`\\Delta t_{\\max}`, in the unit of ``period``.
    """
    return period / math.pi


def temperature_from_kinetic(ke: ArrayLike, dof: int) -> np.ndarray | float:
    """Instantaneous temperature from the kinetic energy (equipartition).

    .. math::

        T = \\frac{2 K}{N_{df}\\, k_B}

    Parameters
    ----------
    ke : array_like
        Kinetic energy (kJ/mol).
    dof : int
        Number of degrees of freedom.

    Returns
    -------
    float or numpy.ndarray
        Temperature (K); a float for a scalar ``ke``.
    """
    result = 2.0 * np.asarray(ke, dtype=float) / (dof * K_B)
    return float(result) if result.ndim == 0 else result
