# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Defines (mathematical) distributions."""

# -------
# Imports
# -------

import numpy as np

from .constants import BOLTZMANN_CONSTANT_EV, PLANCK_CONSTANT_EV
from .utility.numpy_helper import np_check_shape, np_expand_dims

try:
    from numba import njit
except ImportError:
    from .utility.numba_helper import dummy_njit as njit

# ----------
# Lineshapes
# ----------


def gaussian(x, i, mu, sigma):
    r"""Evaluate the Gaussian lineshape with the supplied intensity
    (peak area) `i`, mean `mu` (centre) and standard deviation (width)
    `sigma`.

    Parameters
    ----------
    x : array_like
        Values at which to evaluate the function.
    i, mu, sigma : float
        Intensity, mean and standard deviation.

    Returns
    -------
    g_x : numpy.ndarray
        Gaussian function evaluated at `x`.

    Notes
    -----
    The definition of :math:`G(x, I, \mu, \sigma)` is taken from:
    http://mathworld.wolfram.com/GaussianFunction.html.

    .. math::

        G(x) = \frac{I}{\sqrt{2\pi} \sigma} \exp{-\frac{(x - \mu)^2}{2 \sigma^2}}
    """

    return (i / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(
        -1.0 * (np.asarray(x) - mu) ** 2 / (2 * sigma**2)
    )


def lorentzian(x, i, x0, gamma):
    r"""Evaluate the Lorentzian lineshape with intensity (peak area)
    `i`, central value `x0` and width `gamma`.

    Parameters
    ----------
    x : array_like
        Values at which to evaluate the function.
    i, x0, gamma : float
        Intensity, central value and width.

    Returns
    -------
    l_x : numpy.ndarray
        Lorentzian function evaluated at `x`.

    Notes
    -----
    The definition of :math:`L(x, I, x_0, \Gamma)` is taken from:
    http://mathworld.wolfram.com/LorentzianFunction.html.

    .. math::

        L(x) = \frac{I}{\pi} \frac{\frac{1}{2} \Gamma}{(x - x_0)^2 + (\frac{1}{2} \Gamma)^2}
    """

    return (i / np.pi) * (
        (0.5 * gamma) / ((np.asarray(x) - x0) ** 2 + (0.5 * gamma) ** 2)
    )


def lorentz_oscillator(omega, s, omega_0, eta):
    r"""Evaluate the complex dielectric function for a phonon mode with
    intensity (oscillator strength) `s`, central value `x_0` and width
    ("complex shift") `eta` using the Lorentz oscillator model.

    Parameters
    ----------
    omega : array_like
        Frequencies at which to evaluate the function.
    s : array_like or float
        Tensor (shape: `(3, 3)`) or scalar intensity.
    omega_0, eta : float
        Central value and linewidth.

    Returns
    -------
    dielectric_func : numpy.ndarray
        Complex dielectric function evaluated at `x` (shape `(N,)` for
        scalar `s`, or `(N, 3, 3)` for tensor `s`.

    Notes
    -----
    The definition of the dielectric function is that used in the VASP
    code:
    https://www.vasp.at/wiki/index.php/Category:Dielectric_properties

    .. math::

        f(\omega) = \frac{s}{ \omega_0^2 - (\omega + i\eta)^2 }
    """

    # Denominator.

    denom = omega_0**2 - (np.asarray(omega) + 1.0j * eta) ** 2

    if np.ndim(s) == 0:
        # Scalar oscillator strength.

        return s / denom

    if not np_check_shape(s, (3, 3)):
        # Tensor oscillator strength.

        raise ValueError(
            "s must be a scalar or an array_like with shape (3, 3)."
        )

    dielectric_func = np.zeros((len(omega), 3, 3), dtype=np.complex128)

    dielectric_func += s[np.newaxis, :, :]
    dielectric_func /= denom[:, np.newaxis, np.newaxis]

    return dielectric_func


# -------------
# Bose-Einstein
# -------------


def phonon_occupation_number(nu, t):
    r"""Evaluate the phonon occupation number for frequencies `nu` and
    temperature `t` using the Bose-Einstein distribution.

    Parameters
    ----------
    nu : float or array_like
        Phonon frequency or frequencies in THz.
    t : float
        Temperature in K.

    Returns
    -------
    n : float or array_like
        Phonon occupation number(s) (same shape as `nu`).

    Notes
    -----
    The phonon occupation number :math:`n(\nu, T)` is given by:

    .. math::

        n(\nu, T) = \frac{1}{\exp{[h \nu / k_\mathrm{B} T]} - 1}
    """

    nu, n_dim_add = np_expand_dims(np.asarray(nu), (None,))

    e = PLANCK_CONSTANT_EV * 1.0e12 * nu
    n = 1.0 / (np.exp(e / (BOLTZMANN_CONSTANT_EV * t)) - 1.0)

    return n if n_dim_add == 0 else n[0]


# ---------------------
# Preferred orientation
# ---------------------


@njit
def march_dollase(alpha, r):
    r"""Evaluate the March-Dollase distribution function with the
    supplied angle `alpha` and March parameter `r`.

    Parameters
    ----------
    alpha : float
        Angle in radians.
    r : float
        March parameter.

    Returns
    -------
    v : float
        Value of the March-Dollase function.

    Notes
    -----
    The March-Dollase distribution function is given by:

    .. math::

        f(\alpha, r) = \left[ r^2 \cos^2 \alpha + \frac{1}{r} \sin^2 \alpha \right]^{-3/2}
    """

    return (r**2 * np.cos(alpha) ** 2 + (1.0 / r) * np.sin(alpha) ** 2) ** (
        -3.0 / 2.0
    )


def march_dollase_eta_to_r(eta):
    r"""Convert a crystallte excess fraction `eta` to the corresponding
     March parameter `r` in the March-Dollase distribution function.

    Parameters
    ----------
    eta : float
        Crystallite fraction.

    Returns
    -------
    r : float
        March parameter.

    Notes
    -----
    The conversion between `eta` and `r` is given by:

    ..math::

        \frac{1}{\eta^2 - 1} \left[ -\frac{1}{2} \eta^2 + \sqrt{3} \times \eta \times \sqrt{1 - \frac{1}{4}\eta} - 1 \right]
    """

    if eta < 0.0 or eta >= 1.0:
        raise ValueError(
            "eta must be >= 0 and < 1. eta = 1 causes the "
            "March-Dollase function to diverge."
        )

    return (
        (-0.5 * eta**2)
        + (np.sqrt(3.0) * eta * np.sqrt(1.0 - 0.25 * eta**2))
        - 1.0
    ) / (eta**2 - 1.0)
