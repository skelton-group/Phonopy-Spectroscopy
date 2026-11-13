# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Implements mathematical distributions."""

# -------
# Imports
# -------

import warnings

import numpy as np

from .constants import (
    BOLTZMANN_CONSTANT_EV,
    PLANCK_CONSTANT_EV,
    ZERO_TOLERANCE,
)

try:
    from numba import njit
except ImportError:
    from .utility.numba_helper import dummy_njit as njit

# ------------------
# Spectral functions
# ------------------


def gaussian(x, i, mu, sigma):
    r"""Evaluate the Gaussian lineshape with the supplied intensity
    (peak area) `i`, mean `mu` (centre) and standard deviation (width)
    `sigma`.

    Parameters
    ----------
    x : array_like
        Values at which to evaluate the function (shape: `(O,)`).
    i, mu, sigma : scalar or array_like
        Intensity, mean and standard deviation (shape: `(...)`).

    Returns
    -------
    g_x : numpy.ndarray
        Gaussian function evaluated at `x` (shape: `(..., O)`).

    Notes
    -----
    The definition of :math:`G(x, I, \mu, \sigma)` is taken from:
    http://mathworld.wolfram.com/GaussianFunction.html.

    .. math::

        G(x) = \frac{I}{\sqrt{2\pi} \sigma} \exp{-\frac{(x - \mu)^2}{2 \sigma^2}}
    """

    # Ensure x is a 1D array with shape (O,).

    x = np.atleast_1d(np.asanyarray(x))

    if x.ndim != 1:
        raise ValueError("x must be an array_like with shape (O,).")

    # Broadcast i, mu and sigma to a common shape (...).

    i = np.asanuarray(i)
    mu = np.asanyarray(mu)
    sigma = np.asanyarray(sigma)

    i, mu, sigma = np.broadcast_arrays(i, mu, sigma)

    # Append leading "batch" dimensions to x and trailing "frequency"
    # dimensions to i, mu and sigma.

    x = x.reshape((1,) * i.ndim + (len(x),))

    i = np.expand_dims(i, -1)
    mu = np.expand_dims(mu, -1)
    sigma = np.expand_dims(sigma, -1)

    # x and i/mu/sigma should now broadcast to (..., O).

    return (i / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(
        -1.0 * (x - mu) ** 2 / (2 * sigma**2)
    )


def lorentzian(x, i, x_0, gamma):
    r"""Evaluate the Lorentzian lineshape with intensity (peak area)
    `i`, central value `x0` and width `gamma`.

    Parameters
    ----------
    x : array_like
        Values at which to evaluate the function (shape: `(O,)`).
    i, x_0, gamma : scalar or array_like
        Intensities, central values and widths (shape: `(...)`).

    Returns
    -------
    l_x : numpy.ndarray
        Lorentzian function evaluated at `x` (shape: `(..., O)`).

    Notes
    -----
    The definition of :math:`L(x, I, x_0, \Gamma)` is taken from:
    http://mathworld.wolfram.com/LorentzianFunction.html.

    .. math::

        L(x) = \frac{I}{\pi} \frac{\frac{1}{2} \Gamma}{(x - x_0)^2 + (\frac{1}{2} \Gamma)^2}
    """

    x = np.atleast_1d(np.asanyarray(x))

    if x.ndim != 1:
        raise ValueError("x must be an array_like with shape (O,).")

    i = np.asanyarray(i)
    x_0 = np.asanyarray(x_0)
    gamma = np.asanyarray(gamma)

    i, x_0, gamma = np.broadcast_arrays(i, x_0, gamma)

    x = x.reshape((1,) * x_0.ndim + (len(x),))

    i = np.expand_dims(i, -1)
    x_0 = np.expand_dims(x_0, -1)
    gamma = np.expand_dims(gamma, -1)

    return (i / np.pi) * (
        (0.5 * gamma) / ((x - x_0) ** 2 + (0.5 * gamma) ** 2)
    )


def lorentz_oscillator(x, s, x_0, eta):
    r"""Evaluate the complex dielectric function for a phonon mode with
    intensity (oscillator strength) `s`, central value `x_0` and width
    ("complex shift") `eta` using the Lorentz oscillator model.

    Parameters
    ----------
    x : array_like
        Frequencies at which to evaluate the function (shape: `(O,)`).
    s : array_like
        Oscillator strength (shape: `(..., 3, 3)`)`.
    x_0, eta : scalar or array_like
        Central frequencies/linewidths (shape: `(...)`).

    Returns
    -------
    eps : numpy.ndarray
        Complex dielectric function evaluated at `x` (shape:
        `(..., O, 3, 3)`).

    Notes
    -----
    The definition of the dielectric function is that used in the VASP
    code:
    https://www.vasp.at/wiki/index.php/Category:Dielectric_properties

    .. math::

        f(\omega) = \frac{s}{ \omega_0^2 - (\omega + i\eta)^2 }
    """

    x = np.atleast_1d(np.asanyarray(x))

    if x.ndim != 1:
        raise ValueError("x must be scalar or array_like with shape (O,).")

    # Broadcast x_0 and eta to a common shape (...).

    x_0, eta = np.broadcast_arrays(x_0, eta)

    s = np.asanyarray(s)

    if s.ndim < 2 or s.shape[-2:] != (3, 3):
        raise ValueError("s must be an array_like with shape (..., 3, 3).")

    # Broadcast s to (..., 3, 3).

    s = np.broadcast_to(s, x_0.shape + (3, 3))

    # Append leading batch dimensions to x.

    x = x.reshape((1,) * x_0.ndim + (len(x),))

    # Append trailing frequency dimensions to x_0 and eta.

    x_0 = np.expand_dims(x_0, -1)
    eta = np.expand_dims(eta, -1)

    # Insert a frequency dimension into s.

    s = np.expand_dims(s, -3)

    denom = x_0**2 - (x + 1.0j * eta) ** 2

    return s / denom.reshape(denom.shape + (1, 1))


# --------------------------
# Bose-Einstein distribution
# --------------------------


def phonon_occupation_number(f, t):
    r"""Evaluate the phonon occupation number for frequencies `f` and
    temperature `t` using the Bose-Einstein distribution.

    Parameters
    ----------
    f : scalar or array_like
        Phonon frequency/frequencies in THz (shape: `(...)`).
    t : float
        Temperature in K.

    Returns
    -------
    n : float or array_like
        Phonon occupation number(s) (same shape as `f`).

    Notes
    -----
    The phonon occupation number :math:`n(\nu, T)` is given by:

    .. math::

        n(f, T) = \frac{1}{\exp{[h f / k_\mathrm{B} T]} - 1}
    """

    # The internal convention is that imaginary frequencies are
    # represented by negative (real) values. Complex frequencies
    # could be generated by some parts of the library, but if they
    # are used here it likely indicates a bug elsewhere.

    f = np.real_if_close(np.asanyarray(f), tol=ZERO_TOLERANCE)

    if np.iscomplexobj(f):
        raise ValueError("f must be real.")

    # This is standard practice, but we should warn about it.

    if (f < 0.0).any():
        warnings.warn(
            "Negative frequencies will be converted to absolute values.",
            RuntimeWarning,
        )

        f = np.abs(f)

    e = PLANCK_CONSTANT_EV * 1.0e12 * f
    n = 1.0 / (np.exp(e / (BOLTZMANN_CONSTANT_EV * t)) - 1.0)

    return n


# ---------------------
# Preferred orientation
# ---------------------


@njit(fastmath=True, inline="always")
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
    f : numpy.ndarray
        Value of the March-Dollase function.

    Notes
    -----
    The March-Dollase distribution function is given by:

    .. math::

        f(\alpha, r) = \left[ r^2 \cos^2 \alpha + \frac{1}{r} \sin^2 \alpha \right]^{-3/2}
    """

    # Numba fails to compile if we use convert r to a NumPy array and
    # broadcast the test over elements.

    if r < 0.0:
        raise ValueError("r must be > 0.")

    # Using identity cos^2(x) + sin^2(x) = 1.

    c2_a = np.cos(alpha) ** 2
    s2_a = 1.0 - c2_a

    return (r**2 * c2_a + (1.0 / r) * s2_a) ** (-3.0 / 2.0)


def march_dollase_r_to_eta(r):
    r"""Calculate the excess fraction of crystallites \eta in a
    preferred orientation for a given value of the March parameter r.

    Parameters
    ----------
    r : any
        March parameter.

    Returns
    -------
    eta : numpy.ndarray
        Excess crystallite fraction (same shape as `r`).
    """

    r = np.asanyarray(r)

    if (r < ZERO_TOLERANCE).any():
        raise ValueError("r must be > 0.")

    if np.abs(r - 1.0) < ZERO_TOLERANCE:
        return 0.0

    return np.sqrt((1.0 - r) ** 3 / (1.0 - r**3))


def march_dollase_eta_to_r(eta):
    r"""Calculate the two values of the March parameter r that give the
    fraction of crystallites \eta in a preferred orientation.

    Parameters
    ----------
    eta : any
        Excess crystallite fraction.

    Returns
    -------
    r_vals : tuple of numpy.ndarray
        March parameter r with r < 1 and r > 1 that give \eta (both
        same shape as `eta`).
    """

    eta = np.asanyarray(eta)

    if (eta < 0.0).any() or (eta > 1.0).any():
        raise ValueError("eta must be between 0 and 1.")

    r_1 = (
        (-(eta**2) / 2.0) + (eta * np.sqrt(12.0 - 3.0 * eta**2) / 2.0) - 1.0
    ) / (eta**2 - 1.0)

    r_2 = (
        (-1.0 * eta**2 / 2.0)
        - (eta * np.sqrt(12.0 - 3.0 * eta**2) / 2.0)
        - 1.0
    ) / (eta**2 - 1.0)

    return (r_1, r_2)
