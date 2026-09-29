# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for numerical integration (quadrature)."""

# -------
# Imports
# -------

import numpy as np

from functools import lru_cache

from scipy.integrate import lebedev_rule

from ..constants import ZERO_TOLERANCE

from .geometry import cartesian_to_spherical_polar
from .numpy_helper import np_readonly_view

# -----------
# Unit circle
# -----------


def unit_circle_quad_rule(m, angles=True):
    r"""Generate a list of unit vectors and weights for integrating over
    the unit circle.

    Parameters
    ----------
    m : int
        Number of points for (order of) the quadrature rule.
    angles : bool, optional
        If `True`, return angles \psi in Radians, otherwise return 2D
        (x, y) vectors (default: `True`).

    Returns
    -------
    circle_rule : tuple of numpy.ndarray
        Tuple of `(x, w)` with shape `(N,)` (`angles=True`) or `(N, 2)`
        (`angles=False`) and `(N,)`.

    Notes
    -----
    The integral of :math:`f(\psi)` is given by:

    .. math::

        \int_0^{2\pi} f(\psi) d\psi = 2 \pi \sum_m w_m f(\psi_m)
    """

    x = np.linspace(0.0, 2.0 * np.pi, m, endpoint=False, dtype=np.float64)
    w = np.ones(m, dtype=np.float64) / m

    if angles:
        return (x, w)
    else:
        vecs = np.zeros((m, 2), dtype=np.float64)

        vecs[:, 0] = np.cos(x)
        vecs[:, 1] = np.sin(x)

        return (vecs, w)


# -----------
# Unit sphere
# -----------

_UNIT_SPHERE_LEBEDEV_QUAD_AVAILABLE_PRECS = np.array(
    [
        3,
        5,
        7,
        9,
        11,
        13,
        15,
        17,
        19,
        21,
        23,
        25,
        27,
        29,
        31,
        35,
        41,
        47,
        53,
        59,
        65,
        71,
        77,
        83,
        89,
        95,
        101,
        107,
        113,
        119,
        125,
        131,
    ],
    dtype=int,
)


def unit_sphere_lebedev_quad_rule_available_precs():
    """Return a list available precisions (orders) for Levedev
    quadrature rules.

    Returns
    -------
    p : list of int
        List of available precisions `p` (orders) that can be passed to
        `unit_sphere_lebedev_quad_rule`.

    See Also
    --------
    unit_sphere_lebedev_quad_rule : Return the unit vectors and weights
        for integrating over the unit sphere with Lebedev quadrature.
    """

    return np_readonly_view(_UNIT_SPHERE_LEBEDEV_QUAD_AVAILABLE_PRECS)


@lru_cache(maxsize=64)
def unit_sphere_lebedev_quad_rule(p, angles=True):
    r"""Return a set of polar angles (or vectors) and weights for
    integrating over the unit sphere using Lebedev quadrature.

    Parameters
    ----------
    p : int
        Precision (order) of the Lebedev rule.
    ret_polar : bool, optional
        If `True`, return polar angles instead of vectors (default:
        `False`).

    Returns
    -------
    lebedev_rule : tuple of numpy.ndarray
        Tuple of `(x, w)` (shape: `((N, 2), (N,))` if `ret_polar=True`,
        or `((N, 3), (N,))` if `ret_polar=False`.).

    See Also
    --------
    unit_sphere_lebedev_quad_available_prec :
        List of available precisions `p` for which Lebedev rules are
        available.

    Notes
    -----
    The integral of :math:`f(\phi, \theta)` is given by:

    .. math::

        \int_0^{2\pi} \int_0^{\pi} f(\phi, \theta) \sin \theta \; d\theta \; d\phi = 4 \pi \sum_n w_n f(\phi_n, \theta_n)
    """

    if p not in _UNIT_SPHERE_LEBEDEV_QUAD_AVAILABLE_PRECS:
        raise ValueError(
            "p = {0} is not valid for Lebedev quadrature, or data for "
            "this quadrature rule is not available.".format(p)
        )

    x, w = lebedev_rule(p)

    # (3, N) -> (N, 3).

    x = np.swapaxes(x, 1, 0)

    # The weights from the lebedev_rule are multipled by 4 \pi.

    w /= 4.0 * np.pi

    if angles:
        # Convert vectors to spherical polar coordinates.

        sp = cartesian_to_spherical_polar(x)

        if (np.abs(sp[:, 0] - 1.0) > ZERO_TOLERANCE).any():
            raise RuntimeError(
                "One or more vectors returned by lebedev_rule() has "
                "non-unit length. This is most likely a bug."
            )

        return (sp[:, 1:], w)

    return (x, w)


# -------------
# Product grids
# -------------


@lru_cache(maxsize=32)
def circle_circle_euler_angle_quad_rule(m):
    r"""Return a set of \phi and \psi anfles for integrating over the
    \phi and \psi Euler angles with a circle product rule.

    Parameters
    ----------
    m : int
        Number of points in the circle quadrature rule.

    Returns
    -------
    circ_circ_rule : tuple of numpy.ndarray
        Tuple of `(x, w)` (shape: `((N, 2), (N,))`).
    """

    a, _ = unit_circle_quad_rule(m, angles=True)

    x = np.zeros((m * m, 2), dtype=np.float64)

    x[:, 0] = np.repeat(a, m)
    x[:, 1] = np.tile(a, m)

    w = np.zeros((m * m,), dtype=np.float64) / (m * m)

    return (x, w)


@lru_cache(maxsize=32)
def lebedev_circle_euler_angle_quad_rule(p, m=None):
    r"""Return a set of Euler angles and quadrature weights for
    integrating over Euler angles with a Lebedev + circle product grid.


    Parameters
    ----------
    p : int
        Precision (order) of the Lebedev quadrature rule used to
        generate `\phi` and `\theta`.
    m : int, optional
        Number of points in the circle quadrature rule used to generate
        `\psi` (default: automatically chosen to match the number of
        unique `\phi` in the Lebedev rule).

    Returns
    -------
    leb_circ_rule : tuple of numpy.ndarray
        Tuple of `(x, w)` (shape: `((N, 3), (N,))`).

    See Also
    --------
    unit_sphere_lebedev_quad_available_prec :
        List of available precisions `p` for which Lebedev quadrature
        rules are available.

    Notes
    -----
    The integral of :math:`f(\theta, \phi, \psi)` is computed as:

    .. math::

        \int_0^{2\pi} \int_0^{\pi} \int_0^{2\pi} f(\phi, \theta, \psi) \sin \theta \; d\psi \; d\theta \; d\phi = 8 \pi^2 \sum_n \sum_m w_n w_m f(\phi_n, \theta_n, \psi_m)
    """

    a_n, w_n = unit_sphere_lebedev_quad_rule(p, angles=True)

    n = len(a_n)

    if m is None:
        m = len(np.unique(a_n[:, 0]))

    a_m, w_m = unit_circle_quad_rule(m, angles=True)

    x = np.zeros((n * m, 3), dtype=np.float64)

    x[:, :2] = np.repeat(a_n, m, axis=0)
    x[:, 2] = np.tile(a_m, n)

    w = np.repeat(w_n, m) * np.tile(w_m, n)

    return (x, w)
