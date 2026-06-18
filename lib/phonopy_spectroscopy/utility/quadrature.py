# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for numerical integration (quadrature)."""

# -------
# Imports
# -------

import glob
import os

import numpy as np

from itertools import product

from ..constants import ZERO_TOLERANCE

# ---------
# Constants
# ---------

_UNIT_SPHERE_LEBEDEV_QUAD_DATA_DIR = r"_data/lebedev"

_UNIT_SPHERE_LEBEDEV_QUAD_DATA = None

# -----------
# Unit Circle
# -----------


def unit_circle_quad_rule(n, ret="angles"):
    r"""Generate a list of unit vectors and weights for integrating over
    the unit circle.

    Parameters
    ----------
    n : int
        Number of points for the quadrature rule.
    ret : {"angles", "vectors", "vectors_2d", "vectors_3d"}, optional
        Return angles, 2D vectors (`ret="vectors_2d"`) or 3D vectors
        (`ret="vectors"`, `ret="vectors_3d"`) (default: `"angles"`).

    Returns
    -------
    circle_rule : tuple of numpy.ndarray
        Tuple of `(psi, w)` if `ret="angles"` or `(vecs, w)` if
        `ret="vectors"`/`ret="vectors_3d"` (`vecs` shape: `(n, 3)`) or
        `ret="vectors_2d"` (shape: `(n, 2)`).

    Notes
    -----
    The integral of :math:`f(\psi)` is given by:

    .. math::

        \int_0^{2\pi} f(\psi) d\psi = 2 \pi \sum_i w_i f(\psi_i)
    """

    psi = np.linspace(0.0, 2.0 * np.pi, n + 1)[:-1]
    w = np.ones(n) / n

    if ret == "angles":
        return (psi, w)

    if ret.startswith("vectors"):
        vecs = np.zeros((n, 2 if ret == "vectors_2d" else 3), dtype=np.float64)

        vecs[:, 0] = np.cos(psi)
        vecs[:, 1] = np.sin(psi)

        return (vecs, w)

    raise ValueError('Unknown return format "{0}".'.format(ret))


# -----------
# Unit Sphere
# -----------


def _unit_sphere_lebedev_quad_load_data():
    """Lazily initialise `_UNIT_SPHERE_LEBDEV_QUAD_DATA`.

    The data is read from a set of `lebedev_ppp.txt` files stored in
    `_UNIT_SPHERE_LEBDEV_QUAD_DATA_DIR`. Each file specifies the surface
    normals, in spherical polar coordinates, and integration weights for
    a Lebedev quadrature scheme with precision `p`.

    These files are used under the GNU LGPL license from:
        https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html

    (See `LICENSE` file in `_UNIT_SPHERE_LEBDEV_QUAD_DATA_DIR` for
    details.)
    """

    abs_path = os.path.join(
        os.path.split(__file__)[0], _UNIT_SPHERE_LEBEDEV_QUAD_DATA_DIR
    )

    pattern = r"{0}/*.txt".format(abs_path)

    quad_data = {}

    for f in glob.glob(pattern):
        _, tail = os.path.split(f)
        root, _ = os.path.splitext(tail)

        p = int(root.split("_")[-1])

        if p <= 0:
            raise RuntimeError(
                "One or more stored Lebedev quadrature rules have an "
                "invalid precision (this is likely a bug)."
            )

        temp = [[], [], []]

        with open(f, "r") as input_reader:
            for line in input_reader:
                vals = [float(val) for val in line.strip().split()]

                for col, val in zip(temp, vals):
                    col.append(val)

        phi, theta, w = [np.array(item) for item in temp]

        # Quadrature weights should sum to 1.

        if np.abs(1.0 - np.sum(w)) > ZERO_TOLERANCE:
            raise RuntimeError(
                "One or more stored Lebedev quadrature rules have "
                "invalid integration weights (this is likely a bug)."
            )

        # Convert theta and phi to radians.

        phi, theta = np.radians(phi), np.radians(theta)

        # The "physics" convention for spherical polar coordinates
        # differs from what is used in the data files.

        quad_data[p] = (phi, theta, w)

    if len(quad_data) == 0:
        raise RuntimeError(
            "No stored Lebedev quadrature rules found (this is likely a bug)."
        )

    global _UNIT_SPHERE_LEBEDEV_QUAD_DATA
    _UNIT_SPHERE_LEBEDEV_QUAD_DATA = quad_data


def unit_sphere_lebedev_quad_available_precs():
    """Return a list available precisions for Levedev quadrature rules.

    Returns
    -------
    precs : list of int
        List of available precisions `p` that can be passed to
        `unit_sphere_lebedev_quad_rule`.

    See Also
    --------
    unit_sphere_lebedev_quad_rule : Return the unit vectors and weights
        for integrating over the unit sphere with Lebedev quadrature.
    """

    if _UNIT_SPHERE_LEBEDEV_QUAD_DATA is None:
        _unit_sphere_lebedev_quad_load_data()

    return sorted(_UNIT_SPHERE_LEBEDEV_QUAD_DATA.keys())


def unit_sphere_lebedev_quad_rule(p, ret="angles"):
    r"""Return a set of unit vectors and weights for integrating over
    the unit sphere using Lebedev quadrature.


    Parameters
    ----------
    p : int
        Precision of the quadrature rule.
    ret : {"angles", "vectors"}, optional
        Return angles or vectors (default `"angles"`).

    Returns
    -------
    sphere_rule : tuple of numpy.ndarray
        Tuple of `(phi, theta, w)` if `ret="angles"` or `(vecs, w)` if
        `ret="vectors"` (`vecs` shape: `(N, 3)`).

    See Also
    --------
    unit_sphere_lebedev_quad_available_prec :
        List of available precisions `p` for which Lebedev quadrature
        rules are available.

    Notes
    -----
    The integral of :math:`f(\phi, \theta)` is given by:

    .. math::

        \int_0^{2\pi} \int_0^{\pi} f(\phi, \theta) \sin \theta \; d\theta \; d\phi = 4 \pi \sum_i w_i f(\phi_i, \theta_i)
    """

    if _UNIT_SPHERE_LEBEDEV_QUAD_DATA is None:
        _unit_sphere_lebedev_quad_load_data()

    if p not in _UNIT_SPHERE_LEBEDEV_QUAD_DATA:
        raise ValueError(
            "p = {0} is not valid for Lebedev quadrature, or data for "
            "this quadrature rule is not available.".format(p)
        )

    phi, theta, w = _UNIT_SPHERE_LEBEDEV_QUAD_DATA[p]

    if ret == "angles":
        # Return a copy of the internal data.

        phi, theta, w = _UNIT_SPHERE_LEBEDEV_QUAD_DATA[p]
        return (np.copy(phi), np.copy(theta), np.copy(w))

    if ret == "vectors":
        # Convert angles to vectors in Cartesian coordinates.

        vecs = np.zeros((len(w), 3))

        sin_theta = np.sin(theta)

        vecs[:, 0] = np.cos(phi) * sin_theta
        vecs[:, 1] = np.sin(phi) * sin_theta
        vecs[:, 2] = np.cos(theta)

        return (vecs, np.copy(w))

    raise ValueError('Unknown return format "{0}".'.format(ret))


# ----------------------------------
# Lebedev + circle quadrature scheme
# ----------------------------------


def lebedev_circle_quad(f, p, n=None, args=None):
    r"""Perform a numerical integration of a function over the Euler
    angles `phi`, `\theta` and `\psi` using a combined Lebedev + circle
    quadrature scheme.

    Parameters
    ----------
    f : callable
        Function `f(phi, theta, psi, *args)` to be integrated.
    p : int
        Precision of the Lebedev quadrature rule used to generate `\phi`
        and `\theta`.
    n : int, optional
        Number of points in the circle quadrature rule used to generate
        `\psi` (default: automatically chosen to match the number of
        unique `\phi` in the Lebedev rule).
    args : list or None
        Optional arguments to pass to `f`.

    Returns
    -------
    int : float
        Integral of `f`.

    Notes
    -----
    The integral of :math:`f(\theta, \phi, \psi)` is computed as:

    .. math::

        \int_0^{2\pi} \int_0^{\pi} \int_0^{2\pi} f(\phi, \theta, \psi) \sin \theta \; d\psi \; d\theta \; d\phi = 8 \pi^2 \sum_i \sum_j w_i w_j f(\phi_i, \theta_i, \psi_j)
    """

    phi, theta, w_sph = unit_sphere_lebedev_quad_rule(p, ret="angles")

    # \psi from a circle quadrature rule with n_pts chosen from the
    # number of unique \phi in the Lebedev rule.

    m = len(np.unique(phi))
    psi, w_cir = unit_circle_quad_rule(m, ret="angles")

    prod = product(zip(phi, theta, w_sph), zip(psi, w_cir))

    if args is None:
        args = []

    int_sum = 0.0

    for (phi, theta, w_1), (psi, w_2) in prod:
        int_sum += w_1 * w_2 * f(phi, theta, psi, *args)

    return int_sum
