# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for geometry handling."""

# -------
# Imports
# -------

import numpy as np

from ..constants import ZERO_TOLERANCE

try:
    from numba import njit
except ImportError:
    from .numba_helper import dummy_njit as njit

# ---------
# Functions
# ---------

_DIRECTION_STRING_LUT = {
    "x": [1.0, 0.0, 0.0],
    "+x": [1.0, 0.0, 0.0],
    "-x": [-1.0, 0.0, 0.0],
    "y": [0.0, 1.0, 0.0],
    "+y": [0.0, 1.0, 0.0],
    "-y": [0.0, -1.0, 0.0],
    "z": [0.0, 0.0, 1.0],
    "+z": [0.0, 0.0, 1.0],
    "-z": [0.0, 0.0, -1.0],
}


def parse_direction(x):
    """Parse a direction and return a normalised 3D vector.

    Parameters
    ----------
    dirn : array_like or str
        3D `(x, y, z)` vector, or one of `{"x", "+x", "-x", "y", "+y",
        "-y", "z", "+z", "-z"}`.

    Returns
    -------
    v : numpy.ndarray
        3D direction vector.
    """

    # If x is a string type, assume it specifies one of the hard-coded
    # directions in _DIRECTION_STRING_LUT.

    if isinstance(x, (str, np.str_)):
        dirn_str = str(x).strip().lower()

        if dirn_str in _DIRECTION_STRING_LUT:
            return np.array(_DIRECTION_STRING_LUT[dirn_str], dtype=np.float64)
        else:
            raise ValueError(f'Invalid string specifier "{x}".')

    # Otherwise, assume x is a (potentially) complex number; if it
    # isn't, take the real part to avoid forcing complex arithmetic.

    x = np.real_if_close(
        np.asanyarray(x, dtype=np.complex128), tol=ZERO_TOLERANCE
    )

    if x.ndim != 1:
        raise ValueError("x must be a vector.")

    if len(x) == 2:
        x, y = x
        x = np.array([x, y, 0.0], dtype=x.dtype)

    if len(x) != 3:
        raise ValueError("x must specify a 2D or 3D vector.")

    norm = np.linalg.norm(x)

    if np.isclose(norm, 0.0):
        raise ValueError("x must specify a non-zero vector.")

    return x / norm


def plane_polar_to_cartesian(p, vec_3d=False):
    """Convert plane polar to Cartesian coordinates.

    Parameters
    ----------
    p : array_like
        Polar coordinates `(r, theta)` (shape: `(..., 2)`).
    vec_3d : bool, optional
        If `True`, return 3D `(x, y, z)` vectors, otherwise return 2D
        `(x, y)` vectors (default: `False`).

    Returns
    -------
    v : numpy.ndarray
        Cartesian coordinates (shape: `(..., 2)` or `(..., 3)`).
    """

    p = np.asanyarray(p, dtype=np.float64)

    if p.ndim < 1 or p.shape[-1] != 2:
        raise ValueError("p must be an array_like with shape (..., 2).")

    r, phi = p[..., 0], p[..., 1]

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    if vec_3d:
        z = np.zeros_like(x, dtype=np.float64)
        return np.stack((x, y, z), axis=-1)
    else:
        return np.stack((x, y), axis=-1)


def cartesian_to_spherical_polar(v):
    """Convert vectors in Cartesian coordinates to spherical polar
    coordinates.

    Parameters
    ----------
    v : array_like
        Cartesian coordinates (shape: `(..., 3)`).

    Returns
    -------
    sp : numpy.ndarray
        Spherical polar coordinates `(r, phi, theta)` (same shape as
        `v`).
    """

    v = np.asanyarray(v, dtype=np.float64)

    if v.ndim < 1 or v.shape[-1] != 3:
        raise ValueError("v must be an array_like with shape (..., 3).")

    r = np.linalg.norm(v, axis=-1)
    phi = np.arctan2(v[..., 1], v[..., 0])
    theta = np.arccos(v[..., 2] / r)

    return np.stack((r, phi, theta), axis=-1)


def rotation_matrix_from_vectors(a, b):
    """Compute a 3D rotation matrix that rotates a vector `a` onto
    another vector `b`, given by a 180 degree rotation about the
    midpoint of `a` and `b`.

    Parameters
    ----------
    a : array_like or str
        Initial vector.
    b : array_like or str
        Rotated vector.

    Returns
    -------
    r : numpy.ndarray
        Computed rotation matrix (shape: `(3, 3)`).
    """

    # parse_direction will validate and normalise.

    a = parse_direction(a)
    b = parse_direction(b)

    cos_theta = np.dot(a, b)

    if np.abs(cos_theta - 1.0) < ZERO_TOLERANCE:
        # cos(\theta) = 1 -> vectors parallel.
        return np.identity(3)

    if np.abs(cos_theta + 1.0) < ZERO_TOLERANCE:
        # cos(\theta) = -1 -> vectors antiparallel. A 180 degree
        # rotation about any vector c perpendicular to a is valid. We
        # arbitrarily "guess" c = x and swap to c = y if a is nearly
        # parallel to x.

        c = np.array([1.0, 0.0, 0.0], dtype=np.float64)

        if np.abs(np.dot(a, c)) > 0.9:
            # 0.9 is an arbitrary threshold for "nearly parallel".
            c = np.array([0.0, 1.0, 0.0], dtype=np.float64)

        k = np.cross(a, c)

        return rotation_matrix_from_axis_angle(k, np.pi)

    # Otherwise, a -> b is a 180 degree rotation around the normalised
    # vector bisecting a and b, and the required rotation matrix can be
    # obtained using Rodrigues' formula.

    k = a + b
    k /= np.linalg.norm(k)

    return rotation_matrix_from_axis_angle(k, np.pi)


def rotation_matrix_from_axis_angle(k, theta):
    r"""Compute a 3D rotation matrix for a rotation around an axis `k`
    by angle `theta`.
    
    Parameters
    ----------
    k : array_like or str
        Rotation axis.
    theta : float
        Rotation angle in radians.
    
    Returns
    -------
    r : numpy.ndarray
        Computed rotation matrix (shape: `(3, 3)`).

    Notes
    -----
    The rotation matrix is computed using Rodrigues' rotation formula:

    .. math::

        \textbf{R} = \textbf{I} + \textbf{K} \sin \theta + \textbf{K}^2 [1 - \cos \theta]
        
    where:

    .. math::

        \textbf{K} = \begin{bmatrix}   0   & -k_z &  k_y \\
                                       k_z &  0   & -k_x \\
                                      -k_y &  k_x &  0   \end{bmatrix}
    """

    v_x, v_y, v_z = parse_direction(k)

    w = np.array(
        [[0.0, -v_z, v_y], [v_z, 0.0, -v_x], [-v_y, v_x, 0.0]],
        dtype=np.float64,
    )

    return np.identity(3) + np.sin(theta) * w + (1.0 - np.cos(theta)) * (w @ w)


@njit(fastmath=True, inline="always")
def direction_cosine(phi, theta, psi):
    """Compute the direction cosine for the Euler angles `phi`, `theta`
    and `psi`.

    Parameters
    ----------
    phi, theta, psi : float
        Euler angles in Radians.

    Returns
    -------
    r : numpy.ndarray
        Computed direction cosine (shape: `(3, 3)`).
    """

    s_ph, s_th, s_ps = np.sin(phi), np.sin(theta), np.sin(psi)
    c_ph, c_th, c_ps = np.cos(phi), np.cos(theta), np.cos(psi)

    return np.array(
        [
            [
                c_ph * c_th * c_ps - s_ph * s_ps,
                s_ph * c_th * c_ps + c_ph * s_ps,
                -s_th * c_ps,
            ],
            [
                -c_ph * c_th * s_ps - s_ph * c_ps,
                -s_ph * c_th * s_ps + c_ph * c_ps,
                s_th * s_ps,
            ],
            [c_ph * s_th, s_ph * s_th, c_th],
        ],
        dtype=np.float64,
    )


def rotate_tensors(t, r):
    r"""Rotate tensors `t` by a rotation matrix `r`.

    Parameters
    ----------
    t : array_like
        Tensors (shape: `(..., 3, 3)`).
    r : array_like
        Rotation matrix (shape: `(3, 3)`).

    Returns
    -------
    t_prime : numpy.ndarray
        Rotated tensors (same shape as `t`).

    Notes
    -----
    Tensor rotations are computed as:

    .. math::
        \textbf{T}^\prime = \textbf{R} \textbf{T} \textbf{R}^{-1}
    """

    t = np.asanyarray(t)
    r = np.asanyarray(r)

    if t.ndim < 2 or t.shape[-2:] != (3, 3):
        raise ValueError("t must be an array_like with shape (..., 3, 3).")

    if r.shape != (3, 3):
        raise ValueError("r must be an array_like with shape (3, 3).")

    return r @ t @ r.T
