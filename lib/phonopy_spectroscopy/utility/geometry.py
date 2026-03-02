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
from .numpy_helper import np_check_shape, np_expand_dims

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


def parse_direction(dirn, norm=False):
    """Parse a direction and return a 3D vector.

    Parameters
    ----------
    dirn : array_like or str
        3D `(x, y, z)` vector, or one of `{"x", "+x", "-x", "y", "+y",
        "-y", "z", "+z", "-z"}`.
    norm : bool, optional
        If `True`, normalise the vector (default: `False`).

    Returns
    -------
    v : numpy.ndarray
        3D direction vector.
    """

    dirn_str = str(dirn).strip().lower()

    if dirn_str in _DIRECTION_STRING_LUT:
        return np.array(_DIRECTION_STRING_LUT[dirn_str], dtype=np.float64)

    dirn = np.asarray(dirn)

    if np.ndim(dirn) == 1:
        v = None

        if len(dirn) == 3:
            v = np.array(dirn, dtype=np.float64)

        if len(dirn) == 2:
            x, y = dirn
            v = np.array([x, y, 0.0], dtype=np.float64)

        if norm:
            m = np.linalg.norm(v)

            # Cannot normalise a zero vector.

            if m > 0.0:
                return v / norm

        return v

    raise ValueError('Invalid direction specifier "{0}".'.format(dirn))


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

    a = parse_direction(a, norm=True)
    b = parse_direction(b, norm=True)

    # If the vectors are parallel or orthogonal the solution is trivial.

    cos_theta = np.dot(a, b)

    if np.abs(cos_theta - 1.0) < ZERO_TOLERANCE:
        # cos(\theta) = 1 -> vectors parallel.
        return np.identity(3)

    if np.abs(cos_theta + 1.0) < ZERO_TOLERANCE:
        # cos(\theta) = -1 -> vectors antiparallel.
        return -1.0 * np.identity(3)

    # Otherwise, a -> b is a 180 degree rotation around the normalised
    # vector bisecting a and b, and the required rotation matrix can be
    # obtained using Rodrigues' formula.

    k = a + b
    k /= np.linalg.norm(k)

    return rotation_matrix_from_axis_angle(k, 180.0)


def rotation_matrix_from_axis_angle(k, theta):
    r"""Compute a 3D rotation matrix for a rotation around an axis `k`
    by angle `theta`.
    
    Parameters
    ----------
    k : array_like or str
        Rotation axis.
    theta : float
        Rotation angle.
    
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

    v_x, v_y, v_z = parse_direction(k, norm=True)

    w = np.array(
        [[0.0, -v_z, v_y], [v_z, 0.0, -v_x], [-v_y, v_x, 0.0]],
        dtype=np.float64,
    )

    theta = np.radians(theta)

    return (
        np.identity(3)
        + np.sin(theta) * w
        + (1.0 - np.cos(theta)) * np.matmul(w, w)
    )


@njit
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
    r"""Rotate a tensor or set of tensors `t` by a rotation matrix `r`.

    Parameters
    ----------
    t : array_like
        Single tensor (shape `(3, 3)`) or set of `N` tensors (shape
        `(N, 3, 3)`).
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

    t = np.asarray(t)
    r = np.asarray(r)

    t, n_dim_add = np_expand_dims(t, (None, 3, 3))

    if not np_check_shape(r, (3, 3)):
        raise ValueError("r must be an array_like with shape (3, 3).")

    t_prime = np.matmul(r, np.matmul(t, r.T))

    return t_prime if n_dim_add == 0 else t_prime[0]
