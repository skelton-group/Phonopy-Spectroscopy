# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Miscellaneous routines for implementing the command-line interface."""

# -------
# Imports
# -------

import numpy as np

from ...utility.geometry import rotation_matrix_from_vectors
from ...utility.numpy_helper import np_check_shape

# ---------
# Functions
# ---------


def get_crystal_rotation_matrix(struct, hkl, geom, rot=None):
    """Get the rotation matrix for a single-crystal measurement.

    Parameters
    ----------
    struct : Structure
        Crystal structure.
    hkl : tuple of int
        Miller index of crystal surface to orient along incident
        direction.
    geom : Geometry
        Measurement geometry.
    rot : array_like or None, optional
        Rotation matrix to realign crystal after `hkl` rotation
        (default: `None`).

    Returns
    -------
    rot : numpy.ndarray
        Rotation matrix.
    """

    r = rotation_matrix_from_vectors(
        struct.real_space_normal(hkl, conv=True),
        -1.0 * geom.incident_direction,
    )

    if rot is not None:
        if not np_check_shape(rot, (3, 3)):
            raise ValueError(
                "If supplied, rot must be an array_like with shape (3, 3)."
            )

        r = np.dot(rot, r)

    return r
