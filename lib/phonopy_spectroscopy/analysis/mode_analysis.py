# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines for analysing phonon modes."""


# -------
# Imports
# -------


import warnings

import numpy as np

from ..constants import ZERO_TOLERANCE

from ..utility.numpy_helper import np_check_shape


# ---------
# Functions
# ---------


def participation_ratio(evec):
    """
    Calculate the phonon participation ratio.

    Params
    ------
    evec : array_like
        Phonon eigenvector (shape: `(N, 3)`).

    Returns
    -------
    pr : float
        Participation ratio.

    Notes
    -----
    The phonon participation is given by:

    .. math::
         P = \frac{1}{N \sum_j \left( \sum_\alpha | \boldsymbol{W}_j^\alpha |^2 \right)^2 }

    It can take values between :math:`\frac{1}{N}` and 1, respectively,
    for fully delocalised and fully localised modes.
    """

    evec = np.asarray(evec)

    if not np_check_shape(evec, (None, 3)):
        raise ValueError("evec must be an array_like with shape (N, 3).")

    # Take square modulus (for complex evecs).

    sq_mod = np.abs(evec) ** 2

    # Check eigenvector is normalised.

    norm = sq_mod.sum()

    if np.abs(norm - 1.0) > ZERO_TOLERANCE:
        warnings.warn(
            "evec may not be normalied: norm = {1:.5f}.", UserWarning
        )

    n, _ = sq_mod.shape

    return 1.0 / (n * np.sum(np.sum(sq_mod, axis=1) ** 2))
