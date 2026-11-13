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

# ---------
# Functions
# ---------


def participation_ratio(evecs):
    r"""Calculate the phonon participation ratio.

    Params
    ------
    evecs : array_like
        Phonon eigenvector(s) (shape: `(N, 3)` or `(M, N, 3)`).

    Returns
    -------
    pr : float
        Participation ratio(s) (scalar or shape: `(M,)`).

    Notes
    -----
    The phonon participation is given by:

    .. math::

         P = \frac{1}{N \sum_j \left( \sum_\alpha | \boldsymbol{W}_j^\alpha |^2 \right)^2 }

    :math:`P` can take values between :math:`\frac{1}{N}` and 1,
    respectively, for fully delocalised and fully localised modes.
    """

    evecs = np.asanyarray(evecs, dtype=np.float64)

    ret_single = False

    if evecs.ndim == 2:
        evecs = evecs[np.newaxis, :, :]
        ret_single = True

    if evecs.ndim != 3 or evecs.shape[-1] != 3:
        raise ValueError(
            "evecs must be an array_like with shape (N, 3) or (M, N, 3)."
        )

    # Take square modulus (for complex evecs) and sum over positional
    # components.

    sq_mods = (np.abs(evecs) ** 2).sum(axis=-1)

    # Check eigenvectors are normalised.

    abs_norms = sq_mods.sum(axis=-1)

    if ((abs_norms - 1.0) > ZERO_TOLERANCE).any():
        warnings.warn(
            "One or more eigenvectors may not be normalied: max. norm "
            "= {0:.5f}.".format(abs_norms.max()),
            RuntimeWarning,
        )

    prs = 1.0 / (sq_mods.shape[1] * np.sum(sq_mods**2, axis=-1))

    return prs[0] if ret_single else prs
