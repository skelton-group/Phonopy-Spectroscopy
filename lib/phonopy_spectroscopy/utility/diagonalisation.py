# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for working with eigenvalues and eigenvectors from
diagonalising matrices."""

# -------
# Imports
# -------

import numpy as np

from scipy.optimize import linear_sum_assignment

from ..utility.numpy_helper import atleast_3d

# ---------------
# Branch tracking
# ---------------


def reorder_with_branch_tracking(evals, evecs):
    """Reorder a sequence of complex eigenvalue/eigenvector sets with a
    "velocity-based" branch-tracking algorithm.

    Parameters
    ----------
    evals : array_like
        Eigenvalues (shape: `(O, D)` or `(D,)`).
    evecs : array_like
        Eigenvectors (shape: `(O, D, D)` or `(D, D)`).

    Returns
    -------
    res : tuple of numpy.ndarray
        Tuple of reordered `(evals, evecs)` (shapes: `(O, D)`,
        `(O, D, D)`).

    Notes
    -----
    `evecs` are assumed to be in "colum-major" order, i.e. the `j`th
    eigenvector at the `i`th index in the set is obtained as
    `evecs[i, :, j]`. (This corresponds to the memory layout obtained
    from the diagonalisation routines in `numpy.linalg`.)

    If a single set of eivenvalues/eigenvectors with shapes `(D,)` and
    `(D, D)` is input the output has a trailing initial dimension
    prepended.
    """

    evals = np.atleast_2d(np.asanyarray(evals, dtype=np.complex128))
    evecs = atleast_3d(np.asanyarray(evecs, dtype=np.complex128))

    if evals.ndim != 2:
        raise ValueError("evals must be an array_like with shape (O, D).")

    o, d = evals.shape

    if evecs.shape != (o, d, d):
        raise ValueError("evecs must be an array_like with shape (O, D, D).")

    if d == 1:
        # 1D -> nothing to do.
        return (evals, evecs)

    evals_new = np.zeros_like(evals, dtype=np.complex128)
    evecs_new = np.zeros_like(evecs, dtype=np.complex128)

    evals_new[0] = evals[0]
    evecs_new[0] = evecs[0]

    for i in range(o):
        ev_old, ev_new = evecs_new[i - 1], evecs[i]

        # For the 2nd point, compare to the eigenvectors in the 1st
        # point. For the 3rd point onwards, compare to a predicted
        # eigenvector obtained from the "velocity" between the previous
        # points.

        ev_cmp = ev_old

        if i > 1:
            ev_cmp = ev_cmp + (ev_cmp - evecs_new[i - 2])

        # Reorder the eigenvectors using linear-sum assignment based on
        # the modal assurance criterion (MAC).

        mac = np.abs(np.matmul(ev_new.T.conj(), ev_cmp)) ** 2
        _, col_inds = linear_sum_assignment(-1.0 * mac)

        # col_inds map new eigenvectors to the closest-matching
        # references. We need the indices that reorder the new
        # eigenvectors to match the reference order, which can be
        # obtained using argsort.

        inv = np.argsort(col_inds)

        ev_new = ev_new[:, inv]

        # The phase of the complex vectors is arbitrary -> realign
        # between points.

        for j in range(d):
            dp = np.vdot(ev_old[:, j], ev_new[:, j])
            ev_new[:, j] *= np.exp(-1.0j * np.angle(dp))

        evals_new[i] = evals[i, col_inds]
        evecs_new[i] = ev_new

    return (evals_new, evecs_new)
