# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for working with eigenvalues and eigenvectors from
diagonalisin matrices."""

# -------
# Imports
# -------

import numpy as np

from scipy.optimize import linear_sum_assignment

from ..utility.numpy_helper import np_expand_dims

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

    Input with a single set of eivenvalues/eigenvectors is accepted but
    redundant. In this case, the function returns the input data, with a
    trailing initial dimension prepended if the shapes were `(D,)` and
    `(D, D)`.

    1D input is also considered redundant and the function will again
    return the input data unchanged.
    """

    evals, _ = np_expand_dims(
        np.asarray(evals, dtype=np.complex128), (None, None), expand_order=(0,)
    )

    d_1, d_2 = evals.shape

    evecs, _ = np_expand_dims(
        np.asarray(evecs, dtype=np.complex128),
        (d_1, d_2, d_2),
        expand_order=(0,),
    )

    if d_1 == 1 or d_2 == 1:
        # Nothing to do -> return input.
        return (evals, evecs)

    evals_new = np.zeros_like(evals, dtype=np.complex128)
    evecs_new = np.zeros_like(evecs, dtype=np.complex128)

    evals_new[0] = evals[0]
    evecs_new[0] = evecs[0]

    for i in range(1, d_1):
        ev_old, ev_new = evecs_new[i - 1], evecs[i]

        # For the 2nd point, compare to the eigenvectors in the 1st
        # point. For the 3rd point onwards, compare to a predicted
        # eigenvector obtained from the "velocity" between the previous
        # points.

        ev_cmp = ev_old

        if i > 1:
            ev_cmp = ev_cmp + (ev_cmp - evecs[i - 2])

        # Reorder the eigenvectors using linear-sum assignment based on
        # the modal assurance criterion (MAC).

        mac = np.abs(np.matmul(ev_new.T.conj(), ev_cmp)) ** 2
        _, col_inds = linear_sum_assignment(-1.0 * mac)

        ev_new = ev_new[:, col_inds]

        # The phase of the complex vectors is arbitrary -> realign
        # between points.

        for j in range(d_2):
            dp = np.vdot(ev_old[:, j], ev_new[:, j])
            ev_new[:, j] *= np.exp(-1.0j * np.angle(dp))

        evals_new[i] = evals[i, col_inds]
        evecs_new[i] = ev_new

    return (evals_new, evecs_new)
