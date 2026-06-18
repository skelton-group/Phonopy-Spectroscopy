# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for numerical differentiation."""

# -------
# Imports
# -------

import numpy as np

# --------------
# Data/Constants
# --------------

_CENTRAL_DIFF_COEFFS = {
    1: {
        2: ([-1.0, 1.0], [-0.5, 0.5]),
        4: (
            [-2.0, -1.0, 1.0, 2.0],
            [1.0 / 12.0, -2.0 / 3.0, 2.0 / 3.0, -1.0 / 12.0],
        ),
        6: (
            [-3.0, -2.0, -1.0, 1.0, 2.0, 3.0],
            [
                -1.0 / 60.0,
                3.0 / 20.0,
                -3.0 / 4.0,
                3.0 / 4.0,
                -3.0 / 20.0,
                1.0 / 60.0,
            ],
        ),
        8: (
            [-4.0, -3.0, -2.0, -1.0, 1.0, 2.0, 3.0, 4.0],
            [
                1.0 / 280.0,
                -4.0 / 105.0,
                1.0 / 5.0,
                -4.0 / 5.0,
                4.0 / 5.0,
                -1.0 / 5.0,
                4.0 / 105.0,
                -1.0 / 280.0,
            ],
        ),
    }
}

# -------------------
# Central differences
# -------------------


def central_difference_available_orders():
    """Return a list of orders `n` for which central difference
    coefficients are available.

    Returns
    -------
    orders : list of int
        Available orders `n`.
    """

    return list(_CENTRAL_DIFF_COEFFS.keys())


def central_difference_available_precs(n):
    """Return a list of precisions `p` for which central difference
    coefficients for derivatives of order `n` are available.

    Parameters
    ----------
    n : int
        Order of differentiation `n`.

    Returns
    p : list of int
        Available precisions `p` for order `n`.
    """

    if n in _CENTRAL_DIFF_COEFFS:
        return list(_CENTRAL_DIFF_COEFFS[n].keys())

    raise RuntimeError(
        "Central difference coefficients for derivatives of order "
        "n = {0} are not available.".format(n)
    )


def central_difference_coefficients(n, p):
    """Return the steps and coefficients for calculating derivatives of
    order `n` with a central difference scheme of precision `p`.

    Parameters
    ----------
    n : int
        Order of differentiation `n`.
    p : int
        Precision `p` of central difference scheme.

    Returns
    -------
    steps_coeffs : tuple of numpy.ndarray
        A `(steps, coeffs)` tuple with the steps and coefficients for
        the requested central difference scheme.

    See Also
    --------
    central_difference_available_orders :
        List of available orders `n` for which central difference
        coefficients are available.
    central_difference_available_precs :
        List of available precisions `p` for which central difference
        coefficients are available for order `n`.
    """

    if n not in _CENTRAL_DIFF_COEFFS or p not in _CENTRAL_DIFF_COEFFS[n]:
        raise RuntimeError(
            "Central difference coefficients for order = {0} and "
            "prec = {1} are not available.".format(n, p)
        )

    steps, coeffs = _CENTRAL_DIFF_COEFFS[n][p]

    return (
        np.array(steps, dtype=np.float64),
        np.array(coeffs, dtype=np.float64),
    )
