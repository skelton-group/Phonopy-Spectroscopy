# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Helper routines for working with NumPy."""

# -------
# Imports
# -------

import numpy as np

# ---------
# Functions
# ---------


def atleast_nd(*arys, n):
    """View inputs as arrays with at least `n` dimensions by prepending
    as many singleton dimensions as needed.

    Parameters
    ----------
    arys1, arys2, ... : array_like
        Arrays.
    n : int
        Minimum dimension.

    Returns
    -------
    res, res1, ... : numpy.ndarray or subclass
        Array with at least `n` dimensions.
    """

    res = []

    for a in arys:
        a = np.asanyarray(a)

        while a.ndim < n:
            a = a[np.newaxis, ...]

        res.append(a)

    return res[0] if len(arys) == 1 else tuple(res)


def atleast_3d(*arys):
    """View inputs as arrays with at least three dimensions.

    Parameters
    ----------
    arys1, arys2, ... : array_like
        Arrays.

    Returns
    -------
    res, res1, ... : numpy.ndarray or subclass
        Array with at least three dimensions.

    Notes
    -----
    This differs from NumPy's equivalent function in that any additional
    dimensions are always prepended.
    """

    return atleast_nd(*arys, n=3)
