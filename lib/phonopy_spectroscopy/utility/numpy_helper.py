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


def np_asarray_copy(a, **kwargs):
    """Convert an array_like `a` to a NumPy `ndarray` object using
    `numpy.asarray` with an explicit copy if `a` is already a `ndarray`.

    Parameters
    ----------
    a : array_like
       array_like to convert.

    Returns
    -------
    a_arr : numpy.ndarray
        a as an array.
    **kwargs : any
        Keyword args to `numpy.asarray`.

    Notes
    -----
    This function reproduces some of the functionality of the `copy`
    keyword argument to `numpy.ndarray()` in NumPy >= 2.0 and is for
    compatibility with older versions.
    """

    a_arr = np.asarray(a, **kwargs)
    return a_arr.copy() if a_arr is a else a_arr


def np_readonly_view(a):
    """Return a readonly view of a NumPy ndarray `a`.

    Parameters
    ----------
    input : numpy.ndarray
        NumPy array.

    Returns
    -------
    view : numpy.ndarray
        Readonly view of `a`.
    """

    if not isinstance(a, np.ndarray):
        raise TypeError(
            "a must be a numpy.ndarray (this is most likely a bug)."
        )

    a_v = a.view()
    a_v.flags.writeable = False

    return a_v


def np_check_shape(a, shape):
    """Check a Numpy array `a` has an expected `shape`.

    Parameters
    ----------
    a : numpy.ndarray
        NumPy array to check.
    shape : tuple of (int or None)
        Expected shape of `a`. `None` indicates that a dimension
        is required but can be of variable length.

    Returns
    -------
    valid : bool
        `True` if `a` matches `shape`, otherwise `False`.
    """

    if not isinstance(a, np.ndarray):
        raise TypeError(
            "a must be a numpy.ndarray (this is most likely a bug)."
        )

    if a.ndim != len(shape):
        return False

    for s, e_s in zip(a.shape, shape):
        if e_s is not None and s != e_s:
            return False

    return True


def np_expand_dims(a, shape, expand_order=None):
    """Expand the dimensions of `a` to match a required `shape`.

    Parameters
    ----------
    a : numpy.ndarray
        NumPy array to check and expand.
    shape : tuple of (int or None)
        Expected shape of `a`. `None` specifies variable-length
        dimensions that will be added if not present.
    expand_order : tuple of int or None, optional
        Specifies the order in which to add dimensions (default: add
        starting from innermost dimension).

    Returns
    -------
    res : tuple of (np.ndarray, int)
        Tuple of `(a, added_dims)` containing the `a`, possibly with
        adjusted dimensions, and the number of added dimensions.
    """

    if not isinstance(a, np.ndarray):
        raise TypeError(
            "a must be a numpy.ndarray (this is most likely a bug)."
        )

    if expand_order is None:
        expand_order = [idx for idx, e_s in enumerate(shape) if e_s is None]

    n_dim_add = len(shape) - a.ndim

    if n_dim_add > len(expand_order):
        raise ValueError(
            "numpy.ndarray a with shape {0} has too few dimensions to "
            "match the required shape {1} after expansion."
            "".format(a.shape, shape)
        )

    if a.ndim != len(shape):
        a = np.expand_dims(a, axis=expand_order[:n_dim_add])

    if not np_check_shape(a, shape):
        raise ValueError(
            "numpy.ndarray a with expanded shape {0} does not match "
            "the required shape {1}.".format(a.shape, shape)
        )

    return (a, n_dim_add)


def np_discard_imag_if_real(a):
    """Drop the imaginary part from a NumPy array with
    `dtype=np.complex128` if all elements are real.

    Parameters
    ----------
    a : numpy.ndarray
        NumPy array.

    Returns
    -------
    a : numpy.ndarray
        `a.real` if all elements of `a` are real, otherwise `a`.
    """

    if not isinstance(a, np.ndarray):
        raise TypeError(
            "a must be a numpy.ndarray (this is most likely a bug)."
        )

    if a.dtype != np.complex128:
        raise ValueError(
            "a must have dtype np.complex128 (this is most likely a bug)."
        )

    return a if np.iscomplex(a).any() else a.real
