# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Helper routines for printing to the terminal."""

# -------
# Imports
# -------

import numpy as np

# ---------
# Functions
# ---------


def get_fmt_code(a, dp=3):
    """Get a format code for printing the values in an array with
    the field size set for the largest value with a sign and the
    specified number of decimal places.

    Parameters
    ----------
    a : scalar or array_like
        Values.
    dp : int, optional
        Number of decimal places (default: 3).

    Returns
    -------
    fmt_code : str
        Format code.
    """

    size = len(str(int(np.ceil(a).max())))
    return "{{0: >{0}.{1}f}}".format(size + 1 + dp, dp)


def print_3x3_matrix(m, dp=3):
    """Print a 3x3 matrix with the specified number of decimal places.

    Parameters
    ----------
    m : array_like
        Matrix (shape: `(3, 3)`).
    """

    m = np.asanyarray(m)

    if m.shape != (3, 3):
        raise ValueError("m must be an array_like with shape (3, 3).")

    fmt_code = get_fmt_code(m, dp=dp)

    for r in m:
        print("[  {0}  {1}  {2}  ]".format(*[fmt_code.format(x) for x in r]))
