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


def get_format_string(a, dp=3):
    """Get a format code for printing the values in an array with
    the field size set for the largest value with a sign and the
    specified number of decimal places.

    Params
    ------
    a : scalar or array_like
        Values.
    dp : int, optional
        Number of decimal places (default: 3).

    Returns
    -------
    fmt_code : str
        Format code.
    """

    a = np.asanyarray(a)

    size = len(str(int(np.ceil(np.abs(a)).max())))
    return "{{0: >{0}.{1}f}}".format(size + 2 + dp, dp)


def print_header(h, u="-", newline=False):
    """Print a header with an equal-width underline and optional
    trailing newline.

    Params
    ------
    h : str
        Header.
    u : str, optional
        Underline character (default: "-").
    newline : bool, optional
        Print a trailing newline after the header and underline
        (default: `True`).
    """

    h, u = str(h), str(u)

    if len(u) != 1:
        raise ValueError("u must be a single character.")

    print(h)
    print(u * len(h))

    if newline:
        print("")


def print_3x3_matrix(m, dp=3):
    """Print a 3x3 matrix with the specified number of decimal places.

    Params
    ------
    m : array_like
        Matrix (shape: `(3, 3)`).
    """

    m = np.asanyarray(m)

    if m.shape != (3, 3):
        raise ValueError("m must be an array_like with shape (3, 3).")

    fmt = get_format_string(m, dp=dp)

    for r in m:
        print("[  {0}  {1}  {2}  ]".format(*[fmt.format(x) for x in r]))
