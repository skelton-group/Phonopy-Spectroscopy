# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Helper routines for parsing command-line arguments for implementing
the command-line interface."""

# -------
# Imports
# -------

import numpy as np

# ---------
# Functions
# ---------


def parse_frac(inp):
    """Parse a value specified as a fraction using the `/` character to
    sepearate the numerator and denominator.

    Parameters
    ----------
    inp : str
        Input string to parse.

    Returns
    -------
    val : float
        Parsed value.
    """

    if "/" in inp:
        num, denom = inp.strip().split("/")
        return float(num) / float(denom)

    return float(inp.strip())


def parse_3x3_matrix(inp):
    """Parse a 3x3 matrix from an input string.

    Parameters
    ----------
    inp : str
        Input string to parse.

    Returns
    -------
    mat : numpy.ndarray
        Parsed matrix (shape: (3, 3)).

    See Also
    --------
    parse_frac : Parse a value specified as a fraction.

    Notes
    -----
    inp may contain one, three, six or nine values, any of which may be
    specified as fractions, which are interpreted as follows:
      * 1 value: `xx = yy = zz = val`
      * 3 values: `xx = vals[0]`z`, `yy = vals[1]`, `zz = vals[2]`
      * 6 values: `xx = vals[0]`, `xy = yx = vals[1]`,
        `xz = zx = vals[2]`, `yy = vals[3]`, `yz = zy = vals[4]`,
        `zz = vals[5]`
      * 9 values: `xx, xy, xz, yx, yy, yz, zx, zy, zz = vals[:]`

    If fractional values are used, whitespace between the numerator and
    denominator will result in incorrect parsing and/or an error.
    """

    vals = np.array(
        [parse_frac(v) for v in inp.strip().split()], dtype=np.float64
    )

    if len(vals) == 1:
        # Uniform diagonal: xx = yy = zz.
        return vals[0] * np.eye(3, dtype=np.float64)

    elif len(vals) == 3:
        # Diagonals: xx, yy, zz.
        return np.diag(vals)

    elif len(vals) == 6:
        # Upper triangle: xx, xy, xz, yy, yz, zz.

        xx, xy, xz, yy, yz, zz = vals

        return np.array(
            [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], dtype=np.float64
        )

    elif len(vals) == 9:
        # All elements: xx, xy, xz, yx, yy, yz, zx, zy, zz.
        return vals.reshape((3, 3))

    raise ValueError("inp must specify one, three, six or nine values.")
