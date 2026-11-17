# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Helper routines for handling plot styles."""

# -------
# Imports
# -------

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import hsv_to_rgb

# ---------
# Functions
# ---------


def set_plot_style(style="phonopy_spectroscopy.cli.plotting.default"):
    """Set the Matplotlib plot style.

    Parameters
    ----------
    style : str, optional
        Name, path or URL to style to use (default:
        "phonopy_spectroscopy.cli.plotting.default").
    """

    plt.style.use(style)


def cscale(n, h_s, h_e, s=1.0, v=1.0):
    """Generate a sequence of `n` RGB colours between hues `h_s` and
    `h_e` with uniform saturation and value.

    Params
    ------
    n : int
        Number of colours.
    h_s, h_e : float
        Hues in degrees.
    s, v : float, optional
        Saturation and value (default: 1.0).

    Returns
    -------
    c : numpy.ndarray
        Colours (shape: `(n, 3)`).
    """

    if n < 1:
        raise ValueError("n must be at least 1.")

    if (s < 0.0 or s > 1.0) or (v < 0.0 or v > 1.0):
        raise ValueError("s and v must be between 0 and 1.")

    h = (np.linspace(h_s, h_e, n) / 360.0) % 1.0
    s = np.full((n,), s)
    v = np.full((n,), v)

    return hsv_to_rgb(np.stack((h, s, v), axis=1))


def cscale_fire(n):
    """Generate a sequence of `n` RGB colours between blue and orange in
    HSV colour space.

    Params
    ------
    n : int
        Number of colours.

    Returns
    -------
    c : numpy.ndarray
        Colours (shape: `(n, 3)`).
    """

    return cscale(n, 240.0, 390.0, s=1.0, v=1.0)


def cscale_ice(n):
    """Generate a sequence of `n` RGB colours between blue and cyan in
    HSV colour space.

    Params
    ------
    n : int
        Number of colours.

    Returns
    -------
    c : numpy.ndarray
        Colours (shape: `(n, 3)`).
    """

    return cscale(n, 240.0, 180.0, s=1.0, v=1.0)
