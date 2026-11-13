# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Helper routines for simulating spectra."""

# -------
# Imports
# -------

import numpy as np

from ..constants import ZERO_TOLERANCE

# ---------------
# Automatic range
# ---------------


def auto_range(
    freqs, lws, f_range=None, f_res=None, f_pad_mult=5.0, f_res_min_pts=1000
):
    """Given a set of frequencies and linewidths, determine a suitable
    minimum/maximum frequency and resolution for simulating a spectrum.

    Parameters
    ----------
    freqs : any
        Frequencies.
    lws : any
        Linewidths.
    f_range : tuple of float or None, optional
        `(min, max)` frequency range (default: automatic).
    f_res : float, optional
        Frequency resolution (default: automatically determined).
    f_pad_mult : float, optional
        "Padding" for frequency range as a multiplier of the largest
        linewidth (default: 5.0, ignored when `f_range` is set).
    f_res_min_pts : int, optional
        Minimum number of points between min/max frequency (default:
        1000, ignored when `f_res` is set).

    Returns
    -------
    params : tuple of (tuple of float, int)
        Tuple of `((min, max), res)` for simulating a spectrum.
    """

    freqs, lws = np.asanyarray(freqs), np.asanyarray(lws)

    f_min, f_max = freqs.min(), freqs.max()
    l_max = lws.max()

    if f_range is not None:
        f_min, f_max = f_range

        # The code does not use any supported frequency units for which
        # the frequencies are comparable to or below ZERO_TOLERANCE.

        if np.abs(f_max - f_min) < ZERO_TOLERANCE:
            raise ValueError(
                "If f_range is specified, min and max cannot be the same."
            )

        if f_max < f_min:
            f_min, f_max = f_max, f_min
    else:
        # "Pad" f_min and f_max to allow for line broadening.

        pad = f_pad_mult * l_max

        f_min -= pad
        f_max += pad

    if f_res is None:
        # Choose the resolution to ensure at least f_res_min_pts between
        # f_min and f_max. We estimate it to give f_res_min_pts, then
        # round the result down to the nearest power of 10 to get a
        # "nice" frequency spacing.

        f_res = (f_max - f_min) / f_res_min_pts
        f_res = np.power(10.0, np.floor(np.log10(f_res)))

    # If the range is being set automatically, "align" the minimum and
    # maximum to an order of magnitude.

    if f_range is None:
        align = np.power(10.0, np.floor(np.log10((f_max - f_min) / 10.0)))

        f_min = align * np.floor(f_min / align)
        f_max = align * np.ceil(f_max / align)

    return ((f_min, f_max), f_res)
