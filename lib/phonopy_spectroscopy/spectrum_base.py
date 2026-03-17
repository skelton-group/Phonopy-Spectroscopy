# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Core routines and a base classes for simulating spectra."""


# -------
# Imports
# -------


import math

import numpy as np

from .constants import ZERO_TOLERANCE

from .units import (
    convert_frequency_units,
    get_supported_frequency_units,
    get_frequency_unit_text_label,
    get_frequency_unit_plot_label,
)

from .utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
)


# ----------------------
# Automatic x-axis range
# ----------------------


_AUTO_X_PAD_MULT = 5.0
_AUTO_X_MIN_NUM_PTS = 1000


def auto_x(freqs, lws, x_range=None, x_res=None):
    """Given a set of frequencies and linewidths, determine a suitable
    x-axis minimum, maximum and resolution for simulating a spectrum.

    Parameters
    ----------
    freqs : array_like
        Frequencies.
    lws : array_like
        Linewidths.
    x_range : tuple of float or None. optional
        Specify a (min, max) range (default: automatically determined).
    x_res : float, optional
        Specify a resolution (default: automatically determined).

    Returns
    -------
    params : tuple of (tuple of float, int)
        Tuple of `((min, max), res)` for simulating a spectrum.
    """

    # We only use the minimum and maximum frequencies, and the minimum
    # linewidth, so no need to explicitly convert them to NumPy arrays.

    f_min, f_max = np.min(freqs), np.max(freqs)
    lw_min = np.min(lws)

    auto_range = True

    if x_range is not None:
        f_min, f_max = x_range

        if np.abs(f_max - f_min) < ZERO_TOLERANCE:
            raise ValueError(
                "If a range is specified, min and max cannot be the same."
            )

        if f_max < f_min:
            f_min, f_max = f_max, f_min

        auto_range = False

    # If required, set a default range. To start with, we estimate an
    # approximate range allowing for line broadening.

    if auto_range:
        padding = _AUTO_X_PAD_MULT * np.max(lws)

        f_min -= padding
        f_max += padding

    # If required, automatically determine resolution.

    if x_res is None:
        # Choose the resolution to ensure at least _AUTO_X_MIN_NUM_PTS
        # between f_min and f_max.

        temp = (f_max - f_min) / _AUTO_X_MIN_NUM_PTS
        x_res = np.power(10.0, np.floor(np.log10(temp)))

    # If the range is being set automatically, "align" the minimum and
    # maximum to an order of magnitude.

    if auto_range:
        align = np.power(10.0, np.floor(np.log10((f_max - f_min) / 10.0)))

        f_min = align * np.floor(f_min / align)
        f_max = align * np.ceil(f_max / align)

    return ((f_min, f_max), x_res)


# ------------------
# SpectrumBase class
# ------------------


class SpectrumBase:
    """Base class for simulated spectra."""

    def __init__(self, x_range=None, x_res=None, x=None, x_units="thz"):
        """Create a new instance of the `SpectrumBase` class.

        Parameters
        ----------
        x_range : tuple of float or None, optional
            `(min, max)` range of x-axis of simulated spectrum in
            `x_units` (default: `None`).
        x_res : float or None, optional
            Resolution of x-axis in `x_units` (`None`).
        x : array_like or None, optional
            x-axis values in `x_units` (shape: `(O,)`; overrides
            `x_range` and `x_res`; default: `None`).
        x_units : str, optional
            x-axis units of the simulated spectrum (default: "thz").

        Notes
        -----
        One of `x_range`/`x_res` or `x` must be specified. If `x_range`
        is specified, `x` is generated using the `auto_x()` function
        with `x_res` if supplied.

        See Also
        --------
        auto_x : Function used to set `x` if not specified.
        """

        if x is not None:
            x = np_asarray_copy(x, copy=True, dtype=np.float64)

            if not np_check_shape(x, (None,)):
                raise ValueError("x must be an array_like with shape (O,).")
        else:
            if x_range is None or x_res is None:
                raise ValueError(
                    "x_range and x_res must be specified if x = None."
                )

            x_min, x_max = x_range

            if np.abs(x_max - x_min) < ZERO_TOLERANCE:
                raise ValueError(
                    "The mininum and maximum values specified by "
                    "x_range cannot be equal."
                )

            if x_min > x_max:
                x_min, x_max = x_max, x_min

            x = np.arange(x_min, x_max + x_res / 10.0, x_res, dtype=np.float64)

        x_units = str(x_units).lower()

        if x_units not in get_supported_frequency_units():
            raise ValueError(
                'x_units="{0}" is not a supported frequency/energy '
                "unit.".format(x_units)
            )

        self._x = x
        self._x_units = x_units

    @property
    def x(self):
        """numpy.ndarray : x values for simulated spectrum (shape:
        `(O,)`)."""
        return np_readonly_view(self._x)

    @property
    def x_units(self):
        """str : Code for x units compatible with the unit-handling
        routines in the `units` module."""
        return self._x_units

    @property
    def x_unit_text_label(self):
        """str or None : Label for x-axis units suitable for
        plain-text output."""

        if self._x_units is not None:
            return get_frequency_unit_text_label(self._x_units)

        return None

    @property
    def x_unit_plot_label(self):
        """str or None : Label for x-axis units suitable for plotting
        (may contain TeX strings.)"""

        if self._x_units is not None:
            return get_frequency_unit_plot_label(self._x_units)

        return None


# -----------------------------
# GammaPhononSpectrumBase class
# -----------------------------


class GammaPhononSpectrumBase(SpectrumBase):
    r"""Base class for simulated Gamma-point phonon spectra."""

    def __init__(
        self,
        freqs,
        lws,
        irrep_syms=None,
        x_range=None,
        x_res=None,
        x=None,
        x_units="thz",
        **kwargs
    ):
        """Create a new instance of the SpectrumBase class.

        Parameters
        ----------
        freqs : array_like
            Frequencies in THz (shape: `(N,)`).
        ints : array_like
            Band intensities in Ang^4 / sqrt(amu) (shape: `(N,)`).
        lws : array_like
            Linewidths in THz (shape: `(N,)`).
        irrep_syms : array_like or None, optional
            Irrep symbols of modes.
        x_range : tuple of float or None
            `(min, max)` range of x-axis of simulated spectrum in
            `units` (default: automatically determined).
        x_res : float or None, optional
            Resolution of simulated spectrum in `units` (default:
            automatically determined).
        x : array_like or None, optional
            x-axis values in `units` (overrides `x_range` and `x_res`).
        x_units : str, optional
            x-axis units of the simulated spectrum (default: `"thz"`).

        See Also
        --------
        `units.get_supported_frequency_units` : Get supported
            values of `units`.
        `auto_x` : Algorithm used to automatically determine `x_range`
            and/or `x_res` if not specified.
        """

        freqs = np_asarray_copy(freqs, dtype=np.float64)

        if not np_check_shape(freqs, (None,)):
            raise ValueError("freqs must be an array_like with shape (N,).")

        lws = np_asarray_copy(lws, dtype=np.float64)

        if not np_check_shape(lws, (len(freqs),)):
            raise ValueError("lws must be an array_like with shape (N,).")

        if irrep_syms is not None:
            irrep_syms = np.array(
                [str(sym).capitalize() for sym in irrep_syms], dtype=object
            )

            if not np_check_shape(irrep_syms, (len(freqs),)):
                raise ValueError(
                    "If supplied, irrep_syms must be an array_like "
                    "with shape (N,)."
                )
        else:
            irrep_syms = np.array(["None"] * len(freqs), dtype=object)

        # convert_frequency_units() will raise an error if units are
        # not supported.

        if x_units.lower() != "thz":
            freqs = convert_frequency_units(freqs, "thz", x_units)
            lws = convert_frequency_units(lws, "thz", x_units)

            if x is not None:
                x = convert_frequency_units(x, "thz", x_units)

        # auto_x will do nothing if x_range and x_res are already
        # specified.

        x_range, x_res = auto_x(freqs, lws, x_range=x_range, x_res=x_res)

        self._freqs = freqs
        self._lws = lws

        self._irrep_syms = irrep_syms

        super(GammaPhononSpectrumBase, self).__init__(
            x_range=x_range, x_res=x_res, x=x, x_units=x_units
        )

    @property
    def frequencies(self):
        """numpy.ndarray : Phonon frequencies (shape: `(N,)`)."""
        return np_readonly_view(self._freqs)

    @property
    def linewidths(self):
        """numpy.ndarray : Phonon linewidths (shape: `(N,)`)."""
        return np_readonly_view(self._lws)

    @property
    def irrep_symbols(self):
        """numpy.ndarray : Mode irrep symbols (shape: `(N,)`)."""
        return np_readonly_view(self._irrep_syms)
