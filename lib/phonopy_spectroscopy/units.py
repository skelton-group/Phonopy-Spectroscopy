# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines for unit handling."""


# -------
# Imports
# -------


import numpy as np

from .constants import (
    BOHR_TO_ANG,
    THZ_TO_INV_CM,
    THZ_TO_EV,
    PLANCK_CONSTANT_EV,
    SPEED_OF_LIGHT,
)

from .utility.numpy_helper import np_expand_dims


# --------------
# Distance units
# --------------


_DISTANCE_UNIT_DATA = {
    "ang": {"from_ang": 1.0, "to_ang": 1.0},
    "bohr": {
        "from_ang": 1.0 / BOHR_TO_ANG,
        "to_ang": BOHR_TO_ANG,
    },
}


def get_supported_distance_units():
    """Return a list of supported distance units.

    Returns
    -------
    supported_units : list of str
        List of supported units.
    """

    return [k for k in _DISTANCE_UNIT_DATA.keys()]


def get_distance_unit_conversion_factor(unit_from, unit_to):
    """Return a multiplicative conversion factor to convert distances in
    `units_from` to `units_to`.

    Parameters
    ----------
    unit_from, unit_to : str
        Units to convert from and to.

    Returns
    -------
    conv_fac : float
        Conversion factor.
    """

    unit_from = unit_from.lower()

    if unit_from not in _DISTANCE_UNIT_DATA:
        raise ValueError('Unsupported unit_from="{0}".'.format(unit_from))

    unit_to = unit_to.lower()

    if unit_to not in _DISTANCE_UNIT_DATA:
        raise ValueError('Unsupported unit_to="{0}".'.format(unit_to))

    # Determine conversion factor.

    return (
        _DISTANCE_UNIT_DATA[unit_from]["to_ang"]
        * _DISTANCE_UNIT_DATA[unit_to]["from_ang"]
    )


def convert_distance_units(vals, unit_from, unit_to):
    """Convert distances in `unit_from` to `unit_to`.

    Parameters
    ----------
    vals : array_like or float
        Value(s) to convert.
    unit_from, unit_to : str
        Units to convert from and to.

    Returns
    -------
    conv_vals : array_like or float
        Converted values (shame shape as `vals`).
    """

    vals, n_dim_add = np_expand_dims(np.asarray(vals), (None,))

    vals = vals * get_distance_unit_conversion_factor(unit_from, unit_to)

    return vals if n_dim_add == 0 else vals[0]


# ---------------
# Frequency units
# ---------------


_FREQUENCY_UNIT_DATA = {
    "thz": {
        "text_label": "v / THz",
        "plot_label": r"$\nu$ / THz",
        "from_thz": 1.0,
        "to_thz": 1.0,
    },
    "rad_thz": {
        "text_label": r"\omega / rad THz",
        "plot_label": r"$\omega$ / rad THz",
        "from_thz": 2.0 * np.pi,
        "to_thz": 1.0 / (2.0 * np.pi),
    },
    "inv_cm": {
        "text_label": "v / cm^-1",
        "plot_label": r"$\bar{\nu}$ / cm$^{-1}$",
        "from_thz": THZ_TO_INV_CM,
        "to_thz": 1.0 / THZ_TO_INV_CM,
    },
    "ev": {
        "text_label": "E / eV",
        "plot_label": r"$E$ / eV",
        "from_thz": THZ_TO_EV,
        "to_thz": 1.0 / THZ_TO_EV,
    },
    "mev": {
        "text_label": "E / meV",
        "plot_label": r"$E$ / meV",
        "from_thz": THZ_TO_EV * 1000.0,
        "to_thz": 1.0 / (THZ_TO_EV * 1000.0),
    },
}


def get_supported_frequency_units():
    """Return a list of supported frequency units.

    Returns
    -------
    supported_units : list of str
        List of supported frequency units.
    """

    return [k for k in _FREQUENCY_UNIT_DATA.keys()]


def get_frequency_unit_text_label(unit):
    """Return a label for `unit` suitable for text output.

    Parameters
    ----------
    unit : str
        Frequency unit.

    Returns
    -------
    label : str
        Text label for `unit`.
    """

    if unit is not None:
        k = unit.lower()

        if k in _FREQUENCY_UNIT_DATA:
            return _FREQUENCY_UNIT_DATA[k]["text_label"]

    raise ValueError('Unsupported unit="{0}".'.format(unit))


def get_frequency_unit_plot_label(unit):
    """Return a label for `unit` suitable for plotting, which may
    contain TeX strings.

    Parameters
    ----------
    unit : str
        Frequency unit.

    Returns
    -------
    label : str
        Plot label for `unit`.
    """

    if unit is not None:
        k = unit.lower()

        if k in _FREQUENCY_UNIT_DATA:
            return _FREQUENCY_UNIT_DATA[k]["plot_label"]

    raise ValueError('Unsupported unit="{0}".'.format(unit))


def get_frequency_unit_conversion_factor(unit_from, unit_to):
    """Return a multiplicative convertion factor to convert frequencies
    in `units_from` to `units_to`.

    Parameters
    ----------
    unit_from, unit_to : str
        Units to convert from and to.

    Returns
    -------
    conv_fac : float
        Conversion factor.
    """

    unit_from = unit_from.lower()

    if unit_from not in _FREQUENCY_UNIT_DATA:
        raise ValueError('Unsupported unit_from "{0}".'.format(unit_from))

    unit_to = unit_to.lower()

    if unit_to not in _FREQUENCY_UNIT_DATA:
        raise ValueError('Unsupported unit_to "{0}".'.format(unit_to))

    # Determine conversion factor.

    return (
        _FREQUENCY_UNIT_DATA[unit_from]["to_thz"]
        * _FREQUENCY_UNIT_DATA[unit_to]["from_thz"]
    )


def convert_frequency_units(vals, unit_from, unit_to):
    """Convert frequencies in `unit_from` to `unit_to`.

    Parameters
    ----------
    vals : array_like or float
        Value(s) to convert.
    unit_from, unit_to : str
        Units to convert from and to.

    Returns
    -------
    conv_vals : array_like or float
        Converted values (shame shape as `vals`).
    """

    vals, n_dim_add = np_expand_dims(np.asarray(vals), (None,))

    vals = vals * get_frequency_unit_conversion_factor(unit_from, unit_to)

    return vals if n_dim_add == 0 else vals[0]


# -------------------------
# Miscellaneous conversions
# -------------------------


def nm_to_ev(w):
    """Convert wavelengths in nm to energies in eV.

    Parameters
    ----------
    w : float or array_like
        Wavelengths in nm.

    Returns
    -------
    e : float or array_like
        Energy in eV (same shape as `w`).
    """

    return (PLANCK_CONSTANT_EV * SPEED_OF_LIGHT) / (1.0e-9 * w)


def nm_to_thz(w):
    """Convert wavelengths in nm to frequencies in THz.

    Parameters
    ----------
    w : float or array_like
        Wavelengths in nm.

    Returns
    -------
    e : float or array_like
        Frequencies in THz (same shape as `w`).
    """

    return 1.0e-12 * SPEED_OF_LIGHT / (1.0e-9 * w)
