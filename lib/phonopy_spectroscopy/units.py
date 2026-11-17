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
"""Dictionary of supported distance units and conversion factors
to/from the internal units of Ang."""


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
        raise ValueError(f'Unsupported unit_from "{unit_from}".')

    unit_to = unit_to.lower()

    if unit_to not in _DISTANCE_UNIT_DATA:
        raise ValueError(f'Unsupported unit_to "{unit_to}".')

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
        Converted values (same shape as `vals`).
    """

    factor = get_distance_unit_conversion_factor(unit_from, unit_to)
    return np.asanyarray(vals, dtype=np.float64) * factor


# ----------------------
# Frequency/energy units
# ----------------------

_FREQUENCY_ENERGY_UNIT_DATA = {
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
"""Dictionary of supported frequency/energy units and conversion factors
to/from the internal units of THz."""


def get_supported_frequency_energy_units():
    """Return a list of supported frequency/energy units.

    Returns
    -------
    supported_units : numpy.ndarray
        Supported frequency/energy units.
    """

    return np.array([k for k in _FREQUENCY_ENERGY_UNIT_DATA.keys()], dtype=str)


def get_frequency_energy_unit_text_label(unit):
    """Return a label for `unit` suitable for text output.

    Parameters
    ----------
    unit : str
        Frequency/energy unit.

    Returns
    -------
    label : str
        Text label for `unit`.
    """

    if unit is not None:
        k = unit.lower()

        if k in _FREQUENCY_ENERGY_UNIT_DATA:
            return _FREQUENCY_ENERGY_UNIT_DATA[k]["text_label"]

    raise ValueError(f'Unsupported frequency/energy unit "{unit}".')


def get_frequency_energy_unit_plot_label(unit):
    """Return a label for `unit` suitable for plotting, which may
    contain TeX strings.

    Parameters
    ----------
    unit : str
        Frequency/energy unit.

    Returns
    -------
    label : str
        Plot label for `unit`.
    """

    if unit is not None:
        k = unit.lower()

        if k in _FREQUENCY_ENERGY_UNIT_DATA:
            return _FREQUENCY_ENERGY_UNIT_DATA[k]["plot_label"]

    raise ValueError(f'Unsupported frequency/energy unit "{unit}".')


def get_frequency_energy_unit_conversion_factor(unit_from, unit_to):
    """Return a multiplicative convertion factor to convert
    frequencies/energies in `units_from` to `units_to`.

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

    if unit_from not in _FREQUENCY_ENERGY_UNIT_DATA:
        raise ValueError(f'Unsupported unit_from "{unit_from}".')

    unit_to = unit_to.lower()

    if unit_to not in _FREQUENCY_ENERGY_UNIT_DATA:
        raise ValueError(f'Unsupported unit_to "{unit_to}".')

    # Determine conversion factor.

    return (
        _FREQUENCY_ENERGY_UNIT_DATA[unit_from]["to_thz"]
        * _FREQUENCY_ENERGY_UNIT_DATA[unit_to]["from_thz"]
    )


def convert_frequency_energy_units(vals, unit_from, unit_to):
    """Convert frequencies/energies in `unit_from` to `unit_to`.

    Parameters
    ----------
    vals : array_like or scalar
        Value(s) to convert.
    unit_from, unit_to : str
        Units to convert from and to.

    Returns
    -------
    conv_vals : numpy.ndarray
        Converted values (same shape as `vals`).
    """

    factor = get_frequency_energy_unit_conversion_factor(unit_from, unit_to)
    return factor * np.asanyarray(vals)


# -------------------------
# Miscellaneous conversions
# -------------------------


def ev_to_nm(e):
    """Convert energies in eV to wavelengths in nm.

    Parameters
    ----------
    e : array_like or scalar
        Energies in eV.

    Returns
    -------
    w : numpy.ndarray
        Wavelengths in nm (same shape as `e`).
    """

    return 1.0e9 * (PLANCK_CONSTANT_EV * SPEED_OF_LIGHT) / np.asanyarray(e)


def nm_to_ev(w):
    """Convert wavelengths in nm to energies in eV.

    Parameters
    ----------
    w : array_like or scalar
        Wavelengths in nm.

    Returns
    -------
    e : numpy.ndarray
        Energy in eV (same shape as `w`).
    """

    return (PLANCK_CONSTANT_EV * SPEED_OF_LIGHT) / (1.0e-9 * np.asanyarray(w))


def nm_to_thz(w):
    """Convert wavelengths in nm to frequencies in THz.

    Parameters
    ----------
    w : array_like or scalar
        Wavelengths in nm.

    Returns
    -------
    e : numpy.ndarray
        Frequencies in THz (same shape as `w`).
    """

    return 1.0e-12 * SPEED_OF_LIGHT / (1.0e-9 * np.asanyarray(w))
