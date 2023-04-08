# spectroscopy/units.py


# ---------
# Docstring
# ---------

""" Contains constants and routines for unit handling. """


# -------
# Imports
# -------

import numpy as np


# --------------
# Distance Units
# --------------

""" Conversion factor from Bohr radii to Angstroms. """

BOHR_TO_ANG = 0.529177249

""" Data for supported distance units in the form of dictionaries
with the following keys:
    'text_label' -- a label suitable for text output.
"""

_DISTANCE_UNIT_DATA = {
    'ang': {'text_label': "Ang"},
    'bohr': {'text_label': "Bohr"}
    }

def get_distance_unit_text_label(unit):
    """ Return a label for unit suitable for text outputs. """

    if unit is not None:
        k = unit.lower()
        
        if k in _DISTANCE_UNIT_DATA:
            return _DISTANCE_UNIT_DATA[k]['text_label']
    
    raise Exception("Unsupported unit '{0}'.".format(unit))


# ---------------
# Frequency Units
# ---------------

""" Conversion factor from THz to inverse cm. """

THZ_TO_INV_CM = 33.35641

""" Conversion factor from THz to meV. """

THZ_TO_MEV = 4.13567

""" Data for supported frequency units in the form of dictionaries
with the following keys:
    'text_label' -- a label suitable for text output.
    'plot_label' -- a label suitable for plots, which may  contain TeX
        formatting.
    'from_thz' -- multiplicative conversion factor from SI units of THz.
    'to_thz' -- multiplicative conversion factor to SI units of THz.
"""

_FREQUENCY_UNIT_DATA = {
    'thz': {'text_label': "THz", 'plot_label' : "THz",
        'from_thz': 1.0, 'to_thz': 1.0},
    'inv_cm': {'text_label': "cm^-1", 'plot_label': "cm$^{-1}$",
        'from_thz': THZ_TO_INV_CM, 'to_thz': 1.0 / THZ_TO_INV_CM},
    'mev': {'text_label': "meV", 'plot_label': "meV",
        'from_thz': THZ_TO_MEV, 'to_thz': 1.0 / THZ_TO_MEV}
    }

def get_supported_frequency_units():
    """ Return a list of supported units. """
    
    return [k for k in _FREQUENCY_UNIT_DATA.keys()]

def get_frequency_unit_type(unit):
    """ Return the type of unit ('frequency', 'energy' or
    'wavelength'). """

    if unit is not None:
        k = unit.lower()

        if k in _FREQUENCY_UNIT_DATA:
           return _FREQUENCY_UNIT_DATA[k]['type']
    
    raise Exception("Unsupported unit '{0}'.".format(unit))

def get_frequency_unit_text_label(unit):
    """ Return a label for unit suitable for text outputs. """

    if unit is not None:
        k = unit.lower()
        
        if k in _FREQUENCY_UNIT_DATA:
            return _FREQUENCY_UNIT_DATA[k]['text_label']
    
    raise Exception("Unsupported unit '{0}'.".format(unit))

def get_frequency_unit_plot_label(unit):
    """ Return a label for unit suitable for plots (may contain TeX
    strings). """

    if unit is not None:
        k = unit.lower()
        
        if k in _FREQUENCY_UNIT_DATA:
            return _FREQUENCY_UNIT_DATA[k]['plot_label']
    
    raise Exception("Unsupported unit '{0}'.".format(unit))

def get_frequency_unit_conversion_factor(unit_from, unit_to):
    """ Obtain a multiplicative conversion factor to convert frequencies
    frequencies in units_to to units_from. """

    unit_from = unit_from.lower()
    
    if unit_from not in _FREQUENCY_UNIT_DATA:
        raise Exception("Unsupported unit_from '{0}'.".format(unit_from))
    
    unit_to = unit_to.lower()

    if unit_to not in _FREQUENCY_UNIT_DATA:
        raise Exception("Unsupported unit_to '{0}'.".format(unit_to))

    # Determine conversion factor.

    return (_FREQUENCY_UNIT_DATA[unit_from]['to_thz']
            * _FREQUENCY_UNIT_DATA[unit_to]['from_thz'])

def convert_frequency_units(frequency_or_frequencies, unit_from, unit_to):
    """ Convert a (scalar) frequency or list of frequencies in unit_from
    to unit_to. """

    conversion_factor = get_frequency_unit_conversion_factor(
        unit_from, unit_to)
    
    if np.isscalar(frequency_or_frequencies):
        return frequency_or_frequencies * conversion_factor
    else:
        return [f * conversion_factor for f in frequency_or_frequencies]
