# spectroscopy/constants.py


# ---------
# Docstring
# ---------

""" Defines constants used by other modules. """


# -------
# Imports
# -------

import warnings


# ----------
# Tolerances
# ----------

""" Zero tolerance. """

ZERO_TOLERANCE = 1.0e-8


# ------------------
# Conversion Factors
# ------------------

""" Conversion factor from THz to inverse cm. """

THZ_TO_INV_CM = 33.35641

""" Conversion factor from THz to meV. """

THZ_TO_MEV = 4.13567

""" Conversion factor from THz to _inverse_ microns. """

THZ_TO_INV_UM = 1.0e-4 * THZ_TO_INV_CM


# -------------
# Default Units
# -------------

""" Default frequency units. """

DEFAULT_FREQUENCY_UNITS = "?"

""" Default activity units. """

DEFAULT_ACTIVITY_UNITS = "AU"

""" Default intensity units. """

DEFAULT_INTENSITY_UNITS = "AU"


# -----------
# Unit Labels
# -----------

""" TeX labels for various supported frequency units. """

_FREQUENCY_UNIT_LABELS = {
    'thz': "THz",
    'inv_cm': "cm$^{-1}$",
    'mev': "meV",
    'um': r"$\mu$m"
    }


def get_frequency_unit_label(frequency_units):
    """ Get a label for units frequency_units to be used in plots and
    output files. """

    if frequency_units is None:
        return DEFAULT_FREQUENCY_UNITS
    
    k = frequency_units.lower()

    if k in _FREQUENCY_UNIT_LABELS:
        return _FREQUENCY_UNIT_LABELS[k]
    else:
        return frequency_units


# ------------
# Group Theory
# ------------

_ACTIVE_IRREPS = {
    '1': {},
    '-1': {'IR': ["Au"],
           'Raman': ["Ag"]
           },
    '2': {'IR': ["A", "B"],
          'Raman': ["A", "B"]
          },
    'm': {'IR': ["A'", "A''"],
          'Raman': ["A'", "A''"]
          },
    '2/m': {'IR': ["Au", "Bu"],
            'Raman': ["Ag", "Bg"]
            },
    '222': {'IR': ["B1", "B2", "B3"],
            'Raman': ["A", "B1", "B2", "B3"]
            },
    'mm2': {'IR': ["A1", "B1", "B2"],
            'Raman': ["A1", "A2", "B1", "B2"]
            },
    'mmm': {'IR': ["B1u", "B2u", "B3u"],
            'Raman': ["Ag", "B1g", "B2g", "B3g"]
            },
        
    '4': {'IR': ["A", "1E", "2E"],
          'Raman': ["A", "B", "1E", "2E"]
          },
    '-4': {'IR': ["B", "1E", "2E"],
           'Raman': ["A", "B", "1E", "2E"]
           },
    '4/m': {'IR': ["Au", "1Eu", "2Eu"],
            'Raman': ["Ag", "Bg", "1Eg", "2Eg"]
            },
    '422': {'IR': ["A2", "E"],
            'Raman': ["A1", "B1", "B2", "E"]
            },
    '4mm': {'IR': ["A1", "E"],
            'Raman': ["A1", "B1", "B2", "E"]
            },
    '-42m': {'IR': ["B2", "E"],
             'Raman': ["A1", "B1", "B2", "E"]
             },
    '4/mmm': {'IR': ["A2u", "E"],
              'Raman': ["A1g", "B1g", "B2g", "Eg"]
              },
    '3': {'IR': ["A", "1E", "2E"],
          'Raman': ["A", "1E", "2E"]
          },
    
    '-3': {'IR': ["Au", "1Eu", "2Eu"],
           'Raman': ["Ag", "1Eg", "2Eg"]
           },
    '32': {'IR': ["A2", "E"],
           'Raman': ["A1", "E"]
           },
    '3m': {'IR': ["A1", "E"],
           'Raman': ["A1", "E"]
           },
    '-3m': {'IR': ["A2u", "Eu"],
            'Raman': ["A1g", "Eg"]
            },
    '6': {'IR': ["A", "2E1", "1E1"],
          'Raman': ["A", "1E2", "2E2", "2E1", "1E1"]
          },
    '-6': {'IR': ["A''", "2E'", "1E'"],
           'Raman': ["A'", "2E'", "1E'", "2E''", "1E''"]
           },
    '6/m': {'IR': ["Au", "2E1u", "1E1u"],
            'Raman': ["Ag", "1E2g", "2E2g", "2E1g", "1E1g"]
            },
    '622': {'IR': ["A2", "E1"],
            'Raman': ["A1", "E2", "E1"]
            },
        
    '6mm': {'IR': ["A1", "E1"],
            'Raman': ["A1", "E2", "E1"]
            },
    '-6m2': {'IR': ["A''2", "E'"],
             'Raman': ["A'1", "E'", "E''"]
             },
    '6/mmm': {'IR': ["A2u", "E1u"],
              'Raman': ["A1g", "E2g", "E1g"]
              },
    '23':  {'IR': ["T"],
            'Raman': ["A", "1E", "2E", "T"]
            },
    'm-3': {'IR': ["Tu"],
            'Raman': ["Ag", "1Eg", "2Eg", "Tg"]
            },
    '432': {'IR': ["T1"],
            'Raman': ["A1", "E", "T2"]
            },
    '-43m': {'IR': ["T2"],
             'Raman': ["A1", "E", "T2"]
             },
    'm-3m': {'IR': ["T1u"],
             'Raman': ["A1g", "Eg", "T2g"]
             }
    }


def get_active_irreps(point_group, spectrum_type):
    """ Given a point group, returns a list of irrep symbols that may be
    active for the given spectrum_type (supported values are 'ir' and
    'raman').
    """
    
    point_group = point_group.lower()
    spectrum_type = spectrum_type.lower()

    if point_group not in _ACTIVE_IRREPS:
        warnings.warn(
            "No activity data for point_group '{0}'.".format(point_group),
            RuntimeWarning)
        
        return None
    
    elif spectrum_type in ['ir', 'raman']:
        return [
            symbol for symbol in _ACTIVE_IRREPS[point_group][spectrum_type]
            ]
    
    else:
        raise Exception(
            "Error: Unknown spectrum_type '{0}'.".format(spectrum_type))
