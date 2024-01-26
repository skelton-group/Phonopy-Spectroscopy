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

_IRREP_ACTIVITIES = {
    # Point group C_1.
    '1': {'ir': ["A"],
          'raman': ["A"],
          'all': ["A"]},
    # Point group C_i.
    '-1': {'ir': ["Au"],
           'raman': ["Ag"],
           'all': ["Ag", "Au"]},
    # Point group C_2.
    '2': {'ir': ["A", "B"],
          'raman': ["A", "B"],
          'all': ["A", "B"]},
    # Point group C_s.
    'm': {'ir': ["A'", "A''"],
          'raman': ["A'", "A''"],
          'all': ["A'", "A''"]
          },
    # Point group C_2h.
    '2/m': {'ir': ["Au", "Bu"],
            'raman': ["Ag", "Bg"],
            'all': ["Ag", "Au", "Bg", "Bu"]
            },
    # Point group D_2.
    '222': {'ir': ["B1", "B2", "B3"],
            'raman': ["A", "B1", "B2", "B3"],
            'all': ["A", "B1", "B2", "B3"]
            },
    # Point group C_2v.
    'mm2': {'ir': ["A1", "B1", "B2"],
            'raman': ["A1", "A2", "B1", "B2"],
            'all': ["A1", "A2", "B1", "B2"]
            },
    # Point group D_2h.
    'mmm': {'ir': ["B1u", "B2u", "B3u"],
            'raman': ["Ag", "B1g", "B2g", "B3g"],
            'all': ["Ag", "Au", "B1g", "B1u", "B2g", "B2u", "B3g", "B3u"]
            },
    # Point group C_4.
    '4': {'ir': ["A", "E"],
          'raman': ["A", "B", "E"],
          'all': ["A", "B", "E"]
          },
    # Point group S_4.
    '-4': {'ir': ["B", "E"],
           'raman': ["A", "B", "E"],
           'all': ["A", "B", "E"]
           },
    # Point group C_4h.
    '4/m': {'ir': ["Au", "Eu"],
            'raman': ["Ag", "Bg", "Eg"],
            'all': ["Ag", "Au", "Bg", "Bu", "Eg", "Eu"]
            },
    # Point group D_4.
    '422': {'ir': ["A2", "E"],
            'raman': ["A1", "B1", "B2", "E"],
            'all': ["A1", "A2", "B1", "B2", "E"]
            },
    # Point group C_4v.
    '4mm': {'ir': ["A1", "E"],
            'raman': ["A1", "B1", "B2", "E"],
            'all': ["A1", "A2", "B1", "B2", "E"]
            },
    # Point group D_2d.
    '-42m': {'ir': ["B2", "E"],
             'raman': ["A1", "B1", "B2", "E"],
             'all': ["A1", "A2", "B1", "B2", "E"]
             },
    # Point group D_4h.
    '4/mmm': {'ir': ["A2u", "Eu"],
              'raman': ["A1g", "B1g", "B2g", "Eg"],
              'all': ["A1g", "A1u", "A2g", "A2u", "B1g", "B1u", "B2g", "B2u",
                      "Eg", "Eu"]
              },
    # Point group C_3.
    '3': {'ir': ["A", "E"],
          'raman': ["A", "E"],
          'all': ["A", "E"]
          },
    # Point group C_3i.
    '-3': {'ir': ["Au", "Eu"],
           'raman': ["Ag", "Eg"],
           'all': ["Ag", "Au", "Eg", "Eu"]
           },
    # Point group D_3.
    '32': {'ir': ["A2", "E"],
           'raman': ["A1", "E"],
           'all': ["A1", "A2", "E"]
           },
    # Point group C_3v.
    '3m': {'ir': ["A1", "E"],
           'raman': ["A1", "E"],
           'all': ["A1", "A2", "E"]
           },
    # Point group D_3d.
    '-3m': {'ir': ["A2u", "Eu"],
            'raman': ["A1g", "Eg"],
            'all': ["A1g", "A1u", "A2g", "A2u", "Eg", "Eu"]
            },
    # Point group C_6.
    '6': {'ir': ["A", "E1"],
          'raman': ["A", "E1", "E2"],
          'all': ["A", "B", "E1", "E2"]
          },
    # Point group C_3h.
    '-6': {'ir': ["A''", "E'"],
           'raman': ["A'", "E'", "E''"],
           'all': ["A'", "A''", "E'", "E''"]
           },
    # Point group C_6h.
    '6/m': {'ir': ["Au", "E1u"],
            'raman': ["Ag", "E1g", "E2g"],
            'all': ["Ag", "Au", "Bg", "Bu", "E1g", "E1u", "E2g", "E2u"]
            },
    # Point group D_6.
    '622': {'ir': ["A2", "E1"],
            'raman': ["A1", "E1", "E2"],
            'all': ["A1", "A2", "B1", "B2", "E1", "E2"]
            },
    # Point group C_6v.
    '6mm': {'ir': ["A1", "E1"],
            'raman': ["A1", "E1", "E2"],
            'all': ["A1", "A2", "B1", "B2", "E1", "E2"]
            },
    # Point group D_3h.
    '-6m2': {'ir': ["A2''", "E'"],
             'raman': ["A1'", "E'", "E''"],
             'all': ["A1'", "A1''", "A2'", "A2''", "E'", "E''"]
             },
    # Point group D_6h.
    '6/mmm': {'ir': ["A2u", "E1u"],
              'raman': ["A1g", "E1g", "E2g"],
              'all': ["A1g", "A1u", "A2g", "A2u", "B1g", "B1u", "B2g", "B2u",
                      "E1g", "E1u", "E2g", "E2u"]
              },
    # Point group T.
    '23':  {'ir': ["T"],
            'raman': ["A", "E", "T"],
            'all': ["A", "E", "T"]
            },
    # Point group T_h.
    'm-3': {'ir': ["Tu"],
            'raman': ["Ag", "Eg", "Tg"],
            'all': ["Ag", "Au", "Eg", "Eu", "Tg", "Tu"]
            },
    # Point group O.
    '432': {'ir': ["T1"],
            'raman': ["A1", "E", "T2"],
            'all': ["A1", "A2", "E", "T1", "T2"]
            },
    # Point group T_d.
    '-43m': {'ir': ["T2"],
             'raman': ["A1", "E", "T2"],
             'all': ["A1", "A2", "E", "T1", "T2"]
             },
    # Point group O_h.
    'm-3m': {'ir': ["T1u"],
             'raman': ["A1g", "Eg", "T2g"],
             'all': ["A1g", "A1u", "A2g", "A2u", "Eg", "Eu", "T1g", "T2g",
                     "T1u", "T2u"]
             }
    }

def get_irrep_activities(point_group, spectrum_type):
    """ Given a point group, returns a pair of lists containing the
    active and inactive irreps for a given spectrum_type (supported
    values are 'ir' and 'raman').
    """
    
    point_group = point_group.lower()
    spectrum_type = spectrum_type.lower()

    if point_group not in _IRREP_ACTIVITIES:
        warnings.warn(
            "No activity data for point_group '{0}'.".format(point_group),
            RuntimeWarning)
        
        return None
    
    elif spectrum_type in ['ir', 'raman']:
        active_irreps = [
            symbol for symbol in _IRREP_ACTIVITIES[point_group][spectrum_type]
            ]
        
        inactive_irreps = [
            symbol for symbol in _IRREP_ACTIVITIES[point_group]['all']
                if symbol not in active_irreps
            ]
        
        return (active_irreps, inactive_irreps)
    
    else:
        raise Exception(
            "Error: Unknown spectrum_type '{0}'.".format(spectrum_type))
