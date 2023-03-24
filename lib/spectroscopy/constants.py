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
    }

def get_active_irreps(point_group, spectrum_type):
    """ Given a point group, returns a list of irrep symbols that may be
    active for the given spectrum_type.
    """
    
    if point_group not in _ACTIVE_IRREPS:
        warnings.warn(
            "No activity data for point_group '{0}'.".format(point_group),
            RuntimeWarning)
        
        return None
    
    elif spectrum_type in ['ir', 'raman']:
        return [
            symbol for symbol in _ACTIVE_IRREPS[point_group]
                if symbol[spectrum_type]
            ]
    
    else:
        # Assume all irreps can be active for unknown spectrum_types.

        warnings.warn(
            "Unknown spectrum_type '{0}' -> assuming all irreps can be "
            "active.".format(spectrum_type), RuntimeWarning)

        return [symbol for symbol in _ACTIVE_IRREPS[point_group]]
