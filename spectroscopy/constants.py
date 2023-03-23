# spectroscopy/constants.py


# ---------
# Docstring
# ---------

""" Defines constants used by other modules. """


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

FREQUENCY_UNIT_LABELS = {
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

    if k in FREQUENCY_UNIT_LABELS:
        return FREQUENCY_UNIT_LABELS[k]
    else:
        return frequency_units
