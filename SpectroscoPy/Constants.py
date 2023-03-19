# SpectroscoPy/Constants.py


# ---------
# Docstring
# ---------

""" Defines constants used by other modules. """


# ----------
# Tolerances
# ----------

""" Zero tolerance. """

ZeroTolerance = 1.0e-8


# ------------------
# Conversion Factors
# ------------------

""" Conversion factor from THz to inverse cm. """

THzToInvCm = 33.35641

""" Conversion factor from THz to meV. """

THzToMeV = 4.13567

""" Conversion factor from THz to _inverse_ microns. """

THzToInvUm = 1.0e-4 * THzToInvCm


# -------------
# Default Units
# -------------

""" Default frequency units. """

DefaultFrequencyUnits = "?"

""" Default activity units. """

DefaultActivityUnits = "AU"

""" Default intensity units. """

DefaultIntensityUnits = "AU"


# -----------
# Unit Labels
# -----------

""" TeX labels for various supported frequency units. """

FrequencyUnitLabels = {
    # E203: no whitespace before colon.
    'thz': "THz",
    'inv_cm': "cm$^{-1}$",
    'mev': "meV",
    # W605: invalid escape sequence \m.
    # added r before the quotes to make it a raw string instead of a
    # string literal.
    'um': r"$\mu$m"
    }


# functions name lowercase and two spaces leeway.
def getfrequencyunitlabel(frequencyunits):
    """ Get a label for units frequencyUnits to be used in plots and output
    files. """

    # If frequencyUnits is one of the units supported by
    # ConvertFrequencyUnits(), FrequencyUnitLabels should have a label for it.

    # undefined name 'frequencyUnits' because I changed the argumnent?
    # no longer need key?
    key = frequencyUnits.lower()

    # changed 'if frequencyUnits' to 'if key'.
    if key is None:
        # should DefaultFrequencyunits be DefaultFrequencyUnits as above?
        # Also used in Plotting.py
        return DefaultFrequencyUnits
    elif key in FrequencyUnitLabels:
        return FrequencyUnitLabels[key]
    else:
        return frequencyunits
