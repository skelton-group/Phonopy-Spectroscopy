# SpectroscoPy/Constants.py


# ---------
# Docstring
# ---------

""" Defines constants used by other modules. """


# ----------
# Tolerances
# ----------

""" Zero tolerance. """

ZeroTolerance = 1.0e-8;


# ------------------
# Conversion Factors
# ------------------

""" Conversion factor from THz to inverse cm. """

THzToInvCm = 33.35641;

""" Conversion factor from THz to meV. """

THzToMeV = 4.13567;

""" Conversion factor from THz to _inverse_ microns. """

THzToInvUm = 1.0e-4 * THzToInvCm;


# -------------
# Default Units
# -------------

""" Default frequency units. """

DefaultFrequencyUnits = "?";

""" Default activity units. """

DefaultActivityUnits = "AU";

""" Default intensity units. """

DefaultIntensityUnits = "AU";


# -----------
# Unit Labels
# -----------

""" TeX labels for varuious supported frequency units. """

FrequencyUnitLabels = {
    'thz' : "THz",
    'inv_cm' : "cm$^{-1}$",
    'mev' : "meV",
    'um' : "$\mu$m"
    };

def GetFrequencyUnitLabel(frequencyUnits):
    """ Get a label for units frequencyUnits to be used in plots and output files. """

    # If frequencyUnits is one of the units supported by ConvertFrequencyUnits(), FrequencyUnitLabels should have a label for it.

    key = frequencyUnits.lower();

    if frequencyUnits == None:
        return DefaultFrequencyUnits;
    elif key in FrequencyUnitLabels:
        return FrequencyUnitLabels[key];
    else:
        return frequencyUnits;
