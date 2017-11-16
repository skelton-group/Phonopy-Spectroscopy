# SpectroscoPy/Constants.py


# ---------
# Docstring
# ---------

""" Defines constants used by other modules. """


# ---------
# Constants
# ---------

""" Conversion factor from THz to inverse cm. """

THzToInvCm = 33.35641;

""" Conversion factor from THz to meV. """

THzToMeV = 4.13567;

""" Conversion factor from THz to _inverse_ microns. """

THzToInvUm = 1.0e-4 * THzToInvCm;

""" Labels for varuious supported frequency units. """

FrequencyUnitLabels = {
    'thz' : "THz",
    'inv_cm' : "cm$^{-1}$",
    'mev' : "meV",
    'um' : "$\mu$m"
    };

""" Default frequency units. """

DefaultFrequencyUnits = "?";

""" Default intensity units. """

DefaultIntensityUnits = "AU";
