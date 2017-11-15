# SpectroscoPy/Constants.py


# ---------
# Docstring
# ---------

""" Defines constants used by other modules. """


# ---------
# Constants
# ---------

""" Conversion factors from THz to other supported frequency units. """

FrequencyConversionFactors = {
    # These conversion factors were taken from http://halas.rice.edu/conversions.

    'inv_cm' : 33.35641,
    'mev' : 4.13567
    };

""" Labels for the frequency units supported by the ConvertFrequencyUnits() routine. """

FrequencyUnitLabels = {
    'thz' : "THz",
    'inv_cm' : "cm$^{-1}$",
    'mev' : "meV"
    };

""" Default frequency units. """

DefaultFrequencyUnits = "?";

""" Default intensity units. """

DefaultIntensityUnits = "AU";
