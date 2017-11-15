# SpectroscoPy/Plotting.py


# ---------
# Docstring
# ---------

""" Core routines for plotting spectra with Matplotlib. """


# -------
# Imports
# -------

import matplotlib as mpl;
import matplotlib.pyplot as plt;

from SpectroscoPy import Constants;
from SpectroscoPy import Utilities;


# ---------
# Constants
# ---------

""" Font size for plot text. """

FontSize = 8;

""" Linewidth for lines, axes, etc. """

Linewidth = 0.5;

""" Line colours for different plot types. """

LineColours = {
    'ir' : 'r', 'raman' : 'b'
    };

""" y-axis labels for different plot types. """

YAxisLabels = {
    'ir' : r"$I_{\mathrm{IR}}$",
    'raman' : r"$I_{\mathrm{Raman}}$"
    }


# ---------
# Functions
# ---------

def InitialiseMatplotlib():
    """ Boilerplate initialisation code for setting sensible Matplotlib defaults. """

    mpl.rc('font', **{ 'family' : 'serif', 'size' : FontSize, 'serif' : 'Times New Roman' });
    mpl.rc('mathtext', **{ 'fontset' : 'custom', 'rm' : 'Times New Roman', 'it' : 'Times New Roman:italic', 'bf' : 'Times New Roman:bold' });

    mpl.rc('axes', **{ 'labelsize' : FontSize, 'linewidth' : Linewidth });
    mpl.rc('lines', **{ 'linewidth' : Linewidth });

    tickParams = { 'major.width' : Linewidth, 'minor.width' : Linewidth, 'direction' : 'in' };

    mpl.rc('xtick', **tickParams);
    mpl.rc('ytick', **tickParams);

def PlotSpectrum(spectrum, filePath, frequencyUnits = None, plotType = None):
    """
    Plot a simulated spectrum.

    Arguments:
        spectrum -- a tuple of (frequencies, intensities, intensities_norm) data.
        filePath -- path to save plot to.

    Keyword arguments:
        frequencyUnits -- frequency units (default: Constants.DefaultFrequencyUnits).
        plotType -- adapts plot for specific type of spectra (supported values are 'ir' and 'raman').

    Notes:
        This routine plots the normalised spectrum, and expects that intensities_norm is in the range [0, 1].
    """

    spectrumX, spectrumY, spectrumYNorm = spectrum;

    numDataPoints = len(spectrumX);

    if len(spectrumY) != numDataPoints or len(spectrumYNorm) != numDataPoints:
        raise Exception("Error: The number of frequency and intensity points in spectra are inconsistent.");

    if frequencyUnits == None:
        frequencyUnits = Constants.DefaultFrequencyUnits;

    # Set default y-axis labels and line colours, and override if specified for the supplied plotType.

    plotType = plotType.lower();

    yAxisLabel = "Intensity";

    if plotType in YAxisLabels:
        yAxisLabel = YAxisLabels[plotType];

    lineColour = 'k';

    if plotType in LineColours:
        lineColour = LineColours[plotType];

    # Initialise Matplotlib.

    InitialiseMatplotlib();

    # Generate and save a plot.

    plt.figure(figsize = (8.6 / 2.54, 6.0 / 2.54));

    plt.plot(spectrumX, spectrumYNorm, color = lineColour);

    plt.xlabel(
        r"$\nu$ [{0}]".format(Utilities.GetFrequencyUnitLabel(frequencyUnits))
        );

    plt.ylabel(
        r"{0} [AU]".format(yAxisLabel)
        );

    plt.xlim(min(spectrumX), max(spectrumX));
    plt.ylim(0.0, 1.2);

    # For newer versions of Matplotlib (!).

    axes = plt.gca();

    axes.tick_params(
        direction = 'in', top = True, bottom = True, left = True, right = True
        );

    plt.tight_layout();

    plt.savefig(filePath, format = 'png', dpi = 300);
    plt.close();
