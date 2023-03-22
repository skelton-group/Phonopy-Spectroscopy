# SpectroscoPy/Plotting.py


# ---------
# Docstring
# ---------

""" Routines for plotting spectra with Matplotlib. """

# -------
# Imports
# -------

import matplotlib as mpl
import matplotlib.pyplot as plt

from SpectroscoPy import Constants

# ---------
# Constants
# ---------

""" Font size for plot text. """

FontSize = 8

""" Linewidth for lines, axes, etc. """

Linewidth = 0.5

""" Line colours for different plot types. """

LineColours = {
    'ir': 'r', 'raman': 'b'
}

""" y-axis labels for different plot types. """

YAxisLabels = {
    'ir': r"$I_{\mathrm{IR}}$",
    'raman': r"$I_{\mathrm{Raman}}$"
}


# ---------
# Functions
# ---------

def initialisematplotlib():
    """ Boilerplate initialisation code for setting sensible Matplotlib
    defaults. """

    mpl.rc('font', **{'family': 'serif', 'size': FontSize, 'serif':
                      'Times New Roman'}
           )
    mpl.rc('mathtext', **{'fontset': 'custom', 'rm': 'Times New Roman', 'it':
                          'Times New Roman:italic', 'bf': 'Times New '
                                                          'Roman:bold'}
           )

    mpl.rc('axes', **{'labelsize': FontSize, 'linewidth': Linewidth})
    mpl.rc('lines', **{'linewidth': Linewidth})

    tickparams = {'major.width': Linewidth, 'minor.width': Linewidth,
                  'direction': 'in'
                  }

    mpl.rc('xtick', **tickparams)
    mpl.rc('ytick', **tickparams)


def plotspectrum(spectrum, filepath, frequencyunits=None, plottype=None):
    """
    Plot a simulated spectrum.

    Arguments:
        spectrum -- a tuple of (frequencies, intensities, intensities_norm)
        data.
        filePath -- path to save plot to.

    Keyword arguments:
        frequencyUnits -- frequency units (default:
        Constants.DefaultFrequencyUnits).
        plotType -- adapts plot for specific type of spectra (supported
        values are 'ir' and 'raman').

    Notes:
        This routine plots the normalised spectrum, and expects that
        intensities_norm is in the range [0, 1].
    """

    spectrumx, spectrumy, spectrumynorm = spectrum

    numdatapoints = len(spectrumx)

    if len(spectrumy) != numdatapoints or len(spectrumynorm) != numdatapoints:
        raise Exception("Error: The number of frequency and intensity points "
                        "in spectra are inconsistent."
                        )

    if frequencyunits is None:
        frequencyunits = Constants.DefaultFrequencyUnits

    # Set default y-axis labels and line colours, and override if specified
    # for the supplied plotType.

    yaxislabel = "Intensity"

    if plottype in YAxisLabels:
        yaxislabel = YAxisLabels[plottype]

    linecolour = 'k'

    if plottype in LineColours:
        linecolour = LineColours[plottype]

    # Initialise Matplotlib.

    initialisematplotlib()

    # Generate and save a plot.

    plt.figure(figsize=(8.6 / 2.54, 6.0 / 2.54))

    plt.plot(spectrumx, spectrumynorm, color=linecolour)

    plt.xlabel(
        r"$\nu$ [{0}]".format(Constants.getfrequencyunitlabel(frequencyunits))
    )

    plt.ylabel(
        r"{0} [AU]".format(yaxislabel)
    )

    plt.xlim(min(spectrumx), max(spectrumx))
    plt.ylim(0.0, 1.2)

    # For newer versions of Matplotlib (!).

    axes = plt.gca()

    axes.tick_params(
        direction='in', top=True, bottom=True, left=True, right=True
    )

    plt.tight_layout()

    plt.savefig(filepath, format='png', dpi=300)
    plt.close()
