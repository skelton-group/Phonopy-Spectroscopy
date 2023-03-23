# spectroscopy/plotting.py


# ---------
# Docstring
# ---------

""" Routines for plotting spectra with Matplotlib. """


# -------
# Imports
# -------

import matplotlib as mpl
import matplotlib.pyplot as plt

from spectroscopy import constants


# ---------
# Constants
# ---------

""" Font size for plot text. """

FONT_SIZE = 8

""" Linewidth for lines, axes, etc. """

LINEWIDTH = 0.5

""" Line colours for different plot types. """

LINE_COLOURS = {
    'ir': 'r', 'raman': 'b'
    }

""" y-axis labels for different plot types. """

Y_AXIS_LABELS = {
    'ir': r"$I_{\mathrm{IR}}$",
    'raman': r"$I_{\mathrm{Raman}}$"
    }


# ---------
# Functions
# ---------

def initialise_matplotlib():
    """ Boilerplate initialisation code for setting sensible Matplotlib
    defaults. """

    mpl.rc(
        'font', **{'family': 'serif',
                   'size': FONT_SIZE,
                   'serif': 'Times New Roman'}
        )
    
    mpl.rc(
        'mathtext', **{'fontset': 'custom',
                       'rm': 'Times New Roman',
                       'it': 'Times New Roman:italic',
                       'bf': 'Times New Roman:bold'}
        )

    mpl.rc(
        'axes', **{'labelsize': FONT_SIZE,
                   'linewidth': LINEWIDTH}
        )
    
    mpl.rc('lines', **{'linewidth': LINEWIDTH})

    tickparams = {
        'major.width': LINEWIDTH,
        'minor.width': LINEWIDTH,
        'direction': 'in'
        }

    mpl.rc('xtick', **tickparams)
    mpl.rc('ytick', **tickparams)


def plot_spectrum(
        spectrum, file_path,
        frequency_units=None, plot_type=None
        ):
    """ Plot a simulated spectrum.

    Arguments:
        spectrum -- a tuple of
            (frequencies, intensities, intensities_norm) data.
        file_path -- path to save plot to.

    Keyword arguments:
        frequency_units -- frequency units (default:
        constants.DEFAULT_FREQUENCY_UNITS).
        plot_type -- adapts plot for specific type of spectra (supported
            values are 'ir' and 'raman').

    Notes:
        This routine plots the normalised spectrum, and expects that
        intensities_norm is in the range [0, 1].
    """

    spectrum_x, spectrum_y, spectrum_y_norm = spectrum

    num_data_points = len(spectrum_x)

    if (len(spectrum_y) != num_data_points
            or len(spectrum_y_norm) != num_data_points):
        raise Exception(
            "Error: The number of frequency and intensity points in "
            "spectra are inconsistent.")

    if frequency_units is None:
        frequency_units = constants.DEFAULT_FREQUENCY_UNITS

    # Set default y-axis labels and line colours, and override if
    # specified for the supplied plotType.

    y_axis_label = "Intensity"

    if plot_type in Y_AXIS_LABELS:
        y_axis_label = Y_AXIS_LABELS[plot_type]

    line_colour = 'k'

    if plot_type in LINE_COLOURS:
        line_colour = LINE_COLOURS[plot_type]

    # Initialise Matplotlib.

    initialise_matplotlib()

    # Generate and save a plot.

    plt.figure(
        figsize=(8.6 / 2.54, 6.0 / 2.54)
        )

    plt.plot(
        spectrum_x, spectrum_y_norm,
        color=line_colour
        )

    plt.xlabel(
        r"$\nu$ [{0}]".format(
            constants.get_frequency_unit_label(frequency_units)))

    plt.ylabel(r"{0} [AU]".format(y_axis_label))

    plt.xlim(min(spectrum_x), max(spectrum_x))
    plt.ylim(0.0, 1.2)

    # For newer versions of Matplotlib (!).

    axes = plt.gca()

    axes.tick_params(
        direction='in', top=True, bottom=True, left=True, right=True)

    plt.tight_layout()

    plt.savefig(file_path, format='png', dpi=300)
    
    plt.close()
