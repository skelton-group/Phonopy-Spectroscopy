# spectroscopy/plotting.py


# ---------
# Docstring
# ---------

""" Routines for plotting spectra with Matplotlib. """


# -------
# Imports
# -------

import math

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.cm import hsv
from matplotlib.gridspec import GridSpec

from spectroscopy.constants import ZERO_TOLERANCE


# ---------
# Constants
# ---------

""" Font size for plot text. """

FONT_SIZE = 8

""" Linewidth for lines, axes, etc. """

LINEWIDTH = 0.5


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


def plot_scalar_spectrum_or_spectra(
        frequencies, intensities_or_intensity_sets, frequency_unit_label,
        intensity_unit_label, file_path, normalise=False, legend_labels=None,
        line_colour_or_colours=None):
    """ Plot a simulated spectrum.

    Arguments:
        frequencies -- spectral frequencies.
        intensities_or_intensities_sets -- a single set of spectral
            intensities, or a list of sets of intensities.
        frequency_unit_label, intensity_unit_label -- labels for the
            frequency and intensity units.
        file_path -- output file.
    
    Keyword arguments:
        normalise -- if set, normalise the intensities before plotting
            (default: False).
        legend_labels -- optional legend labels (one per set of intensities).
        line_colour_or_line_colours -- optional line colour(s) to plot with
            (one per set of intensities).
    """

    dim = np.ndim(intensities_or_intensity_sets)

    if dim > 2:
        raise Exception(
            "intensities_or_intensity_sets must be either a single set "
            "of intensities (1D) or list of sets of intensities (2D).")

    num_data_points = len(frequencies)

    if dim == 1:
        intensities_or_intensity_sets = [intensities_or_intensity_sets]

    for intensities in intensities_or_intensity_sets:
        if len(intensities) != num_data_points:
            raise Exception(
                "The number of data points in one or more sets of "
                "intensities are inconsistent with the number of "
                "data points in frequencies.")

    if (legend_labels is not None
            and len(legend_labels) != len(intensities_or_intensity_sets)):
        raise Exception(
            "If supplied, legend_labels must provide a label for every "
            "spectrum.")

    if (line_colour_or_colours is not None
            and len(line_colour_or_colours) != len(
                intensities_or_intensity_sets)):
        raise Exception(
            "If supplied, line_colour_or_colours must provide a label "
            "for every spectrum.")

    # Determine normalisation factor, if required.

    norm = None

    if normalise:
        norm = np.max(intensities_or_intensity_sets)

    # If not specified, automatically generate a gradient of line
    # colours on the HSV scale.

    if line_colour_or_colours is None:
        start = 240.0 / 360.0

        if len(intensities_or_intensity_sets) == 1:
            line_colour_or_colours = [hsv(start % 1.0)]
        else:
            increment = (150.0 / 360.0) / (len(
                intensities_or_intensity_sets) - 1)

            line_colour_or_colours = [
                hsv((start + i * increment) % 1.0)
                for i in range(len(intensities_or_intensity_sets))
                ]

    # Initialise Matplotlib.

    initialise_matplotlib()

    # Generate and save a plot.

    plt.figure(figsize=(8.4 / 2.54, 6.5 / 2.54))

    for i, intensities in enumerate(intensities_or_intensity_sets):
        label, line_colour = None, None
        
        if legend_labels is not None:
            label = legend_labels[i]
        
        if line_colour_or_colours is not None:
            line_colour = line_colour_or_colours[i]
        
        if norm is not None:
            intensities = np.divide(intensities, norm)

        plt.plot(
            frequencies, intensities, label=label, color=line_colour)

    plt.xlabel(r"$\nu$ [{0}]".format(frequency_unit_label))
    
    plt.ylabel(r"$I(\nu)$ [{0}]".format(
        "AU" if normalise else intensity_unit_label))

    plt.xlim(min(frequencies), max(frequencies))

    if normalise:
        plt.ylim(0.0, 1.2)

    # If legend labels are provided, add a legend.

    if legend_labels is not None:
        legend = plt.legend(loc='best', frameon=True)
        legend.get_frame().set_edgecolor('k')
        legend.get_frame().set_linewidth(0.5)

    # For newer versions of Matplotlib (!).

    axes = plt.gca()

    axes.tick_params(
        direction='in', top=True, bottom=True, left=True, right=True)

    plt.tight_layout()

    plt.savefig(file_path, format='png', dpi=300)
    
    plt.close()

def plot_intensity_theta_polar(
        theta_vals, intensities, file_path, plot_label=None):
    """ Plot a set of intensities as a function of polarisation angle
    theta (0 -> 360 degrees) as a polar plot.

    Arguments:
        theta_vals -- polarisation angles (degrees).
        intensities -- band intensity as a function of theta.
        file_path -- file path to save to.
    
    Keyword arguments:
        plot_label -- optional label to be added to plot.
    """
    
    # Initialise Matplotlib.

    initialise_matplotlib()
    
    # Generate and save a polar plot.

    fig, ax = plt.subplots(
        figsize=(8.0 / 2.54, 7.0 / 2.54), subplot_kw={'projection': 'polar'})
    
    # Angles need to be in radians.

    theta_rad = [np.pi * val / 180.0 for val in theta_vals]

    norm = np.max(intensities)

    intensities_norm = intensities

    if norm > ZERO_TOLERANCE:
        intensities_norm = [val / norm for val in intensities]

    # If all the intensities are zero, ax.plot() draws a small circle
    # close to the origin, which looks strange -> better not to plot
    # anything.

    if norm >= ZERO_TOLERANCE:
        ax.plot(
            theta_rad, intensities_norm, label=plot_label, color='b',
            marker='o', mec='b', mfc='none', ms=2.0, markeredgewidth=LINEWIDTH)

    ax.set_rmax(1.2)
    ax.set_rticks([])
    
    if plot_label is not None:
        ax.text(
            0.05, 0.95, plot_label, ha='left', va='top',
            transform=fig.transFigure,
            bbox={'facecolor': 'none', 'edgecolor': 'k',
                  'linewidth': LINEWIDTH})

    ax.grid(color=(0.75, 0.75, 0.75), linewidth=LINEWIDTH)
    
    plt.tight_layout()

    plt.savefig(file_path, format='png', dpi=300)
    
    plt.close()  

def plot_2d_polarised_raman_spectrum(
        frequencies, theta_vals, intensities, frequency_unit_label,
        theta_unit_label, file_path, normalise=False):
    """ Plot a "2D" polarised Raman spectrum as a function of frequency
    and rotation angle theta.

    Arguments:
        frequencies -- spectrum frequencies (x-axis).
        theta_vals -- values of theta (y-axis).
        intensities -- a 2D array of intensities as a function of
            freuqency and theta.
        frequency_unit_label -- x-axis label.
        theta_unit_label -- y-axis label.
        file_path -- file to save to.
    
    Keyword arguments:
        normalise -- if set, normalise the intensities before plotting
            (default: False).
    """

    # Initialise Matplotlib.

    initialise_matplotlib()

    plt.figure(figsize=(8.6 / 2.54, 7.5 / 2.54))

    grid_spec = GridSpec(6, 1)

    cbar_axes = plt.subplot(grid_spec[0, :])
    plot_axes = plt.subplot(grid_spec[1:, :])

    # If required, normalise spectrum.

    if normalise:
        intensities = (np.asarray(intensities, dtype=np.float64)
                       / np.max(intensities))

    mappable = plot_axes.pcolormesh(
        frequencies, theta_vals, intensities, cmap='jet')

    plot_axes.set_xlabel(r"$\nu$ [{0}]".format(frequency_unit_label))
    plot_axes.set_ylabel(r"$\theta$ [{0}]".format(theta_unit_label))

    plt.colorbar(mappable, cax=cbar_axes, orientation='horizontal')

    # This works well for theta from 0 -> 360, which is what this
    # routine will generally be used to do.

    plot_axes.set_yticks(np.linspace(theta_vals[0], theta_vals[-1], 9))

    plt.tight_layout()

    plt.savefig(file_path, dpi=300)

    plt.close()
