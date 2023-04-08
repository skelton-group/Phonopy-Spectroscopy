# spectroscopy/spectrum.py


# ---------
# Docstring
# ---------

""" Core routines for simulating spectra. """


# -------
# Imports
# -------

import math

import numpy as np


# ---------
# Constants
# ---------

""" A "safe" multiplier of the linewidth (e.g. sigma for Gaussian and
gamma for Lorentzian functions) for padding grids used to evaluate peak
functions. """

SPECTRUM_PADDING_MULTIPLIER = 5.0

""" Minimum number of points for automatically determining spectrum 
resolutions. """

SPECTRUM_RESOLUTION_MIN_POINTS = 1000


# ---------
# Functions
# ---------

def gaussian(x, i, mu, sigma):
    """ Return G(x) = (i / (sigma * sqrt(2 * pi))) *
    exp(-1 * (x - mu) ** 2 / (2 * sigma ** 2)). """

    # Definition of the Gaussian function with unit area taken from
    # http://mathworld.wolfram.com/GaussianFunction.html.

    return ((i / (sigma * math.sqrt(2.0 * math.pi)))
            * np.exp(-1.0 * (x - mu) ** 2 / (2 * sigma ** 2)))

def lorentzian(x, i, x_0, gamma):
    """ Return L(x) = (i / pi) *
    ((0.5 * gamma) / ((x - x_0) ** 2) + (0.5 * gamma) ** 2). """

    # Definition of the Lorentzian function with unit area taken from
    # http://mathworld.wolfram.com/LorentzianFunction.html.

    return ((i / math.pi)
            * ((0.5 * gamma) / ((x - x_0) ** 2 + (0.5 * gamma) ** 2)))

def simulate_spectrum(
        frequencies, intensities, linewidths,
        spectrum_range=None, spectrum_resolution=None,
        instrument_broadening_width=None,
        instrument_broadening_shape='gaussian'
        ):

    """ Given a set of mode frequencies and intensities plus a nominal
    linewidth or set of linewidths, simulate a spectrum.
    
    Optionally, apply an instrument broadening by convolving with a
    fixed-width Gaussian or Lorentzian function.

    Arguments:
        frequencies -- mode frequencies.
        intensities -- mode spectroscopic intensities.
        linewidths -- nominal mode linewidth or list of mode linewidths,
            taken to be the full-width at half-maxima (FWHM) of
            Lorentzian peak profiles.

    Keyword arguments:
        spectrum_range -- range of frequencies over which to simulate
            the spectrum (defaults to approx.
            min(frequencies) - SPECTRUM_PADDING_MULTIPLIER * max(linewidths) ->
            max(frequencies) + SPECTRUM_PADDING_MULTIPLIER * max(linewidths)).
        spectrum_resolution -- frequency resolution of the spectrum
            (default: adjust to give at least
            SPECTRUM_RESOLUTION_MIN_POINTS data points, or a minimum
            resolution of 1, whichever is larger).
        instrument_broadening_width -- instrument broadening width
            (default: no instrument broadening).
        instrument_broadening_shape -- shape of instrument broadening
            ('gaussian' or 'lorentzian'; default: 'gaussian').

    Return value:
        A tuple of (spectrum_x, spectrum_y) data.

    Notes:
        If a min and max are specified with spectrum_range, they may be
        modified to "align" with the resolution (i.e. so that the
        frequency axis starts and finishes on an integer multiple of the
        resolution).
    """

    num_modes = len(frequencies)

    if len(intensities) != num_modes:
        raise Exception(
            "Error: The lengths of frequencies and intensities "
            "are inconsistent.")

    if np.isscalar(linewidths):
        linewidths = [float(linewidths)] * num_modes
    else:
        if len(linewidths) != num_modes:
            raise Exception(
                "Error: The lengths of frequencies and linewidths "
                "are inconsistent.")

    # Set a minimum and maximum.

    spectrum_min, spectrum_max = None, None

    if spectrum_range is not None:
        spectrum_min, spectrum_max = spectrum_range
    else:
        max_linewidth = max(linewidths)

        if num_modes == 1 and max_linewidth == 0.0:
            raise Exception(
                "Error: Cannot determine a spectrum range "
                "automatically - please specify manually using the "
                "spectrum_range argument.")

        # +/- 5 sigma/gamma should be sufficient for this.

        spectrum_min = (
            min(frequencies)
            - SPECTRUM_PADDING_MULTIPLIER * max_linewidth)
        
        spectrum_max = (
            max(frequencies)
            + SPECTRUM_PADDING_MULTIPLIER * max_linewidth)

        if spectrum_min == spectrum_max:
            raise Exception(
                "Error: The min and max specified in spectrum_range "
                "are the same.")

        if spectrum_max < spectrum_min:
            spectrum_max, spectrum_min = spectrum_min, spectrum_max

    # Set a resolution if required.

    if spectrum_resolution is None:
        power = math.log10(
            math.ceil(spectrum_max - spectrum_min)
            / SPECTRUM_RESOLUTION_MIN_POINTS)
        
        nominal_resolution = math.pow(10.0, math.floor(power))
        
        spectrum_resolution = min(nominal_resolution, 1.0)

    # If the spectrum range is being automatically determined, make sure
    # it is "aligned" with the resolution.

    if spectrum_range is None:
        spectrum_min = (
            spectrum_resolution
            * math.floor(spectrum_min / spectrum_resolution))
        
        spectrum_max = (
            spectrum_resolution
            * math.ceil(spectrum_max / spectrum_resolution))

    # If applying instrument broadening, calculate the convolution
    # kernel. Also expand the range of the spectrum so the the
    # convolution kernel does not produce boundary effects inside the
    # selected region.

    conv_num_points, conv_kernel = None, None

    if instrument_broadening_width is not None:
        # Calculating the convolution kernel to +/- 5 sigma/gamma should
        # be sufficient. According to
        # https://en.wikipedia.org/wiki/Gaussian_blur, +/- 3 sigma is
        # considered enough for Gaussian blurring kernels.

        conv_num_points = int(
            math.ceil(
                SPECTRUM_PADDING_MULTIPLIER
                * instrument_broadening_width
                / spectrum_resolution))

        conv_x = np.arange(
            -1.0 * conv_num_points * spectrum_resolution,
            conv_num_points + spectrum_resolution / 10.0,
            spectrum_resolution)

        if instrument_broadening_shape == 'gaussian':
            conv_kernel = gaussian(
                conv_x, 1.0, 0.0, instrument_broadening_width)
        elif instrument_broadening_shape == 'lorentzian':
            conv_kernel = lorentzian(
                conv_x, 1.0, 0.0, instrument_broadening_width)
        else:
            raise Exception(
                "Error: Unrecognised instrumentBroadeningShape '{0}'."
                .format(instrument_broadening_shape))

        conv_kernel /= conv_kernel.sum()

        spectrum_min -= conv_num_points * spectrum_resolution
        spectrum_max += conv_num_points * spectrum_resolution

    # Simulate spectrum.

    spectrum_x = np.arange(
        spectrum_min, spectrum_max + spectrum_resolution / 10.0,
        spectrum_resolution, dtype=np.float64)
    
    spectrum_y = np.zeros_like(spectrum_x, dtype=np.float64)

    for frequency, intensity, linewidth in zip(
            frequencies, intensities, linewidths):
        spectrum_y += lorentzian(
            spectrum_x, intensity, frequency, linewidth)

    # Apply instrument broadening if required.

    if conv_kernel is not None:
        # mode = 'valid' will cause np.convolve() to trim the first and
        # last conv_num_points data points from the spectrum.

        spectrum_x = spectrum_x[conv_num_points:-conv_num_points]
        spectrum_y = np.convolve(spectrum_y, conv_kernel, mode='valid')

        # Just in case my maths went wrong somewhere...

        assert len(spectrum_x) == len(spectrum_y)

    # Return spectrum.

    return (spectrum_x, spectrum_y)
