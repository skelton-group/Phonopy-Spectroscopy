# SpectroscoPy/Spectrum.py


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

""" To ensure the maximum is included in grids, np.arange() is called as 
np.arange(min, max + RangeStepMultiplier * step, step). """

RangeStepMultiplier = 1.0e-5

""" A "safe" multiplier of the linewidth (e.g. sigma for Gaussian and gamma 
for Lorentzian functions) for padding grids used to evaluate peak functions. 
"""

SpectrumPaddingMultiplier = 5.0

""" Minimum number of points for automatically determining spectrum 
resolutions. """

SpectrumResolutionMinPoints = 1000


# ---------
# Functions
# ---------

def gaussian(x, i, mu, sigma):
    """ Return G(x) = (i / (sigma * sqrt(2 * pi))) * exp(-1 * (x - mu) ** 2
    / (2 * sigma ** 2)). """

    # Definition of the Gaussian function with unit area taken from
    # http://mathworld.wolfram.com/GaussianFunction.html.

    return (i / (sigma * math.sqrt(2.0 * math.pi))) \
        * np.exp(-1.0 * (x - mu) ** 2 / (2 * sigma ** 2))


def lorentzian(x, i, x0, gamma):
    """ Return L(x) = (i / pi) * ((0.5 * gamma) / ((x - x_0) ** 2)
    + (0.5 * gamma) ** 2). """

    # Definition of the Lorentzian function with unit area taken from
    # http://mathworld.wolfram.com/LorentzianFunction.html.

    return (i / math.pi) * ((0.5 * gamma) / ((x - x0) ** 2 + (0.5 * gamma)
                                             ** 2))


def simulatespectrum(
    frequencies, intensities, linewidths,
    spectrumrange=None, spectrumresolution=None,
    instrumentbroadening=None, instrumentbroadeningshape='gaussian'
):

    """
    Given a set of mode frequencies and intensities plus a nominal linewidth or
     set of linewidths, simulate a spectrum.
    Optionally, apply an instrument broadening by convolving with a fixed-width
     Gaussian or Lorentzian function.

    Arguments:
        frequencies -- mode frequencies.
        intensities -- mode spectroscopic intensities.
        linewidths -- nominal mode linewidth or list of mode linewidths, taken
        to be the full-width at half-maxima (FWHM) of Lorentzian peak profiles.

    Keyword arguments:
        spectrumRange -- range of frequencies over which to simulate the
        spectrum (defaults to approx. min(frequencies)
        - SpectrumPaddingMultiplier * max(linewidths) -> max(frequencies)
        + SpectrumPaddingMultiplier * max(linewidths)).
        spectrumResolution -- frequency resolution of the spectrum (default:
        adjust to give at least SpectrumResolutionMinPoints data points, or a
        minimum resolution of 1, whichever is larger).
        instrumentBroadening -- instrument broadening width (default: no
        instrument broadening).
        instrumentBroadeningShape -- shape of instrument broadening ('gaussian'
        or 'lorentzian'; default: 'gaussian').

    Return value:
        A tuple of (spectrumX, spectrumY, spectrumYNorm) data.

    Notes:
        If a min and max are specified with spectrumRange, they may be modified
        to "align" with the resolution (i.e. so that the frequency axis starts
        and finished on an integer multiple of the resolution).
    """

    nummodes = len(frequencies)

    if len(intensities) != nummodes:
        raise Exception("Error: The lengths of frequencies and intensities "
                        "are inconsistent."
                        )

    try:
        # Assume linewidths is a scalar linewidth and repeat it to a list.

        linewidth = float(linewidths)
        linewidths = [linewidth] * nummodes
    except TypeError:
        # If linewidths is a list or multi-element NumPy array, casting to a
        # float will raise a TypeError.

        if len(linewidths) != nummodes:
            raise Exception("Error: The lengths of frequencies and linewidths "
                            "are inconsistent."
                            )

    # Set a minimum and maximum.

    spectrummin, spectrummax = None, None

    if spectrumrange is None:
        maxlinewidth = max(linewidths)

        if nummodes == 1 and maxlinewidth == 0.0:
            raise Exception("Error: Cannot determine a spectrum range "
                            "automatically - please specify manually using "
                            "the spectrumRange argument."
                            )

        # +/- 5 sigma/gamma should be sufficient for this.

        spectrummin = min(frequencies) - SpectrumPaddingMultiplier \
            * maxlinewidth
        spectrummax = max(frequencies) + SpectrumPaddingMultiplier \
            * maxlinewidth
    else:
        spectrummin, spectrummax = spectrumrange

        if spectrummin == spectrummax:
            raise Exception("Error: The min and max specified in "
                            "spectrumRange are the same."
                            )

        if spectrummax < spectrummin:
            spectrummax, spectrummin = spectrummin, spectrummax

    # Set a resolution if required.

    if spectrumresolution is None:
        nominalresolution = math.pow(
            10.0, math.floor(math.log10(math.ceil(spectrummax - spectrummin)
                                        / SpectrumResolutionMinPoints))
            )

        spectrumresolution = min(nominalresolution, 1.0)

    # If the spectrum range is being automatically determined, make sure it is
    # "aligned" with the resolution.

    if spectrumrange is None:
        spectrummin = spectrumresolution * math.floor(spectrummin /
                                                      spectrumresolution)
        spectrummax = spectrumresolution * math.ceil(spectrummax /
                                                     spectrumresolution)

    # If applying instrument broadening, calculate the convolution kernel.
    # Also expand the range of the spectrum so the the convolution kernel
    # does not produce boundary effects inside the selected region.

    convnumpoints, convkernel = None, None

    if instrumentbroadening is not None:
        # Calculating the convolution kernel to +/- 5 sigma/gamma should be
        # sufficient.
        # According to https://en.wikipedia.org/wiki/Gaussian_blur, +/- 3
        # sigma is considered enough for Gaussian blurring kernels.

        convnumpoints = int(
            math.ceil(SpectrumPaddingMultiplier * instrumentbroadening
                      / spectrumresolution)
            )

        convx = np.arange(
            -1.0 * convnumpoints * spectrumresolution, (convnumpoints +
                                                        1.0e-5) *
            spectrumresolution, spectrumresolution
            )

        if instrumentbroadeningshape == 'gaussian':
            convkernel = gaussian(convx, 1.0, 0.0, instrumentbroadening)
        elif instrumentBroadeningShape == 'lorentzian':
            convkernel = lorentzian(convx, 1.0, 0.0, instrumentbroadening)
        else:
            raise Exception("Error: Unrecognised instrumentBroadeningShape "
                            "'{0}'.".format(instrumentbroadeningshape))

        convkernel /= convkernel.sum()

        spectrummin = spectrummin - convnumpoints * spectrumresolution
        spectrummax = spectrummax + convnumpoints * spectrumresolution

    # Simulate spectrum.

    spectrumx = np.arange(spectrummin, spectrummax + 1.0e-5 *
                          spectrumresolution, spectrumresolution,
                          dtype=np.float64)
    spectrumy = np.zeros_like(spectrumx, dtype=np.float64)

    for frequency, intensity, linewidth in zip(frequencies, intensities,
                                               linewidths):
        spectrumy += lorentzian(spectrumx, intensity, frequency, linewidth)

    # Apply instrument broadening if required.

    if convKernel is not None:
        # mode = 'valid' will cause np.convolve() to trim the first and last
        # convNumPoints data points from the spectrum.

        spectrumx = spectrumx[convnumpoints:-convnumpoints]
        spectrumy = np.convolve(spectrumy, convkernel, mode='valid')

        # Just in case my maths went wrong somewhere...

        assert len(spectrumx) == len(spectrumy)

    # Normalise spectrum.

    spectrumynorm = spectrumy / math.fabs(spectrumy.max())

    # Return simulated spectrum.

    return spectrumx, spectrumy, spectrumynorm
