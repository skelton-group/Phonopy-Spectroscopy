# SpectroscoPy/Spectrum.py


# ---------
# Docstring
# ---------

""" Core routines for simulating spectra. """


# -------
# Imports
# -------

import math;

import numpy as np;


# ---------
# Constants
# ---------

""" To ensure the maximum is included in grids, np.arange() is called as np.arange(min, max + RangeStepMultiplier * step, step). """

RangeStepMultiplier = 1.0e-5;

""" A "safe" multiplier of the linewidth (e.g. sigma for Gaussian and gamma for Lorentzian functions) for padding grids used to evaluate peak functions. """

SpectrumPaddingMultiplier = 5.0;

""" Minimum number of points for automatically determining spectrum resolutions. """

SpectrumResolutionMinPoints = 1000;


# ---------
# Functions
# ---------

def Gaussian(x, i, mu, sigma):
    """ Return G(x) = (i / (sigma * sqrt(2 * pi))) * exp(-1 * (x - mu) ** 2 / (2 * sigma ** 2)). """

    # Definition of the Gaussian function with unit area taken from http://mathworld.wolfram.com/GaussianFunction.html.

    return (i / (sigma * math.sqrt(2.0 * math.pi))) * np.exp(-1.0 * (x - mu) ** 2 / (2 * sigma ** 2));

def Lorentzian(x, i, x0, gamma):
    """ Return L(x) = (i / pi) * ((0.5 * gamma) / ((x - x_0) ** 2) + (0.5 * gamma) ** 2). """

    # Definition of the Lorentzian function with unit area taken from http://mathworld.wolfram.com/LorentzianFunction.html.

    return (i / math.pi) * ((0.5 * gamma) / ((x - x0) ** 2 + (0.5 * gamma) ** 2));

def SimulateSpectrum(
    frequencies, intensities, linewidths,
    spectrumRange = None, spectrumResolution = None,
    instrumentBroadening = None, instrumentBroadeningShape = 'gaussian'
    ):

    """
    Given a set of mode frequencies and intensities plus a nominal linewidth or set of linewidths, simulate a spectrum.
    Optionally, apply an instrument broadening by convolving with a fixed-width Gaussian or Lorentzian function.

    Arguments:
        frequencies -- mode frequencies.
        intensities -- mode spectroscopic intensities.
        linewidths -- nominal mode linewidth or list of mode linewidths, taken to be the full-width at half-maxima (FWHM) of Lorentzian peak profiles.

    Keyword arguments:
        spectrumRange -- range of frequencies over which to simulate the spectrum (defaults to approx. min(frequencies) - SpectrumPaddingMultiplier * max(linewidths) -> max(frequencies) + SpectrumPaddingMultiplier * max(linewidths)).
        spectrumResolution -- frequency resolution of the spectrum (default: adjust to give at least SpectrumResolutionMinPoints data points, or a minimum resolution of 1, whichever is larger).
        instrumentBroadening -- instrumental broadening (default: no instrument broadening).
        instrumentBroadeningShape -- shape of instrumental broadening ('gaussian' or 'lorentzian'; default: 'gaussian').

    Notes:
        If a min and max are specified with spectrumRange, they may be modified to "align" with the resolution (i.e. so that the frequency axis starts and finished on an integer multiple of the resolution).
    """

    numModes = len(frequencies);

    if len(intensities) != numModes:
        raise Exception("Error: The lengths of frequencies and intensities are inconsistent.");

    try:
        # Assume linewidths is a scalar linewidth and repeat it to a list.

        linewidth = float(linewidths);
        linewidths = [linewidth] * numModes;
    except TypeError:
        # If linewidths is a list or multi-element NumPy array, casting to a float will raise a TypeError.

        if len(linewidths) != numModes:
            raise Exception("Error: The lengths of frequencies and linewidths are inconsistent.");

    # Set a minimum and maximum.

    spectrumMin, spectrumMax = None, None;

    if spectrumRange == None:
        maxLinewidth = max(linewidths);

        if numModes == 1 and maxLinewidth == 0.0:
            raise Exception("Error: Cannot determine a spectrum range automatically - please specify manually using the spectrumRange argument.");

        # +/- 5 sigma/gamma should be sufficient for this.

        spectrumMin = min(frequencies) - SpectrumPaddingMultiplier * maxLinewidth;
        spectrumMax = max(frequencies) + SpectrumPaddingMultiplier * maxLinewidth;
    else:
        spectrumMin, spectrumMax = spectrumRange;

        if spectrumMin == spectrumMax:
            raise Exception("Error: The min and max specified in spectrumRange are the same.");

        if spectrumMax < spectrumMin:
            spectrumMax, spectrumMin = spectrumMin, spectrumMax;

    # Set a resolution if required.

    if spectrumResolution == None:
        spectrumResolution = math.pow(
            10.0, math.floor(math.log(math.ceil(spectrumMax - spectrumMin) / SpectrumResolutionMinPoints))
            );

    # Make sure the range is "aligned" with the resolution.

    spectrumMin = spectrumResolution * math.floor(spectrumMin / spectrumResolution);
    spectrumMax = spectrumResolution * math.ceil(spectrumMax / spectrumResolution);

    # If applying instrument broadening, calculate the convolution kernel.
    # Also expand the range of the spectrum so the the convolution kernel does not produce boundary effects inside the selected region.

    convNumPoints, convKernel = None, None;

    if instrumentBroadening != None:
        # Calculating the convolution kernel to +/- 5 sigma/gamma should be sufficient.
        # According to https://en.wikipedia.org/wiki/Gaussian_blur, +/- 3 sigma is considered enough for Gaussian blurring kernels.

        convNumPoints = int(
            math.ceil(SpectrumPaddingMultiplier * instrumentBroadening / spectrumResolution)
            );

        convX = np.arange(
            -1.0 * convNumPoints * spectrumResolution, (convNumPoints + 1.0e-5) * spectrumResolution, spectrumResolution
            );

        if instrumentBroadeningShape == 'gaussian':
            convKernel = Gaussian(convX, 1.0, 0.0, instrumentBroadening);
        elif instrumentBroadeningShape == 'lorentzian':
            convKernel = Lorentzian(convX, 1.0, 0.0, instrumentBroadening);
        else:
            raise Exception("Error: Unrecognised instrumentBroadeningShape '{0}'.".format(instrumentBroadeningShape));

        convKernel /= convKernel.sum();

        spectrumMin = spectrumMin - convNumPoints * spectrumResolution;
        spectrumMax = spectrumMax + convNumPoints * spectrumResolution;

    # Simulate spectrum.

    spectrumX = np.arange(spectrumMin, spectrumMax + 1.0e-5 * spectrumResolution, spectrumResolution, dtype = np.float64);
    spectrumY = np.zeros_like(spectrumX, dtype = np.float64);

    for frequency, intensity, linewidth in zip(frequencies, intensities, linewidths):
        spectrumY += Lorentzian(spectrumX, intensity, frequency, linewidth);

    # Apply instrument broadening if required.

    if convKernel is not None:
        # mode = 'valid' will cause np.convolve() to trim the first and last convNumPoints data points from the spectrum.

        spectrumX = spectrumX[convNumPoints:-convNumPoints];
        spectrumY = np.convolve(spectrumY, convKernel, mode = 'valid');

        # Just in case my maths went wrong somewhere...

        assert len(spectrumX) == len(spectrumY);

    # Return simulated spectrum.

    return (spectrumX, spectrumY);
