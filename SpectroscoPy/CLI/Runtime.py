# SpectroscoPy/CLI/Runtime.py


# ---------
# Docstring
# ---------

""" Implements standard runtime routines. """


# -------
# Imports
# -------

from SpectroscoPy.IR import CalculateIRIntensity;
from SpectroscoPy.Spectrum import SimulateSpectrum;

from SpectroscoPy.DataExport import SavePeakTable, SaveSpectrum;
from SpectroscoPy.Plotting import PlotSpectrum;

from SpectroscoPy.Utilities import ConvertFrequencyUnits, GroupForPeakTable;


# ---------
# Constants
# ---------

""" Units of IR intensity. """

_IRIntensityUnits = "$e^{2}$ amu$^{-1}$";

""" Frequency units used for the peak table if the user requests a wavelength unit for the simulated spectrum. """

_DefaultPeakTableFrequencyUnits = 'inv_cm';


# ---------------
# IR Spectroscopy
# ---------------

def RunMode_IR(frequencies, frequencyUnits, eigendisplacements, becTensors, args, linewidths = None, irRepData = None):
    """
    Implements RunMode = 'ir' and generates a simulated IR spectrum for output.

    Main arguments:
        eigendisplacements -- Gamma-point eigendisplacements.
        becTensors -- Born effective-charge tensors.

    All other arguments are passed through to the RunMode_Output() function.
    """

    # Calculate IR intensities.

    irIntensities = [
        CalculateIRIntensity(eigendisplacement, becTensors)
            for eigendisplacement in eigendisplacements
        ];

    # Hand over to the RunMode_Output() routine.

    RunMode_Output(
        frequencies, frequencyUnits, irIntensities, _IRIntensityUnits, args, linewidths = linewidths, irRepData = irRepData
        );


# -----------
# Data Output
# -----------

def RunMode_Output(frequencies, frequencyUnits, intensities, intensityUnits, args, linewidths = None, irRepData = None):
    """
    Processes input data and arguments and saves a peak table and simulated spectrum.

    Arguments:
        frequencies -- mode frequencies.
        frequencyUnits -- units of frequencies (must be one of the units supported by the ConvertFrequencyUnits() function).
        eigendisplacements -- Gamma-point eigendisplacements.
        becTensors -- Born effective-charge tensors.
        args -- parsed command-line arguments with the options added by the UpdateParser() method in the Parser module.

    Keyword arguments:
        linewidths -- per-mode linewidths.
        irRepData -- set of (symbol, band_indices) tuples assigning ir. rep. symbols to modes.
    """

    # Most of the functions called from the other SpectroscoPy modules have their own error-handling, and we defer to those where possible to reduce "boilerplate" code.

    numModes = len(frequencies);

    if len(intensities) != numModes:
        raise Exception("Error: The numbers of frequencies and intensities are inconsistent.");

    if linewidths == None:
        if args.Linewidth == None:
            raise Exception("Error: A nominal linewidth must be supplied using the --linewidth argument.");

        if args.Linewidth < 0.0:
            raise Exception("Error: Negative mode linewidths are unphysical (!?).");

    if args.IrReps and irRepData == None:
        raise Exception("Error: If the IrReps argument is set, ir. rep. data must be available.");

    # Convert frequency units if required.

    frequencyUnitsConvert = args.FrequencyUnits;

    if frequencyUnitsConvert == 'um':
        # Outputting spectra in wavelength units (e.g. um) is supported, but it does not make sense to output peak parameters like this (central frequencies and linewidths are energies).
        # If the user requests a wavelength unit, we work in a default energy unit, output the peak table, simulate the spectrum, then convert just before output.

        frequencyUnitsConvert = _DefaultPeakTableFrequencyUnits;

        print("INFO: Peak table frequency units reset to '{0}'.".format(frequencyUnitsConvert));

    frequencies = ConvertFrequencyUnits(frequencies, frequencyUnits, frequencyUnitsConvert);

    if linewidths != None:
        linewidths = ConvertFrequencyUnits(linewidths, frequencyUnits, frequencyUnitsConvert);
    else:
        # Convert user-supplied nominal linewidth (in inv. cm) to frequency units.

        args.Linewidth = ConvertFrequencyUnits([args.Linewidth], 'inv_cm', frequencyUnitsConvert)[0];

    frequencyUnits = frequencyUnitsConvert;

    dataFormat = args.DataFormat.lower();

    # Save peak table.

    fileName = "PeakTable.{0}".format(dataFormat);

    if args.OutputPrefix != None:
        fileName = "{0}_{1}".format(args.OutputPrefix, fileName);

    if args.IrReps:
        ptFrequencies, ptIntensities, ptIrRepSymbols, ptLinewidths = GroupForPeakTable(frequencies, intensities, irRepData, linewidths = linewidths);

        SavePeakTable(
            ptFrequencies, ptIntensities, fileName,
            irRepSymbols = ptIrRepSymbols, linewidths = ptLinewidths,
            frequencyUnits = frequencyUnits, intensityUnits = intensityUnits,
            fileFormat = dataFormat
            );
    else:
        SavePeakTable(
            frequencies, intensities, fileName,
            linewidths = linewidths,
            frequencyUnits = frequencyUnits, intensityUnits = intensityUnits,
            fileFormat = dataFormat
            );

    # Simulate spectrum and convert units if required.

    spectrumRange = args.SpectrumRange;

    if spectrumRange != None and args.FrequencyUnits != frequencyUnits:
        spectrumMin, spectrumMax = args.SpectrumRange;

        if args.FrequencyUnits == 'um':
            # If the user specifies a spectrum min and/or max <= 0 when using wavelength units, this will cause problems when converting units.

            if spectrumMin <= 0.0 or spectrumMax <= 0.0:
                raise Exception("Error: If a spectrum range is specified for wavelength units, both min and max must be > 0.");

        spectrumRange = ConvertFrequencyUnits(spectrumRange, args.FrequencyUnits, frequencyUnits);

    spectrum = SimulateSpectrum(
        frequencies, intensities, linewidths if linewidths != None else args.Linewidth,
        spectrumRange = spectrumRange, spectrumResolution = args.SpectrumResolution,
        instrumentBroadeningWidth = args.InstrumentBroadeningWidth, instrumentBroadeningShape = args.InstrumentBroadeningShape
        );

    if args.FrequencyUnits != frequencyUnits:
        spectrumX, spectrumY, spectrumYNorm = spectrum;

        if args.FrequencyUnits == 'um':
            # Reverse arrays.

            spectrumX = spectrumX[::-1];
            spectrumY, spectrumYNorm = spectrumX[::-1], spectrumYNorm[::-1];

            # Frequencies <= 0 cause problems when converting to wavelength units.
            # Since wavelength units are not anticipated to be a routine use of this program, for now we simply mask them out and print a warning message for the user.

            mask = spectrumX > 0.0;

            clip = not mask.all();

            if clip:
                spectrumX = spectrumX[mask];
                spectrumY, spectrumYNorm = spectrumY[mask], spectrumYNorm[mask];

            spectrumX = ConvertFrequencyUnits(spectrumX, frequencyUnits, args.FrequencyUnits);

            if clip:
                print("INFO: Spectrum clipped to ({0:.2f}, {1:.2f}).".format(spectrumX[0], spectrumX[-1]));

        spectrum = (spectrumX, spectrumY, spectrumYNorm);

    # Save spectrum.

    fileName = "IRSpectrum.{0}".format(dataFormat);

    if args.OutputPrefix != None:
        fileName = "{0}_{1}".format(args.OutputPrefix, fileName);

    SaveSpectrum(
        spectrum, fileName,
        frequencyUnits = args.FrequencyUnits, intensityUnits = intensityUnits,
        fileFormat = dataFormat
        );

    # Plot spectrum.

    fileName = "IRSpectrum.png";

    if args.OutputPrefix != None:
        fileName = "{0}_{1}".format(args.OutputPrefix, fileName);

    PlotSpectrum(
        spectrum, fileName,
        frequencyUnits = args.FrequencyUnits,
        plotType = 'ir'
        );
