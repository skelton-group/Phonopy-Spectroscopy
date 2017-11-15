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


# ---------------
# IR Spectroscopy
# ---------------

_IRIntensityUnits = "$e^{2}$ amu$^{-1}$";

def RunMode_IR(frequencies, frequencyUnits, eigendisplacements, becTensors, args, linewidths = None, irRepData = None):
    """
    Implements RunMode = 'ir' and generates and outputs a simulated IR spectrum.

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

    # Most of the functions in the other SpectroscoPy modules have their own error-handling, and we defer to those where possible to reduce "boilerplate" code.

    numModes = len(frequencies);

    if len(eigendisplacements) != numModes:
        raise Exception("Error: The numbers of frequencies and eigendisplacements are inconsistent.");

    if linewidths == None:
        if args.Linewidth == None:
            raise Exception("Error: A nominal linewidth must be supplied using the --linewidth argument.");

        if args.Linewidth < 0.0:
            raise Exception("Error: Negative mode linewidths are unphysical (!?).");

    if args.IrReps and irRepData == None:
        raise Exception("Error: If the IrReps argument is set, ir. rep. data must be available.");

    # Make sure frequencies are in supported units.

    frequencies = ConvertFrequencyUnits(
        frequencies, frequencyUnits, args.FrequencyUnits if args.FrequencyUnits != None else frequencyUnits
        );

    if args.FrequencyUnits != None:
        frequencyUnits = args.FrequencyUnits;

    # Calculate IR intensities.

    irIntensities = [
        CalculateIRIntensity(eigendisplacement, becTensors)
            for eigendisplacement in eigendisplacements
        ];

    dataFormat = args.DataFormat.lower();

    # Save peak table.

    fileName = "PeakTable.{0}".format(dataFormat);

    if args.OutputPrefix != None:
        fileName = "{0}_{1}".format(args.OutputPrefix, fileName);

    if args.IrReps:
        ptFrequencies, ptIntensities, ptIrRepSymbols, ptLinewidths = GroupForPeakTable(frequencies, irIntensities, irRepData, linewidths = linewidths);

        SavePeakTable(
            ptFrequencies, ptIntensities, fileName,
            irRepSymbols = ptIrRepSymbols, linewidths = ptLinewidths,
            frequencyUnits = frequencyUnits, intensityUnits = _IRIntensityUnits,
            fileFormat = dataFormat
            );
    else:
        SavePeakTable(
            frequencies, irIntensities, fileName,
            linewidths = linewidths,
            frequencyUnits = frequencyUnits, intensityUnits = _IRIntensityUnits,
            fileFormat = dataFormat
            );

    # Simulate spectrum.

    spectrum = SimulateSpectrum(
        frequencies, irIntensities, linewidths if linewidths != None else args.Linewidth,
        spectrumRange = args.SpectrumRange, spectrumResolution = args.SpectrumResolution,
        instrumentBroadeningWidth = args.InstrumentBroadeningWidth, instrumentBroadeningShape = args.InstrumentBroadeningShape
        );

    # Save spectrum.

    fileName = "IRSpectrum.{0}".format(dataFormat);

    if args.OutputPrefix != None:
        fileName = "{0}_{1}".format(args.OutputPrefix, fileName);

    SaveSpectrum(
        spectrum, fileName,
        frequencyUnits = frequencyUnits, intensityUnits = _IRIntensityUnits,
        fileFormat = 'dat'
        );

    # Plot spectrum.

    fileName = "IRSpectrum.png";

    if args.OutputPrefix != None:
        fileName = "{0}_{1}".format(args.OutputPrefix, fileName);

    PlotSpectrum(
        spectrum, fileName,
        frequencyUnits = frequencyUnits,
        plotType = 'ir'
        );
