# SpectroscoPy/CLI/Runtime.py


# ---------
# Docstring
# ---------

""" Implements standard runtime routines. """


# -------
# Imports
# -------

import math;
import os;

from SpectroscoPy.IR import CalculateIRIntensity;
from SpectroscoPy.Raman import GenerateDisplacedStructures, CalculateRamanActivityTensor, ScalarAverageRamanActivityTensor;
from SpectroscoPy.Spectrum import SimulateSpectrum;

from SpectroscoPy.TextExport import GroupForPeakTable, SavePeakTable, SaveSpectrum;
from SpectroscoPy.YAMLExport import SaveRamanTensors;
from SpectroscoPy.Plotting import PlotSpectrum;

from SpectroscoPy.Constants import GetFrequencyUnitLabel;
from SpectroscoPy.Utilities import CalculateCellVolume, ConvertFrequencyUnits;

from SpectroscoPy.CLI.IntermediateIO import ReadRamanDataSet, WriteRamanDataSet;


# ---------
# Constants
# ---------

""" Frequency units used for the peak table if the user requests a wavelength unit for the simulated spectrum. """

_DefaultPeakTableFrequencyUnits = 'inv_cm';


""" Units of Raman displacement step (normal-mode coordinate) for plain-text intermediate YAML data sets. """

_Raman_StepUnits = "sqrt(amu) * Ang";

""" Units of maximum displacement distance in Raman calculations for plain-text intermediate YAML data sets. """

_Raman_DistanceUnits = "Ang";

""" Units of the cell volume for plain-text intermediate YAML data sets. """

_Raman_VolumeUnits = "Ang ^ 3";


""" Units of (tensor) Raman activity for plain-text YAML output file. """

_Raman_ActivityUnits = "sqrt(amu) * Ang ^ 2";


""" Units of IR intensity for plots/text output files. """

_IR_IntensityUnits = "$e^{2}$ amu$^{-1}$";

""" Units of (scalar) Raman intensity for plots/text output files. """

_Raman_IntensityUnits = "$\AA^{4}$ amu$^{-1}$";


# ---------------
# IR Spectroscopy
# ---------------

def RunMode_IR(frequencies, frequencyUnits, eigendisplacements, becTensors, args, linewidths = None, irRepData = None):
    """
    Calculates IR intensities for spectrum simulation.

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
        'ir', frequencies, frequencyUnits, irIntensities, _IR_IntensityUnits, args, linewidths = linewidths, irRepData = irRepData
        );


# ------------------
# Raman Spectroscopy
# ------------------

def RunMode_Raman_Disp(structure, frequencies, frequencyUnits, eigendisplacements, args):
    """
    Prepare and returns a sequence of displaced structures for a Raman activity calculation.
    Metadata about the displacements performed is written to "[<output_prefix>_]Raman.yaml".

    Arguments:
        structure -- structure to displace as a tuple of (lattice_vectors, atomic_symbols, atom_positions_frac).
        frequencies -- Gamma-point mode frequencies.
        frequencyUnits -- units of mode frequencies.
        eigendisplacements -- Gamma-point eigendisplacements.
        args -- parsed command-line arguments with the RamanStep, RamanStepType, RamanBandIndices and, optionally, OutputPrefix attributes.

    Return value:
        A list of tuples of (band_index, band_frequency, structures, disp_steps, max_disps) data.
        Band indices are zero-based.
    """

    numModes = len(frequencies);

    bandIndices = args.RamanBandIndices;

    if bandIndices == None:
        # Exclude the acoustic modes by identifying the three modes whose absolute frequency is closest to zero.
        # This is fairly crude compared to e.g. checking the eigenvector coefficients, but the user can always override it using command-line arguments.

        temp = [
            (i, math.fabs(frequency)) for i, frequency in enumerate(frequencies)
            ];

        temp.sort(key = lambda item : item[1]);

        excludeIndices = [i for i, _ in temp[:3]];

        print("INFO: Automatically excluding bands {0}, {1}, {2} from displacement set (to override use --raman_bands).".format(*[index + 1 for index in excludeIndices]));

        bandIndices = [
            i for i in range(0, numModes)
                if i not in excludeIndices
            ];

    for index in bandIndices:
        if index < 0 or index >= numModes:
            raise Exception("Error: Band index {0} out of range for # modes = {1}.".format(index + 1, numModes));

    # Build a dictionary of displacement sets: key -> band index, value -> (frequency, [(step, structure) tuples]).
    # At present, we use a default of two finite-difference steps per mode.

    dispSets = { };

    for index in bandIndices:
        dispStructures = GenerateDisplacedStructures(
            structure, eigendisplacements[index], args.RamanStep, 2, stepType = args.RamanStepType, omitZero = True
            );

        dispSets[index] = (frequencies[index], dispStructures);

    # Write data set to a YAML file.

    latticeVectors, _, _ = structure;
    cellVolume = CalculateCellVolume(latticeVectors);

    frequencies, modeDispSteps, modeMaxDisps = [], [], [];

    for index in bandIndices:
        frequency, (dispSteps, _, maxDisps) =  dispSets[index];

        frequencies.append(frequency);
        modeDispSteps.append(dispSteps);
        modeMaxDisps.append(maxDisps);

    fileName = "Raman.yaml";

    if args.OutputPrefix != None:
        fileName = "{0}_{1}".format(args.OutputPrefix, fileName);

    WriteRamanDataSet(
        fileName,
        cellVolume, _Raman_VolumeUnits,
        bandIndices, frequencies, frequencyUnits,
        modeDispSteps, _Raman_StepUnits, modeMaxDisps, _Raman_DistanceUnits,
        epsTensorSets = None
        );

    # Return the band indices, frequencies and displacement sets.

    dispSetsList = [];

    for index in bandIndices:
        frequency, (dispSteps, structures, maxDisps) = dispSets[index];

        dispSetsList.append(
            (index, frequency, structures, dispSteps, maxDisps)
            );

    return dispSetsList;

def RunMode_Raman_Read(epsTensors, args):
    """
    Collect dielectric tensors calculated for displaced structures prepared for a Raman calculation.

    Metadata on the displaced structures is read from "[<output_prefix>_]Raman.yaml".
    An updated file is then written which contains the displacement metadata and calculated dielectric tensors.

    Arguments:
        epsTensors -- list of dielectric tensors calculated for displaced structures.
        args -- parsed command-line arguments with the optional OutputPrefix attribute.

    Notes:
        The dielectric tensors are assumed to be in the same order as the displacements in the YAML data set written out by RunMode_Raman_Disp().
    """

    # Get the name of the intermediate YAML file.

    dataSetFile = "Raman.yaml";

    if args.OutputPrefix != None:
        dataSetFile = "{0}_{1}".format(args.OutputPrefix, dataSetFile);

    if not os.path.isfile(dataSetFile):
        raise Exception("Error: Displacement data set file \"{0}\" not found.".format(dataSetFile));

    # Read displacement data set: tuple of (cell_volume, volume_units, band_indices, frequencies, frequency_units, disp_step_sets, step_units, max_disps_sets, distance_units, eps_tensor_sets).

    cellVolume, volumeUnits, bandIndices, frequencies, frequencyUnits, dispStepSets, stepUnits, maxDispsSets, distanceUnits, _ = ReadRamanDataSet(dataSetFile);

    # Verify units.

    if volumeUnits != _Raman_VolumeUnits:
        raise Exception("Error: Unexpected volume units '{0}' read from displacement data set file \"{1}\".".format(volumeUnits, dataSetFile));

    if stepUnits != _Raman_StepUnits:
        raise Exception("Error: Unexpected step units '{0}' read from displacement data set file \"{1}\".".format(stepUnits, dataSetFile));

    if distanceUnits != _Raman_DistanceUnits:
        raise Exception("Error: Unexpected distance units '{0}' read from displacement data set file \"{1}\".".format(ditanceUnits, dataSetFile));

    # Check the expected number of dielectric tensors has been supplied.

    numEpsTensorsExpected = sum(
        len(dispSteps) for dispSteps in dispStepSets
        );

    if len(epsTensors) != numEpsTensorsExpected:
        raise Exception("Error: Incorrect number of dielectric tensors supplied - {0} expected, {1} supplied.".format(numFilesExpected, len(epsTensors)));

    # Sort the dielectric tensors into sets.

    epsTensorSets = [];

    inputPointer = 0;

    for dispStepSet in dispStepSets:
        epsTensorSet = [];

        for i in range(0, len(dispStepSet)):
            epsTensorSet.append(epsTensors[inputPointer]);

            inputPointer += 1;

        epsTensorSets.append(epsTensorSet);

    # Write an updated data set containing the extracted dielectric tensors.

    WriteRamanDataSet(
        dataSetFile,
        cellVolume, volumeUnits,
        bandIndices, frequencies, frequencyUnits,
        dispStepSets, stepUnits, maxDispsSets, distanceUnits,
        epsTensorSets
        );

def RunMode_Raman_PostProc(args, linewidths = None, irRepData = None):
    """
    Post process a data set prepared for a Raman calculation.

    A data set containing frequencies, displacement steps and dielectric tensors is read from "[<output_prefix>_]Raman.yaml".
    The Raman activity tensors are calculated from the derivatives of the static dielectric tensor with respect to the normal-mode coordinates and written to "[<output_prefix>_]RamanActivity.yaml".
    The activity tensors are averaged to obtai scalar intensities that are then used to plot a spectrum.

    Arguments:
        args -- parsed command-line arguments with the optional OutputPrefix attribute, plus settings to be passed to the RunMode_Output() routine.

    Keyword arguments:
        TODO!
    """

    # Work out the name of the intermediate YAML file and check the file exists.

    dataSetFile = "Raman.yaml";

    if args.OutputPrefix != None:
        dataSetFile = "{0}_{1}".format(args.OutputPrefix, dataSetFile);

    if not os.path.isfile(dataSetFile):
        raise Exception("Error: Displacement data set file \"{0}\" not found.".format(dataSetFile));

    # Read in the Raman data set and verify hard-coded units.

    cellVolume, volumeUnits, bandIndices, frequencies, frequencyUnits, dispStepSets, stepUnits, maxDispsSets, distanceUnits, epsTensorSets = ReadRamanDataSet(dataSetFile);

    if volumeUnits != _Raman_VolumeUnits:
        raise Exception("Error: Unexpected volume units '{0}' read from Raman data set file \"{1}\".".format(volumeUnits, dataSetFile));

    if stepUnits != _Raman_StepUnits:
        raise Exception("Error: Unexpected step units '{0}' read from Raman data set file \"{1}\".".format(stepUnits, dataSetFile));

    if distanceUnits != _Raman_DistanceUnits:
        raise Exception("Error: Unexpected distance units '{0}' read from Raman data set file \"{1}\".".format(ditanceUnits, dataSetFile));

    # Calculate Raman activity tensors and averaged scalar intensities.

    ramanTensors, ramanIntensities = [], [];

    for dispSteps, epsTensors in zip(dispStepSets, epsTensorSets):
        ramanTensor = CalculateRamanActivityTensor(epsTensors, dispSteps, cellVolume, 'finite_difference');

        ramanTensors.append(ramanTensor);

        ramanIntensities.append(
            ScalarAverageRamanActivityTensor(ramanTensor)
            );

    # Write out Raman tensors.

    fileName = "Raman-Tensors.yaml";

    if args.OutputPrefix != None:
        fileName = "{0}_{1}".format(args.OutputPrefix, fileName);

    # Convert frequency units, if required.

    frequenciesOutput = None;

    if frequencyUnits != args.FrequencyUnits:
        frequenciesOutput = ConvertFrequencyUnits(frequencies, frequencyUnits, args.FrequencyUnits);
    else:
        frequenciesOutput = frequencies;

    SaveRamanTensors(
        bandIndices, frequencies, ramanTensors, fileName,
        frequencyUnits = frequencyUnits, activityUnits = _Raman_ActivityUnits
        );

    # If linewidths are supplied, discard data for modes for which Raman activities were not calculated.

    if linewidths != None:
        linewidths = [
            linewidths[index] for index in bandIndices
            ];

    # If ir. rep. data is supplied, we need to:
    # (1) re-map the sets of band indices to account for modes that have been discarded; and
    # (2) remove sets of data for which all modes have been discarded.

    indexMapping = {
        index : i for i, index in enumerate(bandIndices)
        };

    if irRepData != None:
        irRepDataNew = [];

        for symbol, indices in irRepData:
            indicesNew = [];

            for index in indices:
                if index in indexMapping:
                    indicesNew.append(
                        indexMapping[index]
                        );

            if len(indicesNew) > 0:
                irRepDataNew.append(
                    (symbol, indicesNew)
                    );

        irRepData = irRepDataNew;

    # Hand over to RunMode_Output() routine with the scalar Raman intensities.

    RunMode_Output(
        'raman', frequencies, frequencyUnits, ramanIntensities, _Raman_IntensityUnits, args, linewidths = linewidths, irRepData = irRepData
        );

# -----------
# Data Output
# -----------

def RunMode_Output(spectrumType, frequencies, frequencyUnits, intensities, intensityUnits, args, linewidths = None, irRepData = None):
    """
    Processes input data and arguments and saves a peak table and simulated spectrum.

    Arguments:
        spectrumType -- type of spectrum being output ('ir' or 'raman').
        frequencies -- mode frequencies.
        frequencyUnits -- units of frequencies (must be one of the units supported by the ConvertFrequencyUnits() function).
        intensities -- band intensities.
        intensityUnits -- units of intensities.
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

    # Set prefix for output files.

    outputPrefix = None;

    if spectrumType == 'ir':
        outputPrefix = "IR";
    elif spectrumType == 'raman':
        outputPrefix = "Raman";
    else:
        # We _could_ work around this, but if it happens it's almost certainly a programming error.

        raise Exception("Error: Unknown spectrumType '{0}' - if you are running one of the command-line scripts this is probably a bug.".format(spectrumType));

    if args.OutputPrefix != None:
        outputPrefix = "{0}_{1}".format(args.OutputPrefix, outputPrefix);

    # Save peak table.

    fileName = "{0}-PeakTable.{1}".format(outputPrefix, dataFormat);

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
        instrumentBroadening = args.InstrumentBroadening, instrumentBroadeningShape = args.InstrumentBroadeningShape
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

    fileName = "{0}-Spectrum.{1}".format(outputPrefix, dataFormat);

    SaveSpectrum(
        spectrum, fileName,
        frequencyUnits = args.FrequencyUnits, intensityUnits = intensityUnits,
        fileFormat = dataFormat
        );

    # Plot spectrum.

    fileName = "{0}-Spectrum.png".format(outputPrefix);

    PlotSpectrum(
        spectrum, fileName,
        frequencyUnits = args.FrequencyUnits,
        plotType = spectrumType
        );
