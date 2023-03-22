# SpectroscoPy/CLI/Runtime.py


# ---------
# Docstring
# ---------

""" Implements standard runtime routines. """


# -------
# Imports
# -------

import math
import os

from SpectroscoPy.IR import calculateirintensity

from SpectroscoPy.Raman import (generatedisplacedstructures,
                                calculateramanactivitytensor,
                                scalaraverageramanactivitytensor
                                )

from SpectroscoPy.Spectrum import simulatespectrum

from SpectroscoPy.TextExport import (groupforpeaktable, savepeaktable,
                                     savespectrum
                                     )

from SpectroscoPy.YAMLExport import saveramantensors

from SpectroscoPy.Plotting import plotspectrum

from SpectroscoPy.Constants import getfrequencyunitlabel

from SpectroscoPy.Utilities import calculatecellvolume, convertfrequencyunits

from SpectroscoPy.CLI.IntermediateIO import (readramandataset,
                                             writeramandataset
                                             )


# ---------
# Constants
# ---------

""" Frequency units used for the peak table if the user requests a 
wavelength unit for the simulated spectrum. """

_DefaultPeakTableFrequencyUnits = 'inv_cm'


""" Units of Raman displacement step (normal-mode coordinate) for plain-text 
intermediate YAML data sets. """

_Raman_StepUnits = "sqrt(amu) * Ang"

""" Units of maximum displacement distance in Raman calculations for 
plain-text intermediate YAML data sets. """

_Raman_DistanceUnits = "Ang"

""" Units of the cell volume for plain-text intermediate YAML data sets. """

_Raman_VolumeUnits = "Ang ^ 3"


""" Units of (tensor) Raman activity for plain-text YAML output file. """

_Raman_ActivityUnits = "sqrt(amu) * Ang ^ 2"


""" Units of IR intensity for plots/text output files. """

_IR_IntensityUnits = "$e^{2}$ amu$^{-1}$"

""" Units of (scalar) Raman intensity for plots/text output files. """

_Raman_IntensityUnits = r"$\AA^{4}$ amu$^{-1}$"


# ---------------
# IR Spectroscopy
# ---------------

def runmode_ir(frequencies, frequencyunits, eigendisplacements, bectensors,
               args, linewidths=None, irrepdata=None):
    """
    Calculates IR intensities for spectrum simulation.

    Main arguments:
        eigendisplacements -- Gamma-point eigendisplacements.
        becTensors -- Born effective-charge tensors.

    All other arguments are passed through to the RunMode_Output() function.
    """

    # Calculate IR intensities.

    irintensities = [
        calculateirintensity(eigendisplacement, bectensors)
        for eigendisplacement in eigendisplacements
        ]

    # Hand over to the RunMode_Output() routine.

    runmode_output(
        'ir', frequencies, frequencyunits, irintensities,
        _IR_IntensityUnits, args, linewidths=linewidths,
        irrepdata=irrepdata
        )


# ------------------
# Raman Spectroscopy
# ------------------

def runmode_raman_disp(structure, frequencies, frequencyunits,
                       eigendisplacements, args):
    """
    Prepare and returns a sequence of displaced structures for a Raman
    activity calculation.
    Metadata about the displacements performed is written to "[
    <output_prefix>_]Raman.yaml".

    Arguments:
        structure -- structure to displace as a tuple of (lattice_vectors,
        atomic_symbols, atom_positions_frac).
        frequencies -- Gamma-point mode frequencies.
        frequencyUnits -- units of mode frequencies.
        eigendisplacements -- Gamma-point eigendisplacements.
        args -- parsed command-line arguments with the RamanStep,
        RamanStepType, RamanBandIndices and, optionally, OutputPrefix
        attributes.

    Return value:
        A list of tuples of (band_index, band_frequency, structures,
        disp_steps, max_disps) data.
        Band indices are zero-based.
    """

    nummodes = len(frequencies)

    bandindices = args.RamanBandIndices

    if bandindices is None:
        # Exclude the acoustic modes by identifying the three modes whose
        # absolute frequency is closest to zero.
        # This is fairly crude compared to e.g. checking the eigenvector
        # coefficients, but the user can always override it using
        # command-line arguments.

        temp = [
            (i, math.fabs(frequency)) for i, frequency in enumerate(
                frequencies)
            ]

        temp.sort(key=lambda item: item[1])

        excludeindices = [i for i, _ in temp[:3]]

        print("INFO: Automatically excluding bands {0}, {1}, {2} from "
              "displacement set (to override use --raman_bands).".format(*[
                index + 1 for index in excludeindices])
              )

        bandindices = [
            i for i in range(0, nummodes)
            if i not in excludeindices
            ]

    for index in bandindices:
        if index < 0 or index >= nummodes:
            raise Exception("Error: Band index {0} out of range for # modes "
                            "= {1}.".format(index + 1, nummodes)
                            )

    # Build a dictionary of displacement sets: key -> band index, value -> (
    # frequency, [(step, structure) tuples]).
    # At present, we use a default of two finite-difference steps per mode.

    dispsets = {}

    for index in bandindices:
        dispstructures = generatedisplacedstructures(
            structure, eigendisplacements[index], args.RamanStep, 2,
            steptype=args.RamanStepType, omitzero=True
            )

        dispsets[index] = (frequencies[index], dispstructures)

    # Write data set to a YAML file.

    latticevectors, _, _ = structure
    cellvolume = calculatecellvolume(latticevectors)

    frequencies, modedispsteps, modemaxdisps = [], [], []

    for index in bandindices:
        frequency, (dispsteps, _, maxdisps) = dispsets[index]

        frequencies.append(frequency)
        modedispsteps.append(dispsteps)
        modemaxdisps.append(maxdisps)

    filename = "Raman.yaml"

    if args.OutputPrefix is not None:
        filename = "{0}_{1}".format(args.OutputPrefix, filename)

    writeramandataset(filename, cellvolume, _Raman_VolumeUnits, bandindices,
                      frequencies, frequencyunits, modedispsteps,
                      _Raman_StepUnits, modemaxdisps, _Raman_DistanceUnits,
                      epstensorsets=None
                      )

    # Return the band indices, frequencies and displacement sets.

    dispsetslist = []

    for index in bandindices:
        frequency, (dispsteps, structures, maxdisps) = dispsets[index]

        dispsetslist.append(
            (index, frequency, structures, dispsteps, maxdisps)
            )

    return dispsetslist


def runmode_raman_read(epstensors, args):
    """
    Collect dielectric tensors calculated for displaced structures prepared
    for a Raman calculation.

    Metadata on the displaced structures is read from "[
    <output_prefix>_]Raman.yaml".
    An updated file is then written which contains the displacement metadata
    and calculated dielectric tensors.

    Arguments:
        epsTensors -- list of dielectric tensors calculated for displaced
        structures.
        args -- parsed command-line arguments with the optional OutputPrefix
        attribute.

    Notes:
        The dielectric tensors are assumed to be in the same order as the
        displacements in the YAML data set written out by RunMode_Raman_Disp().
    """

    # Get the name of the intermediate YAML file.

    datasetfile = "Raman.yaml"

    if args.OutputPrefix is not None:
        datasetfile = "{0}_{1}".format(args.OutputPrefix, datasetfile)

    if not os.path.isfile(datasetfile):
        raise Exception("Error: Displacement data set file \"{0}\" not "
                        "found.".format(datasetfile)
                        )

    # Read displacement data set: tuple of (cell_volume, volume_units,
    # band_indices, frequencies, frequency_units, disp_step_sets,
    # step_units, max_disps_sets, distance_units, eps_tensor_sets).

    cellvolume, volumeunits, bandindices, frequencies, frequencyunits, \
        dispstepsets, stepunits, maxdispssets, distanceunits, \
        _ = readramandataset(datasetfile)

    # Verify units.

    if volumeunits != _Raman_VolumeUnits:
        raise Exception("Error: Unexpected volume units '{0}' read from "
                        "displacement data set file \"{1}\".".format(
                                        volumeunits, datasetfile)
                        )

    if stepunits != _Raman_StepUnits:
        raise Exception("Error: Unexpected step units '{0}' read from "
                        "displacement data set file \"{1}\".".format(
                                        stepunits, datasetfile)
                        )

    if distanceunits != _Raman_DistanceUnits:
        raise Exception("Error: Unexpected distance units '{0}' read from "
                        "displacement data set file \"{1}\".".format(
                                        distanceunits, datasetfile)
                        )

    # Check the expected number of dielectric tensors has been supplied.

    numepstensorsexpected = sum(
        len(dispsteps) for dispsteps in dispstepsets
        )

    if len(epstensors) != numepstensorsexpected:
        raise Exception("Error: Incorrect number of dielectric tensors "
                        "supplied - {0} expected, {1} supplied.".format(
                                        numfilesexpected, len(epstensors))
                        )

    # Sort the dielectric tensors into sets.

    epstensorsets = []

    inputpointer = 0

    for dispstepset in dispstepsets:
        epstensorset = []

        for i in range(0, len(dispstepset)):
            epstensorset.append(epstensors[inputpointer])

            inputpointer += 1

        epstensorsets.append(epstensorset)

    # Write an updated data set containing the extracted dielectric tensors.

    writeramandataset(datasetfile, cellvolume, volumeunits, bandindices,
                      frequencies, frequencyunits, dispstepsets, stepunits,
                      maxdispssets, distanceunits, epstensorsets
                      )


def runmode_raman_postproc(args, linewidths=None, irrepdata=None):
    """
    Post process a data set prepared for a Raman calculation.

    A data set containing frequencies, displacement steps and dielectric
    tensors is read from "[<output_prefix>_]Raman.yaml".
    The Raman activity tensors are calculated from the derivatives of the
    static dielectric tensor with respect to the normal-mode coordinates and
    written to "[<output_prefix>_]RamanActivity.yaml".
    The activity tensors are averaged to obtai scalar intensities that are
    then used to plot a spectrum.

    Arguments:
        args -- parsed command-line arguments with the optional OutputPrefix
        attribute, plus settings to be passed to the RunMode_Output() routine.

    Keyword arguments:
        TODO!
    """

    # Work out the name of the intermediate YAML file and check the file
    # exists.

    datasetfile = "Raman.yaml"

    if args.OutputPrefix is not None:
        datasetfile = "{0}_{1}".format(args.OutputPrefix, datasetfile)

    if not os.path.isfile(datasetfile):
        raise Exception("Error: Displacement data set file \"{0}\" not "
                        "found.".format(datasetfile)
                        )

    # Read in the Raman data set and verify hard-coded units.

    cellvolume, volumeunits, bandindices, frequencies, frequencyunits, \
        dispstepsets, stepunits, maxdispssets, distanceunits, epstensorsets \
        = readramandataset(datasetfile)

    if volumeunits != _Raman_VolumeUnits:
        raise Exception("Error: Unexpected volume units '{0}' read from "
                        "Raman data set file \"{1}\".".format(volumeunits,
                                                              datasetfile)
                        )

    if stepunits != _Raman_StepUnits:
        raise Exception("Error: Unexpected step units '{0}' read from Raman "
                        "data set file \"{1}\".".format(stepunits,
                                                        datasetfile)
                        )

    if distanceunits != _Raman_DistanceUnits:
        raise Exception("Error: Unexpected distance units '{0}' read from "
                        "Raman data set file \"{1}\".".format(distanceunits,
                                                              datasetfile)
                        )

    # Calculate Raman activity tensors and averaged scalar intensities.

    ramantensors, ramanintensities = [], []

    for dispsteps, epstensors in zip(dispstepsets, epstensorsets):
        ramantensor = calculateramanactivitytensor(epstensors, dispsteps,
                                                   cellvolume,
                                                   'finite_difference')

        ramantensors.append(ramantensor)

        ramanintensities.append(
            scalaraverageramanactivitytensor(ramantensor)
            )

    # Write out Raman tensors.

    filename = "Raman-Tensors.yaml"

    if args.OutputPrefix is not None:
        filename = "{0}_{1}".format(args.OutputPrefix, filename)

    # Convert frequency units, if required.

    frequenciesoutput = None

    if frequencyUnits != args.FrequencyUnits:
        frequenciesoutput = convertfrequencyunits(frequencies,
                                                  frequencyunits,
                                                  args.FrequencyUnits)
    else:
        frequenciesoutput = frequencies

    saveramantensors(
        bandindices, frequencies, ramantensors, filename,
        frequencyunits=frequencyunits, activityunits=_Raman_ActivityUnits
        )

    # If linewidths are supplied, discard data for modes for which Raman
    # activities were not calculated.

    if linewidths is not None:
        linewidths = [
            linewidths[index] for index in bandindices
            ]

    # If ir. rep. data is supplied, we need to:
    # (1) re-map the sets of band indices to account for modes that have
    # been discarded; and (2) remove sets of data for which all modes have
    # been discarded.

    indexmapping = {
        index: i for i, index in enumerate(bandindices)
        }

    if irrepdata is not None:
        irrepdatanew = []

        for symbol, indices in irrepdata:
            indicesnew = []

            for index in indices:
                if index in indexmapping:
                    indicesnew.append(
                        indexmapping[index]
                        )

            if len(indicesnew) > 0:
                irrepdatanew.append(
                    (symbol, indicesnew)
                    )

        irrepdata = irrepdatanew

    # Hand over to RunMode_Output() routine with the scalar Raman intensities.

    Runmode_output(
        'raman', frequencies, frequencyunits, ramanintensities,
        _Raman_IntensityUnits, args, linewidths=linewidths, irrepdata=irrepdata
        )

# -----------
# Data Output
# -----------


def runmode_output(spectrumtype, frequencies, frequencyunits, intensities,
                   intensityunits, args, linewidths=None, irrepdata=None):
    """
    Processes input data and arguments and saves a peak table and simulated
    spectrum.

    Arguments:
        spectrumType -- type of spectrum being output ('ir' or 'raman').
        frequencies -- mode frequencies.
        frequencyUnits -- units of frequencies (must be one of the units
        supported by the ConvertFrequencyUnits() function).
        intensities -- band intensities.
        intensityUnits -- units of intensities.
        args -- parsed command-line arguments with the options added by the
        update_parser() method in the Parser module.

    Keyword arguments:
        linewidths -- per-mode linewidths.
        irRepData -- set of (symbol, band_indices) tuples assigning ir. rep.
        symbols to modes.
    """

    # Most of the functions called from the other SpectroscoPy modules have
    # their own error-handling, and we defer to those where possible to
    # reduce "boilerplate" code.

    nummodes = len(frequencies)

    if len(intensities) != nummodes:
        raise Exception("Error: The numbers of frequencies and intensities "
                        "are inconsistent.")

    if linewidths is None:
        if args.Linewidth < 0.0:
            raise Exception("Error: Negative mode linewidths are unphysical "
                            "(!?).")

    if args.IrReps and irrepdata is None:
        raise Exception("Error: If the IrReps argument is set, ir. rep. data "
                        "must be available.")

    # Convert frequency units if required.

    frequencyunitsconvert = args.FrequencyUnits

    if frequencyunitsconvert == 'um':
        # Outputting spectra in wavelength units (e.g. um) is supported,
        # but it does not make sense to output peak parameters like this (
        # central frequencies and linewidths are energies).
        # If the user requests a wavelength unit, we work in a default
        # energy unit, output the peak table, simulate the spectrum,
        # then convert just before output.

        frequencyunitsconvert = _DefaultPeakTableFrequencyUnits

        print("INFO: Peak table frequency units reset to '{0}'.".format(
            frequencyunitsconvert)
            )

    frequencies = convertfrequencyunits(frequencies, frequencyunits,
                                        frequencyunitsconvert)

    if linewidths is not None:
        linewidths = convertfrequencyunits(linewidths, frequencyunits,
                                           frequencyunitsconvert)
    else:
        # Convert user-supplied nominal linewidth (in inv. cm) to frequency
        # units.

        args.Linewidth = convertfrequencyunits([args.Linewidth], 'inv_cm',
                                               frequencyunitsconvert)[0]

    frequencyunits = frequencyunitsconvert

    dataformat = args.DataFormat.lower()

    # Set prefix for output files.

    outputprefix = None

    if spectrumtype == 'ir':
        outputprefix = "IR"
    elif spectrumtype == 'raman':
        outputprefix = "Raman"
    else:
        # We _could_ work around this, but if it happens it's almost
        # certainly a programming error.

        raise Exception("Error: Unknown spectrumType '{0}' - if you are "
                        "running one of the command-line scripts this is "
                        "probably a bug.".format(spectrumtype)
                        )

    if args.OutputPrefix is not None:
        outputprefix = "{0}_{1}".format(args.OutputPrefix, outputprefix)

    # Save peak table.

    filename = "{0}-PeakTable.{1}".format(outputprefix, dataformat)

    if args.IrReps:
        ptfrequencies, ptintensities, ptirrepsymbols, ptlinewidths = \
            groupforpeaktable(frequencies, intensities, irrepdata,
                              linewidths=linewidths)

        savepeaktable(
            ptfrequencies, ptintensities, filename,
            irrepsymbols=ptirrepsymbols, linewidths=ptlinewidths,
            frequencyunits=frequencyunits, intensityunits=intensityunits,
            fileformat=dataformat
            )
    else:
        savepeaktable(
            frequencies, intensities, filename,
            linewidths=linewidths,
            frequencyunits=frequencyunits, intensityunits=intensityunits,
            fileformat=dataformat
            )

    # Simulate spectrum and convert units if required.

    spectrumrange = args.SpectrumRange

    if spectrumrange is not None and args.FrequencyUnits != frequencyunits:
        spectrummin, spectrummax = args.SpectrumRange

        if args.FrequencyUnits == 'um':
            # If the user specifies a spectrum min and/or max <= 0 when
            # using wavelength units, this will cause problems when
            # converting units.

            if spectrummin <= 0.0 or spectrummax <= 0.0:
                raise Exception("Error: If a spectrum range is specified for "
                                "wavelength units, both min and max must be "
                                "> 0.")

        spectrumrange = convertfrequencyunits(spectrumrange,
                                              args.FrequencyUnits,
                                              frequencyunits)

    spectrum = simulatespectrum(
        frequencies, intensities, linewidths if linewidths is not None else
        args.Linewidth, spectrumrange=spectrumrange,
        spectrumresolution=args.SpectrumResolution,
        instrumentbroadening=args.InstrumentBroadening,
        instrumentbroadeningshape=args.InstrumentBroadeningShape
        )

    if args.FrequencyUnits != frequencyunits:
        spectrumx, spectrumy, spectrumynorm = spectrum

        if args.FrequencyUnits == 'um':
            # Reverse arrays.

            spectrumx = spectrumx[::-1]
            spectrumy, spectrumynorm = spectrumx[::-1], spectrumynorm[::-1]

            # Frequencies <= 0 cause problems when converting to wavelength
            # units.
            # Since wavelength units are not anticipated to be a routine use
            # of this program, for now we simply mask them out and print a
            # warning message for the user.

            mask = spectrumx > 0.0

            clip = not mask.all()

            if clip:
                spectrumx = spectrumx[mask]
                spectrumy, spectrumynorm = spectrumy[mask], spectrumynorm[mask]

            spectrumx = convertfrequencyunits(spectrumx, frequencyunits,
                                              args.FrequencyUnits)

            if clip:
                print("INFO: Spectrum clipped to ({0:.2f}, "
                      "{1:.2f}).".format(spectrumx[0], spectrumx[-1])
                      )

        spectrum = (spectrumx, spectrumy, spectrumynorm)

    # Save spectrum.

    filename = "{0}-Spectrum.{1}".format(outputprefix, dataformat)

    savespectrum(spectrum, filename, frequencyunits=args.FrequencyUnits,
                 intensityunits=intensityunits, fileformat=dataformat
                 )

    # Plot spectrum.

    filename = "{0}-Spectrum.png".format(outputprefix)

    plotspectrum(spectrum, filename, frequencyunits=args.FrequencyUnits,
                 plottype=spectrumtype
                 )
