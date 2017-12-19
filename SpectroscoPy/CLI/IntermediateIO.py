# SpectroscoPy/CLI/IntermediateIO.py


# ---------
# Docstring
# ---------

""" Routines for storing and retrieving intermediate data in multi-step calculations. """


# ---------
# Constants
# ---------

""" Default step unit written to intermediate YAML files. """

_DefaultStepUnit = "Unspecified";

""" Default distance unit written to intermediate YAML files. """

_DefaultDistanceUnit = "Unspecified";


# -------
# Imports
# -------

import yaml;

# Try to use the C YAML Loader if possible; if not, fall back to the standard routines.

try:
    from yaml import CLoader as Loader;
except ImportError:
    from yaml import Loader;

import numpy as np;


# ------------------
# Raman Calculations
# ------------------

""" Default string written to intermediate YAML files when a unit is unspecified. """

_RamanDataSet_UnspecifiedUnitsDescription = "Unspecified";

def ReadRamanDataSet(filePath):
    """
    Read data for a Raman calculation from a YAML-format file produced by WriteRamanDataSet().

    Return value:
        A tuple of (cell_volume, volume_units, band_indices, frequencies, frequency_units, disp_step_sets, step_units, max_disps_sets, distance_units, eps_tensor_sets) data read from the file.
        (The order of the elements matches the arguments passed to WriteRamanDataSet()).
    """

    # Read and parse data set.

    dataSet = None;

    with open(filePath, 'r') as inputReader:
        dataSet = yaml.load(inputReader, Loader = Loader);

    # Validate.

    for requiredKey in ['frequency_units', 'step_units', 'distance_units', 'volume_units', 'cell_volume', 'displacement_sets']:
        if requiredKey not in dataSet:
            raise Exception("Error: Missing key '{0}' in displacement data set.".format(requiredKey));

    if len(dataSet['displacement_sets']) == 0:
        raise Exception("Error: Displacement dataset contains no displacement sets.");

    hasEpsTensors = False;

    for item1 in dataSet['displacement_sets']:
        for requiredKey in ['band_index', 'frequency', 'displacements']:
            if requiredKey not in item1:
                raise Exception("Error: Missing key '{0}' in one or more displacement sets.".format(requiredKey));

        for item2 in item1['displacements']:
            for requiredKey in ['displacement_step', 'max_cartesian_displacement']:
                if requiredKey not in item2:
                    raise Exception("Error: Missing key '{0}' in one or more displacement entries.".format(requiredKey));

            if 'epsilon_static' in item2:
                hasEpsTensors = True;

    if hasEpsTensors:
        for item1 in dataSet['displacement_sets']:
            for item2 in item1['displacements']:
                if 'epsilon_static' not in item2:
                    raise Exception("Error: If dielectric tensors are supplied, they must be present for all displacements.");
                else:
                    dimOK = True;

                    if len(item2['epsilon_static']) != 3:
                        dimOK = False;
                    else:
                        for row in item2['epsilon_static']:
                            if len(row) != 3:
                                dimOK = False;

                    if not dimOK:
                        raise Exception("Error: One or more dielectric tensors in the data set does not have the correct dimension.");

    # Extract, reformat and return data.

    bandIndices = [];
    frequencies = [];
    dispStepSets, maxDispsSets = [], [];

    epsTensorSets = [] if hasEpsTensors else None;

    for item1 in dataSet['displacement_sets']:
        bandIndices.append(item1['band_index'] - 1);

        frequencies.append(item1['frequency']);

        dispStepSets.append(
            [item2['displacement_step'] for item2 in item1['displacements']]
            );

        maxDispsSets.append(
            [item2['max_cartesian_displacement'] for item2 in item1['displacements']]
            );

        if hasEpsTensors:
            epsTensorSet = [
                np.array(item2['epsilon_static'], dtype = np.float64)
                    for item2 in item1['displacements']
                ];

            epsTensorSets.append(epsTensorSet);

    return (
        dataSet['cell_volume'], dataSet['volume_units'],
        bandIndices, frequencies, dataSet['frequency_units'],
        dispStepSets, dataSet['step_units'], maxDispsSets, dataSet['distance_units'],
        epsTensorSets
        );

def WriteRamanDataSet(
    filePath,
    cellVolume, volumeUnits,
    bandIndices, frequencies, frequencyUnits,
    dispStepSets, stepUnits, maxDispsSets, distanceUnits,
    epsTensorSets = None
    ):
    """
    Write band indices, frequencies, displacement data and, optionally, calculated dielectric tensors to a YAML-format file.

    Arguments:
        filePath -- path to the output file.
        cellVolume -- volume of the unit cell (stored for normalising the Raman activity during post processing).
        bandIndices -- indices of the bands (modes) for which data is to be output.
        frequencies -- mode frequencies.
        frequencyUnits -- units of the frequencies.
        dispStepSets -- sets of displacements performed along each mode in normal-mode coordinates.
        stepUnits -- units of the normal-mode coordinate.
        maxDispsSets -- sets of maximum Cartesian atomic displacements across all atoms for each displacement.
        distanceUnits -- units of the maximum distance.

    Keyword arguments:
        epsTensorSets -- sets of 3x3 dielectric tensors calculated at each displacement.
    """

    with open(filePath, 'w') as outputWriter:
        # Write units for physical quantities.

        if frequencyUnits == None:
            frequencyUnits = _RamanDataSet_UnspecifiedUnitsDescription;

        outputWriter.write(
            "frequency_units: {0}\n".format(frequencyUnits)
            );

        if stepUnits == None:
            stepUnits = _RamanDataSet_UnspecifiedUnitsDescription;

        outputWriter.write(
            "step_units: {0}\n".format(stepUnits)
            );

        if distanceUnits == None:
            distanceUnits = _RamanDataSet_UnspecifiedUnitsDescription;

        outputWriter.write(
            "distance_units: {0}\n".format(distanceUnits)
            );

        if volumeUnits == None:
            volumeUnits = _RamanDataSet_UnspecifiedUnitsDescription;

        outputWriter.write(
            "volume_units: {0}\n".format(volumeUnits)
            );

        outputWriter.write("\n");

        # Write cell volume.

        outputWriter.write(
            "cell_volume: {0: 12.4f}\n".format(cellVolume)
            );

        outputWriter.write("\n");

        # Write displacement data for each mode.

        outputWriter.write("displacement_sets:\n");

        for i, bandIndex in enumerate(bandIndices):
            frequency = frequencies[i];
            dispSteps, maxDisps = dispStepSets[i], maxDispsSets[i];

            epsTensors = epsTensorSets[i] if epsTensorSets != None else None;

            outputWriter.write(
                "- # {0}\n".format(i + 1)
                );

            # Convert band indices to "one-based".

            outputWriter.write(
                "  band_index: {0}\n".format(bandIndex + 1)
                );

            outputWriter.write(
                "  frequency: {0}\n".format(frequency)
                );

            outputWriter.write("  displacements:\n");

            # Write displacement steps, max. Cartesian displacements, and calculated dielectric tensors, if available.

            for j in range(0, len(dispSteps)):
                outputWriter.write(
                    "   - # Step {0}\n".format(j + 1)
                    );

                outputWriter.write(
                    "     displacement_step: {0: 12.8f}\n".format(dispSteps[j])
                    );

                outputWriter.write(
                    "     max_cartesian_displacement: {0: 12.8f}\n".format(maxDisps[j])
                    );

                if epsTensors != None:
                    outputWriter.write("     epsilon_static:\n");

                    for k in range(0, 3):
                        outputWriter.write(
                            "     - [  {0: 15.8f},  {1: 15.8f},  {2: 15.8f}  ]\n".format(*[item for item in epsTensors[j][k]])
                            );
