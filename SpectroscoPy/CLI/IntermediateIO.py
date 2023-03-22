# SpectroscoPy/CLI/IntermediateIO.py


# ---------
# Docstring
# ---------

""" Routines for storing and retrieving intermediate data in multi-step
calculations. """


# -------
# Imports
# -------

import yaml

# Try to use the C YAML Loader if possible; if not, fall back to the
# standard routines.

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

import numpy as np


# ---------
# Constants
# ---------

""" Default step unit written to intermediate YAML files. """

_DefaultStepUnit = "Unspecified"

""" Default distance unit written to intermediate YAML files. """

_DefaultDistanceUnit = "Unspecified"


# ------------------
# Raman Calculations
# ------------------

""" Default string written to intermediate YAML files when a unit is 
unspecified. """

_RamanDataSet_UnspecifiedUnitsDescription = "Unspecified"


def readramandataset(filepath):
    """
    Read data for a Raman calculation from a YAML-format file produced by
    WriteRamanDataSet().

    Return value:
        A tuple of (cell_volume, volume_units, band_indices, frequencies,
        frequency_units, disp_step_sets, step_units, max_disps_sets,
        distance_units, eps_tensor_sets) data read from the file.
        (The order of the elements matches the arguments passed to
        WriteRamanDataSet()).
    """

    # Read and parse data set.

    dataset = None

    with open(filepath, 'r') as inputReader:
        dataset = yaml.load(inputReader, Loader=Loader)

    # Validate.

    for requiredkey in ['frequency_units', 'step_units', 'distance_units',
                        'volume_units', 'cell_volume', 'displacement_sets']:
        if requiredkey not in dataset:
            raise Exception("Error: Missing key '{0}' in displacement data "
                            "set.".format(requiredkey)
                            )

    if len(dataset['displacement_sets']) == 0:
        raise Exception("Error: Displacement dataset contains no "
                        "displacement sets.")

    hasepstensors = False

    for item1 in dataset['displacement_sets']:
        for requiredkey in ['band_index', 'frequency', 'displacements']:
            if requiredkey not in item1:
                raise Exception("Error: Missing key '{0}' in one or more "
                                "displacement sets.".format(requiredkey)
                                )

        for item2 in item1['displacements']:
            for requiredkey in ['displacement_step',
                                'max_cartesian_displacement']:
                if requiredkey not in item2:
                    raise Exception("Error: Missing key '{0}' in one or "
                                    "more displacement entries.".format(
                                                            requiredkey)
                                    )

            if 'epsilon_static' in item2:
                hasepstensors = True

    if hasepstensors:
        for item1 in dataset['displacement_sets']:
            for item2 in item1['displacements']:
                if 'epsilon_static' not in item2:
                    raise Exception("Error: If dielectric tensors are "
                                    "supplied, they must be present for all "
                                    "displacements.")
                else:
                    dimok = True

                    if len(item2['epsilon_static']) != 3:
                        dimok = False
                    else:
                        for row in item2['epsilon_static']:
                            if len(row) != 3:
                                dimok = False

                    if not dimok:
                        raise Exception("Error: One or more dielectric "
                                        "tensors in the data set does not "
                                        "have the correct dimension.")

    # Extract, reformat and return data.

    bandindices = []
    frequencies = []
    dispstepsets, maxdispssets = [], []

    epstensorsets = [] if hasepstensors else None

    for item1 in dataset['displacement_sets']:
        bandindices.append(item1['band_index'] - 1)

        frequencies.append(item1['frequency'])

        dispstepsets.append(
            [item2['displacement_step'] for item2 in item1['displacements']]
            )

        maxdispssets.append(
            [item2['max_cartesian_displacement'] for item2 in item1[
                'displacements']]
            )

        if hasepstensors:
            epstensorset = [
                np.array(item2['epsilon_static'], dtype=np.float64)
                for item2 in item1['displacements']
                ]

            epstensorsets.append(epstensorset)

    return (
        dataset['cell_volume'], dataset['volume_units'],
        bandindices, frequencies, dataset['frequency_units'],
        dispstepsets, dataset['step_units'], maxdispssets, dataset[
            'distance_units'],
        epstensorsets
        )


def writeramandataset(filepath, cellvolume, volumeunits, bandindices,
                      frequencies, frequencyunits, dispstepsets, stepunits,
                      maxdispssets, distanceunits, epstensorsets=None):
    """
    Write band indices, frequencies, displacement data and, optionally,
    calculated dielectric tensors to a YAML-format file.

    Arguments:
        filePath -- path to the output file.
        cellVolume -- volume of the unit cell (stored for normalising the
        Raman activity during post processing).
        bandIndices -- indices of the bands (modes) for which data is to be
        output.
        frequencies -- mode frequencies.
        frequencyUnits -- units of the frequencies.
        dispStepSets -- sets of displacements performed along each mode in
        normal-mode coordinates.
        stepUnits -- units of the normal-mode coordinate.
        maxDispsSets -- sets of maximum Cartesian atomic displacements
        across all atoms for each displacement.
        distanceUnits -- units of the maximum distance.

    Keyword arguments:
        epsTensorSets -- sets of 3x3 dielectric tensors calculated at each
        displacement.
    """

    with open(filepath, 'w') as outputWriter:
        # Write units for physical quantities.

        if frequencyunits is None:
            frequencyunits = _RamanDataSet_UnspecifiedUnitsDescription

        outputWriter.write(
            "frequency_units: {0}\n".format(frequencyunits)
            )

        if stepunits is None:
            stepunits = _RamanDataSet_UnspecifiedUnitsDescription

        outputWriter.write(
            "step_units: {0}\n".format(stepunits)
            )

        if distanceunits is None:
            distanceunits = _RamanDataSet_UnspecifiedUnitsDescription

        outputWriter.write(
            "distance_units: {0}\n".format(distanceunits)
            )

        if volumeunits is None:
            volumeunits = _RamanDataSet_UnspecifiedUnitsDescription

        outputWriter.write(
            "volume_units: {0}\n".format(volumeunits)
            )

        outputWriter.write("\n")

        # Write cell volume.

        outputWriter.write(
            "cell_volume: {0: 12.4f}\n".format(cellvolume)
            )

        outputWriter.write("\n")

        # Write displacement data for each mode.

        outputWriter.write("displacement_sets:\n")

        for i, bandindex in enumerate(bandindices):
            frequency = frequencies[i]
            dispsteps, maxdisps = dispstepsets[i], maxdispssets[i]

            epstensors = epstensorsets[i] if epstensorsets is not None else \
                None

            outputWriter.write(
                "- # {0}\n".format(i + 1)
                )

            # Convert band indices to "one-based".

            outputWriter.write(
                "  band_index: {0}\n".format(bandindex + 1)
                )

            outputWriter.write(
                "  frequency: {0}\n".format(frequency)
                )

            outputWriter.write("  displacements:\n")

            # Write displacement steps, max. Cartesian displacements,
            # and calculated dielectric tensors, if available.

            for j in range(0, len(dispsteps)):
                outputWriter.write(
                    "   - # Step {0}\n".format(j + 1)
                    )

                outputWriter.write(
                    "     displacement_step: {0: 12.8f}\n".format(dispsteps[j])
                    )

                outputWriter.write(
                    "     max_cartesian_displacement: {0: 12.8f}\n".format(
                        maxdisps[j])
                    )

                if epstensors is not None:
                    outputWriter.write("     epsilon_static:\n")

                    for k in range(0, 3):
                        outputWriter.write(
                            "     - [  {0: 15.8f},  {1: 15.8f},  {2: 15.8f}  "
                            "]\n".format(*[item for item in epstensors[j][k]])
                            )
