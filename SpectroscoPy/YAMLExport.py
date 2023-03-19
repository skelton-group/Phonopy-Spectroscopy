# SpectroscoPy/YAMLExport.py


# ---------
# Docstring
# ---------

""" Routines for exporting data in YAML format. """


# -------
# Imports
# -------

import numpy as np


# ---------
# Functions
# ---------

def saveramantensors(bandindices, bandfrequencies, ramantensors, filepath,
                     frequencyunits=None, activityunits=None):
    """
    Write Raman tensors to a YAML-format file.

    Arguments:
        bandIndices -- a list of band (mode) indices for which the Raman
        tensors were calculated.
        bandFrequencies -- a list of band (mode) frequencies.
        ramanTensors -- a list of Raman tensors (3x3 matrices).
        filePath -- file to write the data to.

    Keyword arguments:
        frequencyUnits -- units of frequencies (defaults to
        Constants.DefaultFrequencyUnits).
        activityUnits -- units of the Raman activity (defaults to
        Constants.DefaultActivityUnits).

    Notes:
        Indices in bandIndices should be zero-based and will be incremented by
        1 automatically when written to the
        output file.
    """

    if len(bandindices) == 0:
        raise Exception("Error: bandIndices must contain at least one "
                        "element."
                        )

    numbands = len(bandindices)

    if len(bandfrequencies) != numbands:
        raise Exception("Error: The number of band frequencies is "
                        "inconsistent with the number of band indices."
                        )

    if len(ramantensors) != numbands:
        raise Exception("Error: The number of Raman tensors is inconsistent "
                        "with the number of band indices."
                        )

    for ramantensor in ramantensors:
        dim1, dim2 = np.shape(ramantensor)

        if dim1 != 3 or dim2 != 3:
            raise Exception("Error: Raman tensors should be 3x3 matrices.")

    # Set default values for frequency and activity units if required.

    if frequencyunits is None:
        frequencyunits = Constants.DefaultFrequencyUnits

    if activityunits is None:
        activityunits = Constants.DefaultActivityUnits

    # Write output file.

    with open(filepath, 'w') as outputWriter:
        # In principle, we could build a dictionary structure and use the yaml
        # module to dump the data.
        # However, the format of the dumped files is often not terribly
        # "human-readable", and it's not hard to write short YAML files "by
        # hand".

        outputWriter.write(
            "frequency_units: {0}\n".format(frequencyunits)
            )

        outputWriter.write(
            "activity_units: {0}\n".format(activityunits)
            )

        outputWriter.write("\n")

        outputWriter.write("raman_activities:\n")

        for i, bandindex in enumerate(bandindices):
            frequency, ramantensor = bandfrequencies[i], ramantensors[i]

            outputWriter.write(
                "- #{0}\n".format(i)
                )

            outputWriter.write(
                "  band_index: {0}\n".format(bandindex + 1)
                )

            outputWriter.write(
                "  frequency: {0: 16.8f}\n".format(frequency)
                )

            outputWriter.write("  raman_tensor:\n")

            for j in range(0, 3):
                outputWriter.write(
                    "  - [  {0: 15.8f},  {1: 15.8f},  {2: 15.8f}  ]\n"
                    .format(*[item for item in ramantensor[j]])
                    )
