# spectroscopy/yaml_export.py


# ---------
# Docstring
# ---------

""" Routines for exporting data in YAML format. """


# -------
# Imports
# -------

import numpy as np

from spectroscopy import constants


# ---------
# Functions
# ---------

def save_raman_tensors(
        band_indices, band_frequencies, raman_tensors, file_path,
        frequency_units=constants.DEFAULT_FREQUENCY_UNITS,
        activity_units=constants.DEFAULT_ACTIVITY_UNITS
        ):
    """ Write Raman tensors to a YAML-format file.

    Arguments:
        band_indices -- a list of band (mode) indices for which the
            Raman tensors were calculated.
        band_frequencies -- a list of band (mode) frequencies.
        raman_tensors -- a list of Raman tensors (3x3 matrices).
        file_path -- file to write the data to.

    Keyword arguments:
        frequency_units -- units of frequencies (defaults to
            constants.DEFAULT_FREQUENCY_UNITS).
        activity_units -- units of the Raman activity (defaults to
            constants.DEFAULT_ACTIVITY_UNITS).

    Notes:
        Indices in band_indices should be zero-based and will be
        incremented by 1 automatically when written to the output file.
    """

    if len(band_indices) == 0:
        raise Exception(
            "Error: band_indices must contain at least one element.")

    num_bands = len(band_indices)

    if len(band_frequencies) != num_bands:
        raise Exception(
            "Error: The number of band frequencies is inconsistent "
            "with the number of band indices.")

    if len(raman_tensors) != num_bands:
        raise Exception(
            "Error: The number of Raman tensors is inconsistent with "
            "the number of band indices.")

    for raman_tensor in raman_tensors:
        dim_1, dim_2 = np.shape(raman_tensor)

        if dim_1 != 3 or dim_2 != 3:
            raise Exception(
                "Error: Raman tensors should be 3x3 matrices.")

    # Write output file.

    with open(file_path, 'w') as output_writer:
        # In principle, we could build a dictionary structure and use
        # the yaml module to dump the data. However, the format of the
        # dumped files is often not terribly "human-readable", and it's
        # not hard to write short YAML files "by hand".

        output_writer.write("frequency_units: {0}\n".format(frequency_units))
        output_writer.write("activity_units: {0}\n".format(activity_units))
        output_writer.write("\n")

        output_writer.write("raman_activities:\n")

        for i, band_index in enumerate(band_indices):
            frequency, raman_tensor = band_frequencies[i], raman_tensors[i]

            output_writer.write("- #{0}\n".format(i))

            output_writer.write("  band_index: {0}\n".format(band_index + 1))
            output_writer.write("  frequency: {0: 16.8f}\n".format(frequency))

            output_writer.write("  raman_tensor:\n")

            for j in range(0, 3):
                output_writer.write(
                    "  - [  {0: 15.8f},  {1: 15.8f},  {2: 15.8f}  ]\n"
                    .format(*[item for item in raman_tensor[j]]))
