# spectroscopy/yaml_export.py


# ---------
# Docstring
# ---------

""" Routines for exporting data in YAML format. """


# -------
# Imports
# -------

import numpy as np

from spectroscopy.utilities import numpy_check_shape


# ---------
# Functions
# ---------

def save_raman_tensors(
        band_indices, frequencies, raman_tensors, frequency_unit,
        activity_unit, file_path,
        surface_hkl=None, theta=None, surface_rotation=None,
        theta_rotation=None, combined_rotation=None
        ):
    """ Write Raman tensors to a YAML-format file.

    Arguments:
        band_indices -- a list of band (mode) indices for which the
            Raman tensors were calculated.
        band_frequencies -- a list of band (mode) frequencies.
        raman_tensors -- a list of Raman tensors (3x3 matrices).
        frequency_unit, activity_unit -- units of frequency and Raman
            activity.
        file_path -- file to save to.
    
    Keyword arguments (polarised Raman calculations):
        surface_hkl -- (h, k, l) values specifying the Miller indices
            of the surface aligned to the z direction.
        theta -- rotation angle in the xy plane (degrees).
        surface_rotation, theta_rotation, total_rotation -- 3x3 rotation
            matrices specifying (1) the rotation required to align the
            specified surface to the z direction; (2) the rotation for
            the specified theta; and (3) the combined rotation applied
            to the Raman tensors.
    
    Notes:
        * Indices in band_indices should be zero-based and will be
            incremented by 1 automatically when written to the output
            file.
        * If surface_hkl is specified, surface_rotation must also be
            specified; if theta is also supplied, then theta_rotation
            and total_rotation must also be supplied. 
    """
    
    if len(band_indices) == 0:
        raise Exception(
            "band_indices must contain at least one element.")

    num_bands = len(band_indices)

    if len(frequencies) != num_bands:
        raise Exception(
            "The number of band frequencies is inconsistent with the "
            "wnumber of band indices.")

    if not numpy_check_shape(raman_tensors, (num_bands, 3, 3)):
            raise Exception(
                "Raman tensors should specify one 3x3 matrix for each "
                "mode.")
    
    if surface_hkl is not None:
        if not numpy_check_shape(surface_hkl, (3, )):
            raise Exception(
                "If supplied, surface_hkl must specify three values.")
        
        h, k, l = surface_hkl

        try:
            int(h)
            int(k)
            int(l)
        except ValueError:
            raise Exception(
                "If supplied, the (h, k, l) in surface_hkl must be "
                "integers.")

        if surface_rotation is None:
            raise Exception(
                "If surface_hkl is supplied, surface_rotation must "
                "also be supplied.")

        if theta is not None:
            if theta_rotation is None or combined_rotation is None:
                raise Exception(
                    "If surface_hkl and theta are supplied, "
                    "theta_rotation and combined_rotation must also be "
                    "supplied.")
        
        for m in surface_rotation, theta_rotation, combined_rotation:
            if m is not None and not numpy_check_shape(m, (3, 3)):
                raise Exception(
                    "Rotation matrices must be 3x3 matrices.")
    
    # Write output file.

    with open(file_path, 'w') as output_writer:
        # In principle, we could build a dictionary structure and use
        # the yaml module to dump the data. However, the format of the
        # dumped files is often not terribly "human-readable", and it's
        # not hard to write short YAML files "by hand".

        output_writer.write("frequency_units: {0}\n".format(frequency_unit))
        output_writer.write("activity_units: {0}\n".format(activity_unit))

        if surface_hkl is not None and theta is not None:
            output_writer.write("angle_units: deg\n")

        output_writer.write("\n")

        if surface_hkl is not None:
            output_writer.write(
                "surface_hkl: [ {0}, {1}, {2} ]\n".format(*surface_hkl))
            
            if theta is not None:
                output_writer.write("theta: {0:.3f}\n".format(theta))
            
            output_writer.write("\n")

            matrix_output_list = [('surface_rotation', surface_rotation)]

            if theta is not None:
                matrix_output_list.append(('theta_rotation', theta_rotation))
                matrix_output_list.append(
                    ('comvined_rotation', combined_rotation))
            
            for tag, matrix in matrix_output_list:
                output_writer.write("{0}:\n".format(tag))

                for row in matrix:
                    output_writer.write(
                        "  - [  {0: 12.8f},  {1: 12.8f},  {2: 12.8f}  ]\n"
                        .format(*row))
            
            output_writer.write("\n")

        output_writer.write("raman_activities:\n")

        for i, band_index in enumerate(band_indices):
            frequency, raman_tensor = frequencies[i], raman_tensors[i]

            output_writer.write("- #{0}\n".format(i))

            output_writer.write("  band_index: {0}\n".format(band_index + 1))
            output_writer.write("  frequency: {0: 16.8f}\n".format(frequency))

            output_writer.write("  raman_tensor:\n")

            for row in raman_tensor:
                output_writer.write(
                    "  - [  {0: 15.8f},  {1: 15.8f},  {2: 15.8f}  ]\n"
                    .format(*row))
