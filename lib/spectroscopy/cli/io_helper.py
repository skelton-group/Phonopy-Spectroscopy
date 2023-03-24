# Spectroscopy/cli/io_helper.py


# TODO: Only print message "spacer" if required.


# ---------
# Docstring
# ---------

""" Contains helper routines for importing and exporting data in a subset of
formats supported by the Interfaces module. """


# -------
# Imports
# -------

from spectroscopy.interfaces.vasp_interface import write_poscar, parse_outcar


# ----------------------
# Batch Structure Export
# ----------------------

def write_structures_raman(disp_sets, file_format, output_prefix=None):
    """ Write out a set of displaced structures for a Raman activity
    calculation.

    Arguments:
        disp_sets -- a list of (band_index, frequency, structures,
            disp_steps, max_disps) tuples providing sets of displaced
            structures for each mode to be calculated.
        file_format -- file format to write displaced structures into
            (currently only 'vasp_poscar' is supported).

    Keyword arguments:
        output_prefix -- prefix to prepend to output files.

    Notes:
        Band indices should be zero based.
    """

    # Check format support.

    if file_format not in ['vasp_poscar']:
        raise Exception(
            "Error: The file format '{0}' is not currently supported "
            "by this function.".format(file_format))

    # Set a prefix for the output files.

    if output_prefix is not None:
        output_prefix = "{0}_Raman".format(output_prefix)
    else:
        output_prefix = "Raman"

    if file_format == 'vasp_poscar':
        # Output structures in the VASP POSCAR format.

        for index, _, structures, disp_steps, max_disps in disp_sets:
            for i, (structure, disp_step, max_disp) in \
                    enumerate(zip(structures, disp_steps, max_disps)):
                file_name = "{0}-POSCAR.{1:0>4}.{2:0>3}.vasp".format(
                    output_prefix, index + 1, i + 1)

                title_line = (
                    "Mode = {0:0>4}, Disp. # = {1:0>3}, Step = "
                    "{2: 8.4f} (max. Cart. disp. = {3: 7.4f})"
                    .format(index + 1, i + 1, disp_step, max_disp))

                write_poscar(structure, file_name, title_line)

# ------------------------------
# Batch Dielectric Tensor Import
# ------------------------------

def read_dielectric_tensors(input_files, file_format):
    """ Read and return a list of dielectric tensors from a set of input
    files in the specified file format (currently only 'vasp_outcar' is
    supported). """

    # Check file format.

    if file_format not in ['vasp_outcar']:
        raise Exception(
            "Error: The file format '{0}' is currently not supported "
            "by this function.".format(file_format))

    # Read input files.

    eps_tensors = []

    if file_format == 'vasp_outcar':
        # Parse OUTCAR file and obtain static dielectric constant.

        for input_file in input_files:
            print("Reading input file: {0}".format(input_file))

            outcar_data = parse_outcar(
                input_file, extract_list=['epsilon_static'])

            eps_tensors.append(outcar_data['epsilon_static'])

        if len(eps_tensors) > 0:
            print("")

    return eps_tensors
