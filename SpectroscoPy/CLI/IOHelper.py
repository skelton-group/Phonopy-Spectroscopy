# SpectroscoPy/CLI/IOHelper.py


# TODO: Only print message "spacer" if required.


# ---------
# Docstring
# ---------

""" Contains helper routines for importing and exporting data in a subset of
formats supported by the Interfaces module. """


# -------
# Imports
# -------

from SpectroscoPy.Interfaces.VASP import writeposcar, parseoutcar


# ----------------------
# Batch Structure Export
# ----------------------

def writestructuresraman(dispsets, fileformat, outputprefix=None):
    """
    Write out a set of displaced structures for a Raman activity calculation.

    Arguments:
        dispSets -- a list of (band_index, frequency, structures,
        disp_steps, max_disps) tuples providing sets of displaced structures
        for each mode to be calculated.
        fileFormat -- file format to write displaced structures into (
        currently only 'vasp_poscar' is supported).

    Keyword arguments:
        outputPrefix -- prefix to prepend to output files.

    Notes:
        Band indices should be zero based.
    """

    # Check format support.

    if fileformat not in ['vasp_poscar']:
        raise Exception("Error: The file format '{0}' is not currently "
                        "supported by this function.".format(fileformat)
                        )

    # Set a prefix for the output files.

    if outputprefix is not None:
        outputprefix = "{0}_Raman".format(outputprefix)
    else:
        outputprefix = "Raman"

    if fileformat == 'vasp_poscar':
        # Output structures in the VASP POSCAR format.

        for index, _, structures, dispsteps, maxdisps in dispsets:
            for i, (structure, dispstep, maxdisp) in enumerate(zip(
                    structures, dispsteps, maxdisps)):
                filename = "{0}-POSCAR.{1:0>4}.{2:0>3}.vasp".format(
                    outputprefix, index + 1, i + 1
                )

                titleline = "Mode = {0:0>4}, Disp. # = {1:0>3}, " \
                            "Step = {2: 8.4f} (max. Cart. disp. " \
                            "= {3: 7.4f})".format(
                             index + 1, i + 1, dispstep, maxdisp
                            )

                writeposcar(structure, filename, titleline)


# ------------------------------
# Batch Dielectric Tensor Import
# ------------------------------

def read_dielectric_tensors(inputfiles, fileformat):
    """ Read and return a list of dielectric tensors from a set of input
    files in the specified file format (currently only 'vasp_outcar' is
    supported). """

    # Check file format.

    if fileformat not in ['vasp_outcar']:
        raise Exception("Error: The file format '{0}' is currently not "
                        "supported by this function.".format(fileformat)
                        )

    # Read input files.

    epstensors = []

    if fileformat == 'vasp_outcar':
        # Parse OUTCAR file and obtain static dielectric constant.

        for inputfile in inputfiles:
            print("Reading input file: {0}".format(inputfile))

            outcardata = parseoutcar(
                inputfile, extractlist=['epsilon_static']
                )

            epstensors.append(
                outcardata['epsilon_static']
                )

        if len(epstensors) > 0:
            print("")

    return epstensors
