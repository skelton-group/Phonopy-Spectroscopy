# SpectroscoPy/CLI/IOHelper.py


# TODO: Only print message "spacer" if required.


# ---------
# Docstring
# ---------

""" Contains helper routines for importing and exporting data in a subset of formats supported by the Interfaces module. """


# -------
# Imports
# -------

from SpectroscoPy.Interfaces.VASP import WritePOSCAR, ParseOUTCAR;


# ----------------------
# Batch Structure Export
# ----------------------

def WriteStructuresRaman(dispSets, fileFormat, outputPrefix = None):
    """
    Write out a set of displaced structures for a Raman activity calculation.

    Arguments:
        dispSets -- a list of (band_index, frequency, structures, disp_steps, max_disps) tuples providing sets of displaced structures for each mode to be calculated.
        fileFormat -- file format to write displaced structures into (currently only 'vasp_poscar' is supported).

    Keyword arguments:
        outputPrefix -- prefix to prepend to output files.

    Notes:
        Band indices should be zero based.
    """

    # Check format support.

    if fileFormat not in ['vasp_poscar']:
        raise Exception("Error: The file format '{0}' is not currently supported by this function.".format(fileFormat));

    # Set a prefix for the output files.

    if outputPrefix != None:
        outputPrefix = "{0}_Raman".format(outputPrefix);
    else:
        outputPrefix = "Raman";

    if fileFormat == 'vasp_poscar':
        # Output structures in the VASP POSCAR format.

        for index, _, structures, dispSteps, maxDisps in dispSets:
            for i, (structure, dispStep, maxDisp) in enumerate(zip(structures, dispSteps, maxDisps)):
                fileName = "{0}-POSCAR.{1:0>4}.{2:0>3}.vasp".format(outputPrefix, index + 1, i + 1);

                titleLine = "Mode = {0:0>4}, Disp. # = {1:0>3}, Step = {2: 8.4f} (max. Cart. disp. = {3: 7.4f})".format(index + 1, i + 1, dispStep, maxDisp);

                WritePOSCAR(structure, fileName, titleLine);


# ------------------------------
# Batch Dielectric Tensor Import
# ------------------------------

def ReadDielectricTensors(inputFiles, fileFormat):
    """ Read and return a list of dielectric tensors from a set of input files in the specified file format (currently only 'vasp_outcar' is supported). """

    # Check file format.

    if fileFormat not in ['vasp_outcar']:
        raise Exception("Error: The file format '{0}' is currently not supported by this function.".format(fileFormat));

    # Read input files.

    epsTensors = [];

    if fileFormat == 'vasp_outcar':
        # Parse OUTCAR file and obtain static dielectric constant.

        for inputFile in inputFiles:
            print("Reading input file: {0}".format(inputFile));

            outcarData = ParseOUTCAR(
                inputFile, extractList = ['epsilon_static']
                );

            epsTensors.append(
                outcarData['epsilon_static']
                );

        if len(epsTensors) > 0:
            print("");

    return epsTensors;
