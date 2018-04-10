# SpectroscoPy/CLI/Phonopy.py


# ---------
# Docstring
# ---------

""" Shared routines for the phonopy-* scripts, for handling multiple input sources. """


# -------
# Imports
# -------

import os;

import yaml;

from phonopy.structure.atoms import atom_data;

from SpectroscoPy.Interfaces.Phonopy import ReadPhonopyYAML, ReadMeshYAML, ReadIrRepsYAML, ReadMeshHDF5, ReadPhono3pyHDF5;
from SpectroscoPy.Interfaces.VASP import ReadPOSCAR;


# ------------------------------
# Command-Line Argument Handling
# ------------------------------

def Phonopy_UpdateParser(parser, spectrumType):
    """ Update parser with Phonopy-specific command-line options for spectrumType. """

    # Options for input of structure, atomic masses, frequencies/eigenvectors, and Born effective-charge tensors (for IR spectra).

    group = parser.add_argument_group("Input files");

    group.add_argument(
        "--mesh_yaml",
        metavar = "<file_path>",
        type = str, dest = "MeshYAML",
        default = "mesh.yaml",
        help = "mesh.yaml file to read Gamma-point frequencies/eigenvectors and structure, if present, from (default: mesh.yaml)"
        );

    group.add_argument(
        "--mesh_hdf5",
        metavar = "<file_path>",
        type = str, dest = "MeshHDF5",
        default = "mesh.hdf5",
        help = "mesh.hdf5 file to read Gamma-point frequencies and eigenvectors from (default: mesh.hdf5)"
        );

    group.add_argument(
        "--phonopy_yaml",
        metavar = "<file_path>",
        type = str, dest = "PhonopyYAML",
        default = "phonopy.yaml",
        help = "phonopy.yaml file to read structure and atomic masses from (optional; default: phonopy.yaml)"
        );

    group.add_argument(
        "-c", "--cell",
        metavar = "<file_path>",
        type = str, dest = "CellFile",
        default = "POSCAR",
        help = "Input file to read structure from (optional; VASP 5.x POSCAR files only; default: POSCAR)"
        );

    group.add_argument(
        "--ir_reps_yaml",
        metavar = "<file_path>",
        type = str, dest = "IrRepsYAML",
        default = "irreps.yaml",
        help = "irreps.yaml file to read if using the --ir_reps argument (default: irreps.yaml)"
        );

    if spectrumType == 'ir':
        group.add_argument(
            "--born_file",
            metavar = "<file_path>",
            type = str, dest = "BORNFile",
            default = "BORN",
            help = "BORN file to read Born charges from (default: BORN)"
            );

    # Options for importing Phono3py linewidths.

    group = parser.add_argument_group("Linewidths");

    group.add_argument(
        "--linewidth_hdf5",
        metavar = "<file_path>",
        type = str, dest = "LinewidthHDF5",
        default = None,
        help = "Phono3py HDF5 file to read mode linewidths from"
        );

    group.add_argument(
        "--linewidth_temperature",
        metavar = "<temperature>",
        type = float, dest = "LinewidthTemperature",
        default = None,
        help = "Temperature to read linewidths at"
        );


# ------------
# I/O Handling
# ------------

def Phonopy_LoadData_Core(args, extractList):
    """
    Reads input files specified in the parsed command-line arguments args and loads items specified in extractList.

    exctractList may contain the following keys:
        'phonon_modes' -- phonon frequencies (in THz, rad. THz, inv. cm and meV) and eigenvectors.
        'structure' -- structure specified as a tuple of (lattice_vectors, atomic_symbols, atom_positions_frac) data.
        'atomic_masses' -- list of atomic masses.

    If items in the extract list cannot be taken from the supplied input file(s), an error is raised.

    Return value:
        A dictionary containing the requested data, keyed with the entries in extractList.
    """

    # We have several options for reading in a structure, a list of atomic masses and a set of Gamma-point frequencies and eigenvectors from Phonopy output.
    # If the user specifies input files, we use those preferentially; if not, we check for standard files in order of ease of reading.

    # Data to capture.

    structure, atomicMasses, phononModes = None, None, None;

    # Variable to keep track of whether we've attempted to read the mesh.yaml file.
    # If mesh.yaml is used as the source of frequencies and eigenvectors, but is an old-format file without structural data, we want to avoid reading it a second time to try and load a structure.

    readMeshYAML = False;

    # Variable to keep track of whether status messages have been printed (used to control printing a "spacer" for output formatting).

    printSpacer = False;

    # Read frequencies and eigenvectors if required.

    if 'phonon_modes' in extractList:
        if os.path.isfile(args.MeshHDF5):
            phononModes = ReadMeshHDF5(args.MeshHDF5);

            print("Frequencies and eigenvectors read from \"{0}\"".format(args.MeshHDF5));

            printSpacer = True;
        elif os.path.isfile(args.MeshYAML):
            phononModes, structure, atomicMasses = ReadMeshYAML(args.MeshYAML);

            print("Frequencies and eigenvectors read from \"{0}\"".format(args.MeshYAML));

            if structure != None:
                print("Structure read from \"{0}\"".format(args.MeshYAML));

            if atomicMasses != None:
                print("Atomic masses read from \"{0}\"".format(args.MeshYAML));

            readMeshYAML = True;

            printSpacer = True;

        # Check we were able to load them.

        if 'phonon_modes' in extractList and phononModes == None:
            raise Exception("Error: Unable to load frequencies and eigenvectors - please check input files.");

    # Read structure if required.

    if 'structure' in extractList and structure == None:
        # If available, read the structure and atomic masses from a phonopy.yaml file.

        if os.path.isfile(args.PhonopyYAML):
            # In some cases, spglib returns quotes in the Hall symbol, and Phonopy writes the symbol into the Phonopy YAML file without escaping, which prevents the YAML being parsed correctly.
            # As a workaround, we catch the yaml.scanner.ScannerError, print a warning message and continue.

            try:
                structure, atomicMasses = ReadPhonopyYAML(args.PhonopyYAML);

                # The ReadPhonopyYAML() routine fails if the correct data is not present.

                print("Structure read from \"{0}\"".format(args.PhonopyYAML));
                print("Atomic masses read from \"{0}\"".format(args.PhonopyYAML));

            except yaml.scanner.ScannerError:
                print("WARNING: Error parsing \"{0}\" -> try to read structure from other input files.".format(args.PhonopyYAML));

                printSpacer = True;

        if structure == None:
            # Next, try to read the structure from a cell file (currently only VASP 5.x POSCAR files are supported).

            if os.path.isfile(args.CellFile):
                _, latticeVectors, atomicSymbols, atomPositions = ReadPOSCAR(args.CellFile);

                structure = (latticeVectors, atomicSymbols, atomPositions);

                print("Structure read from \"{0}\"".format(args.CellFile));
                printSpacer = True;
            else:
                # Finally, try and read from a mesh.yaml file if we haven't already done so.
                # If the file contains data for a dense q-point mesh and/or has eigenvectors, this could be very slow, so we only do it as a last resort.

                if not readMeshYAML and os.path.isfile(args.MeshYAML):
                    _, structure, atomicMasses = ReadMeshYAML(args.MeshYAML);

                    if structure != None:
                        print("Structure read from \"{0}\"".format(args.MeshYAML));
                        printSpacer = True;

                    if atomicMasses != None:
                        print("Atomic masses read from \"{0}\"".format(args.MeshYAML));
                        printSpacer = True;

        # Check we have a structure if we need one.

        if 'structure' in extractList and structure == None:
            raise Exception("Error: Unable to read a structure - please check input files.");

    # Finally, if we need atomic masses and we haven't read them from an input file, get them from the phonopy.structure.atoms database.

    if 'atomic_masses' in extractList and atomicMasses == None:
        _, atomicSymbols, _ = structure;

        atomicMasses = [];

        for symbol in atomicSymbols:
            symbol = symbol.title();

            mass = None;

            for _, dbSymbol, _, dbMass in atom_data:
                if symbol == dbSymbol:
                    mass = dbMass;
                    break;

            if mass == None:
                raise Exception("Error: Atomic mass for symbol '{0}' not found in Phonopy database.".format(symbol));

            atomicMasses.append(mass);

        print("Atomic masses taken from Phonopy database");
        printSpacer = True;

    # If required, print a blank line after status messages.

    if printSpacer:
        print("");

    # Return a dictonary containing the requested data.

    outputData = { };

    if 'phonon_modes' in extractList:
        outputData['phonon_modes'] = phononModes;

    if 'structure' in extractList:
        outputData['structure'] = structure;

    if 'atomic_masses' in extractList:
        outputData['atomic_masses'] = atomicMasses;

    return outputData;

def Phonopy_LoadData_Optional(args):
    """
    Reads and returns optional irreducible representation and linewidth data if the relevant arguments are set in the supplied command-line arguments.

    Return value:
        A tuple of (linewidths, ir_rep_data) containing the linewidths and/or irreducible representations.
        If the associated optional parameters are not set in args, either or both will be set to None.
    """

    printSpacer = False;

    # Read linewidths if required.

    linewidths = None;

    if args.LinewidthHDF5:
        if args.LinewidthTemperature == None:
            raise Exception("Error: If reading linewidths from a Phono3py HDF5 file, a temperature must be supplied using the --linewidth_temperature argument.");

        linewidths = ReadPhono3pyHDF5(args.LinewidthHDF5, args.LinewidthTemperature);

        print("Linewidths read from \"{0}\" (T = {1:.2f} K)".format(args.LinewidthHDF5, args.LinewidthTemperature));
        printSpacer = True;

    # Read ir. rep. data if required.

    irRepData = None;

    if args.IrReps:
        irRepData = ReadIrRepsYAML(filePath = args.IrRepsYAML);

        print("Irreducible representations read from \"{0}\"".format(args.IrRepsYAML));
        printSpacer = True;

    # Spacer after status messages if required.

    if printSpacer:
        print("");

    # Return data.

    return (linewidths, irRepData);
