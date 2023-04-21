# spectroscopy/cli/phonopy_helper.py


# ---------
# Docstring
# ---------

""" Shared routines for the phonopy-* scripts, for handling multiple input
sources. """


# -------
# Imports
# -------

import os

import yaml

from phonopy.structure.atoms import atom_data

from spectroscopy.interfaces.phonopy_interface import (
    read_phonopy_yaml, read_mesh_yaml, read_irreps_yaml, read_mesh_hdf5,
    read_phono3py_hdf5)

from spectroscopy.interfaces.vasp_interface import read_poscar

from spectroscopy.cli.parser import update_parser


# ------------------------------
# Command-Line Argument Handling
# ------------------------------

def phonopy_update_parser(parser, spectrum_type):
    """ Update parser with Phonopy-specific command-line options for
    spectrum_type. """

    # Options for input of structure, atomic masses,
    # frequencies/eigenvectors, and Born effective-charge tensors (for
    # IR spectra).

    group = parser.add_argument_group("Input files")

    group.add_argument(
        "--mesh-yaml",
        metavar="<file_path>", type=str, dest="MeshYAML", default="mesh.yaml",
        help="mesh.yaml file to read Gamma-point "
             "frequencies/eigenvectors and structure, if present, from "
             "(default: mesh.yaml)")

    group.add_argument(
        "--mesh-hdf5",
        metavar="<file_path>", type=str, dest="MeshHDF5", default="mesh.hdf5",
        help="mesh.hdf5 file to read Gamma-point frequencies and "
             "eigenvectors from (default: mesh.hdf5)")

    group.add_argument(
        "--phonopy-yaml",
        metavar="<file_path>",
        type=str, dest="PhonopyYAML", default="phonopy.yaml",
        help="phonopy.yaml file to read structure and atomic masses "
             "from (optional; default: phonopy.yaml)")

    group.add_argument(
        "-c", "--cell",
        metavar="<file_path>",
        type=str, dest="CellFile", default="POSCAR",
        help="Input file to read structure from (optional; VASP 5.x "
             "POSCAR files only; default: POSCAR)")

    group.add_argument(
        "--irreps-yaml",
        metavar="<file_path>",
        type=str, dest="IrrepsYAML", default="irreps.yaml",
        help="irreps.yaml file to read if using the --irreps argument "
             "(default: irreps.yaml)")

    if spectrum_type == 'ir':
        group.add_argument(
            "--born-file",
            metavar="<file_path>", type=str, dest="BORNFile", default="BORN",
            help="BORN file to read Born charges from (default: BORN)")

    # Options for importing Phono3py linewidths.

    group = parser.add_argument_group("Linewidths")

    group.add_argument(
        "--linewidth-hdf5",
        metavar="<file_path>", type=str, dest="LinewidthHDF5", default=None,
        help="Phono3py HDF5 file to read mode linewidths from")

    group.add_argument(
        "--linewidth-temperature",
        metavar="<temperature>",
        type=float, dest="LinewidthTemperature", default=300,
        help="Temperature to read linewidths at (default: 300 K)")
    
    # Hand over to the "general" update_paraer() routine to add other
    # arguments.

    update_parser(parser, spectrum_type=spectrum_type)


# ------------
# I/O Handling
# ------------

def phonopy_load_core(args, extract_list):
    """ Reads input files specified in the parsed command-line arguments
    args and loads items specified in extract_list.

    exctract_list may contain the following keys:
        'phonon_modes' -- phonon frequencies (in THz, rad. THz, inv. cm
            and meV) and eigenvectors.
        'structure' -- structure specified as a tuple of
            (lattice_vectors, atomic_symbols, atom_positions_frac) data.
        'atomic_masses' -- list of atomic masses.

    If items in the extract list cannot be taken from the supplied input
    file(s), an error is raised.

    Return value:
        A dictionary containing the requested data, keyed with the
        entries in extract_list.
    """

    # We have several options for reading in a structure, a list of
    # atomic masses and a set of Gamma-point frequencies and
    # eigenvectors from Phonopy output. If the user specifies input
    # files, we use those preferentially; if not, we check for standard
    # files in order of "ease of reading".

    # Data to capture.

    structure, atomic_masses, phonon_modes = None, None, None

    # Variable to keep track of whether we've attempted to read the
    # mesh.yaml file. If mesh.yaml is used as the source of frequencies
    # and eigenvectors, but is an old-format file without structural
    # data, we want to avoid reading it a second time to try and load a
    # structure.

    mesh_yaml_read = False

    # Variable to keep track of whether status messages have been
    # printed (used to control printing a "spacer" for output
    # formatting).

    print_spacer = False

    # Read frequencies and eigenvectors if required.

    if 'phonon_modes' in extract_list:
        if os.path.isfile(args.MeshHDF5):
            phonon_modes = read_mesh_hdf5(args.MeshHDF5)

            print(
                "Frequencies and eigenvectors read from \"{0}\""
                .format(args.MeshHDF5))

            print_spacer = True
        elif os.path.isfile(args.MeshYAML):
            phonon_modes, structure, atomic_masses = read_mesh_yaml(
                args.MeshYAML)

            print("Frequencies and eigenvectors read from \"{0}\""
                .format(args.MeshYAML))

            if structure is not None:
                print("Structure read from \"{0}\"".format(args.MeshYAML))

            if atomic_masses is not None:
                print("Atomic masses read from \"{0}\"".format(args.MeshYAML))

            mesh_yaml_read = True

            print_spacer = True

        # Check we were able to load them.

        if 'phonon_modes' in extract_list and phonon_modes is None:
            raise Exception(
                "Error: Unable to load frequencies and eigenvectors - "
                "please check input files.")

    # Read structure if required.

    if 'structure' in extract_list and structure is None:
        # If available, read the structure and atomic masses from a
        # phonopy.yaml file.

        if os.path.isfile(args.PhonopyYAML):
            # In some cases, spglib returns quotes in the Hall symbol,
            # and Phonopy writes the symbol into the Phonopy YAML file
            # without escaping, which prevents the YAML being parsed
            # correctly. As a workaround, we catch the
            # yaml.scanner.ScannerError, print a warning message and
            # continue.

            try:
                structure, atomic_masses = read_phonopy_yaml(args.PhonopyYAML)

                # The read_phonopy_yaml() routine fails if the correct 
                # data is not present.

                print("Structure read from \"{0}\"".format(args.PhonopyYAML))
                
                print(
                    "Atomic masses read from \"{0}\""
                    .format(args.PhonopyYAML))

            except yaml.scanner.ScannerError:
                print(
                    "WARNING: Error parsing \"{0}\" -> try to read "
                    "structure from other input files."
                    .format(args.PhonopyYAML))

                print_spacer = True

        if structure is None:
            # Next, try to read the structure from a cell file
            # (currently only VASP 5.x POSCAR files are supported).

            if os.path.isfile(args.CellFile):
                _, lattice_vectors, atomic_symbols, atom_positions = \
                    read_poscar(args.CellFile)

                structure = (lattice_vectors, atomic_symbols, atom_positions)

                print("Structure read from \"{0}\"".format(args.CellFile))
                
                print_spacer = True
            else:
                # Finally, try and read from a mesh.yaml file if we
                # haven't already done so. If the file contains data for
                # a dense q-point mesh and/or has eigenvectors, this
                # could be very slow, so we only do it as a last resort.

                if not mesh_yaml_read and os.path.isfile(args.MeshYAML):
                    _, structure, atomic_masses = read_mesh_yaml(args.MeshYAML)

                    if structure is not None:
                        print(
                            "Structure read from \"{0}\""
                            .format(args.MeshYAML))
                        
                        print_spacer = True

                    if atomic_masses is not None:
                        print(
                            "Atomic masses read from \"{0}\""
                            .format(args.MeshYAML))
                        
                        print_spacer = True

        # Check we have a structure if we need one.

        if 'structure' in extract_list and structure is None:
            raise Exception(
                "Error: Unable to read a structure - please check "
                " input files.")

    # Finally, if we need atomic masses and we haven't read them from an
    # input file, get them from the phonopy.structure.atoms database.

    if 'atomic_masses' in extract_list and atomic_masses is None:
        _, atomic_symbols, _ = structure

        atomic_masses = []

        for symbol in atomic_symbols:
            symbol = symbol.title()

            mass = None

            for _, db_symbol, _, db_mass in atom_data:
                if symbol == db_symbol:
                    mass = db_mass
                    break

            if mass is None:
                raise Exception(
                    "Error: Atomic mass for symbol '{0}' not found in "
                    "Phonopy database.".format(symbol))

            atomic_masses.append(mass)

        print("Atomic masses taken from Phonopy database")
        
        print_spacer = True

    # If required, print a blank line after status messages.

    if print_spacer:
        print("")

    # Return a dictonary containing the requested data.

    output_data = {}

    if 'phonon_modes' in extract_list:
        output_data['phonon_modes'] = phonon_modes

    if 'structure' in extract_list:
        output_data['structure'] = structure

    if 'atomic_masses' in extract_list:
        output_data['atomic_masses'] = atomic_masses

    return output_data

def phonopy_load_optional(args):
    """ Reads and returns optional irreducible representation and
    linewidth data.

    Return value:
        A dictionary containing some or all of the following keys:
            'point_group' : point group of the structure.
            'irrep_data' : list of (irrep_symbol, band_indices) tuples
                assigning irreps to modes.
            'linewidths' : list of mode linewidths.
    """

    print_spacer = False

    # Read irrep data if available.

    point_group, irrep_data = None, None

    if os.path.isfile(args.IrrepsYAML):
        try:
            point_group, irrep_data = read_irreps_yaml(args.IrrepsYAML)

            print(
                "Irreducible representations read from \"{0}\""
                .format(args.IrrepsYAML))
        except:
            print(
                "WARNING: Failed to read irreducible representations "
                "from \"{0}\". Irrep-related functionality will not be "
                "available.".format(args.IrrepsYAML))
        
        print_spacer = True

    # Read linewidths if required.

    linewidths = None

    if args.LinewidthHDF5:
        if args.LinewidthTemperature is None:
            raise Exception(
                "Error: If reading linewidths from a Phono3py HDF5 "
                "file, a temperature must be supplied using the "
                "--linewidth-temperature argument.")

        linewidths = read_phono3py_hdf5(
            args.LinewidthHDF5, args.LinewidthTemperature)

        print(
            "Linewidths read from \"{0}\" (T = {1:.2f} K)"
            .format(args.LinewidthHDF5, args.LinewidthTemperature))
        
        print_spacer = True
            
    # Spacer after status messages if required.

    if print_spacer:
        print("")

    # Return data.

    output_data = {}

    if point_group is not None and irrep_data is not None:
        output_data['point_group'] = point_group
        output_data['irrep_data'] = irrep_data
    
    if linewidths is not None:
        output_data['linewidths'] = linewidths

    return output_data
