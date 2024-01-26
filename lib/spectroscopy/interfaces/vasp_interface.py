# spectroscopy/interfaces/vasp_interface.py


# ---------
# Docstring
# ---------

""" Routines for interfacing with the Vienna Ab-initio Simulation
Package (VASP) code. """


# -------
# Imports
# -------

import math
import re

import numpy as np

from spectroscopy.constants import ZERO_TOLERANCE
from spectroscopy import utilities


# ------------
# POSCAR Files
# ------------

""" Default system name for POSCAR files written by the write_poscar() 
function. """

_POSCAR_DEFAULT_SYSTEM_NAME = "Unknown Structure"


def read_poscar(file_path):
    """ Read a crystal structure from a VASP POSCAR file.

    Return value:
        A tuple of (system_name, lattice_vectors, atomic_symbols,
        atom_positions_frac) data.

    Notes:
        If the title line is blank, it will be set to the
        _POSCAR_DEFAULT_SYSTEM_NAME constant.
        This function only supports VASP 5.x-format POSCAR files (i.e.
        with atomic symbols as well as counts given in the header).
        Negative scale factors on line 2, interpreted by VASP as the
        volume of the cell, are not supported.
    """
    system_name = None
    lattice_vectors = None
    atomic_symbols, atom_positions = None, None

    with open(file_path, 'r') as input_reader:
        # Read title from line 1.

        system_name = next(input_reader).strip()

        if system_name == "":
            system_name = _POSCAR_DEFAULT_SYSTEM_NAME

        # Read scale factor.

        scale_factor = float(next(input_reader).strip())

        if scale_factor < 0.0:
            raise Exception(
                "Error: Unsupported negative scale factor in input "
                "file \"{0}\".".format(file_path))

        # Read lattice vectors.

        lattice_vectors = []

        for i in range(0, 3):
            lattice_vectors.append(
                [float(item)
                    for item in next(input_reader).strip().split()[:3]])

        lattice_vectors = [
            scale_factor * np.array(v, dtype=np.float64)
                for v in lattice_vectors
            ]

        # Read atom types and positions.

        atom_types = next(input_reader).strip().split()

        atom_counts = None

        try:
            atom_counts = [
                int(item) for item in next(input_reader).strip().split()
                ]
        except ValueError:
            raise Exception(
                "Error: Failed to read atom counts from input file "
                "\"{0}\".".format(file_path))

        if len(atom_types) != len(atom_counts):
            raise Exception(
                "Error: Inconsistent number of atom types and atom "
                "counts in input file \"{0}\".".format(file_path))

        atomic_symbols = []

        for atom_type, atom_count in zip(atom_types, atom_counts):
            atomic_symbols = atomic_symbols + [atom_type] * atom_count

        # Check for direct or Cartesian coordinates.

        key = next(input_reader).strip()[0].lower()

        if key == 's':
            # Selective dynamics keyword -> read first letter of the
            # following line.

            key = next(input_reader).strip()[0].lower()

        coordinate_type = None

        if key == 'c' or key == 'k':
            coordinate_type = 'cartesian'
        elif key == 'd':
            coordinate_type = 'fractional'
        else:
            raise Exception(
                "Error: Unknown coordinate type in input file \"{0}\"."
                .format(file_path))

        # Read atomic positions.

        atom_positions = []

        for i in range(0, sum(atom_counts)):
            atom_positions.append(
                [float(item)
                    for item in next(input_reader).strip().split()[:3]])

        atom_positions = [
            np.array(v, dtype=np.float64) for v in atom_positions
            ]

        if coordinate_type == 'cartesian':
            # Scale and convert to fractional coordinates.

            atom_positions = utilities.cartesian_to_fractional_coordinates(
                [scale_factor * v for v in atom_positions], lattice_vectors)

    return system_name, lattice_vectors, atomic_symbols, atom_positions


def write_poscar(
        structure, file_path, system_name=_POSCAR_DEFAULT_SYSTEM_NAME
        ):
    """ Write a structure to a VASP POSCAR file (VASP 5 format).

    Arguments:
        structure -- a tuple of (lattice_vectors, atomic_symbols,
            atom_positions_frac) data.
        file_path -- file path to save structure to.

    Keyword arguments:
        system_name -- system name to include on line 1 of the POSCAR 
        file (default: _POSCAR_DEFAULT_SYSTEM_NAME constant).

    Notes:
        The atom positions should be in fractional coordinates.
    """
    lattice_vectors, atomic_symbols, atom_positions = structure

    # Group atom positions by symbol. To avoid reordering atom types, we
    # keep track of the order in which new symbols were added.

    atom_groups_symbols, atom_groups = [], {}

    for symbol, position in zip(atomic_symbols, atom_positions):
        if symbol not in atom_groups_symbols:
            atom_groups[symbol] = []
            atom_groups_symbols.append(symbol)

        atom_groups[symbol].append(position)

    with open(file_path, 'w') as output_writer:
        output_writer.write("{0}\n".format(system_name))

        # Write scale factor.

        output_writer.write("  {0: >19.16f}\n".format(1.0))

        # Write lattice vectors.

        for i in range(0, 3):
            v_x, v_y, v_z = lattice_vectors[i]
            output_writer.write(
                "  {0: >21.16f}  {1: >21.16f}  {2: >21.16f}\n"
                .format(v_x, v_y, v_z))

        # Write atomic symbols and atom counts.

        symbols_line = "  ".join(
            ["{0: >3}".format(symbol) for symbol in atom_groups_symbols])

        output_writer.write("  {0}\n".format(symbols_line))

        atom_counts_line = "  ".join(
            ["{0: >3}".format(len(atom_groups[symbol]))
                for symbol in atom_groups_symbols])

        output_writer.write("  {0}\n".format(atom_counts_line))

        # Write positions.

        output_writer.write("Direct\n")

        for symbol in atom_groups_symbols:
            for f_1, f_2, f_3 in atom_groups[symbol]:
                output_writer.write(
                    "  {0: >21.16f}  {1: >21.16f}  {2: >21.16f}\n"
                    .format(f_1, f_2, f_3))


# ------------
# OUTCAR files
# ------------

""" Generic regex for capturing integer values. """

_OUTCAR_GENERIC_INTEGER_REGEX = re.compile(r"(?P<value>\d+)")

""" Regex for extrating atomic masses from the POMASS tags in re-printed
POTCAR headers at the top of the OUTCAR file. """

_OUTCAR_POMASS_REGEX = re.compile(r"POMASS =\s*(?P<pomass>\d+\.\d+)")

""" Regex for exctracting the number of atoms in the cell. """

_OUTCAR_NIONS_REGEX = re.compile(r"NIONS\s*=\s*(?P<nions>\d+)")

""" Regex for extracting phonon frequencies in THz, rad. THz, inv. Cm
and meV. """

_OUTCAR_MODE_FREQUENCY_REGEX = re.compile(
    r"^\s*(?P<mode_index>\d+) (?P<freq_sign>f|f/i)\s*=\s*"
    r"(?P<freq_thz>\d+\.\d+) THz\s*(?P<freq_radthz>\d+\.\d+) 2PiTHz\s*"
    r"(?P<freq_invcm>\d+\.\d+) cm-1\s*(?P<freq_mev>\d+\.\d+) meV\s*$"
    )


def parse_outcar(file_path, extract_list):
    """Parse a VASP OUTCAR file and extract the data specified by
    extract_list.

    extract_list may contain any of the following items:
        'atomic_masses' -- atomic masses (taken from the POMASS tags in
            the re-printed POTCAR headers at the top of the OUTCAR
            file).
        'born_charges' -- Born effective-charge tensors.
        'epsilon_static' -- macroscopic dielectric tensor (_not_ the
            ionic contribution).
        'phonon_modes' -- phonon frequencies (in THz, rad. THz, inv. cm
            and meV) and eigenvectors.

    If one or more requested data items cannot be extracted, an error is
    raised.

    Return value:
        A dictionary with the requested data stored under the keys
        specified in extract_list.

    Notes:
        If the atomic masses are overridden in the INCAR file, this will
        not be picked up by this function; this is because the fixed
        column width used to write out the POMASS INCAR tag makes it
        awkward to parse.
    """

    extract_list = [tag.lower() for tag in extract_list]

    for tag in extract_list:
        if tag not in ['atomic_masses', 'born_charges', 'epsilon_static',
                'phonon_modes']:
            raise Exception(
                "Error: Unrecognised extract tag '{0}'.".format(tag))

    output_data = {}

    with open(file_path, 'r') as input_reader:
        # Regardless of what we're returning, we parse the first part of
        # the OUTCAR file to obtain:
        # (1) the POMASS values written in the pseudopotential files;
        # (2) the number of atoms; and
        # (3) the number of atoms of each type.

        num_atoms = None
        type_masses, type_counts = None, None

        for line in input_reader:
            if "POMASS" in line:
                match = _OUTCAR_POMASS_REGEX.search(line)

                if match:
                    if type_masses is None:
                        type_masses = []

                    type_masses.append(float(match.group('pomass')))

            elif "NIONS" in line:
                match = _OUTCAR_NIONS_REGEX.search(line)
                num_atoms = int(match.group('nions'))

            elif "ions per type =" in line:
                matches = _OUTCAR_GENERIC_INTEGER_REGEX.findall(line)

                type_counts = [int(match) for match in matches]

            if (num_atoms is not None
                    and type_masses is not None
                    and type_counts is not None):
                # Sanity check.

                if (sum(type_counts) != num_atoms
                        or len(type_masses) != len(type_counts)):
                    raise Exception(
                        "Error: Failed to read atom information from "
                        "OUTCAR file \"{0}\".".format(file_path))

                break

        if 'atomic_masses' in extract_list:
            atomic_masses = []

            for mass, count in zip(type_masses, type_counts):
                atomic_masses = atomic_masses + [mass] * count

            output_data['atomic_masses'] = atomic_masses

        # Loop over the remaining lines and check for other things to extract.

        for line in input_reader:
            line = line.strip()

            if 'born_charges' in extract_list:
                if line == "BORN EFFECTIVE CHARGES (in e, cummulative output)":
                    next(input_reader)

                    bec_tensors = []

                    for _ in range(num_atoms):
                        next(input_reader)

                        bec_tensor = []

                        for i in range(3):
                            _, i_1, i_2, i_3 = \
                                next(input_reader).strip().split()

                            bec_tensor.append(
                                [float(i_1), float(i_2), float(i_3)])

                        bec_tensors.append(
                            np.array(bec_tensor, dtype=np.float64))

                    output_data['born_charges'] = bec_tensors

            if 'epsilon_static' in extract_list:
                # Make sure we don't capture the ionic part computed
                # with certain combinations of tags -- this is not what
                # we want for Raman calculations.

                if (line.startswith("MACROSCOPIC STATIC DIELECTRIC TENSOR")
                        and "IONIC CONTRIBUTION" not in line):
                    next(input_reader)

                    dielectric_tensor = []

                    for i in range(3):
                        i_1, i_2, i_3 = next(input_reader).strip().split()

                        dielectric_tensor.append(
                            [float(i_1), float(i_2), float(i_3)])

                    output_data['epsilon_static'] = np.array(
                        dielectric_tensor, dtype=np.float64)
                
                # Make sure we capture the "density-density" calculation of the
                # real dielectric function -- again, this is what we want for
                # the Raman calculations.
                
                elif ("REAL DIELECTRIC FUNCTION" in line
                      and "density-density" in line):
                    
                    for _ in range(2):
                        next(input_reader)
                    
                    eps_e, eps_r_e = [], []
                    
                    for line in input_reader:
                        vals = line.strip().split()
                        
                        if len(vals) == 7:
                            e, xx, yy, zz, xy, yz, xz = [
                                float(val) for val in vals]
                            
                            eps_r = np.array(
                                [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]],
                                dtype=np.float64)
                            
                            eps_e.append(e)
                            eps_r_e.append(eps_r)
                        else:
                            break
                    
                    if len(eps_e) == 0:
                        raise Exception("Error: Failed to parse real "
                                        "dielectric function.")
                    
                    if math.fabs(eps_e[0] > ZERO_TOLERANCE):
                        raise Exception("Error: Real dielectric function does "
                                        "does not start with E = 0.")
                    
                    output_data['epsilon_static'] = eps_r_e[0]
            
            if 'phonon_modes' in extract_list:
                if line == "Eigenvectors and eigenvalues of the dynamical " \
                        "matrix":
                    if 'phonon_modes' in output_data:
                        # If NWRITE = 3 was set in the INCAR file, a
                        # second set of frequencies and eigenvectors
                        # after division by sqrt(mass) are written to
                        # the OUTCAR file. Therefore, if we encounter
                        # multiple sets of phonon data, we skip it.

                        continue

                    frequencies = [[] for _ in range(4)]
                    eigenvectors = []

                    for i in range(3):
                        next(input_reader)

                    for i in range(3 * num_atoms):
                        # Read frequencies in various units.

                        match = _OUTCAR_MODE_FREQUENCY_REGEX.match(
                            next(input_reader))

                        if int(match.group('mode_index')) != i + 1:
                            raise Exception(
                                "Error: Unexpected mode ordering in a "
                                "frequencies/eigenvectors block in "
                                "OUTCAR file \"{0}\".".format(file_path))

                        freq_sign = -1.0 if match.group('freq_sign') == 'f/i' \
                            else 1.0

                        for e, group_name in enumerate(
                                ['freq_thz', 'freq_radthz', 'freq_invcm',
                                'freq_mev']):
                            
                            freq = math.copysign(
                                    float(match.group(group_name)), freq_sign)

                            frequencies[e].append(freq)

                        next(input_reader)

                        # Read eigenvector.

                        eigenvector = []

                        for j in range(0, num_atoms):
                            vals = next(input_reader).strip().split()[-3:]
                              
                            eigenvector.append([float(val) for val in vals])

                        eigenvectors.append(
                            np.array(eigenvector, dtype=np.float64))

                        next(input_reader)

                    output_data['phonon_modes'] = (frequencies, eigenvectors)

    # Check we were able to extract everything that was asked for.

    for tag in extract_list:
        if tag not in output_data:
            raise Exception(
                "Error: Failed to extract data for tag '{0}' from "
                "OUTCAR file \"{1}\".".format(tag, file_path))

    return output_data
