# spectroscopy/interfaces/phonopy_interface.py


# ---------
# Docstring
# ---------

""" Routines for interfacing with the Phono(3)py packages. """


# -------
# Imports
# -------

import math

import h5py

import yaml

# This code tries to use the C version of the YAML Loader and falls back
# to the standard Loader if this is not available. The code was taken
# from https://pyyaml.org/wiki/PyYAMLDocumentation.

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

import numpy as np

# Needed for parsing Phonopy BORN files.

from phonopy.file_IO import parse_BORN
from phonopy.structure.atoms import Atoms

from spectroscopy.constants import ZERO_TOLERANCE


# ----------
# YAML Files
# ----------

def read_phonopy_yaml(file_path=r"phonopy.yaml"):
    """ Read a Phonopy phonopy.yaml file and extract the structure and
    atomic masses.

    Return value:
        A (structure, atomic_masses) tuple containing the extracted
        data, where structure is a tuple of (lattice_vectors,
        atomic_symbols, atom_positions_frac) data.
    """

    # Load and parse YAML file.

    input_yaml = None

    with open(file_path, 'r') as input_reader:
        input_yaml = yaml.load(input_reader, Loader=Loader)

    # Get primitive cell used in the calculation.

    primitive_cell = input_yaml['primitive_cell']

    # Get structure.

    lattice_vectors = [
        np.array(v, dtype=np.float64)
            for v in primitive_cell['lattice']
        ]

    atomic_symbols = [
        atom['symbol'] for atom in primitive_cell['points']
        ]

    atom_positions = [
        np.array(atom['coordinates'], dtype=np.float64)
            for atom in primitive_cell['points']
        ]

    structure = (lattice_vectors, atomic_symbols, atom_positions)

    # Get atomic masses.

    atomic_masses = [
        atom['mass'] for atom in primitive_cell['points']
        ]

    # Return data.

    return (structure, atomic_masses)


def read_mesh_yaml(file_path=r"mesh.yaml"):
    """ Read a Phonopy mesh.yaml file and extract the Gamma-point
    eigenvectors and a structure and list of atomic masses if available.

    Return value:
        A tuple of (phonon_modes, structure, atomic_masses) containing
        the Gamma-point frequencies/eigenvectors, and the structure and
        a list of atomic masses if available.
        The phonon modes are returned as a tuple of (frequencies,
        eigenvectors) lists with 3N elements; the eigenvectors are Nx3
        NumPy arrays.
        The structure is returned as a tuple of (lattice_vectors,
        atomic_symbols, atom_positions_frac).

    Notes:
        If the input file does not contain frequencies and eigenvectors
        at q_Gamma = (0, 0, 0), errors are raised.
        YAML files exported by older versions of Phonopy do not contain
        the structure and atomic masses - in this case, structure and
        atomic_masses are set to None.
    """
    # Load and parse YAML file.

    input_yaml = None

    with open(file_path, 'r') as input_reader:
        input_yaml = yaml.load(input_reader, Loader=Loader)

    # Get structure and atomic masses if available.

    structure, atomic_masses = None, None

    lattice_vectors, atom_positions, atomic_symbols = None, None, None

    if 'lattice' in input_yaml:
        lattice_vectors = [
            np.array(v, dtype=np.float64)
                for v in input_yaml['lattice']
            ]

    # Get atom positions, atomic symbols and atomic masses.

    if 'atoms' or 'points' in input_yaml:
        # Workaround for differences between older and newer versions of
        # Phonopy...

        atoms_key = 'atoms' if 'atoms' in input_yaml else 'points'

        atomic_symbols = [
            atom['symbol'] for atom in input_yaml[atoms_key]
            ]

        # ... and again!

        atom = input_yaml[atoms_key][0]

        pos_key = 'position' if 'position' in atom else 'coordinates'

        atom_positions = [
            np.array(atom[pos_key], dtype=np.float64)
                for atom in input_yaml[atoms_key]
            ]

        atomic_masses = [
            atom['mass'] for atom in input_yaml[atoms_key]
            ]

    # If structural information was read, create a (lattice_vectors,
    # atomic_symbols, atom_positions) tuple.

    if (lattice_vectors is not None
            and atom_positions is not None
            and atomic_symbols is not None):
        structure = (lattice_vectors, atomic_symbols, atom_positions)

    # Get Gamma-point frequencies and eigenvectors.

    phonon_modes = None

    frequencies, eigenvectors = None, None

    for q_point in input_yaml['phonon']:
        q_x, q_y, q_z = q_point['q-position']

        if (math.fabs(q_x) < ZERO_TOLERANCE
                and math.fabs(q_y) < ZERO_TOLERANCE
                and math.fabs(q_z) < ZERO_TOLERANCE):
            frequencies, eigenvectors = [], []

            for mode in q_point['band']:
                frequencies.append(mode['frequency'])

                if 'eigenvector' not in mode:
                    raise Exception(
                        "Error: Eigenvectors not found in YAML file "
                        "\"{0}\".".format(file_path))

                eigenvector = [
                    [atom_disp[i][0] for i in range(0, 3)]
                        for atom_disp in mode['eigenvector']
                    ]

                eigenvectors.append(
                    np.array(eigenvector, dtype=np.float64)
                    )

            break

    # Check Gamma-point frequencies and eigenvectors were found.

    if frequencies is None:
        raise Exception(
            "Error: Data for q_Gamma = (0, 0, 0) not found in YAML "
            "file \"{0}\".".format(file_path))

    phonon_modes = (frequencies, eigenvectors)

    # Return Gamma-point frequencies and eigenvectors, structure and
    # atomic masses.

    return phonon_modes, structure, atomic_masses

def read_irreps_yaml(file_path=r"irreps.yaml"):
    """
    Read a Phonopy irreps.yaml file and extract the irreps and
    associated band indices. The irreps must have been calculated at the
    Gamma point.

    Return value:
        A list of (irrep_symbol, band_indices) tuples mapping irreps
        to band indices.
    """

    input_yaml = None

    with open(file_path, 'r') as input_reader:
        input_yaml = yaml.load(input_reader, Loader=Loader)

    # Check the irreps are given for the Gamma point modes.

    q_x, q_y, q_z = input_yaml['q-position']

    if (math.fabs(q_x) >= ZERO_TOLERANCE
            or math.fabs(q_y) >= ZERO_TOLERANCE
            or math.fabs(q_z) >= ZERO_TOLERANCE):
        raise Exception(
            "Error: Irreps. are required for the Gamma point modes.")

    # Integer point group symbols are parsed as integers.

    point_group = str(input_yaml['point_group'])

    irrep_labels, band_indices = [], []

    for mode in input_yaml['normal_modes']:
        # Unassigned irreps are assigned as (text) "None".
        
        symbol = mode['ir_label']

        irrep_labels.append(symbol if symbol.lower() != "none" else None)

        # Band indices in the irreps.yaml files start from one -> adjust
        # to zero-based.

        band_indices.append([index - 1 for index in mode['band_indices']])

    return (point_group, [item for item in zip(irrep_labels, band_indices)])


# ----------
# HDF5 Files
# ----------

def read_mesh_hdf5(file_path):
    """ Read a Phonopy mesh.hdf5 file and return the Gamma-point
    eigenvectors.

    Return value:
        A tuple of (frequencies, eigenvectors) lists containing the 3N
        Gamma-point frequencies and eigenvectors. The eigenvectors are
        returned as Nx3 NumPy arrays.

    Notes:
        Errors are raised if the file does not contain an entry for
        q_Gamma = (0, 0, 0) or does not contain eigenvectors.
    """

    frequencies, eigenvectors = None, None

    mesh_hdf5 = h5py.File(file_path, 'r')

    try:
        # Locate the Gamma point in the q-point sampling mesh.

        q_index = None

        for i, (q_x, q_y, q_z) in enumerate(mesh_hdf5['qpoint']):
            if (math.fabs(q_x) < ZERO_TOLERANCE
                    and math.fabs(q_y) < ZERO_TOLERANCE
                    and math.fabs(q_z) < ZERO_TOLERANCE):
                q_index = i
                break

        if q_index is None:
            raise Exception(
                "Error: q_Gamma = (0, 0, 0) not found in the q-point "
                "mesh specified in \"{0}\".".format(file_path))

        # Extract a list of Gamma-point frequencies.

        frequencies = [
            frequency for frequency in mesh_hdf5['frequency'][q_index]
            ]

        if 'eigenvector' not in mesh_hdf5:
            raise Exception(
                "Error: Eigenvectors not found in \"{0}\"."
                .format(file_path))

        # Sanity check.

        assert len(frequencies) % 3 == 0

        n_atoms = len(frequencies) // 3

        # Extract the Gamma-point eigenvectors.

        eigenvectors = None
        
        temp = mesh_hdf5['eigenvector'][q_index].real

        eigenvectors = [
            temp[:, i].reshape(n_atoms, 3) for i in range(len(frequencies))]

    finally:
        # Make sure the input file gets closed.

        mesh_hdf5.close()

    # Return the Gamma-point frequencies and eigenvectors.

    return frequencies, eigenvectors


def read_phono3py_hdf5(file_path, linewidth_temperature):
    """ Read a Phono3py HDF5 file and return the Gamma-point mode
    linewidths at the specified temperature. """

    linewidths = None

    phono3py_hdf5 = h5py.File(file_path, 'r')

    try:
        # Locate the temperature index corresponding to
        # linewidth_temperature.

        t_index = None

        for i, t in enumerate(phono3py_hdf5['temperature']):
            if math.fabs(t - linewidth_temperature) < ZERO_TOLERANCE:
                t_index = i
                break

        if t_index is None:
            raise Exception(
                "Error: Linewidth temperature T = {0} not found under "
                "the 'temperature' key.".format(linewidth_temperature))

        q_index = None

        if 'qpoint' in phono3py_hdf5:
            # If the HDF5 file has a 'qpoint' key, search the list for
            # q_Gamma = (0, 0, 0). (If it doesn't, assume the HDF5 file
            # is an "intermediate" file containing the linewidths at
            # q_Gamma.

            for i, (q_x, q_y, q_z) in enumerate(
                    phono3py_hdf5['qpoint']):
                if (math.fabs(q_x) < ZERO_TOLERANCE
                        and math.fabs(q_y) < ZERO_TOLERANCE
                        and math.fabs(q_z) < ZERO_TOLERANCE):
                    q_index = i
                    break

                if q_index is None:
                    raise Exception(
                        "Error: q = (0, 0, 0) not found under the "
                        "'qpoint' key.")
        
        # Extract the linewidths.

        if q_index is not None:
            # If the 'qpoint' key is present in the file, the linewidths
            # dataset has shape (num_temperatures, num_q_points, num_bands).

            linewidths = [
                linewidth
                    for linewidth in phono3py_hdf5['gamma'][t_index][q_index]
                ]
        else:
            # If not, the dataset has shape (num_temperatures, num_bands).

            linewidths = [
                linewidth for linewidth in phono3py_hdf5['gamma'][t_index]
                ]
    
    finally:
        # Make sure we close the input file.

        phono3py_hdf5.close()

    return linewidths


# ----------
# Misc Files
# ----------

def read_born(structure, file_path="BORN"):
    """ Read a Phonopy BORN file, expand the charges for the
    supplied structure, and return a list of Born effective-charge
    tensors for each atom in the structure.

    Arguments:
        structure -- a tuple of (lattice_vectors, atomic_symbols,
            atom_positions) specifying the structure to be used to
            expand the Born charges to the full set of atoms.

    Keyword arguments:
        file_path -- path to the Phonopy BORN-format file to read
            (default: BORN).

    Return value:
        A list of 3x3 Born effective-charge tensors for each atom in the
        supplied structure.

    Notes:
        The atom positions supplied with structure are assumed to be in
        fractional coordinates.
    """

    # Convert the supplied structure to a Phonopy Atoms object.

    lattice_vectors, atomic_symbols, atom_positions = structure

    cell = Atoms(
        symbols=atomic_symbols, cell=lattice_vectors,
        scaled_positions=atom_positions)

    # Read the BORN file.

    born_data = parse_BORN(cell, filename=file_path)

    # Return the Born effective-charge tensors from the dictionary
    # returned by parse_BORN.

    return [bec_tensor for bec_tensor in born_data['born']]
