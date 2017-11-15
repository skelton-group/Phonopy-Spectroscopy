# SpectroscoPy/Interfaces/Phonopy.py


# ---------
# Docstring
# ---------

""" Routines for interfacing with the Phono(3)py packages. """


# -------
# Imports
# -------

import yaml;

import numpy as np;

# Needed for parsing Phonopy BORN files.

from phonopy.file_IO import parse_BORN;
from phonopy.structure.atoms import Atoms


# ---------
# Functions
# ---------

def ReadMeshYAML(filePath = r"mesh.yaml"):
    """
    Read a Phonopy mesh.yaml file and extract the structure, atomic masses and Gamma-point eigenvectors.

    If the specified file is missing one or more of these items, an error is raised.

    Return value:
        A dictionary of containing the extracted data under the following keys:
            'structure' : a tuple of (lattice_vectors, atomic_symbols, atom_positions) data.
            'atomic_masses' : a list of atomic masses.
            'gamma_modes' : lists of 3N Gamma-point frequencies and eigenvectors as Nx3 NumPy arrays.

    Notes:
        The atom positions under the 'structure' key in the output dictionary are in fractional coordinates.
    """

    # Load and parse YAML file.

    inputYAML = None;

    with open(filePath, 'r') as inputReader:
        inputYAML = yaml.load(inputReader);

    outputData = { };

    # Get structure.

    latticeVectors = [
        np.array(vector, dtype = np.float64)
            for vector in inputYAML['lattice']
        ];

    atomicSymbols = [
        atom['symbol'] for atom in inputYAML['atoms']
        ];

    atomPositions = [
        np.array(atom['position'], dtype = np.float64)
            for atom in inputYAML['atoms']
        ];

    outputData['structure'] = (latticeVectors, atomicSymbols, atomPositions);

    # Get atomic masses.

    atomicMasses = [
        atom['mass'] for atom in inputYAML['atoms']
        ];

    outputData['atomic_masses'] = atomicMasses;

    # Try to get Gamma-point eigenvectors.

    for qPoint in inputYAML['phonon']:
        qx, qy, qz = qPoint['q-position'];

        if qx == 0.0 and qy == 0.0 and qz == 0.0:
            frequencies, eigenvectors = [], [];

            for mode in qPoint['band']:
                frequencies.append(mode['frequency']);

                eigenvector = [
                     [atomDisp[i][0] for i in range(0, 3)]
                        for atomDisp in mode['eigenvector']
                    ];

                eigenvectors.append(
                    np.array(eigenvector, dtype = np.float64)
                    );

            outputData['gamma_modes'] = (frequencies, eigenvectors);

            break;

    if 'gamma_modes' not in outputData:
        raise Exception("Error: Gamma-point phonon eigenvectors could not be read from YAML file \"{0}\".".format(filePath));

    return outputData;

def ReadIrRepsYAML(filePath = r"irreps.yaml"):
    """
    Read a Phonopy irreps.yaml file and extract the ir. reps. and associated band indices.

    The ir. reps. must have been calculated at the Gamma point.

    Return value:
        A list of (irrep_symbol, band_indices) tuples mapping ir. reps. to band indices.
    """

    inputYAML = None;

    with open(filePath, 'r') as inputReader:
        inputYAML = yaml.load(inputReader);

    # Check the ir. reps. are given for the Gamma point modes.

    qx, qy, qz = inputYAML['q-position'];

    if qx != 0.0 or qy != 0.0 or qz != 0.0:
        raise Exception("Error: Ir. reps. are required for the Gamma point modes.");

    irRepLabels, bandIndices = [], [];

    for mode in inputYAML['normal_modes']:
        irRepLabels.append(mode['ir_label']);
        bandIndices.append(mode['band_indices']);

    return [item for item in zip(irRepLabels, bandIndices)];

def ReadBORN(structure, filePath = "BORN"):
    """
    Read a Phonopy BORN file, and expand the charges for the supplied structure, and return a list of Born effective-charge tensors for each atom in the structure.

    Arguments:
        structure -- a tuple of (lattice_vectors, atomic_symbols, atom_positions) specifying the structure to be used to expand the Born charges to the full set of atoms.

    Keyword arguments:
        filePath -- path to the Phonopy BORN-format file to read (default: BORN).

    Return value:
        A list of 3x3 Born effective-charge tensors for each atom in the supplied structure.

    Notes:
        The atom positions supplied with structure are assumed to be in fractional coordinates.
    """

    # Convert the supplied structure to a Phonopy Atoms object.

    latticeVectors, atomicSymbols, atomPositions = structure;

    cell = Atoms(
        symbols = atomicSymbols, cell = latticeVectors, scaled_positions = atomPositions
        );

    # Read the BORN file.

    bornData = parse_BORN(cell, filename = filePath);

    # Return the Born effective-charge tensors from the dictionary returned by parse_BORN.

    return [becTensor for becTensor in bornData['born']];
