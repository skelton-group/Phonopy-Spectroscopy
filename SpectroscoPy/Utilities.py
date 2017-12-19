# SpectroscoPy/Utilities.py


# ---------
# Docstring
# ---------

""" Utility functions and constants. """


# -------
# Imports
# -------

import numpy as np;

from SpectroscoPy import Constants;


# ---------------
# Unit Conversion
# ---------------

def ConvertFrequencyUnits(frequencies, unitsFrom, unitsTo):
    """
    Convert frequencies in unitsFrom to unitsTo.
    Supported values of unitsFrom/unitsTo are: 'thz', 'inv_cm', 'mev' and 'um'.
    """

    unitsFrom, unitsTo = unitsFrom.lower(), unitsTo.lower();

    # No point in doing any work if we don't have to...

    if unitsFrom == unitsTo:
        return frequencies;

    # First convert frequencies to THz...

    if unitsFrom != 'thz':
        if unitsFrom == 'inv_cm':
            frequencies = [frequency / Constants.THzToInvCm for frequency in frequencies];
        elif unitsFrom == 'mev':
            frequencies = [frequency / Constants.THzToMeV for frequency in frequencies];
        elif unitsFrom == 'um':
            frequencies = [1.0 / (frequency * Constants.THzToInvUm) for frequency in frequencies];
        else:
            raise Exception("Error: Unsupported units '{0}'.".format(unitsFrom));

    # ... and now convert to the desired unit.

    if unitsTo != 'thz':
        if unitsTo == 'inv_cm':
            frequencies = [frequency * Constants.THzToInvCm for frequency in frequencies];
        elif unitsTo == 'mev':
            frequencies = [frequency * Constants.THzToMeV for frequency in frequencies];
        elif unitsTo == 'um':
            frequencies = [1.0 / (frequency * Constants.THzToInvUm) for frequency in frequencies];
        else:
            raise Exception("Error: Unsupported units '{0}'.".format(unitsTo));

    return frequencies;


# ------------------
# Structure Handling
# ------------------

def CalculateCellVolume(latticeVectors):
    """ Calculate the cell volume from a set of lattice vectors by computing the scalar triple product. """

    dim1, dim2 = np.shape(latticeVectors);

    if dim1 != 3 or dim2 != 3:
        raise Exception("Error: Lattice vectors must be a 3x3 matrix.");

    v1, v2, v3 = latticeVectors;

    return np.dot(v1, np.cross(v2, v3));

def CartesianToFractionalCoordinates(positions, latticeVectors):
    """
    Convert positions from Cartesian to fractional coordinates.

    Arguments:
        positions -- a list of N three-component vectors in Cartesian coordinates.
        latticeVectors -- must be convertible to a 3x3 NumPy matrix.
    """

    dim1, dim2 = np.shape(latticeVectors);

    if dim1 != 3 or dim2 != 3:
        raise Exception("Error: latticeVectors must be a 3x3 matrix.");

    # The transformation matrix from cartesian to fractional coordinates is the inverse of a matrix built from the lattice vectors.

    transformationMatrix = np.linalg.inv(latticeVectors);

    # If the positions are not three-component vectors, np.dot() raises a readable error message.

    return [
        np.dot(position, transformationMatrix) % 1.0
            for position in positions
        ];

def FractionalToCartesianCoordinates(positions, latticeVectors):
    """
    Convert positions from fractional to Cartesian coordinates.

    Arguments:
        positions -- a list of N three-component vectors in Cartesian coordinates.
        latticeVectors -- must be convertible to a 3x3 NumPy matrix.
    """

    dim1, dim2 = np.shape(latticeVectors);

    if dim1 != 3 or dim2 != 3:
        raise Exception("Error: latticeVectors must be a 3x3 matrix.");

    v1, v2, v3 = latticeVectors[0], latticeVectors[1], latticeVectors[2];

    return [
        np.multiply(f1, v1) + np.multiply(f2, v2) + np.multiply(f3, v3)
            for f1, f2, f3 in positions
        ];


# --------------------
# Eigenvector handling
# --------------------

def EigenvectorsToEigendisplacements(eigenvectors, atomicMasses):
    """
    Return the eigenvectors after division of each component by sqrt(mass).

    Arguments:
        eigenvectors -- eigenvectors as 3N x N three-component vectors.
        atomicMasses -- set of N atomic masses.
    """

    numModes, numAtoms, eigDim3 = len(eigenvectors), len(eigenvectors[0]), len(eigenvectors[0][0]);

    if numModes != 3 * numAtoms or eigDim3 != 3:
        raise Exception("Error: eigenvectors should be a 3N x N x 3 matrix.");

    massesDim1 = len(atomicMasses);

    if massesDim1 != numAtoms:
        raise Exception("Error: The number of supplied atomic masses is inconsistent with the number of displacements in the eigenvectors.");

    sqrtMasses = np.sqrt(atomicMasses);

    eigendisplacements = np.zeros(
        (numModes, numAtoms, 3), dtype = np.float64
        );

    for i in range(0, numModes):
        # This should avoid requiring external code to supply eigenvectors as NumPy arrays.

        eigenvector = eigenvectors[i];

        for j in range(0, numAtoms):
            eigendisplacements[i, j, :] = np.divide(eigenvector[j], sqrtMasses[j]);

    return eigendisplacements;
