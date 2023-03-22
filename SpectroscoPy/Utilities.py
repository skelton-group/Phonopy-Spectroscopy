# SpectroscoPy/Utilities.py


# ---------
# Docstring
# ---------

""" Utility functions and constants. """


# -------
# Imports
# -------

import numpy as np

from SpectroscoPy import Constants


# ---------------
# Unit Conversion
# ---------------

def convertfrequencyunits(frequencies, unitsfrom, unitsto):
    """
    Convert frequencies in unitsfrom to unitsto.
    Supported values of unitsfrom/unitsto are: 'thz', 'inv_cm', 'mev' and 'um'.
    """

    # No point in doing any work if we don't have to...

    if unitsfrom == unitsto:
        return frequencies

    # First convert frequencies to THz...

    if unitsfrom != 'thz':
        if unitsfrom == 'inv_cm':
            frequencies = [frequency / Constants.THzToInvCm for frequency in
                           frequencies
                           ]
        elif unitsfrom == 'mev':
            frequencies = [frequency / Constants.THzToMeV for frequency in
                           frequencies
                           ]
        elif unitsfrom == 'um':
            frequencies = [1.0 / (frequency * Constants.THzToInvUm) for
                           frequency in frequencies
                           ]
        else:
            raise Exception("Error: Unsupported units '{0}'.".format(
                unitsfrom))

    # ... and now convert to the desired unit.

    if unitsto != 'thz':
        if unitsto == 'inv_cm':
            frequencies = [frequency * Constants.THzToInvCm for frequency in
                           frequencies
                           ]
        elif unitsto == 'mev':
            frequencies = [frequency * Constants.THzToMeV for frequency in
                           frequencies
                           ]
        elif unitsto == 'um':
            frequencies = [1.0 / (frequency * Constants.THzToInvUm) for
                           frequency in frequencies
                           ]
        else:
            raise Exception("Error: Unsupported units '{0}'.".format(unitsto))

    return frequencies


# ------------------
# Structure Handling
# ------------------

def calculatecellvolume(latticevectors):
    """ Calculate the cell volume from a set of lattice vectors by computing
    the scalar triple product. """

    dim1, dim2 = np.shape(latticevectors)

    if dim1 != 3 or dim2 != 3:
        raise Exception("Error: Lattice vectors must be a 3x3 matrix.")

    v1, v2, v3 = latticevectors

    return np.dot(v1, np.cross(v2, v3))


def cartesiantofractionalcoordinates(positions, latticevectors):
    """
    Convert positions from Cartesian to fractional coordinates.

    Arguments:
        positions -- a list of N three-component vectors in Cartesian
        coordinates.
        latticeVectors -- must be convertible to a 3x3 NumPy matrix.
    """

    dim1, dim2 = np.shape(latticevectors)

    if dim1 != 3 or dim2 != 3:
        raise Exception("Error: latticeVectors must be a 3x3 matrix.")

    # The transformation matrix from cartesian to fractional coordinates is
    # the inverse of a matrix built from the lattice vectors.

    transformationmatrix = np.linalg.inv(latticevectors)

    # If the positions are not three-component vectors, np.dot() raises a
    # readable error message.

    return [
        np.dot(position, transformationmatrix) % 1.0
        for position in positions
        ]


def fractionaltocartesiancoordinates(positions, latticevectors):
    """
    Convert positions from fractional to Cartesian coordinates.

    Arguments:
        positions -- a list of N three-component vectors in Cartesian
        coordinates.
        latticeVectors -- must be convertible to a 3x3 NumPy matrix.
    """

    dim1, dim2 = np.shape(latticevectors)

    if dim1 != 3 or dim2 != 3:
        raise Exception("Error: latticeVectors must be a 3x3 matrix.")

    v1, v2, v3 = latticevectors[0], latticevectors[1], latticevectors[2]

    return [
        np.multiply(f1, v1) + np.multiply(f2, v2) + np.multiply(f3, v3)
        for f1, f2, f3 in positions
        ]


# --------------------
# Eigenvector handling
# --------------------

def eigenvectors_to_eigendisplacements(eigenvectors, atomicmasses):
    """
    Return the eigenvectors after division of each component by sqrt(mass).

    Arguments:
        eigenvectors -- eigenvectors as 3N x N three-component vectors.
        atomicMasses -- set of N atomic masses.
    """

    nummodes, numatoms, eigdim3 = len(eigenvectors), len(eigenvectors[0]), \
        len(eigenvectors[0][0])

    if nummodes != 3 * numatoms or eigdim3 != 3:
        raise Exception("Error: eigenvectors should be a 3N x N x 3 matrix.")

    massesdim1 = len(atomicmasses)

    if massesdim1 != numatoms:
        raise Exception("Error: The number of supplied atomic masses is "
                        "inconsistent with the number of displacements in the "
                        "eigenvectors."
                        )

    sqrtmasses = np.sqrt(atomicmasses)

    eigendisplacements = np.zeros(
        (nummodes, numatoms, 3), dtype=np.float64
        )

    for i in range(0, nummodes):
        # This should avoid requiring external code to supply eigenvectors as
        # NumPy arrays.

        eigenvector = eigenvectors[i]

        for j in range(0, numatoms):
            eigendisplacements[i, j, :] = np.divide(eigenvector[j],
                                                    sqrtmasses[j])

    return eigendisplacements
