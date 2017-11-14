# SpectroscoPy/Utilities.py


# ---------
# Docstring
# ---------

""" Utility functions. """


# -------
# Imports
# -------

import numpy as np;


# ---------
# Functions
# ---------

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
