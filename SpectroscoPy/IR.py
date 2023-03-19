# SpectroscoPy/IR.py


# ---------
# Docstring
# ---------

""" Core routines for calculating infrared (IR) intensities. """


# -------
# Imports
# -------

import numpy as np


# ---------
# Functions
# ---------

# lowercase functions and arguments.
def calculateirintensity(eigendisplacement, bectensors):
    """
    Calculate the infrared intensity of a mode with eigendisplacement vector
    X_j,b from the Born effective-charge tensor Z_j,ab.

    Arguments:
        eigendisplacement -- mode eigendisplacement (N x three-component
        vectors).
        becTensors -- Born effective-charge tensors (N x 3x3 matrices).
    """

    # This should in principle allow external code to supply the parameters
    # as (correctly-dimensioned) lists.

    eigdim1, eigdim2 = np.shape(eigendisplacement)

    if eigdim2 != 3:
        raise Exception("Error: eigendisplacement must be a set of N "
                        "three-component vectors."
                        )

    becdim1, becdim2, becdim3 = np.shape(bectensors)

    if becdim2 != 3 or becdim3 != 3:
        raise Exception("Error: becTensors must be a set of N 3x3 matrices.")

    if eigdim1 != becdim1:
        raise Exception("Error: The first dimensions of eigendisplacement and "
                        "becTensors (= the number of atoms) must be the same."
                        )

    # This could be written more efficiently using e.g. np.einsum.
    # However, this routine is unlikely to be called inside a tight loop, and
    # having the formula laid out clearly makes it a lot more readable.

    irintensity = 0.0

    for a in range(0, 3):
        # Sum(a).

        sumtemp1 = 0.0

        # eigdim1 = becdim1 is the number of atoms.

        for j in range(0, eigdim1):
            # Sum(j).

            sumtemp2 = 0.0

            for b in range(0, 3):
                # Sum(b) Z_j,ab x X_j,b

                sumtemp2 += bectensors[j][a][b] * eigendisplacement[j][b]

            sumtemp1 += sumtemp2

        irintensity += sumtemp1 ** 2

    return irintensity
