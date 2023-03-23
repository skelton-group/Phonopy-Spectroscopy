# spectroscopy/ir.py


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

def calculate_ir_intensity(eigendisplacement, bec_tensors):
    """ Calculate the infrared intensity of a mode with
    eigendisplacement vector X_j,b from the Born effective-charge tensor
    Z_j,ab.

    Arguments:
        eigendisplacement -- mode eigendisplacement (N x three-component
            vectors).
        bec_tensors -- Born effective-charge tensors (N x 3x3 matrices).
    """

    # This should in principle allow external code to supply the
    # parameters as (correctly-dimensioned) lists.

    eig_dim_1, eig_dim_2 = np.shape(eigendisplacement)

    if eig_dim_2 != 3:
        raise Exception(
            "Error: eigendisplacement must be a set of N "
            "three-component vectors.")

    bec_dim_1, bec_dim_2, bec_dim_3 = np.shape(bec_tensors)

    if bec_dim_2 != 3 or bec_dim_3 != 3:
        raise Exception(
            "Error: becTensors must be a set of N 3x3 matrices.")

    if eig_dim_1 != bec_dim_1:
        raise Exception(
            "Error: The first dimensions of eigendisplacement and "
            "bec_tensors (= the number of atoms) must be the same.")

    # This could be written more efficiently using e.g. np.einsum.
    # However, this routine is unlikely to be called inside a tight
    # loop, and having the formula laid out clearly makes it a lot more
    # readable.

    ir_intensity = 0.0

    for a in range(0, 3):
        # Sum(a).

        sum_temp_1 = 0.0

        # eig_dim_1 = bec_dim_1 is the number of atoms.

        for j in range(0, eig_dim_1):
            # Sum(j).

            sum_temp_2 = 0.0

            for b in range(0, 3):
                # Sum(b) Z_j,ab x X_j,b

                sum_temp_2 += bec_tensors[j][a][b] * eigendisplacement[j][b]

            sum_temp_1 += sum_temp_2

        ir_intensity += sum_temp_1 ** 2

    return ir_intensity
