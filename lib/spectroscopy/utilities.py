# spectroscopy/utilities.py


# ---------
# Docstring
# ---------

""" Utility functions and constants. """


# -------
# Imports
# -------

import numpy as np

from spectroscopy import constants


# ---------------
# Unit Conversion
# ---------------

def convert_frequency_units(frequencies, units_from, units_to):
    """ Convert frequencies in units_from to units_to.
    Supported values of units_from/units_to are: 'thz', 'inv_cm', 'mev'
    and 'um'.
    """

    # No point in doing any work if we don't have to...

    if units_from == units_to:
        return frequencies

    # First convert frequencies to THz...

    if units_from != 'thz':
        if units_from == 'inv_cm':
            frequencies = [
                frequency / constants.THZ_TO_INV_CM
                    for frequency in frequencies
                ]
        elif units_from == 'mev':
            frequencies = [
                frequency / constants.THZ_TO_MEV 
                    for frequency in frequencies
                ]
        elif units_from == 'um':
            frequencies = [
                1.0 / (frequency * constants.THZ_TO_INV_UM)
                    for frequency in frequencies
                ]
        else:
            raise Exception(
                "Error: Unsupported units '{0}'.".format(units_from)
                )

    # ... and now convert to the desired unit.

    if units_to != 'thz':
        if units_to == 'inv_cm':
            frequencies = [
                frequency * constants.THZ_TO_INV_CM
                    for frequency in frequencies
                ]
        elif units_to == 'mev':
            frequencies = [
                frequency * constants.THZ_TO_MEV
                    for frequency in frequencies
                ]
        elif units_to == 'um':
            frequencies = [
                1.0 / (frequency * constants.THZ_TO_INV_UM)
                    for frequency in frequencies
                ]
        else:
            raise Exception(
                "Error: Unsupported units '{0}'.".format(units_to)
                )

    return frequencies


# ------------------
# Structure Handling
# ------------------

def calculate_cell_volume(lattice_vectors):
    """ Calculate the cell volume from a set of lattice vectors by
    computing the scalar triple product. """

    dim_1, dim_2 = np.shape(lattice_vectors)

    if dim_1 != 3 or dim_2 != 3:
        raise Exception("Error: Lattice vectors must be a 3x3 matrix.")

    v_1, v_2, v_3 = lattice_vectors

    return np.dot(v_1, np.cross(v_2, v_3))


def cartesian_to_fractional_coordinates(positions, lattice_vectors):
    """ Convert positions from Cartesian to fractional coordinates.

    Arguments:
        positions -- a list of N three-component vectors in Cartesian
            coordinates.
        lattice_vectors -- must be convertible to a 3x3 NumPy matrix.
    """

    dim_1, dim_2 = np.shape(lattice_vectors)

    if dim_1 != 3 or dim_2 != 3:
        raise Exception("Error: lattice_vectors must be a 3x3 matrix.")

    # The transformation matrix from cartesian to fractional coordinates
    # is the inverse of a matrix built from the lattice vectors.

    trans_mat = np.linalg.inv(lattice_vectors)

    # If the positions are not three-component vectors, np.dot() raises
    # a readable error message.

    return [
        np.dot(pos, trans_mat) % 1.0
            for pos in positions
        ]

def fractional_to_cartesian_coordinates(positions, lattice_vectors):
    """ Convert positions from fractional to Cartesian coordinates.

    Arguments:
        positions -- a list of N three-component vectors in Cartesian
            coordinates.
        lattice_vectors -- must be convertible to a 3x3 NumPy matrix.
    """

    dim_1, dim_2 = np.shape(lattice_vectors)

    if dim_1 != 3 or dim_2 != 3:
        raise Exception("Error: lattice_vectors must be a 3x3 matrix.")

    v_1, v_2, v_3 = lattice_vectors[0], lattice_vectors[1], lattice_vectors[2]

    return [
        np.multiply(f_1, v_1) + np.multiply(f_2, v_2) + np.multiply(f_3, v_3)
            for f_1, f_2, f_3 in positions
        ]

def hkl_to_real_space_normal(hkl, lattice_vectors):
    """ Given a vector of (h, k, l) integers specifying a plane, return
    the real-space vector normal to the plane. """

    raise NotImplementedError()


# --------------------
# Eigenvector handling
# --------------------

def eigenvectors_to_eigendisplacements(eigenvectors, atomic_masses):
    """ Return the eigenvectors after division of each component by
    sqrt(mass).

    Arguments:
        eigenvectors -- eigenvectors as 3N x N three-component vectors.
        atomic_masses -- set of N atomic masses.
    """

    num_modes, num_atoms = None, None

    error = False

    if np.ndim(eigenvectors) == 3:
        num_modes, num_atoms, eig_dim_3 = np.shape(eigenvectors)

        if num_modes != 3 * num_atoms or eig_dim_3 != 3:
            raise Exception(
                "Error: eigenvectors should be a 3N x N x 3 matrix.")

    if len(atomic_masses) != num_atoms:
        raise Exception(
            "Error: The number of supplied atomic masses is "
            "inconsistent with the number of displacements in the "
            "eigenvectors.")

    sqrt_masses = np.sqrt(atomic_masses)

    return np.divide(eigenvectors, sqrt_masses[np.newaxis, :, np.newaxis])


# ----
# Misc
# ----

def rotation_matrix_xy(theta):
    """ Given an angle theta, return the rotation matrix for a 2D
    rotation about theta in the xy plane. """

    raise NotImplementedError()

def rotation_matrix_from_vectors(v_1, v_2):
    """ Given a pair of 3D vectors v_1 and v_2, return the 3D rotation
    matrix that rotates v_1 to v_2. """

    raise NotImplementedError()
