# spectroscopy/utilities.py


# ---------
# Docstring
# ---------

""" Utility functions and constants. """


# -------
# Imports
# -------

import numpy as np
from numpy.linalg import norm, inv

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


def flatten(l):
    """ Flatten a list containing nested lists as elements.
    """
    
    result = [] 
    
    # for each list in an nested list, append each element within a sublist.
    
    for neslist in l: 
        for i in neslist: 
            result.append(i) 
    
    return result


def hkl_to_real_space_normal(hkl, lattice_vectors):
    """ Given a vector of (h, k, l) integers specifying a plane, return
    the real-space vector normal to the plane.
    
    Arguments:
        hkl --a 3x1 NumPy array.
        lattixe_vectors -- must be convertible to a 3x3 NumPy matrix.
        """
    
    a1 = lattice_vectors[0, :]
    a2 = lattice_vectors[1, :]
    a3 = lattice_vectors[2, :]
    
    h, k, l = hkl[0], hkl[1], hkl[2]
    
    a = norm(a1)
    b = norm(a2)
    c = norm(a3)
    
    # generate the reciprocal lattice vector.
    
    a2a3 = np.cross(a2, a3)
    V = np.dot(a1, a2a3)
    a3a1 = np.cross(a3, a1)
    a1a2 = np.cross(a1, a2)

    b1 = a2a3 / V
    b2 = a3a1 / V
    b3 = a1a2 / V

    G = (h*b1) + (k*b2) + (l*b3)
    
    # convert reciprocal lattice vector into real space.
    
    recip_metric = np.matrix([[np.dot(b1, b1), np.dot(b1, b2), np.dot(b1, b3)],
                              [np.dot(b2, b1), np.dot(b2, b2), np.dot(b2, b3)],
                              [np.dot(b3, b1), np.dot(b3, b2), np.dot(b3, b3)]]
                             )
    
    real_metric = inv(recip_metric)
    
    component = np.dot(recip_metric, hkl)
    component = np.array(flatten(component.tolist()))
    
    temp = (component[0] * a1) + (component[1] * a2) + (component[2] * a3)
    normal = temp / norm(temp)
    
    const = 2 / (3 * np.square(a))
    primitive_norm = normal / const
    
    return normal


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
    
    cth = np.cos(theta)
    sth = np.sin(theta)
    
    z_rot = np.matrix([[cth, -sth, 0],
                       [sth, cth, 0],
                       [0, 0, 1]]
                      )
    
    return z_rot


def rotation_matrix_from_vectors(a, b):
    """ Given a pair of 3D vectors a and b, return the 3D rotation
    matrix that rotates a into b. """
    
    # obtain the axis of rotation normal to the ab plane and the angle of
    # rotation.
    
    v1 = a / norm(a)
    v2 = b / norm(b)
    x = np.cross(v1, v2)
    k = x / (norm(v1) * norm(v2))
    
    c = np.dot(v1, v2)
    s = norm(x)
    
    # obtain the skew-symmetric cross product matrix and find the rotation
    # matrix using Rodrigues' formula.
    
    K = np.matrix([[0, -k[2], k[1]],
                   [k[2], 0, -k[0]],
                   [-k[1], k[0], 0]]
                  )
    
    size = len(v1)
    I = np.identity(size)
    
    # check if the vectors are parallel or orthogonal to each other.
    
    min = 1.0e-10
    
    if np.abs(c - 1.0) < min:
        R = I
    elif np.abs(c + 1.0) < min:
        R = -I
    else:
        temp = (1 - c) / np.square(s)
        R = I + K + ((K @ K) * temp)
    
    return R


def rotate_tensors(rotation_matrix, tensors):
    """ Rotate a 3x3 tensor or list of tensors by the supplied
    3D rotation matrix.
    
    Notes:
        tensors may be a single tensor or list of tensors. If a single
        tensor is supplied, a single tensor will be returned; if a list
        is supplied, then a list will be returned.
    """

    tensors_ok = True

    n_dim = np.ndim(tensors)

    if n_dim < 2 or n_dim > 3:
        tensors_ok = False
    else:
        dim_1, dim_2 = np.shape(tensors)[-2:]

        if dim_1 != 3 or dim_2 != 3:
            tensors_ok = False
    
    if not tensors_ok:
        raise Exception(
            "Error: tensors must be a 3x3 matrix or set of 3x3 "
            "matrices.")
    
    if n_dim == 2:
        tensors = [tensors]
    
    rotated_tensors = []
    
    for tensor in tensors:
        temp = rotation_matrix * tensor * inv(rotation_matrix)
        rotated_tensors.append(temp)

    if len(tensors) == 1:
        return tensors[0]
    else:
        return tensors
