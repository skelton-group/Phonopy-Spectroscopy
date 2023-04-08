# spectroscopy/utilities.py


# ---------
# Docstring
# ---------

""" Utility functions and constants. """


# -------
# Imports
# -------

import warnings

import numpy as np

from spectroscopy.constants import (
    ZERO_TOLERANCE, THZ_TO_INV_CM, THZ_TO_MEV, THZ_TO_INV_UM)


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

    return [np.dot(pos, trans_mat) % 1.0 for pos in positions]

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
    the real-space vector normal to the plane.
    
    Arguments:
        hkl -- three integers specifying the Miller indices of a
            surface.
        lattixe_vectors -- must be convertible to a 3x3 NumPy matrix.
        """
    
    h, k, l = hkl
    a_1, a_2, a_3 = lattice_vectors
    
    a, b, c = np.linalg.norm(lattice_vectors, axis = 1)
    
    # Generate reciprocal lattice vectors.

    v = calculate_cell_volume(lattice_vectors)

    b_1 = np.cross(a_2, a_3) / v
    b_2 = np.cross(a_3, a_1) / v
    b_3 = np.cross(a_1, a_2) / v

    g = h * b_1 + k * b_2 + l * b_3
    
    # Metric tensors.
    
    recip_metric = np.array(
        [[np.dot(b_1, b_1), np.dot(b_1, b_2), np.dot(b_1, b_3)],
         [np.dot(b_2, b_1), np.dot(b_2, b_2), np.dot(b_2, b_3)],
         [np.dot(b_3, b_1), np.dot(b_3, b_2), np.dot(b_3, b_3)]])
    
    f_1, f_2, f_3 = np.dot(recip_metric, hkl)

    normal = (np.multiply(f_1, a_1)
              + np.multiply(f_2, a_2)
              + np.multiply(f_3, a_3))
    
    return normal / np.linalg.norm(normal)


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
    """ Given an angle theta in degrees, return the rotation matrix for
    a 2D rotation about theta in the xy plane. """
    
    theta_rad = np.pi * theta / 180.0

    cos_theta = np.cos(theta_rad)
    sin_theta = np.sin(theta_rad)
    
    return np.array(
        [[cos_theta, -sin_theta, 0.0], [sin_theta, cos_theta, 0.0],
        [0.0, 0.0, 1.0]], dtype = np.float64)

def rotation_matrix_from_vectors(a, b):
    """ Given a pair of 3D vectors a and b, return the 3D rotation
    matrix that rotates a into b. """
    
    # Obtain the axis of rotation normal to the ab plane and the angle
    # of rotation.
    
    a = np.array(a, dtype = np.float64) / np.linalg.norm(a)
    b = np.array(b, dtype = np.float64) / np.linalg.norm(b)
    
    c = np.cross(a, b)
    
    cos_theta = np.dot(a, b)
    sin_theta = np.linalg.norm(c)
    
    # Generate the skew-symmetric cross product matrix and find the
    # rotation matrix using Rodrigues' formula.
    
    k = np.array(
        [[0, -c[2], c[1]], [c[2], 0, -c[0]], [-c[1], c[0], 0]],
        dtype = np.float64)
    
    i = np.identity(3)
    
    # The rotation matrix depends on whether the vectors are parallel
    # or orthogonal.
    
    if np.abs(cos_theta - 1.0) < ZERO_TOLERANCE:
        return i
    elif np.abs(cos_theta + 1.0) < ZERO_TOLERANCE:
        return -i
    else:
        temp = (1 - cos_theta) / (sin_theta ** 2)
        return i + k + temp * np.matmul(k, k)

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
    
    inv_rotation_matrix = np.linalg.inv(rotation_matrix)
    
    tensors = [
        np.matmul(rotation_matrix, np.matmul(tensor, rotation_matrix))
            for tensor in tensors
        ]

    if n_dim == 2:
        return tensors[0]
    else:
        return tensors


# ---------------------
# NumPy Helper Routines
# ---------------------

def numpy_check_shape(a, expected_shape):
    """ Verifies that the array-like object a has the shape specified by
    the tuple of integers expected_shape.
    """

    if np.ndim(a) != len(expected_shape):
        return False
    
    shape = np.shape(a)

    for s, e_s in zip(shape, expected_shape):
        if s != e_s:
            return False
    
    return True

def numpy_create_or_copy(a):
    """ If a is a NumPy, return a copy; otherwise, create and return
    a new NumPy array. """

    b = np.asarray(a)
    
    if b is a:
        return b.copy()
    else:
        return b

def numpy_as_readonly(a):
    """ Returns a view of the NumPy array a with the writeable flag set
    to False. """

    b = np.asarray(a)

    if not b is a:
        warnings.warn(
            "a was converted to a NumPy array before returning a "
            "readonly view", UserWarning)
    
    v = b.view()
    v.flags.writeable = False

    return v
