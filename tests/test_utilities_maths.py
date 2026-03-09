# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Tests for the the maths routines in phonopy_spectroscopy.utilities."""


# -------
# Imports
# -------


import itertools
import unittest

import numpy as np

from phonopy_spectroscopy.utility.geometry import (
    parse_direction,
    rotation_matrix_from_axis_angle,
    rotation_matrix_from_vectors,
    rotate_tensors,
)

from phonopy_spectroscopy.utility.differentiation import (
    central_difference_available_precs,
    central_difference_coefficients,
)

from phonopy_spectroscopy.utility.quadrature import (
    unit_circle_quad_rule,
    unit_sphere_lebedev_quad_available_precs,
    unit_sphere_lebedev_quad_rule,
)


# ----------------
# Helper functions
# ----------------


def _rotation_matrix(axis, theta):
    """Return the rotation matrix for rotation about axis by the angle
    theta using the analytical formulae."""

    theta = np.radians(theta)

    c = np.cos(theta)
    s = np.sin(theta)

    if axis == "x":
        return [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]]

    if axis == "y":
        return [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]]

    if axis == "z":
        return [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]

    raise Exception("Unknown axis '{0}'.".format(axis))


def _circle(x, y):
    """Equation for a circle with unit radius."""
    return x**2 + y**2


def _sphere(x, y, z):
    """Equation for a sphere with unit radius."""
    return x**2 + y**2 + z**2


# ------------------------------------
# Tests for geometry-handling routines
# ------------------------------------


class TestGeometry(unittest.TestCase):
    """Class implementing unit tests for geometry-handling routines."""

    def test_rotation_matrix_from_vectors_1(self):
        """Test `rotation_matrix_from_vectors` with all combinations
        of +/-{x, y, z}."""

        cart_dirn = ["+x", "-x", "+y", "-y", "+z", "-z"]

        for d_1, d_2 in itertools.product(cart_dirn, cart_dirn):
            v_1, v_2 = parse_direction(d_1), parse_direction(d_2)

            r = rotation_matrix_from_vectors(v_1, v_2)

            self.assertTrue(np.allclose(np.dot(r, v_1), v_2))

    def test_rotation_matrix_from_vectors_2(self):
        """Test `rotation_matrix_from_vectors` with all combinations
        of 25 random vectors."""

        num_vecs = 25

        vecs = (2.0 * np.random.random(num_vecs * 3) - 0.5).reshape(
            num_vecs, 3
        )
        vecs /= np.linalg.norm(vecs, axis=-1)[:, np.newaxis]

        for v_1, v_2 in itertools.product(vecs, vecs):
            v_1 /= np.linalg.norm(v_1)
            v_2 /= np.linalg.norm(v_2)

            r = rotation_matrix_from_vectors(v_1, v_2)

            self.assertTrue(np.allclose(np.dot(r, v_1), v_2))

    def test_rotation_matrix_from_axis_angle(self):
        """Test matrices from `rotation_matrix_from_axis_angle` for
        rotations about the three Cartesian axes, with 25 random angles,
        against analytical formulae."""

        num_angs = 25

        angles = np.random.random(num_angs) * 360.0

        for cart_dir in "x", "y", "z":
            for a in angles:
                r = rotation_matrix_from_axis_angle(
                    parse_direction(cart_dir), a
                )

                r_cmp = _rotation_matrix(cart_dir, a)

                self.assertTrue(np.allclose(r, r_cmp))

    def test_rotate_tensors(self):
        """Test `rotate_tensors` with the six worked examples from:
        https://empossible.net/wp-content/uploads/2020/06/Lecture-Tensor-Math.pdf
        """

        t = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])

        r_x_20 = rotation_matrix_from_axis_angle([1.0, 0.0, 0.0], 20.0)
        r_y_45 = rotation_matrix_from_axis_angle([0.0, 1.0, 0.0], 45.0)
        r_z_60 = rotation_matrix_from_axis_angle([0.0, 0.0, 1.0], 60.0)

        rots = [
            r_x_20,
            r_y_45,
            r_z_60,
            np.dot(r_y_45, r_x_20),
            np.dot(r_y_45, r_z_60),
            np.dot(r_z_60, np.dot(r_y_45, r_x_20)),
        ]

        ref_rot_t = [
            np.array(
                [
                    [1.0000, 0.0000, 0.0000],
                    [0.0000, 2.1170, -0.3214],
                    [0.0000, -0.3214, 2.8830],
                ],
                dtype=np.float64,
            ),
            np.array(
                [
                    [2.0000, 0.0000, 1.0000],
                    [0.0000, 2.0000, 0.0000],
                    [1.0000, 0.0000, 2.0000],
                ],
                dtype=np.float64,
            ),
            np.array(
                [
                    [1.7500, -0.4330, 0.0000],
                    [-0.4330, 1.2500, 0.0000],
                    [0.0000, 0.0000, 3.0000],
                ],
                dtype=np.float64,
            ),
            np.array(
                [
                    [1.9415, -0.2273, 0.9415],
                    [-0.2273, 2.1170, -0.2273],
                    [0.9415, -0.2273, 1.9415],
                ],
                dtype=np.float64,
            ),
            np.array(
                [
                    [2.3750, -0.3062, 0.6250],
                    [-0.3062, 1.2500, 0.3062],
                    [0.6250, 0.3062, 2.3750],
                ],
                dtype=np.float64,
            ),
            np.array(
                [
                    [2.2699, 0.0377, 0.6676],
                    [0.0377, 1.7886, 0.7017],
                    [0.6676, 0.7017, 1.9415],
                ],
                dtype=np.float64,
            ),
        ]

        for rot, ref in zip(rots, ref_rot_t):
            rot_t = rotate_tensors(t, rot)
            self.assertTrue(np.allclose(np.round(rot_t, 4), ref))


# ----------------------------------
# Tests for differentiation routines
# ----------------------------------


class TestDifferentiation(unittest.TestCase):
    """Class implementing unit tests for numerical differentiation
    routines."""

    def test_central_difference_order_1(self):
        """Test central difference schemes for calculating the
        derivatives of some polynomials."""

        precs = central_difference_available_precs(1)

        scale = 1.0e-2

        for p in precs:
            steps, coeffs = central_difference_coefficients(1, p)

            # Set absolute tolerance based on precision of central
            # difference scheme.

            abs_tol = np.power(10.0, -p)

            for f, f_prime in [
                (lambda x: x**2, lambda x: 2.0 * x),
                (lambda x: x**3, lambda x: 3.0 * x**2),
                (lambda x: x**4, lambda x: 4.0 * x**3),
            ]:
                x_vals = np.random.random((10,))

                for x in x_vals:
                    der = (coeffs * f(x + (scale * steps))).sum() / scale
                    self.assertTrue(np.isclose(der, f_prime(x), atol=abs_tol))


# -----------------------------
# Tests for quadrature routines
# -----------------------------


class TestQuadrature(unittest.TestCase):
    """Class implementing unit tests for numerical integration
    routines."""

    def test_unit_circle_rule(self):
        """Test vectors produced by `unit_circle_test_rule`."""

        num_pts = 16

        vecs, w = unit_circle_quad_rule(num_pts, ret="vectors_2d")

        ang_inc = (2.0 * np.pi) / num_pts

        for i in range(len(vecs)):
            v_1, v_2 = vecs[i], vecs[(i + 1) % num_pts]

            a = np.arccos(np.dot(v_1, v_2))

            self.assertTrue(np.isclose(a, ang_inc))

        p = 2.0 * np.pi * sum(w * _circle(*v) for v, w in zip(vecs, w))

        self.assertTrue(np.isclose(p, 2.0 * np.pi))

    def test_unit_sphere_lebedev_rule(self):
        """Test integration over the surface of the unit sphere using
        `unit_sphere_levedev_quad_rule`."""

        for p in unit_sphere_lebedev_quad_available_precs():
            vecs, w = unit_sphere_lebedev_quad_rule(p, ret="vectors")

            a = 4.0 * np.pi * sum(w * _sphere(*v) for v, w in zip(vecs, w))

            self.assertTrue(np.isclose(a, 4.0 * np.pi))


# ----
# Main
# ----


if __name__ == "__main__":
    unittest.main()
