# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Test routines for parts of the phonopy_spectroscopy.raman package
using a set of model Raman tensors for Si and (in part) results from the
literature.

The model tensors are taken from:

    U. Ramabadran and B. Roughani, Mater. Sci. Eng. B 230, 31 (2018),
    DOI: 10.1016/j.mseb.2017.12.040

The single-crystal tests replicate a subset of the simulations in the
paper and verify the results against the analytical formulae derived
therein.

The powder Raman tests use the model tensors to perform a simple chi
rotation and compare results from an analytical formula with numerical
integration.

The tests also compare the depolarisation ratio to the analytical
formula in the original work in:

    D. Porezag and M. R. Pederson, Phys. Rev. B 54, 7830 (1996), DOI:
    10.1103/PhysRevB.54.7830
"""


# -------
# Imports
# -------


import unittest

import numpy as np

from phonopy_spectroscopy.instrument import Geometry, Polarisation

from phonopy_spectroscopy.raman.intensity import (
    calculate_single_crystal_raman_intensities,
    calculate_powder_raman_intensities,
)

from phonopy_spectroscopy.structure import Structure

from phonopy_spectroscopy.utility.geometry import (
    rotation_matrix_from_vectors,
    rotation_matrix_from_axis_angle,
)


# ----------------
# Helper functions
# ----------------


def _get_pol_vecs_for_psi_rot(theta, phi, psi, psi_0):
    r"""Generate incident polarisation vectors e_i for a psi rotation
    using Eq. 7a-7c in the paper linked in the top-level docstring."""

    theta, phi = np.radians(theta), np.radians(phi)
    psi = np.radians(psi + psi_0)

    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    cos_phi, sin_phi = np.cos(phi), np.sin(phi)
    cos_psi, sin_psi = np.cos(psi), np.sin(psi)

    e_is = np.array(
        [
            cos_theta * cos_phi * cos_psi - sin_phi * sin_psi,
            cos_theta * sin_phi * cos_psi + cos_phi * sin_psi,
            -1.0 * sin_theta * cos_psi,
        ],
        dtype=np.float64,
    ).T

    return e_is


def _calculate_powder_raman_intensity_legacy(r_t):
    r"""Given a Raman tensor `r_t`, calculate the powder intensities
    I_\par/I_\per parallel and perpendicular to the laser polarisation.
    """

    a_p_2 = ((1.0 / 3.0) * np.trace(r_t)) ** 2

    b_p_2 = 0.5 * (
        (r_t[0, 0] - r_t[1, 1]) ** 2
        + (r_t[0, 0] - r_t[2, 2]) ** 2
        + (r_t[1, 1] - r_t[2, 2]) ** 2
        + 6.0 * (r_t[0, 1] ** 2 + r_t[0, 2] ** 2 + r_t[1, 2] ** 2)
    )

    return (a_p_2 + (4.0 / 45.0) * b_p_2, (3.0 / 45.0) * b_p_2)


# ---------------------------
# Tests for Raman simulations
# ---------------------------


class TestRamanSimulation(unittest.TestCase):
    """Class implementing unit tests for Raman simulations."""

    def setUp(self):
        """Perform setup."""

        # "Empty" structure for obtaining surface normals.

        a = 5.431020511

        self._struct = Structure(
            [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]], [], [], []
        )

        # Model Raman tensors.

        b = 1.0

        r_x = [[0.0, 0.0, 0.0], [0.0, 0.0, b], [0.0, b, 0.0]]
        r_y = [[0.0, 0.0, b], [0.0, 0.0, 0.0], [b, 0.0, 0.0]]
        r_z = [[0.0, b, 0.0], [b, 0.0, 0.0], [0.0, 0.0, 0.0]]

        self._r_t = np.array([r_x, r_y, r_z], dtype=np.float64)

        # Measurement geometry.

        self._geom = Geometry("z", "-z")

    def test_single_crystal(self):
        """Test the single-crystal Raman simulations in Figs 3, 5, 8 and
        10 in the paper linked in the top-level docstring."""

        # (theta, phi) rotations used in different simulations.

        rotations = [
            (0.0, 0.0),
            (45.0, 0.0),
            (90.0, 0.0),
            (0.0, 10.0),
            (45.0, 15.0),
            (45.0, 30.0),
            (45.0, 45.0),
        ]

        # Generate rotation matrix to align the crystal (0, 0, 1)
        # direction with -z.

        rot_001 = rotation_matrix_from_vectors(
            self._struct.real_space_normal((0, 0, 1)),
            -1.0 * self._geom.collection_direction,
        )

        # Generate polarisation vectors for the \psi rotations performed
        # in the paper.

        pols = Polarisation.from_rotation("z")

        # Loop over simulation parameters.

        for theta, phi in rotations:
            # \theta and \phi rotations correspond to rotations about
            # the y and z axes, respectively.

            r = np.dot(
                rotation_matrix_from_axis_angle("z", phi),
                rotation_matrix_from_axis_angle("y", theta),
            )

            # Reference polarisation vectors.

            pol_vecs_ref = _get_pol_vecs_for_psi_rot(
                theta, phi, np.arange(0.0, 360.0 + 0.25, 2.5), 0.0
            )

            # Rotating the polarisation vectors defined in the (x, y)
            # plane should give the same results as the analytical
            # formulae.

            for p, p_ref in zip(pols, pol_vecs_ref):
                self.assertTrue(np.allclose(np.dot(r, p.vectors[0]), p_ref))

            # To emulate the experiments in the paper we need to rotate
            # the crystal by the inverse of the polarisation rotation.

            rot = np.dot(rot_001, np.linalg.inv(r))

            ints = [
                calculate_single_crystal_raman_intensities(
                    self._r_t,
                    self._geom,
                    p,
                    p,
                    rot=rot,
                )
                for p in pols
            ]

            # Reference intensnties.

            ints_ref = []

            for e_i, e_s in zip(pol_vecs_ref, pol_vecs_ref):
                ints_ref.append(
                    [np.dot(e_s, np.dot(t, e_i)) ** 2 for t in self._r_t]
                )

            self.assertTrue(np.allclose(ints, ints_ref))

    def test_powder_1(self):
        """Test a powder chi rotation with the analytical formula and
        numerical integration."""

        # Generate a list of Polarisation objects for a rotation around
        # the z-axis with the default parameters.

        i_pol = Polarisation.from_direction("x")
        s_pols = Polarisation.from_rotation("z", step=90.0)

        ints_1 = [
            calculate_powder_raman_intensities(
                self._r_t,
                self._geom,
                i_pol,
                s_pol,
            )
            for s_pol in s_pols
        ]

        ints_2 = [
            calculate_powder_raman_intensities(
                self._r_t, self._geom, i_pol, s_pol, method="quad"
            )
            for s_pol in s_pols
        ]

        ints_3 = [
            calculate_powder_raman_intensities(
                self._r_t,
                self._geom,
                i_pol,
                s_pol,
                method="leb+circ",
                lc_prec=5,
            )
            for s_pol in s_pols
        ]

        self.assertTrue(np.allclose(ints_2, ints_1))
        self.assertTrue(np.allclose(ints_3, ints_1))

    def test_powder_2(self):
        """Test calculations of parallel and perpendicular powder
        intensities with the analytical formula against the original
        paper linked in the top-level docstring."""

        powder_ref = [
            _calculate_powder_raman_intensity_legacy(t) for t in self._r_t
        ]

        pol_par = Polarisation.from_direction("x")
        pol_per = Polarisation.from_direction("y")

        ints_par = calculate_powder_raman_intensities(
            self._r_t, self._geom, pol_par, pol_par
        )

        ints_per = calculate_powder_raman_intensities(
            self._r_t, self._geom, pol_par, pol_per
        )

        for (i_par_ref, i_per_ref), i_par, i_per in zip(
            powder_ref, ints_par, ints_per
        ):
            self.assertTrue(np.isclose(i_par_ref, i_par))
            self.assertTrue(np.isclose(i_per_ref, i_per))

    def test_powder_with_po_1(self):
        """Test calculations of the depolarisation ratio of an ideal
        powder and a powder with a small preferred orientation."""

        pol_par = Polarisation.from_direction("x")
        pol_per = Polarisation.from_direction("y")

        ints_par = calculate_powder_raman_intensities(
            self._r_t, self._geom, pol_par, pol_par
        )

        ints_per = calculate_powder_raman_intensities(
            self._r_t, self._geom, pol_par, pol_per
        )

        rho = ints_per / ints_par

        # The three Si Raman tensors only have off-diagonal components
        # and should all have \rho = 0.75.

        self.assertTrue(np.allclose(rho, 0.75))

        rho_pos = []

        for hkl in (1, 0, 0), (0, 1, 0), (0, 0, 1):
            po_surf_norm = self._struct.real_space_normal(hkl)

            ints_par_po = calculate_powder_raman_intensities(
                self._r_t,
                self._geom,
                pol_par,
                pol_par,
                po_surf_norm=po_surf_norm,
                po_eta=0.1,
            )

            ints_per_po = calculate_powder_raman_intensities(
                self._r_t,
                self._geom,
                pol_par,
                pol_per,
                po_surf_norm=po_surf_norm,
                po_eta=0.1,
            )

            rho_po = ints_per_po / ints_par_po

            # Preferred orientation should result in \rho != 0.75.

            self.assertTrue((rho_po != 0.75).all())

        # Different preferred orientation should not affect the three
        # Raman tensors equally.

        for i in range(len(rho_pos)):
            for j in range(i + 1, len(rho_pos)):
                self.assertTrue((rho_pos[i] != rho_pos[j]).all())

    def test_powder_with_po_2(self):
        """Test calculations of the parallel and perpenducular powder
        intensnties for a powder with a small preferred orientation
        along (0, 0, 1) with the two numerical integration approaches.
        """

        # This may need adjusting depending on the value of \eta.

        lc_prec = 21

        pol_par = Polarisation.from_direction("x")
        pol_per = Polarisation.from_direction("y")

        po_surf_norm = self._struct.real_space_normal((0, 0, 1))

        for s_pol in pol_par, pol_per:
            ints_po_nquad = calculate_powder_raman_intensities(
                self._r_t,
                self._geom,
                pol_par,
                s_pol,
                po_eta=0.1,
                po_surf_norm=po_surf_norm,
                method="quad",
            )

            ints_po_lebedev = calculate_powder_raman_intensities(
                self._r_t,
                self._geom,
                pol_par,
                s_pol,
                po_eta=0.1,
                po_surf_norm=po_surf_norm,
                method="leb+circ",
                lc_prec=lc_prec,
            )

            self.assertTrue(np.allclose(ints_po_lebedev, ints_po_nquad))


# ----
# Main
# ----


if __name__ == "__main__":
    unittest.main()
