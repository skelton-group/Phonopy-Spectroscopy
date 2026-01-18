# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Test routines for infrared simulations with a complete calculation on
Pnma SnSe."""


# -------
# Imports
# -------


import os
import unittest

import xml.etree.ElementTree as ET

import numpy as np

from phonopy_spectroscopy.constants import VASP_TO_THZ, ZERO_TOLERANCE

from phonopy_spectroscopy.interfaces.phonopy_interface import (
    gamma_phonons_from_phono3py,
    hf_dielectric_and_born_from_born,
)

from phonopy_spectroscopy.interfaces.vasp_interface import (
    _parse_dielectric_constant,
    _parse_dielectric_function,
)

from phonopy_spectroscopy.ir.calculation import InfraredCalculation

from phonopy_spectroscopy.phonon import GammaPhonons

import matplotlib.pyplot as plt


# ---------
# Constants
# ---------


_EXAMPLE_BASE_DIR = r"../example/snse-pnma"


# ------------------------------
# Tests for infrared simulations
# ------------------------------


class TestInfratedSimulations(unittest.TestCase):
    def setUp(self):
        """Perform setup."""

        # Set up and store an InfraredCalculator object for testing.

        lws_file = os.path.join(
            _EXAMPLE_BASE_DIR, r"kappa-m16816.Prim.FullPP.hdf5"
        )

        gamma_ph = gamma_phonons_from_phono3py(
            os.path.join(_EXAMPLE_BASE_DIR, r"POSCAR.Opt"),
            os.path.join(_EXAMPLE_BASE_DIR, r"mesh.yaml"),
            lws_file=lws_file,
            irreps_file=os.path.join(_EXAMPLE_BASE_DIR, r"irreps.yaml"),
        )

        eps_inf, born_charges = hf_dielectric_and_born_from_born(
            os.path.join(_EXAMPLE_BASE_DIR, r"BORN"),
            gamma_ph.structure,
        )

        self._calc = InfraredCalculation(
            gamma_ph, born_charges, eps_inf=eps_inf
        )

    def test_hessian(self):
        """Test calculation of the Hessian matrix."""

        gamma_ph = self._calc.gamma_phonons

        h = gamma_ph.hessian()

        # Reconstruct the dynamical matrix, diagonalise, and ensure we
        # get the same frequencies and eigenvectors.

        sqrt_m_rep = np.zeros((gamma_ph.num_modes,), dtype=np.float64)

        for i in range(3):
            sqrt_m_rep[i::3] = gamma_ph.structure.atomic_masses

        sqrt_m_rep = np.sqrt(sqrt_m_rep)

        # Construct the dynamical matrix from the Hessian (= 2nd-order
        # force constsnta, FC2) by dividing through by sqrt(m_i * mj).

        d = h.copy()
        d /= sqrt_m_rep[:, np.newaxis]
        d /= sqrt_m_rep[np.newaxis, :]

        # Diagonalising the dynamical matrix yields \omega^2 as the
        # eigenvalues.

        freqs_sq, evecs = np.linalg.eigh(d)

        # We can't sqrt() a -ve number, so we follow convention and take
        # the square root of the absolute value and show imaginary
        # frequencies as -ve numbers.

        freqs = np.copysign(np.sqrt(np.abs(freqs_sq)), freqs_sq) * VASP_TO_THZ

        # Rearrange the eigenvectors to the same data layout as the
        # GammaPhonons object.

        evecs_temp = np.zeros(
            (gamma_ph.num_modes, gamma_ph.structure.num_atoms, 3),
            dtype=np.float64,
        )

        for i in range(gamma_ph.num_modes):
            evecs_temp[i, :, :] = evecs[:, i].reshape(
                gamma_ph.structure.num_atoms, 3
            )

        evecs = evecs_temp

        # The forward/reverse transformation can produce small changes
        # to the frequencies of the acoustic modes and rotations of the
        # eigenvectors. We therefore exclude them from the comparison.

        # (This problem could also affect degenerate optic modes, but
        # the optic modes in Pnma SnSe are al singly degenerate.)

        excl_inds = gamma_ph.get_acoustic_mode_indices()

        band_inds = [
            i for i in range(gamma_ph.num_modes) if i not in excl_inds
        ]

        self.assertTrue(
            np.allclose(freqs[band_inds], gamma_ph.frequencies[band_inds])
        )

        for e_cmp, e_ref in zip(
            evecs[band_inds], gamma_ph.eigenvectors[band_inds]
        ):
            # The forward/reverse transformation can invert the sign of
            # the eigenvectors.

            e_cmp = np.copysign(e_cmp, e_ref)

            self.assertTrue(np.allclose(e_cmp, e_ref))

    def test_pop_freq(self):
        """Test calculation of the polar-optic phonon (POP) frequency
        against a reference value."""

        w_po = self._calc.pop_frequency()

        # Reference data from CalcPOP.py, a standalone script that
        # implements the same algorithm as used in AMSET, including
        # using the same source of Lebedev quadrature weights.

        w_po_ref = 3.265018552508473

        self.assertTrue(np.allclose(w_po, w_po_ref))

    def test_eps_ionic(self):
        """Test calculation of the ionic contribution to the static
        dielectric constant (epsilon_ionic) against a reference produced
        by the Vienna Ab initio Simulation Package (VASP) code."""

        # Load reference data.

        vasprun_xml = os.path.join(
            _EXAMPLE_BASE_DIR, r"ir_ref/vasprun-PBEsol+D3-Epsilon.xml"
        )

        tree = ET.parse(vasprun_xml)
        root = tree.getroot()

        eps_ionic_ref = _parse_dielectric_constant(
            vasprun_xml,
            root.findall('./calculation/varray[@name="epsilon_ion"]')[0],
        )

        eps_ionic = self._calc.epsilon_ionic

        # Check the the calculated and reference \epsilon_ionic have
        # zero elements in the same place.

        mask = np.abs(eps_ionic) < ZERO_TOLERANCE
        mask_ref = np.abs(eps_ionic_ref) < ZERO_TOLERANCE

        self.assertTrue((mask == mask_ref).all())

        # Check the non-zero elements are the same to within a
        # tolerance. Calculating \epsilon_ionic involves inverting the
        # Hessian matrix, which is generally an ill-conditioned problem.
        # In this case, some variation between the calculated and
        # reference results is inevitable.

        mask_nonzero = np.logical_not(mask)

        tolerance = 0.2

        diff = (
            eps_ionic[mask_nonzero] - eps_ionic_ref[mask_nonzero]
        ) / eps_ionic_ref[mask_nonzero]

        self.assertTrue((np.abs(diff) < tolerance).all())

    def test_tensor_dielectric_func_vasp(self):
        """Compare the simulated dielectric function to a reference
        produced by the Vienna Ab initio Simulation Package (VASP) code.
        """

        # Load reference data.

        vasprun_xml = os.path.join(
            _EXAMPLE_BASE_DIR, r"ir_ref/vasprun-PBEsol+D3-Epsilon.xml"
        )

        tree = ET.parse(vasprun_xml)
        root = tree.getroot()

        e_ref, eps_e_ref = _parse_dielectric_function(
            vasprun_xml,
            root.findall("./calculation/dielectricfunction")[0],
        )

        e_ref /= 2.0 * np.pi

        # To compare to the reference calculation, we need to generate
        # new GammaPhonons and InfraredCalculator objects without the
        # calculated linewidths.

        calc = self._calc
        gamma_ph = calc.gamma_phonons

        gamma_ph_new = GammaPhonons(
            gamma_ph.structure,
            gamma_ph.frequencies,
            gamma_ph.eigenvectors,
            irreps=gamma_ph.irreps,
        )

        calc_new = InfraredCalculation(
            gamma_ph_new, calc.born_effective_charges, eps_inf=calc.epsilon_inf
        )

        dielectric_func = calc_new.tensor_dielectric_function(
            lw=2.0 * np.mean(e_ref[1:] - e_ref[:-1]),
            add_eps_inf=False,
            x=e_ref,
        )

        e, eps_e = dielectric_func.x, dielectric_func.epsilon

        # The off-diagonal components should be close to zero and
        # therefore equal between the two sets of data.

        inds = [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]]

        for i1, i2 in inds:
            self.assertTrue(
                np.allclose(eps_e[:, i1, i2], eps_e_ref[:, i1, i2])
            )

        # The diagonal components have some slight frequency shifts
        # between the two sets of data. These are to be expected, and
        # are visually insignificant, but make direct comparison of the
        # two sets of data difficult. As a workaround, we apply a small
        # boxcar averaging, then check the two functions are "almost
        # equal", with only the largest five differences >1%.

        conv_kernel = np.ones((5,)) / 5.0

        ave_eps_e = np.zeros((3, len(e)), dtype=np.complex128)

        for i in range(3):
            ave_eps_e[i, :] = np.convolve(
                eps_e[:, i, i], conv_kernel, mode="same"
            )

        ave_eps_e_ref = np.zeros((3, len(e)), dtype=np.complex128)

        for i in range(3):
            ave_eps_e_ref[i, :] = np.convolve(
                eps_e_ref[:, i, i], conv_kernel, mode="same"
            )

        diff = np.sort(np.abs((ave_eps_e - ave_eps_e_ref) / ave_eps_e_ref))

        s_tolerance = 1.0e-2
        l_tolerance = 0.25

        self.assertTrue((diff[:-5] < s_tolerance).all())
        self.assertTrue((diff[-5:] < l_tolerance).all())


# ----
# Main
# ----


if __name__ == "__main__":
    unittest.main()
