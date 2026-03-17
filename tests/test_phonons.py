# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Test routines for core functionality of the `PolarGammaPhonons`
object."""


# -------
# Imports
# -------


import os
import unittest

import xml.etree.ElementTree as ET

import numpy as np

from phonopy_spectroscopy.constants import ZERO_TOLERANCE

from phonopy_spectroscopy.interfaces.phonopy_interface import (
    gamma_phonons_from_phono3py,
    gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml,
)

from phonopy_spectroscopy.interfaces.vasp_interface import (
    _parse_dielectric_constant,
)


# ---------
# Constants
# ---------


_EXAMPLE_BASE_DIR_SNSE = r"../example/snse-pnma"


# -----------------------------------
# Tests for  PolarGammaPhonons object
# -----------------------------------


class TestPolarGammaPhonons(unittest.TestCase):
    """Class implementing unit tests for the `PolarGammaPhonons` object
    functionality."""

    def setUp(self):
        """Perform setup."""

        # Set up and store a PolarGammaPhonons object for testing.

        self._gamma_ph = gamma_phonons_from_phono3py(
            os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"POSCAR.Opt"),
            os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"mesh.yaml"),
            lws_file=os.path.join(
                _EXAMPLE_BASE_DIR_SNSE, r"kappa-m323216-g0.hdf5"
            ),
            irreps_file=os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"irreps.yaml"),
            born_file=os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"BORN"),
        )

    def test_eps_ionic(self):
        r"""Test calculation of the ionic contribution to the static
        dielectric constant (\eps_ionic) against a reference produced
        by the Vienna Ab initio Simulation Package (VASP) code."""

        # Load reference data.

        vasprun_xml = os.path.join(
            _EXAMPLE_BASE_DIR_SNSE, r"ir_ref/vasprun-Epsilon-FD+DFPT.xml"
        )

        tree = ET.parse(vasprun_xml)
        root = tree.getroot()

        eps_ionic_ref = _parse_dielectric_constant(
            vasprun_xml,
            root.findall('./calculation/varray[@name="epsilon_ion"]')[0],
        )

        eps_ionic = self._gamma_ph.epsilon_ionic

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

        tolerance = 0.15

        diff = (
            eps_ionic[mask_nonzero] - eps_ionic_ref[mask_nonzero]
        ) / eps_ionic_ref[mask_nonzero]

        self.assertTrue((np.abs(diff) < tolerance).all())

    def test_nac(self):
        """Test non-analytical correction (NAC, LO/TO splitting)."""

        # "Sanity check" frequencies of original calculation.

        freqs, _ = gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml(
            os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"nac_ref/qpoints.yaml")
        )

        # Very small differences in frequency, likely because the "base"
        # calculation was run with a different installation/version of
        # Phonopy to that used to generate the reference data for
        # testing the NAC implementation.

        self.assertTrue(
            np.allclose(freqs, self._gamma_ph.frequencies, atol=1.0e-5)
        )

        # Recalculate frequencies with NAC correction for q along the
        # three principal axes and compare frequencies to reference
        # values from Phonopy.

        for q in (1, 0, 0), (0, 1, 0), (0, 0, 1):
            f = os.path.join(
                _EXAMPLE_BASE_DIR_SNSE,
                r"nac_ref/qpoints-{0}{1}{2}.yaml".format(*q),
            )

            freqs, _ = gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml(f)
            gamma_ph = self._gamma_ph.gamma_phonons_with_nac(q)

            self.assertTrue(
                np.allclose(np.sort(gamma_ph.frequencies), freqs, atol=1.0e-4)
            )

    def test_pop_freq(self):
        """Test calculation of the polar-optic phonon (POP) frequency
        against a reference value."""

        w_po = self._gamma_ph.pop_frequency()

        # Reference data from CalcPOP.py, a standalone script that
        # implements the same algorithm as used in AMSET, including
        # using the same source of Lebedev quadrature weights.

        w_po_ref = 3.284157105842913

        self.assertTrue(np.isclose(w_po, w_po_ref))


# ----
# Main
# ----


if __name__ == "__main__":
    unittest.main()
