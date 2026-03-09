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

from comparison_helper import compare_irreps


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
            gamma_ph, evec_proj = self._gamma_ph.gamma_phonons_with_nac(q)

            self.assertTrue(
                np.allclose(np.sort(gamma_ph.frequencies), freqs, atol=1.0e-5)
            )

            self.assertTrue(np.allclose(np.diag(evec_proj), 1.0, atol=0.1))

            # Check individual eignvector projections/sorting. SnSe does
            # not have degenerate modes, so the projection coefficients
            # should be high (>0.9).

            for idx, (evec_old, evec_new) in enumerate(
                zip(self._gamma_ph.eigenvectors, gamma_ph.eigenvectors)
            ):
                proj = np.abs(np.dot(evec_new.flat, evec_old.flat))

                self.assertTrue(np.isclose(proj, evec_proj[idx, idx]))
                self.assertTrue(proj > 0.9)

            # Linewidths should be passed through.

            self.assertTrue(
                np.allclose(gamma_ph.linewidths, self._gamma_ph.linewidths)
            )

            # Since SnSe does not have degenerate modes, the Irreps
            # objects from the original and corrected calculations
            # should be equivalent.

            self.assertTrue(
                compare_irreps(gamma_ph.irreps, self._gamma_ph.irreps)
            )


# ----
# Main
# ----


if __name__ == "__main__":
    unittest.main()
