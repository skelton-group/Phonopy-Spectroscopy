# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Test routines for infrared simulations on Pnma SnSe."""

# -------
# Imports
# -------

import os
import unittest

import xml.etree.ElementTree as ET

import numpy as np

from scipy.signal import find_peaks

from phonopy_spectroscopy.interfaces.phonopy_interface import (
    gamma_phonons_from_phono3py,
)

from phonopy_spectroscopy.interfaces.vasp_interface import (
    _parse_dielectric_function,
)

from phonopy_spectroscopy.ir.calculation import InfraredCalculation

from phonopy_spectroscopy.optics.optical_eigenmodes import (
    optical_spectra_from_optical_properties,
)

from phonopy_spectroscopy.gamma_phonons import get_gamma_phonons
from phonopy_spectroscopy.instrument import SPPolarisation

from phonopy_spectroscopy.utility.geometry import parse_direction

# ---------
# Constants
# ---------

_EXAMPLE_BASE_DIR_SNSE = r"../example/snse-pnma"

"""Path to SnSe (Pnma) example directory."""

# ----------------
# Helper functions
# ----------------


def verify_bruggeman_equation(eps, fracs, eps_eff):
    """Test whether an effective dielectric constant/function for a
    multiphase mixture satisfies the Bruggeman equation.

    Parameters
    ----------
    eps : list of (complex or array_like)
        Dielectric constants/functions of the components.
    fracs : list of float
        Volume fractions of the components.
    eps_eff : array_like
        Effective dielectric constant/function/

    Returns
    -------
    test : bool
        `True` if `eps_eff` satisfies the Bruggeman equation, otherwise
        `False`.
    """

    res = np.zeros_like(eps_eff, dtype=np.complex128)

    for e, f in zip(eps, fracs):
        res += f * (e - eps_eff) / (e + 2.0 * eps_eff)

    return np.allclose(res, 0.0)


# ------------------------------
# Tests for infrared simulations
# ------------------------------


class TestInfraredSimulation(unittest.TestCase):
    """Class implementing unit tests for infrared simulations."""

    def setUp(self):
        """Perform setup."""

        # Set up and store an InfraredCalculator object for testing.

        gamma_ph = gamma_phonons_from_phono3py(
            os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"POSCAR.Opt"),
            os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"mesh.yaml"),
            lws_file=os.path.join(
                _EXAMPLE_BASE_DIR_SNSE, r"kappa-m323216-g0.hdf5"
            ),
            irreps_file=os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"irreps.yaml"),
            born_file=os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"BORN"),
        )

        self._calc = InfraredCalculation(gamma_ph)

    def test_dielectric_function_vasp(self):
        """Compare the simulated dielectric function to a reference
        produced by the Vienna Ab initio Simulation Package (VASP) code.
        """

        # Load reference data.

        vasprun_xml = os.path.join(
            _EXAMPLE_BASE_DIR_SNSE, r"ir_ref/vasprun-Epsilon-FD+DFPT.xml"
        )

        tree = ET.parse(vasprun_xml)
        root = tree.getroot()

        eps_vasp = _parse_dielectric_function(
            vasprun_xml,
            root.findall("./calculation/dielectricfunction")[0],
        )

        e = eps_vasp.energy / (2.0 * np.pi)
        eps_e_ref = eps_vasp.epsilon

        # To compare to the reference calculation, we need to generate
        # a new PolarGammaPhonons object without calculated linewidths.

        calc = self._calc
        gamma_ph = calc.phonon_calculation

        gamma_ph_new = get_gamma_phonons(
            gamma_ph.structure,
            gamma_ph.frequencies,
            gamma_ph.eigenvectors,
            gamma_ph.epsilon_inf,
            gamma_ph.born_effective_charges,
            irreps=gamma_ph.irreps,
        )

        calc_new = InfraredCalculation(gamma_ph_new)

        e_step = np.mean(e[1:] - e[:-1])

        eps_ir = calc_new.dielectric_function(
            lw=(2.0 * e_step), e_range=(e[0], e[-1]), e_res=e_step, nac=False
        )

        # VASP does not include \eps_inf in its calculated dielectric
        # function.

        eps_f = eps_ir.epsilon - eps_ir.epsilon_inf[np.newaxis, :, :]

        # The diagonal components have some slight energy shifts and
        # differences in magnitude, both of which are expected, but are
        # otherwise visually very similar. We compare the real and
        # imaginary parts separately by using scipy.signal.find_peaks to
        # compare: (1) the number of peaks; (2) the absolute energy
        # shifts; and (3) and the differences in magnitude relative to
        # the highest feature in the reference.

        freq_tol = 0.1
        ints_tol = 5.0e-2

        for i in range(3):
            for eps_ii, eps_ii_ref in [
                (eps_f.real[:, i, i], eps_e_ref.real[:, i, i]),
                (eps_f.imag[:, i, i], eps_e_ref.imag[:, i, i]),
            ]:
                inds, _ = find_peaks(eps_ii)
                inds_ref, _ = find_peaks(eps_ii_ref)

                self.assertTrue(len(inds) == len(inds_ref))

                for idx, idx_ref in zip(inds, inds_ref):
                    self.assertTrue(np.abs(e[idx] - e[idx_ref]) < freq_tol)

                    scale = np.abs(
                        (eps_ii[idx] - eps_ii_ref[idx_ref]) / eps_ii_ref.max()
                    )

                    self.assertTrue(scale < ints_tol)

        # The off-diagonal components should be close to zero and
        # therefore equal between the two sets of data.

        inds = [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]]

        for i1, i2 in inds:
            self.assertTrue(
                np.allclose(eps_f[:, i1, i2], eps_e_ref[:, i1, i2])
            )

    def test_single_crystal_eigenmode_projection(self):
        """Check consistency of the optical spectra obtained with
        the single-crystal eigenmode projection method."""

        calc = self._calc

        eps_ir = calc.dielectric_function(nac=False)

        # Ensure the energy axes of subsequent simulations match the
        # "base" dielectric function.

        e_range = (eps_ir.energy[0], eps_ir.energy[-1])
        e_res = np.mean(eps_ir.energy[1:] - eps_ir.energy[:-1])

        # The three optical eigenmodes should be orthogonal, each should
        # be aligned to one of the Cartesian directions, and the
        # eigenvalues should be equal to the diagonal components of the
        # "source" dielectric function.

        oes = calc.optical_eigenmodes(e_range=e_range, e_res=e_res)

        oe_evecs_row = oes.eigenvectors.swapaxes(1, 2)

        for idx, v in enumerate(oe_evecs_row[0]):
            self.assertTrue(np.allclose(oe_evecs_row[:, idx], v))

        mapping = []

        for v in oe_evecs_row[0]:
            for idx, dirn in enumerate(["x", "y", "z"]):
                if np.allclose(v, parse_direction(dirn)):
                    mapping.append(idx)

        self.assertTrue(len(mapping) == len(set(mapping)))

        for i, idx in enumerate(mapping):
            assert np.allclose(
                oes.eigenvalues[:, i], eps_ir.epsilon[:, idx, idx]
            )

        # Eigenmode projections with suitable crystal faces and incident
        # polarisations should yield the same optical spectra as
        # calculated for the three diagonal components of \eps_ir.

        for i, (hkl, i_pol) in enumerate(
            [
                ((0, 0, 1), SPPolarisation([[0.0, -1.0]])),
                ((0, 0, 1), SPPolarisation([[-1.0, 0.0]])),
                ((1, 0, 0), SPPolarisation([[0.0, -1.0]])),
            ]
        ):
            idx = mapping.index(i)

            t_i, r_s, r_t, t_t = optical_spectra_from_optical_properties(
                oes.refractive_index[:, idx],
                oes.absorption_coefficient[:, idx],
                t=1.0e-3,
            )

            sp = calc.incoherent_single_crystal_on_axis_spectrum(
                hkl, i_pol=i_pol, t=1.0e-3, e_range=e_range, e_res=e_res
            )

            self.assertTrue(
                np.allclose(sp.intrinsic_absorbance, -1.0 * np.log10(t_i))
            )

            self.assertTrue(np.allclose(sp.single_reflectivity, r_s))
            self.assertTrue(np.allclose(sp.total_reflectivity, r_t))
            self.assertTrue(np.allclose(sp.transmission, t_t))


# ----
# Main
# ----

if __name__ == "__main__":
    unittest.main()
