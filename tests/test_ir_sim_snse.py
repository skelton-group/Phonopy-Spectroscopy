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

from phonopy_spectroscopy.ir.multiphase_mixture import (
    bruggeman_two_phase_scalar,
    bruggeman_three_phase_scalar,
    bruggeman_multiphase_scalar,
)

from phonopy_spectroscopy.ir.optical_eigenmodes import (
    optical_spectra_from_optical_properties,
)

from phonopy_spectroscopy.gamma_phonons import PolarGammaPhonons
from phonopy_spectroscopy.instrument import Polarisation

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

        e, eps_e_ref = _parse_dielectric_function(
            vasprun_xml,
            root.findall("./calculation/dielectricfunction")[0],
        )

        e /= 2.0 * np.pi

        # To compare to the reference calculation, we need to generate
        # a new PolarGammaPhonons object without calculated linewidths.

        calc = self._calc
        gamma_ph = calc.phonon_calculation

        gamma_ph_new = PolarGammaPhonons(
            gamma_ph.structure,
            gamma_ph.frequencies,
            gamma_ph.eigenvectors,
            gamma_ph.epsilon_inf,
            gamma_ph.born_effective_charges,
            irreps=gamma_ph.irreps,
        )

        calc_new = InfraredCalculation(gamma_ph_new)

        eps_ir = calc_new.dielectric_function(
            lw=2.0 * np.mean(e[1:] - e[:-1]),
            x=e,
        )

        # VASP does not include \eps_inf in its calculated dielectric
        # function.

        eps_e = eps_ir.epsilon - eps_ir.epsilon_inf[np.newaxis, :, :]

        # The diagonal components have some slight frequency shifts and
        # differences in magnitude, both of which are expected, but are
        # otherwise visually very similar. We compare the real and
        # imaginary parts separately by using scipy.signal.find_peaks to
        # compare: (1) the number of peaks; (2) the absolute frequency
        # shifts; and (3) and the differences in magnitude relative to
        # the highest feature in the reference.

        freq_tol = 0.1
        ints_tol = 5.0e-2

        for i in range(3):
            for f, f_ref in [
                (eps_e.real[:, i, i], eps_e_ref.real[:, i, i]),
                (eps_e.imag[:, i, i], eps_e_ref.imag[:, i, i]),
            ]:
                inds, _ = find_peaks(f)
                inds_ref, _ = find_peaks(f_ref)

                self.assertTrue(len(inds) == len(inds_ref))

                for idx, idx_ref in zip(inds, inds_ref):
                    self.assertTrue(np.abs(e[idx] - e[idx_ref]) < freq_tol)

                    scale = np.abs((f[idx] - f_ref[idx_ref]) / f_ref.max())
                    self.assertTrue(scale < ints_tol)

        # The off-diagonal components should be close to zero and
        # therefore equal between the two sets of data.

        inds = [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]]

        for i1, i2 in inds:
            self.assertTrue(
                np.allclose(eps_e[:, i1, i2], eps_e_ref[:, i1, i2])
            )

    def test_powder(self):
        r"""Check consistency of the effective dielectric functions
        \eps_eff obtained with different powder methods to the "source"
        bulk infrared dielectric function."""

        calc = self._calc

        # Reference IR dielectric function.

        eps_ir = calc.dielectric_function().epsilon

        # SnSe has orthorhombic symmetry, so eps_ir should be diagonal
        # and the eps_eff for the EMA and the average of the optical
        # eigenmode eigenvalues should both be equivalent to averaging
        # the trace.

        eps_eff = np.trace(eps_ir, axis1=1, axis2=2) / 3.0

        sp_ema = calc.powder_optical_spectrum_ema()

        self.assertTrue(
            np.allclose(sp_ema.optical_eigenmodes.eigenvalues[:, 0], eps_eff)
        )

        sp_oea = calc.powder_optical_spectrum_eigenmode_average()

        eps_eff_oea = np.mean(sp_oea.optical_eigenmodes.eigenvalues, axis=1)

        self.assertTrue(np.allclose(eps_eff_oea, eps_eff))

    def test_single_crystal_eigenmode_projection(self):
        """Check consistency of the optical spectra obtained with
        the single-crystal eigenmode projection method."""

        calc = self._calc

        eps_ir = calc.dielectric_function().epsilon

        # The three optical eigenmodes should be orthogonal, each should
        # be aligned to one of teh Cartesian directions, and the
        # eigenvalues should be equal to the diagonal components of the
        # "source" dielectric function.

        oes = calc.optical_eigenmodes()

        vecs = oes.eigenvectors_row[0]

        for idx, v in enumerate(vecs):
            self.assertTrue(np.allclose(oes.eigenvectors[:, idx], v))

        mapping = []

        for v in vecs:
            for idx, dirn in enumerate(["x", "y", "z"]):
                if np.allclose(v, parse_direction(dirn)):
                    mapping.append(idx)

        self.assertTrue(len(mapping) == len(set(mapping)))

        for i, idx in enumerate(mapping):
            assert np.allclose(oes.eigenvalues[:, i], eps_ir[:, idx, idx])

        # Eigenmode projections with suitable crystal faces and incident
        # polarisations should yield the same optical spectra as
        # calculated for the three diagonal components of \eps_ir.

        for i, (hkl, i_pol) in enumerate(
            [
                ((0, 0, 1), Polarisation.from_direction("-x")),
                ((0, 0, 1), Polarisation.from_direction("-y")),
                ((1, 0, 0), Polarisation.from_direction("-x")),
            ]
        ):
            idx = mapping.index(i)

            a_int, r_s, r_t, t = optical_spectra_from_optical_properties(
                oes.refractive_index[:, idx],
                oes.absorption_coefficient[:, idx],
                t=1.0e-3,
            )

            sp = calc.single_crystal_eigenmode_projection(
                hkl, i_pol=i_pol, t=1.0e-3
            )

            self.assertTrue(np.allclose(sp.intrinsic_absorbance, a_int))
            self.assertTrue(np.allclose(sp.single_reflectivity, r_s))
            self.assertTrue(np.allclose(sp.total_reflectivity, r_t))
            self.assertTrue(np.allclose(sp.transmission, t))

    def test_bruggeman_multiphase(self):
        """Test the implementation of the Bruggeman multiphase models."""

        calc = self._calc

        eps_kbr = 4.9
        eps_air = 1.0

        sp_ema = calc.powder_optical_spectrum_ema()

        eps_snse = sp_ema.optical_eigenmodes.eigenvalues[:, 0]

        # Pure powder with 90% theoretical denisty.

        eps_eff_ld = bruggeman_two_phase_scalar(eps_snse, eps_air, 0.9, 0.1)

        self.assertTrue(
            verify_bruggeman_equation(
                [eps_snse, eps_air], [0.9, 0.1], eps_eff_ld
            )
        )

        # 5% KBr pellet.

        eps_eff_kbr = bruggeman_two_phase_scalar(eps_snse, eps_kbr, 0.05, 0.95)

        self.assertTrue(
            verify_bruggeman_equation(
                [eps_snse, eps_kbr], [0.05, 0.95], eps_eff_kbr
            )
        )

        # 5% KBr pellet with 90% denisity.

        eps_eff_kbr_ld = bruggeman_three_phase_scalar(
            eps_snse, eps_kbr, eps_air, 0.9 * 0.05, 0.9 * 0.95, 0.1
        )

        self.assertTrue(
            verify_bruggeman_equation(
                [eps_snse, eps_kbr, eps_air],
                [0.9 * 0.05, 0.9 * 0.95, 0.1],
                eps_eff_kbr_ld,
            )
        )

        # Test the implementation via keywords to
        # powder_optical_spectrum_ema().

        sp_ema_ld = calc.powder_optical_spectrum_ema(m_rho=0.9)

        self.assertTrue(
            np.allclose(
                sp_ema_ld.optical_eigenmodes.eigenvalues[:, 0], eps_eff_ld
            )
        )

        sp_ema_kbr = calc.powder_optical_spectrum_ema(m_f=0.05, m_eps=eps_kbr)

        self.assertTrue(
            np.allclose(
                sp_ema_kbr.optical_eigenmodes.eigenvalues[:, 0],
                eps_eff_kbr,
            )
        )

        sp_ema_kbr_ld = calc.powder_optical_spectrum_ema(
            m_f=0.05, m_eps=eps_kbr, m_rho=0.9
        )

        self.assertTrue(
            np.allclose(
                sp_ema_kbr_ld.optical_eigenmodes.eigenvalues[:, 0],
                eps_eff_kbr_ld,
            )
        )

        # General multiphase solver.

        for eps, fracs, eps_eff_ref in [
            ([eps_snse, eps_air], [0.9, 0.1], eps_eff_ld),
            ([eps_snse, eps_kbr], [0.05, 0.95], eps_eff_kbr),
            (
                [eps_snse, eps_kbr, eps_air],
                [0.9 * 0.05, 0.9 * 0.95, 0.1],
                eps_eff_kbr_ld,
            ),
        ]:
            eps_eff = bruggeman_multiphase_scalar(eps, fracs)

            self.assertTrue(verify_bruggeman_equation(eps, fracs, eps_eff))

            # Verify that the general solution matches the analytical
            # formulae.

            self.assertTrue(np.allclose(eps_eff, eps_eff_ref))


# ----
# Main
# ----

if __name__ == "__main__":
    unittest.main()
