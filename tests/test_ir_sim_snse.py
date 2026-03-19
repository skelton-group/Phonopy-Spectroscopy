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

from phonopy_spectroscopy.instrument import Polarisation
from phonopy_spectroscopy.ir.calculation import InfraredCalculation

from phonopy_spectroscopy.ir.spectrum_funcs import (
    bruggeman_two_phase_mixture,
    bruggeman_three_phase_mixture,
    bruggeman_multiphase_mixture,
)

from phonopy_spectroscopy.phonon import PolarGammaPhonons


# ---------
# Constants
# ---------


_EXAMPLE_BASE_DIR_SNSE = r"../example/snse-pnma"


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

    def test_optical_spectra_1(self):
        r"""Compare the effective dielectric functions \eps_eff for
        optical spectra obtained with different methods to the "source"
        bulk infrared dielectric function."""

        calc = self._calc

        # Reference IR dielectric function.

        eps_ir = calc.dielectric_function().epsilon

        # SnSe has orthorhombic symmetry, so eps_ir should be diagonal
        # and both powder models should give the same eps_eff as
        # averaging the trace.

        eps_eff = np.trace(eps_ir, axis1=1, axis2=2) / 3.0

        sp_ema = calc.powder_optical_spectrum_ema(t=1.0e-6)
        self.assertTrue(np.allclose(sp_ema.epsilon, eps_eff))

        sp_ave_oe = calc.powder_optical_spectrum_eigenmode_average(t=1.0e-6)
        self.assertTrue(np.allclose(sp_ave_oe.epsilon, eps_eff))

        # Polarised single-crystal measurements along the (001)
        # should select the xx and yy diagonal components of eps_ir.

        sp_sc_pol_x = calc.single_crystal_input_polarised_optical_spectrum(
            (0, 0, 1), Polarisation.from_direction("x"), t=1.0e-6
        )

        sp_sc_pol_y = calc.single_crystal_input_polarised_optical_spectrum(
            (0, 0, 1), Polarisation.from_direction("y"), t=1.0e-6
        )

        self.assertTrue(np.allclose(sp_sc_pol_x.epsilon, eps_ir[:, 0, 0]))
        self.assertTrue(np.allclose(sp_sc_pol_y.epsilon, eps_ir[:, 1, 1]))

        # An unpolarised spectrum along the (001) should give the
        # average of the x- and y-polarised spectra

        sp_sc_unpol = calc.single_crystal_unpolarised_optical_spectrum(
            (0, 0, 1), t=1.0e-6
        )

        sp_sc_pol_ave = (sp_sc_pol_x.epsilon + sp_sc_pol_y.epsilon) / 2.0

        self.assertTrue(np.allclose(sp_sc_unpol.epsilon, sp_sc_pol_ave))

    def test_optical_spectra_2(self):
        """Test the "synchronisation" of the sample thickness between
        the `EigenmodeAverageOpticalSpectrum` and the underlying
        `OpticalEigenmodeSpectrum`."""

        sp_ema = self._calc.powder_optical_spectrum_ema(t=1.0e-6)

        trans_int_1u = sp_ema.intrinsic_transmission
        trans_norm_1u = sp_ema.normal_transmission
        trans_incoh_1u = sp_ema.incoherent_transmission

        sp_ema.sample_thickness = 1.0e-7

        self.assertEqual(
            sp_ema.sample_thickness,
            sp_ema.optical_eigenmode_spectrum.sample_thickness,
        )

        self.assertFalse(
            np.equal(sp_ema.intrinsic_transmission, trans_int_1u).all()
        )

        self.assertFalse(
            np.equal(sp_ema.normal_transmission, trans_norm_1u).all()
        )

        self.assertFalse(
            np.equal(sp_ema.incoherent_transmission, trans_incoh_1u).all()
        )

        sp_ema.optical_eigenmode_spectrum.sample_thickness = 1.0e-6

        self.assertEqual(
            sp_ema.sample_thickness,
            sp_ema.optical_eigenmode_spectrum.sample_thickness,
        )

        self.assertTrue(
            np.equal(sp_ema.intrinsic_transmission, trans_int_1u).all()
        )

        self.assertTrue(
            np.equal(sp_ema.normal_transmission, trans_norm_1u).all()
        )

        self.assertTrue(
            np.equal(sp_ema.incoherent_transmission, trans_incoh_1u).all()
        )

    def test_bruggeman_multiphase(self):
        """Test the implementation of the Bruggeman multiphase models."""

        calc = self._calc

        eps_kbr = 4.9
        eps_air = 1.0

        sp_ema = calc.powder_optical_spectrum_ema(t=1.0e-6)

        eps_snse = sp_ema.epsilon

        # Pure powder with 90% theoretical denisty.

        eps_eff_ld = bruggeman_two_phase_mixture(eps_snse, eps_air, 0.9, 0.1)

        self.assertTrue(
            verify_bruggeman_equation(
                [eps_snse, eps_air], [0.9, 0.1], eps_eff_ld
            )
        )

        # 5% KBr pellet.

        eps_eff_kbr = bruggeman_two_phase_mixture(
            eps_snse, eps_kbr, 0.05, 0.95
        )

        self.assertTrue(
            verify_bruggeman_equation(
                [eps_snse, eps_kbr], [0.05, 0.95], eps_eff_kbr
            )
        )

        # 5% KBr pellet with 90% denisity.

        eps_eff_kbr_ld = bruggeman_three_phase_mixture(
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

        sp_ema_ld = calc.powder_optical_spectrum_ema(t=1.0e-6, p_den=0.9)
        self.assertTrue(np.allclose(sp_ema_ld.epsilon, eps_eff_ld))

        sp_ema_kbr = calc.powder_optical_spectrum_ema(
            p_vol_frac=0.05, p_binder_eps=eps_kbr
        )

        self.assertTrue(np.allclose(sp_ema_kbr.epsilon, eps_eff_kbr))

        sp_ema_kbr_ld = calc.powder_optical_spectrum_ema(
            p_vol_frac=0.05, p_binder_eps=eps_kbr, p_den=0.9
        )

        self.assertTrue(np.allclose(sp_ema_kbr_ld.epsilon, eps_eff_kbr_ld))

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
            eps_eff = bruggeman_multiphase_mixture(eps, fracs)

            self.assertTrue(verify_bruggeman_equation(eps, fracs, eps_eff))

            # Verify that the general solution matches the analytical
            # formulae.

            self.assertTrue(np.allclose(eps_eff, eps_eff_ref))


# ----
# Main
# ----


if __name__ == "__main__":
    unittest.main()
