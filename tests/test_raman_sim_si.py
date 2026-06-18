# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Test routines for Raman simulations with a complete calculation on
Si.

Some of the results are tested against simulations with the model Raman
tensors from:

    U. Ramabadran and B. Roughani, Mater. Sci. Eng. B 230, 31 (2018),
    DOI: 10.1016/j.mseb.2017.12.040

(The unit tests in `test_raman_sim_si_model.py` verify single-crystal
and powder results obtained with these tensors.)
"""

# -------
# Imports
# -------

import os
import unittest

import numpy as np

from phonopy_spectroscopy.cli.utility.raman_io import (
    fd_read_dielectrics_vasp,
)

from phonopy_spectroscopy.instrument import Geometry, Polarisation

from phonopy_spectroscopy.interfaces.phonopy_interface import (
    gamma_phonons_from_phono3py,
)

from phonopy_spectroscopy.raman.calculators import (
    FiniteDisplacementRamanTensorCalculator,
)

from phonopy_spectroscopy.raman.intensity import (
    calculate_single_crystal_raman_intensities,
    calculate_powder_raman_intensities,
)

from phonopy_spectroscopy.structure import Structure

from phonopy_spectroscopy.utility.geometry import (
    rotation_matrix_from_vectors,
    rotation_matrix_from_axis_angle,
)

from io_helper import generate_fd_raman_dielectric_input_file_list

# ---------
# Constants
# ---------

_EXAMPLE_BASE_DIR_SI = r"../example/si"

"""Path to Si example directory."""


# ---------------------------
# Tests for Raman simulations
# ---------------------------


class TestRamanSimulation(unittest.TestCase):
    """Class implementing unit tests for Raman simulations."""

    def setUp(self):
        """Perform setup."""

        # Load a Gamma-point phonon calculation using the high-level
        # Phono(3)py "loader" function.

        conv_trans = np.linalg.inv(
            [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        )

        gamma_ph = gamma_phonons_from_phono3py(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"POSCAR.Opt.Prim"),
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"mesh.yaml"),
            lws_file=os.path.join(
                _EXAMPLE_BASE_DIR_SI, r"kappa-m646464-g0.hdf5"
            ),
            irreps_file=os.path.join(_EXAMPLE_BASE_DIR_SI, r"irreps.yaml"),
            conv_trans=conv_trans,
        )

        # Construct a FiniteDisplacementRamanTensorCalculator object
        # with the default parameters.

        fd_calc = FiniteDisplacementRamanTensorCalculator(gamma_ph)

        # Load energy-dependent dielectric functions for displaced
        # structutes.

        file_path_template = os.path.join(
            _EXAMPLE_BASE_DIR_SI,
            r"raman_ref",
            r"vasprun-HSE06-LinearOptics-350eV-{0:0>4}-{1:0>2}.xml",
        )

        file_list = generate_fd_raman_dielectric_input_file_list(
            fd_calc.band_indices, fd_calc.num_steps, file_path_template
        )

        e, eps_e = fd_read_dielectrics_vasp(
            file_list, fd_calc.num_bands, fd_calc.num_steps
        )

        # Obtain a RamanCalculation object for performing test Raman
        # simulations.

        self._calc = fd_calc.get_raman_calculation(eps_e, e)

        # Define an "empty" structure and set of model Raman tensors.

        a = 5.431020511

        self._model_struct = Structure(
            [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]], [], [], []
        )

        b = 1.0

        r_x = [[0.0, 0.0, 0.0], [0.0, 0.0, b], [0.0, b, 0.0]]
        r_y = [[0.0, 0.0, b], [0.0, 0.0, 0.0], [b, 0.0, 0.0]]
        r_z = [[0.0, b, 0.0], [b, 0.0, 0.0], [0.0, 0.0, 0.0]]

        self._model_r_t = np.array([r_x, r_y, r_z], dtype=np.float64)

        # Define a measurement geometry.

        self._geom = Geometry("z", "-z")

    def test_single_crystal_xtal_rot(self):
        """Test a simulation of a single-crystal crystal rotation Raman
        experiment."""

        pol = Polarisation.from_direction("x")

        # Perform simulation with calculated Raman tensors.

        sp = self._calc.single_crystal_crystal_rotation(
            (0, 0, 1),
            geom=self._geom,
            i_pol=pol,
            s_pol=pol,
            phi_start=0.0,
            phi_end=360.0,
            phi_step=2.5,
        )

        # Perform simulation "manually" with model Raman tensors.

        r_hkl = rotation_matrix_from_vectors(
            self._model_struct.real_space_normal((0, 0, 1)),
            -1.0 * self._geom.incident_direction,
        )

        angles = np.linspace(0.0, 360.0, sp.intensities.shape[1])

        ints = np.zeros(
            (self._model_r_t.shape[0], len(angles)), dtype=np.float64
        )

        for i, a in enumerate(angles):
            r = rotation_matrix_from_axis_angle(
                self._geom.incident_direction, a
            )

            ints[:, i] = calculate_single_crystal_raman_intensities(
                self._model_r_t, self._geom, pol, pol, rot=np.dot(r, r_hkl)
            )

        ints = ints.sum(axis=0)

        # Normalise both sets of intensities for comparison.

        ints_ref = ints / ints.max()
        ints_cmp = sp.intensities[0, :] / sp.intensities[0, :].max()

        # An absolute tolerance of 1.0e-2 requires that the two sets of
        # relative intensities are within 1% of each other.

        self.assertTrue(np.allclose(ints_cmp, ints_ref, atol=1.0e-2))

    def test_powder_pol_rot(self):
        """Test a simulation of a polarisation rotation powder Raman
        experiment."""

        i_pol = Polarisation.from_direction("x")

        sp = self._calc.powder_polarisation_rotation(
            geom=self._geom,
            i_pol=i_pol,
            s_pol="rot",
            chi_start=0.0,
            chi_end=360.0,
            chi_step=2.5,
            e_rt=0.0,
        )

        s_pols = Polarisation.from_rotation(
            self._geom.collection_direction, start=0.0, end=360.0, step=2.5
        )

        ints = np.zeros(
            (self._model_r_t.shape[0], len(s_pols)), dtype=np.float64
        )

        for i, s_pol in enumerate(s_pols):
            ints[:, i] = calculate_powder_raman_intensities(
                self._model_r_t, self._geom, i_pol, s_pol
            )

        ints = ints.sum(axis=0)

        ints_ref = ints / ints.max()
        ints_cmp = sp.intensities[0, :] / sp.intensities[0, :].max()

        # Powder averaging effectively removes the variation between the
        # model and calculated Raman tensors, so the relative
        # intensities can be compared with the standard settings.

        self.assertTrue(np.allclose(ints_cmp, ints_ref))


# ----
# Main
# ----

if __name__ == "__main__":
    unittest.main()
