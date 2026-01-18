# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Test routines for the infrared calculation workflow."""


# -------
# Imports
# -------

import os
import unittest

from phonopy_spectroscopy.interfaces.phonopy_interface import (
    gamma_phonons_from_phono3py,
    hf_dielectric_and_born_from_born,
)

from phonopy_spectroscopy.ir.calculation import InfraredCalculation

from phonopy_spectroscopy.utility.io_helper import load_json, save_json

from comparison_helper import compare_infrared_calculators


# ---------
# Constants
# ---------


_EXAMPLE_BASE_DIR = r"../example/snse-pnma"


# -------------------------------
# Tests for infrared calculations
# -------------------------------


class TestInfraredCalculations(unittest.TestCase):
    def test_ir_calc(self):
        # Load a Gamma-point phonon calculation using the high-level
        # Phono(3)py "loader" function.

        lws_file = os.path.join(
            _EXAMPLE_BASE_DIR, r"kappa-m16816.Prim.FullPP.hdf5"
        )

        gamma_ph = gamma_phonons_from_phono3py(
            os.path.join(_EXAMPLE_BASE_DIR, r"POSCAR.Opt"),
            os.path.join(_EXAMPLE_BASE_DIR, r"mesh.yaml"),
            lws_file=lws_file,
            irreps_file=os.path.join(_EXAMPLE_BASE_DIR, r"irreps.yaml"),
        )

        # Load high-frequency dielectric constant (eps_inf) and Born
        # charges from a Phonopy BORN file.

        eps_inf, born_charges = hf_dielectric_and_born_from_born(
            os.path.join(_EXAMPLE_BASE_DIR, r"BORN"),
            gamma_ph.structure,
        )

        # Construct an InfraredCalculator object.

        calc = InfraredCalculation(gamma_ph, born_charges, eps_inf=eps_inf)

        # Test serialisation/deserialisation of the InfraredCalculator
        # class.

        save_json(calc.to_dict(), r"ir_calculator.json.tmp")

        calc_cmp = InfraredCalculation.from_dict(
            load_json(r"ir_calculator.json.tmp")
        )

        os.remove(r"ir_calculator.json.tmp")

        self.assertTrue(compare_infrared_calculators(calc, calc_cmp))


# ----
# Main
# ----


if __name__ == "__main__":
    unittest.main()
