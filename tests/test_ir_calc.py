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
)

from phonopy_spectroscopy.ir.calculation import InfraredCalculation
from phonopy_spectroscopy.utility.io_helper import load_json, save_json

from comparison_helper import compare_infrared_calculations

# ---------
# Constants
# ---------

_EXAMPLE_BASE_DIR_SNSE = r"../example/snse-pnma"

"""Path to SnSe (Pnma) example directory."""


# ---------------------------------------
# Tests for infrared calculation workflow
# ---------------------------------------


class TestInfraredCalculationWorkflow(unittest.TestCase):
    """Class implementing unit tests for the infrared calculation
    workflow."""

    def test_ir_calc(self):
        """Test the construction and serialisation/deserialisation
        of the `InfraredCalculation` class."""

        # Load a Gamma-point phonon calculation using the high-level
        # Phono(3)py "loader" function.

        gamma_ph = gamma_phonons_from_phono3py(
            os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"POSCAR.Opt"),
            os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"mesh.yaml"),
            lws_file=os.path.join(
                _EXAMPLE_BASE_DIR_SNSE, r"kappa-m323216-g0.hdf5"
            ),
            irreps_file=os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"irreps.yaml"),
            born_file=os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"BORN"),
        )

        # Construct an InfraredCalculator object.

        calc = InfraredCalculation(gamma_ph)

        # Test serialisation/deserialisation of the InfraredCalculator
        # class.

        save_json(calc.to_dict(), r"ir_calculator.json.tmp")

        calc_cmp = InfraredCalculation.from_dict(
            load_json(r"ir_calculator.json.tmp")
        )

        os.remove(r"ir_calculator.json.tmp")

        self.assertTrue(compare_infrared_calculations(calc, calc_cmp))


# ----
# Main
# ----

if __name__ == "__main__":
    unittest.main()
