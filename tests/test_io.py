# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Test routines for the core I/O."""

# -------
# Imports
# -------

import os
import unittest

import numpy as np

from phonopy_spectroscopy.interfaces.phonopy_interface import (
    gamma_phonons_from_phono3py,
    structure_from_phonopy_yaml,
    gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml,
    gamma_freqs_evecs_from_mesh_qpoints_or_band_hdf5,
)

from phonopy_spectroscopy.interfaces.vasp_interface import (
    structure_from_poscar,
    structure_to_poscar,
    dielectric_from_vasprun_xml,
)

from phonopy_spectroscopy.gamma_phonons import GammaPhonons, PolarGammaPhonons
from phonopy_spectroscopy.structure import Structure

from phonopy_spectroscopy.utility.io_helper import load_json, save_json

from comparison_helper import (
    compare_structures,
    compare_gamma_phonons,
    compare_polar_gamma_phonons,
)

# ---------
# Constants
# ---------

_EXAMPLE_BASE_DIR_SI = r"../example/si"

"""Path to Si example directory."""

_EXAMPLE_BASE_DIR_SNSE = r"../example/snse-pnma"

"""Path to SnSe (Pnma) example directory."""

# -------------
# Tests for I/O
# -------------


class TestIO(unittest.TestCase):
    """Class implementing unit tests for I/O routines."""

    def test_structure_io(self):
        """Test structure input/output routines."""

        # Load reference structures from VASP POSCAR and Phonopy
        # phonopy.yaml files and compare.

        struct_ref_1 = structure_from_poscar(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"POSCAR.Opt.Prim")
        )

        struct_ref_2 = structure_from_phonopy_yaml(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"struct_ref/phonopy.yaml")
        )

        self.assertTrue(compare_structures(struct_ref_2, struct_ref_1))

        # Write a structure to a VASP POSCAR file, read it back, and
        # ensure the two Structure objects are equivalent.

        structure_to_poscar(struct_ref_1, r"POSCAR.vasp.tmp")
        struct_cmp = structure_from_poscar(r"POSCAR.vasp.tmp")

        self.assertTrue(compare_structures(struct_cmp, struct_ref_1))

        os.remove(r"POSCAR.vasp.tmp")

        # Serialise the structure to a dictionary, write it to a JSON
        # file, reload and recreate it, and check the Structure objects
        # are equivalent.

        save_json(struct_ref_1.to_dict(), r"structure.json.tmp")
        struct_cmp = Structure.from_dict(load_json(r"structure.json.tmp"))

        self.assertTrue(compare_structures(struct_cmp, struct_ref_1))

        os.remove(r"structure.json.tmp")

    def test_freqs_evecs_io(self):
        """Test routines for reading phonon frequencies and eigenvectors
        from Phonopy calculations."""

        # Load frequencies and eigenvectors from mesh/band YAML and HDF5
        # files and check the data are eqivalent.

        freqs_evecs_1 = gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"freqs_evecs_ref/mesh.yaml")
        )

        freqs_evecs_2 = gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"freqs_evecs_ref/qpoints.yaml")
        )

        freqs_evecs_3 = gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"freqs_evecs_ref/band.yaml")
        )

        freqs_evecs_4 = gamma_freqs_evecs_from_mesh_qpoints_or_band_hdf5(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"freqs_evecs_ref/mesh.hdf5")
        )

        freqs_evecs_5 = gamma_freqs_evecs_from_mesh_qpoints_or_band_hdf5(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"freqs_evecs_ref/qpoints.hdf5")
        )

        freqs_evecs_6 = gamma_freqs_evecs_from_mesh_qpoints_or_band_hdf5(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"freqs_evecs_ref/band.hdf5")
        )

        freqs_ref, evecs_ref = freqs_evecs_1

        for freqs_cmp, evecs_cmp in (
            freqs_evecs_2,
            freqs_evecs_3,
            freqs_evecs_4,
            freqs_evecs_5,
            freqs_evecs_6,
        ):
            self.assertTrue(np.allclose(freqs_cmp, freqs_ref))
            self.assertTrue(np.allclose(evecs_cmp, evecs_ref))

    def test_gamma_phonons_io_1(self):
        """Test high-level Phono(3)py "loader" and
        serialisation/deserialisation of `GammaPhonons` objects."""

        # Test the construction of a complete GammaPhonons object
        # including a structure, frequencies/eigenvectors, linewidths
        # and irreps.

        gamma_ph = gamma_phonons_from_phono3py(
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"POSCAR.Opt.Prim"),
            os.path.join(_EXAMPLE_BASE_DIR_SI, r"mesh.yaml"),
            lws_file=os.path.join(
                _EXAMPLE_BASE_DIR_SI, r"kappa-m646464-g0.hdf5"
            ),
            irreps_file=os.path.join(_EXAMPLE_BASE_DIR_SI, r"irreps.yaml"),
        )

        # The Gamma-point eigenvectors should be real and should be
        # stored as np.float64.

        self.assertEqual(gamma_ph.eigenvectors.dtype, np.float64)

        # Serialise the GammaPhonons to a dictionary, write it to a JSON
        # file, reload and recreate it, and check the two objects
        # are equivalent.

        save_json(gamma_ph.to_dict(), r"gamma_phonons.json.tmp")

        gamma_ph_cmp = GammaPhonons.from_dict(
            load_json(r"gamma_phonons.json.tmp")
        )

        self.assertTrue(compare_gamma_phonons(gamma_ph_cmp, gamma_ph))

        os.remove(r"gamma_phonons.json.tmp")

    def test_gamma_phonons_io_2(self):
        """Test high-level Phono(3)py "loader" and
        serialisation/deserialisation of `PolarGammaPhonons` objects."""

        cell_file = os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"POSCAR.Opt")
        freqs_evecs_file = os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"mesh.yaml")

        lws_file = os.path.join(
            _EXAMPLE_BASE_DIR_SNSE, r"kappa-m323216-g0.hdf5"
        )

        irreps_file = os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"irreps.yaml")
        born_file = os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"BORN")

        gamma_ph_1 = gamma_phonons_from_phono3py(
            cell_file,
            freqs_evecs_file,
            lws_file=lws_file,
            irreps_file=irreps_file,
        )

        gamma_ph_2 = gamma_phonons_from_phono3py(
            cell_file,
            freqs_evecs_file,
            lws_file=lws_file,
            irreps_file=irreps_file,
            born_file=born_file,
        )

        # Check gamma_phonons_from_phono3py() returns a
        # PolarGammaPhonons object when passed the optional born_file
        # argument, and a GammaPhonons object otherwise.

        self.assertTrue(isinstance(gamma_ph_1, GammaPhonons))
        self.assertFalse(isinstance(gamma_ph_1, PolarGammaPhonons))

        self.assertTrue(isinstance(gamma_ph_2, GammaPhonons))
        self.assertTrue(isinstance(gamma_ph_2, PolarGammaPhonons))

        # Check the PolarGammaPhonons object has the same "base" data
        # as the GammaPhonons initialised from the same files.

        self.assertTrue(compare_gamma_phonons(gamma_ph_2, gamma_ph_1))

        # Test the serialisation/deserialisation of the
        # PolarGammaPhonons object.

        save_json(gamma_ph_2.to_dict(), r"polar_gamma_phonons.json.tmp")

        gamma_ph_2_cmp = PolarGammaPhonons.from_dict(
            load_json(r"polar_gamma_phonons.json.tmp")
        )

        self.assertTrue(
            compare_polar_gamma_phonons(gamma_ph_2_cmp, gamma_ph_2)
        )

        os.remove(r"polar_gamma_phonons.json.tmp")

    def test_gamma_phonons_io_3(self):
        """Test high-level Phono(3)py "loader" and
        serialisation/deserialisation of `GammaPhonons` objects with
        complex eigenvectors."""

        cell_file = os.path.join(_EXAMPLE_BASE_DIR_SNSE, r"POSCAR.Opt")

        # Reference calculation with non-analytical correction has
        # complex eigenvectors.

        freqs_evecs_file = os.path.join(
            _EXAMPLE_BASE_DIR_SNSE, r"nac_ref/qpoints-001.yaml"
        )

        # Test the option to discard the imaginary parts of complex
        # eigenvectors from the gamma_phonons_from_phono3py() interface.

        gamma_ph = gamma_phonons_from_phono3py(
            cell_file,
            freqs_evecs_file,
            discard_imag=False,
        )

        # Loading complex eigenvectors with discard_imag=True should
        # result in a UserWarning.

        with self.assertRaises(RuntimeError):
            gamma_phonons_from_phono3py(
                cell_file,
                freqs_evecs_file,
                discard_imag=True,
            )

        self.assertTrue(np.iscomplexobj(gamma_ph.eigenvectors))

        # Test the serialisation/deserialisation of the GammaPhonons
        # object with complex eigenvectors.

        save_json(gamma_ph.to_dict(), r"gamma_phonons.json.tmp")

        gamma_ph_cmp = GammaPhonons.from_dict(
            load_json(r"gamma_phonons.json.tmp")
        )

        self.assertTrue(compare_gamma_phonons(gamma_ph_cmp, gamma_ph))

        os.remove(r"gamma_phonons.json.tmp")

    def test_dielectric_io(self):
        """Test routines for reading dielectric data from vasprun.xml
        files."""

        input_files = [
            r"vasprun-PBEsol-DFPT.xml",
            r"vasprun-PBEsol-FiniteField.xml",
            r"vasprun-PBEsol-LinearOptics.xml",
            r"vasprun-r2SCAN-FiniteField.xml",
            r"vasprun-r2SCAN-LinearOptics.xml",
            r"vasprun-mBJ-LinearOptics.xml",
            r"vasprun-HSE06-FiniteField.xml",
            r"vasprun-HSE06-LinearOptics.xml",
        ]

        ref_eps_hf = [
            13.49657154,
            13.27554306,
            14.2169,
            11.78495469,
            11.8552,
            10.3019,
            11.09544159,
            10.3564,
        ]

        for f, ref_eps in zip(input_files, ref_eps_hf):
            e, eps_e = dielectric_from_vasprun_xml(
                os.path.join(_EXAMPLE_BASE_DIR_SI, r"raman_ref", f)
            )

            # First energy should be E = 0.

            self.assertTrue(np.isclose(e[0], 0.0))

            # Dielectric constant at E = 0 should be real.

            self.assertFalse(np.iscomplex(eps_e[0]).any())

            # Si is isotropic so the three diagonal elements should be
            # equal by symmetry.

            self.assertTrue((np.diag(eps_e[0]) == ref_eps).all())


# ----
# Main
# ----

if __name__ == "__main__":
    unittest.main()
