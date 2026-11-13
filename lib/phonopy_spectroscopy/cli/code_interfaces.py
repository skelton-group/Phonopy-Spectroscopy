# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Interfaces for reading code-specific input files."""

# -------
# Imports
# -------

import os

from ..interfaces.phonopy_interface import (
    structure_from_phonopy_yaml,
    gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml,
    irreps_from_irreps_yaml,
    gamma_freqs_evecs_from_mesh_qpoints_or_band_hdf5,
    gamma_linewidths_from_kappa_hdf5,
    eps_inf_and_born_from_born,
)

from ..interfaces.vasp_interface import (
    structure_from_poscar,
    structure_to_poscar,
)

from .io.raman import fd_read_dielectrics_vasp

# ---------
# Constants
# ---------

_SUPPORTED_INTERFACES = ["vasp"]

"""List of supported code interfaces."""

# ---------
# Functions
# ---------


def get_supported_code_interfaces():
    """Get a list of supported code interfaces.

    Returns
    -------
    interfaces : list of str
        List of interfaces.
    """

    return list(_SUPPORTED_INTERFACES)


def get_code_interface(code):
    """Get a loader object for interfacing with a supported code.

    Parameters
    ----------
    code : str
        Code to get loader object for.

    Returns
    -------
    loader : LoaderBase
        Loader object for `code`.
    """

    code = code.lower()

    match code:
        case "vasp":
            return VASPLoader()
        case _:
            if code in _SUPPORTED_INTERFACES:
                raise RuntimeError("")

            raise RuntimeError("")


# ----------------
# LoaderBase class
# ----------------


class LoaderBase:
    """Base class for code-specific Loader objects.

    Notes
    -----
    "Loader" objects are used by the CLI to read input files.

    This class defines standard methods for reading from Phonopy
    input/output files, which dervied classes override to support
    reading from other input files (e.g. code-specific structure
    formats) with the default Phonopy files as a fallback.

    The class also exposes properties indicating support for Raman
    calculations, and defined methods for writing displaced structures
    and reading back data from calculations on these structures. These
    are skeleton implementations that raise a `NotImplementedError` and
    must be overridden and implemented in derived classes.
    """

    def __init__(self):
        """Create a new instance of the `LoaderBase` class."""
        pass

    def try_load_structure(self, file_path=None, **kwargs):
        """Read a structure from a phonopy.yaml file and return a
        `Structure` object.

        Parameters
        ----------
        file_path : str, optional
            File path (default: None -> "phonopy.yaml").
        **kwargs : any, optional
            Keyword arguments to `structure_from_phonopy_yaml()`.

        Returns
        -------
        res : tuple or None
            Tuple of `(struct, file_path)` with a `Structure` object if
            successfully read, otherwise `None`.

        See Also
        --------
        interfaces.phonopy_interface.structure_from_phonopy_yaml :
            Read a structure from a phonopy.yaml file.
        """

        if file_path is None and os.path.isfile(r"phonopy.yaml"):
            file_path = r"phonopy.yaml"

        if file_path is not None:
            return (
                structure_from_phonopy_yaml(file_path=file_path, **kwargs),
                file_path,
            )

        return None

    def try_load_gamma_freqs_evecs(self, file_path=None):
        r"""Read \Gamma-point frequencies and eigenvectors from a
        Phonopy mesh.hdf5/mesh.yaml or band.hdf5/band.yaml file.

        Parameters
        ----------
        file_path : str, optional
            File path (default: None -> "mesh.hdf5", "mesh.yaml",
            "band.hdf5", "band.yaml").

        Returns
        -------
        res : tuple or None
            Tuple of `((freqs, evecs), file_path)` if successfully read,
            otherwise `None` (shapes: `(3N,)`, `(3N, N, 3)`).

        See Also
        --------
        interfaces.phonopy_interface.gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml :
            Read \Gamma-point frequencies and eigenvectors from a
            mesh.yaml or band.yaml file.
        interfaces.phonopy_interface.gamma_freqs_evecs_from_mesh_qpoints_or_band_hdf5 :
            Read \Gamma-point frequencies and eigenvectors from a
            mesh.hdf5 or band.hdf5 file.
        """

        if file_path is None:
            for f in ["mesh.hdf5", "mesh.yaml", "band.hdf5", "band.yaml"]:
                if os.path.isfile(f):
                    file_path = f
                    break

        if file_path is not None:
            _, ext = os.path.splitext(file_path)

            match ext:
                case ".hdf5":
                    return (
                        gamma_freqs_evecs_from_mesh_qpoints_or_band_hdf5(
                            file_path
                        ),
                        file_path,
                    )
                case ".yaml":
                    return (
                        gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml(
                            file_path
                        ),
                        file_path,
                    )

        return None

    def try_load_linewidths(self, file_path=None, t=300.0):
        r"""Read \Gamma-point linewidths at a specified temperature
        from a Phonopy kappa-m*.hdf5 file.

        Parameters
        ----------
        file_path : str, optional
            File path (default: None)
        t : float, optional
            Temperature (default: 300.0)

        Returns
        -------
        res : tuple or None
            Tuple of `(lws, file_path)` if `file_path` is set, otherwise
            otherwise `None` (shape: `(3N,)`).

        See Also
        --------
        interfaces.phonopy_interface.gamma_linewidths_from_kappa_hdf5 :
            Read \Gamma-point linewidths from a kappa-m*.hdf5 file.
        """

        if file_path is not None:
            return (
                gamma_linewidths_from_kappa_hdf5(file_path, t=t),
                file_path,
            )

        return None

    def try_load_irreps(self, file_path=None):
        r"""Read irreps from a Phonopy irreps.yaml file.

        Parameters
        ----------
        file_path : str, optional
            File path (default: None -> "irreps.yaml").

        Returns
        -------
        res : tuple or None
            Tuple of `((irreps, born), file_path)` with an `Irreps`
            object if successful, otherwise `None`.

        See Also
        --------
        interfaces.phonopy_interface.irreps_from_irreps_yaml :
            Read irreps from an irreps.yaml file.
        """

        if file_path is None and os.path.isfile(r"irreps.yaml"):
            file_path = r"irreps.yaml"

        if file_path is not None:
            return (irreps_from_irreps_yaml(file_path), file_path)

        return None

    def try_load_eps_inf_and_born(self, struct, file_path=None):
        r"""Read high-frequency dielectric constant \eps_inf and Born
        effective charges Z* from a Phonopy BORN file.

        Parameters
        ----------
        struct : Structure
            Crystal structure.
        file_path : str, optional
            File path (default: None -> "BORN").

        Returns
        -------
        res : tuple or None
            Tuple of `((eps_inf, born), file_path)` if successful,
            otherwise `None` (shapes: `(3, 3)`, `(N, 3, 3)`).

        See Also
        --------
        interfaces.phonopy_interface.eps_inf_and_born_from_born :
            Read high-frequency dielectric constant and Born charges
            from a BORN file.
        """

        if file_path is None and os.path.isfile(r"BORN"):
            file_path = r"BORN"

        if file_path is not None:
            return (eps_inf_and_born_from_born(file_path, struct), file_path)

        return None

    @property
    def SupportsFDRamanCalculations(self):
        """bool : `True` if the loader supports finite-displacement
        calculations, otherwise `False`."""

        return False

    def save_fd_raman_structures(self, structs, band_inds, num_disps):
        """Save a set of `Structure` objects generated by performing M
        displacements along Raman-active N modes to code-specific
        structure files.

        Parameters
        ----------
        structs : numpy.ndarray
            Structures (shape: `(N, M)`).
        band_inds : numpy.ndarray
            Band indices (shape: `(N,)`).
        num_disps : int
            Number of displacements (M).
        """

        raise NotImplementedError(
            "The base loader does not support Raman calculations."
        )

    def read_fd_raman_dielectrics(self, input_files, band_inds, num_disps):
        """Read high-frequency dielectric constants or energy-dependent
        dielectric functions calculated for the structures generated by
        `save_fd_raman_structures()` from a set of input files.

        Parameters
        ----------
        input_files : list of str
            Input files (length: N * M).
        band_inds : numpy.ndarray
            Band indices (shape: `(N,)`).
        num_disps : int
            Number of displacements (M).

        Returns
        -------
        dielectrics : numpy.ndarray or tuple of numpy.ndarray
            High-frequency dielectric constants (shape: `(N, M, 3, 3)`)
            or `(e, eps_e)` tuple of energy-dependent dielectric
            functions (shape: `(O,)` and `(N, M, O, 3, 3)`).
        """

        raise NotImplementedError(
            "The base loader does not support Raman calculations."
        )


# ----------------
# VASPLoader class
# ----------------


class VASPLoader(LoaderBase):
    pass

    def __init__(self):
        """Create a new instance of the VASPLoader class."""
        pass

    def try_load_structure(self, file_path=None, **kwargs):
        """Override of base method to support reading structures from
        VASP POSCAR files."""

        if file_path is None and os.path.isfile(r"POSCAR"):
            file_path = r"POSCAR"

        if file_path is not None:
            _, ext = os.path.splitext(file_path)

            if ext.lower() != ".yaml":
                # Assume a VASP POSCAR file.

                return (structure_from_poscar(file_path, **kwargs), file_path)

        return super().try_load_structure(file_path, **kwargs)

    @property
    def SupportsFDRamanCalculations(self):
        return True

    def save_fd_raman_structures(self, structs, band_inds, num_disps):
        """Save a set of `Structure` objects with finite displacements
        along the Raman-active normal modes to VASP POSCAR files."""

        for i, idx in enumerate(band_inds):
            for j in range(num_disps):
                structure_to_poscar(
                    structs[i, j],
                    "POSCAR-{0:0>4}-{1:0>2}".format(idx + 1, j + 1),
                )

    def read_fd_raman_dielectrics(self, input_files, band_inds, num_disps):
        """Read high-frequency calculated dielectric constants or
        energy-dependent dielectric functions from a set of vasprun.xml
        files."""

        # Default energy tolerance of 0.1 eV for matching
        # energy-dependent dielectric functions.

        return fd_read_dielectrics_vasp(
            input_files, len(band_inds), num_disps, e_tol=0.1
        )
