# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for interfacing with the Phono(3)py code."""

# -------
# Imports
# -------

import os
import warnings

import h5py

import numpy as np

from ..constants import ZERO_TOLERANCE
from ..irreps import Irreps
from ..gamma_phonons import get_gamma_phonons
from ..structure import Structure
from ..utility.io_helper import load_yaml

from .vasp_interface import structure_from_poscar

_PHONOPY_AVAILABLE = False

try:
    from phonopy.file_IO import parse_BORN

    _PHONOPY_AVAILABLE = True
except ImportError:
    warnings.warn(
        "Imports from Phonopy failed - some functions require Phonopy "
        "and will raise exceptions if it is not installed.",
        RuntimeWarning,
    )

# -----------------------
# Code interface handling
# -----------------------

_INTERFACE_UNITS = {
    "abacus": {"distance": "bohr"},
    "abinit": {"distance": "bohr"},
    "aims": {"distance": "ang"},
    "castep": {"distance": "ang"},
    "cp2k": {"distance": "ang"},
    "crystal": {"distance": "ang"},
    "dftbp": {"distance": "bohr"},
    "elk": {"distance": "bohr"},
    "fleur": {"distance": "bohr"},
    "lammps": {"distance": "ang"},
    "qlm": {"distance": "bohr"},
    "qe": {"distance": "bohr"},
    "siesta": {"distance": "bohr"},
    "turbomole": {"distance": "bohr"},
    "vasp": {"distance": "ang"},
    "wien2k": {"distance": "bohr"},
    "pwmat": {"distance": "ang"},
}

"""Unit systems used in codes supported by Phonopy."""


def get_distance_unit_for_interface(calculator):
    """Return the distance units used by a named calculator.

    Parameters
    ----------
    calculator : str
        Calculator name.

    Returns
    -------
    units : str
        Distance units.
    """

    k = calculator.lower()

    if k in _INTERFACE_UNITS:
        return _INTERFACE_UNITS[k]["distance"]

    raise ValueError(f'Unknown interface "{calculator}".')


# -------------------
# High-level "loader"
# -------------------


def gamma_phonons_from_phono3py(
    struct,
    freqs_evecs_file,
    lws_file=None,
    lws_t=300.0,
    irreps_file=None,
    born_file=None,
    allow_imag_evecs=False,
    **kwargs,
):
    r"""Read a complete Phono(3)py calculation and return a
    `GammaPhonons` or `PolarGammaPhonons` object.

    Parameters
    ----------
    struct : Structure or str
        `Structure` object or path to a POSCCAR or phonopy.yaml file to
        read a structure from.
    freqs_evecs_file : str
        mesh.yaml, band.yaml, mesh.hdf5 or band.hdf5 file to read phonon
        frequencies/eigenvectors from.
    lws_file : str, optional
        kappa-m*.hdf5 file to read phonon linewidths from (default:
        `None`).
    lws_t : float, optional
        Temperature to read linewidths at (default: 300 K)
    irreps_file : str, optional
        irreps.yaml file to read irreps from (default: `None`).
    born_file : str, optional
        BORN file to read high-frequency dielectric constant \eps_inf
        and Born effective charges from (default: `None`).
    allow_imag_evecs : bool, optional
        If `True`, allow imaginary eigenvectors (default: `False`).
    **kwargs : any, optional
        Keyword arguments for reading structure (ignored if `struct` is
        a `Structure` object).

    Returns
    -------
    gamma_ph : `GammaPhonons`
        `GammaPhonons` object containing the calculation data.

    Notes
    -----
    The phonon eigenvectors at q = \Gamma should be real. If this is not
    the case (to within a tight numerical tolerance), a `RuntimeError`
    is raised. Setting `allow_imag_evecs` downgrades this to a
    `RuntimeWarning`.
    """

    # If struct is not a Structure object, assume it is a file path. If
    # if has a .yaml extension, assume it is a phonopy.yaml file;
    # otherwise, assume it is a VASP POSCAR file.

    if not isinstance(struct, Structure):
        _, ext = os.path.splitext(struct)

        if ext.lower() == ".yaml":
            struct = structure_from_phonopy_yaml(struct, **kwargs)
        else:
            struct = structure_from_poscar(struct, **kwargs)

    # Read a set of frequencies/eigenvectors. If freqs_evecs_file has a
    # .yaml extension, assume it is a mesh.yaml/band.yaml file. If the
    # file has a .hdf5 extension, assume it is a mesh.hdf5/band.hdf5
    # file. Otherwise, raise an error and defer responsibility to the
    # calling code.

    freqs, evecs = None, None

    _, ext = os.path.splitext(freqs_evecs_file)

    if ext.lower() == ".yaml":
        freqs, evecs = gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml(
            freqs_evecs_file
        )
    elif ext.lower() == ".hdf5":
        freqs, evecs = gamma_freqs_evecs_from_mesh_qpoints_or_band_hdf5(
            freqs_evecs_file
        )
    else:
        raise RuntimeError(
            f"Frequencies/eigenvectors file {freqs_evecs_file}: "
            "unknown format."
        )

    # The functions for reading eigenvectors from YAML/HDF5 files both
    # strip small imaginary parts -> a complex dtype indicates
    # non-negligible imaginary parts.

    if np.iscomplexobj(evecs):
        # Gamma-point eigenvectors should be real unless a
        # non-analytical correction is applied to the dynamical
        # matrix. If this is the case, it will conflict with the
        # internal NAC implementation and is very likely a user
        # error.

        max_abs_imag = np.abs(evecs.imag).max()

        msg = (
            "Discarding imaginary parts of eigenvectors "
            f"with maximum absolute value {max_abs_imag:.3e} > "
            f"ZERO_TOLERANCE = {ZERO_TOLERANCE:.3e}. This may "
            "indicate a calculation performed with a "
            "non-analytical correction to the dynamical matrix."
        )

        if allow_imag_evecs:
            warnings.warn(msg, RuntimeWarning)
        else:
            raise RuntimeError(msg)

    # If a lws_file is specified, read linewidths; otherwise, set a
    # uniform linewidth of lw.

    lws = None

    if lws_file is not None:
        lws = gamma_linewidths_from_kappa_hdf5(lws_file, lws_t)

    # If irreps_file is specified, read irrep data.

    irreps = None

    if irreps_file is not None:
        irreps = irreps_from_irreps_yaml(irreps_file)

    t = lws_t if lws is not None else None

    eps_inf, born_charges = None, None

    if born_file is not None:
        # If born_file is specified, read high-frequency dielectric
        # constant and Born charges.

        eps_inf, born_charges = eps_inf_and_born_from_born(born_file, struct)

    # Automatically return "best" GammaPhonons object depending on data.

    return get_gamma_phonons(
        struct,
        freqs,
        evecs,
        eps_inf=eps_inf,
        born_charges=born_charges,
        lws=lws,
        irreps=irreps,
        t=t,
    )


# ----------
# YAML files
# ----------


def structure_from_phonopy_yaml(file_path, conv_trans=None):
    """Read a structure from a phonopy.yaml file and return a
    `Structure` object.

    Parameters
    ----------
    file_path : str
        Input file.
    conv_trans : array_like, optional
        Transformation matrix to convert the structure to its
        conventional cell (shape: `(3, 3)`, default: `None`).

    Returns
    -------
    struct : Structure
        `Structure` object containing the structure.
    """

    data = load_yaml(file_path)

    cell = data["primitive_cell"]

    return Structure(
        cell["lattice"],
        [atom["coordinates"] for atom in cell["points"]],
        [atom["symbol"] for atom in cell["points"]],
        at_m=[atom["mass"] for atom in cell["points"]],
        conv_trans=conv_trans,
    )


def gamma_freqs_evecs_from_mesh_qpoints_or_band_yaml(file_path):
    """Read Gamma-point phonon frequencies and eigenvectors
    from a mesh.yaml or band.yaml file.

    Parameters
    ----------
    file_path : str
        Input file.

    Returns
    -------
    freqs_evecs : tuple of numpy.ndarray
        A `(freqs, evecs)` tuple of arrays with shapes `(3N,)` and
        `(3N, N, 3)`.
    """

    data = load_yaml(file_path)

    # Get Gamma-point frequencies and eigenvectors.

    for qpt in data["phonon"]:
        if np.allclose(qpt["q-position"], 0.0, atol=ZERO_TOLERANCE):
            # q = (0, 0, 0) = \Gamma.

            if "eigenvector" not in qpt["band"][0]:
                raise RuntimeError(
                    f"mesh.yaml/band.yaml file {file_path}: "
                    "Eigenvectors not found."
                )

            freqs = np.array(
                [mode["frequency"] for mode in qpt["band"]], dtype=np.float64
            )

            evecs = np.array(
                [mode["eigenvector"] for mode in qpt["band"]], dtype=np.float64
            )

            evecs = evecs[:, :, :, 0] + 1.0j * evecs[:, :, :, 1]
            evecs = evecs.reshape(-1, len(freqs) // 3, 3)

            evecs = np.real_if_close(evecs, tol=ZERO_TOLERANCE)

            return (freqs, evecs)

    raise RuntimeError(
        f"mesh.yaml/band.yaml file {file_path}: Gamma-point "
        "frequencies/eigenvectors not found."
    )


def irreps_from_irreps_yaml(file_path):
    """Read Gamma-point mode irreducible representations
    (irreps) from an irreps.yaml file and return an `Irreps` object.

    Parameters
    ----------
    file_path : str
        File path.

    Returns
    -------
    irrep_data : Irreps
        `Irreps` object containing the irreps.
    """

    data = load_yaml(file_path)

    if not np.allclose(data["q-position"], 0.0, atol=ZERO_TOLERANCE):
        raise RuntimeError(
            f"irreps.yaml file {file_path}: Irreps are for a non-Gamma q."
        )

    num_bands = sum(len(mode["band_indices"]) for mode in data["normal_modes"])

    ir_syms = np.array(
        [mode["ir_label"] for mode in data["normal_modes"]], dtype=np.str_
    )

    ir_band_inds = np.zeros((num_bands,), dtype=int)

    for i, mode in enumerate(data["normal_modes"]):
        for j in mode["band_indices"]:
            ir_band_inds[j - 1] = i

    return Irreps(data["point_group"], ir_syms, ir_band_inds)


# ----------
# HDF5 files
# ----------


def gamma_freqs_evecs_from_mesh_qpoints_or_band_hdf5(file_path):
    """Read Gamma-point phonon frequencies and eigenvectors from a
    mesh.hdf5 or band.hdf5 file.

    Parameters
    ----------
    file_path : str
        File path.

    Returns
    -------
    freqs_evecs : tuple of numpy.ndarray
        A `(freqs, evecs)` tuple of arrays with shapes `(3N,)` and
        `(3N, N, 3)`.
    """

    with h5py.File(file_path, "r") as f:
        if "eigenvector" not in f:
            raise RuntimeError(
                f"mesh.hdf5/band.hdf5 file {file_path}: Eigenvectors "
                "not found."
            )

        # mesh.hdf5 and band.hdf5 files have slightly different layouts.

        q_pts, freqs, evecs = None, None, None

        if "qpoint" in f:
            # mesh.hdf5 file.

            q_pts = f["qpoint"][:]
            freqs, evecs = f["frequency"][:], f["eigenvector"][:]
        elif "nqpoint" in f:
            path = f["path"][:]
            frequency = f["frequency"][:]

            n_seg, n_qpts, _ = path.shape
            _, _, n_bnd = frequency.shape

            q_pts = path.reshape((n_seg * n_qpts, 3))
            freqs = frequency.reshape((n_seg * n_qpts, n_bnd))

            evecs = f["eigenvector"][:].reshape((n_seg * n_qpts, n_bnd, n_bnd))
        else:
            raise RuntimeError(
                f"mesh.hdf5/band.hdf5 file {file_path}: Unknown data "
                "format."
            )

        for idx, q_pos in enumerate(q_pts):
            if np.allclose(q_pos, 0.0, atol=ZERO_TOLERANCE):
                freqs = freqs[idx]
                evecs = evecs[idx].T.reshape(-1, len(freqs) // 3, 3)

                evecs = np.real_if_close(evecs, tol=ZERO_TOLERANCE)

                return (freqs, evecs)

    raise RuntimeError(
        f"mesh.hdf5/band.hdf5 file {file_path}: Gamma-point "
        "frequencies/eigenvectors not found."
    )


def gamma_linewidths_from_kappa_hdf5(file_path, t=300.0):
    """Read Gamma-point linewidths at the specified temperature from a
    kappa-m*.hdf5 file.

    Parameters
    ----------
    file_path : str
        File path.

    Returns
    -------
    lws : numpy.ndarray
        Linewidths (shape: `(3N,)`).
    """

    with h5py.File(file_path, "r") as f:
        cond = np.isclose(f["temperature"][:], t, atol=ZERO_TOLERANCE)

        if cond.sum() != 1:
            raise RuntimeError(
                f"kappa-m*.hdf5 file {file_path}: Requested t = "
                f"{t:.2f} not found."
            )

        (t_idx,) = np.flatnonzero(cond)

        lws = None

        if "qpoint" in f:
            # gamma has shape (n_t, n_q, 3 n_a).

            for q_idx, q_pos in enumerate(f["qpoint"]):
                if np.allclose(q_pos, 0.0, atol=ZERO_TOLERANCE):
                    lws = f["gamma"][t_idx, q_idx]
        else:
            # gamma has shape (n_t, 3 n_a).
            lws = f["gamma"][t_idx]

        if lws is not None:
            # The "gamma" key in the kappa-m*.hdf5 files is defined such
            # that the phonon linewidths are 2 \Gamma.

            return 2.0 * lws

        raise RuntimeError(
            f"kappa-m*.hdf5 file {file_path}: Gamma-point linewidths "
            "not found."
        )


# ---------
# BORN file
# ---------


def eps_inf_and_born_from_born(file_path, struct):
    """Read the high-frequency dielectric constant and Born effective
    charges from a Phonopy BORN file and expand the charges for the
    supplied structure.

    Parameters
    ----------
    file_path : str
        File path.
    struct : Structure
        Crystal structure as a `Structure` object.

    Returns
    -------
    eps_born : tuple of numpt.ndarray
        A `(eps, born)` tuple with the high-frequency dielectric
        constant (shape: `(3, 3)`) and Born charges (shape:
        `(N, 3, 3)`).
    """

    if not _PHONOPY_AVAILABLE:
        raise RuntimeError(
            "eps_inf_and_born_from_born() requires the "
            "phonopy.file.IO.parse_born function."
        )

    born_data = parse_BORN(struct.to_phonopy_atoms(), filename=file_path)

    return (
        np.asanyarray(born_data["dielectric"]),
        np.asanyarray(born_data["born"]),
    )
