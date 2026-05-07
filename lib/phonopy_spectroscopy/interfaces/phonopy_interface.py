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


from ..constants import ZERO_TOLERANCE
from ..irreps import Irreps
from ..gamma_phonons import GammaPhonons, PolarGammaPhonons
from ..structure import Structure
from ..utility.io_helper import load_yaml

from .vasp_interface import structure_from_poscar


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

    raise ValueError('Unknown interface "{0}".'.format(calculator))


# -------------------
# High-level "loader"
# -------------------


def gamma_phonons_from_phono3py(
    cell_file,
    freqs_evecs_file,
    lws_file=None,
    lws_t=300.0,
    irreps_file=None,
    born_file=None,
    at_m=None,
    conv_trans=None,
    discard_imag=True,
):
    r"""Read a complete Phono(3)py calculation and return a
    `GammaPhonons` or `PolarGammaPhonons` object.

    Parameters
    ----------
    cell_file : str
        VASP POSCAR or phonopy.yaml file to read structure from.
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
    at_m : array_like, optional
        Atomic masses (optional, default: `None`; overridden if
        `cell_file` is a phonopy.yaml file).
    conv_trans : array_like, optional
        Transformation matrix to convert the structure to its
        conventional cell (shape: `(3, 3)`, default: `None`).
    discard_imag : bool, optional
        Discard the imaginary part of complex eigenvectors (default:
        `True`).

    Returns
    -------
    gamma_ph : `GammaPhonons`
        `GammaPhonons` object containing the calculation data.
    """

    # Read a structure. If cell_file has a .yaml extension, assume it
    # is a phonopy.yaml file; otherwise, assume it is a VASP POSCAR
    # file.

    struct = None

    _, ext = os.path.splitext(cell_file)

    if ext.lower() == ".yaml":
        struct = structure_from_phonopy_yaml(cell_file, conv_trans=conv_trans)
    else:
        struct = structure_from_poscar(
            cell_file, at_m=at_m, conv_trans=conv_trans
        )

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
            "Frequencies/eigenvectors file {0}: unknown format.".format(
                freqs_evecs_file
            )
        )

    if discard_imag:
        if np.iscomplex(evecs).any():
            # Gamma-point eigenvectors should be real unless a
            # non-analytical correction is applied to the dynamical
            # matrix. If this is the case, it will conflict with the
            # internal NAC implementation and is very likely a user
            # error.

            max_abs_imag = np.abs(evecs.imag).max()

            if max_abs_imag > ZERO_TOLERANCE:
                raise RuntimeError(
                    "Discarding imaginary parts of eigenvectors "
                    "with maximum absolute value {0:.3e} > "
                    "ZERO_TOLERANCE = {1:.3e}. This may indicate "
                    "a calculation performed with a non-analytical "
                    "correction to the dynamical matrix."
                    "".format(max_abs_imag, ZERO_TOLERANCE)
                )

        evecs = evecs.real

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

    if born_file is not None:
        # Read BORN file and return a PolarGammaPhonons object.

        eps_inf, born_charges = hf_dielectric_and_born_from_born(
            born_file, struct
        )

        return PolarGammaPhonons(
            struct,
            freqs,
            evecs,
            eps_inf,
            born_charges,
            lws=lws,
            irreps=irreps,
            t=t,
        )

    # Return a GammaPhonons object.

    return GammaPhonons(struct, freqs, evecs, lws=lws, irreps=irreps, t=t)


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
    r"""Read Gamma-point phonon frequencies and eigenvectors
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
                    "mesh.yaml/band.yaml file {0}: Eigenvectors not found."
                    "".format(file_path)
                )

            freqs = np.array(
                [mode["frequency"] for mode in qpt["band"]], dtype=np.float64
            )

            evecs = np.array(
                [mode["eigenvector"] for mode in qpt["band"]], dtype=np.float64
            )

            evecs = evecs[:, :, :, 0] + 1.0j * evecs[:, :, :, 1]
            evecs = evecs.reshape(-1, len(freqs) // 3, 3)

            if not np.iscomplex(evecs).any():
                evecs = evecs.real

            return (freqs, evecs)

    raise RuntimeError(
        "mesh.yaml/band.yaml file {0}: Gamma-point "
        "frequencies/eigenvectors not found.".format(file_path)
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
            "irreps.yaml file {0}: Irreps are for a non-Gamma q."
            "".format(file_path)
        )

    return Irreps(
        str(data["point_group"]),
        [mode["ir_label"] for mode in data["normal_modes"]],
        [
            [idx - 1 for idx in mode["band_indices"]]
            for mode in data["normal_modes"]
        ],
    )


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
                "mesh.hdf5/band.hdf5 file {0}: Eigenvectors not found."
                "".format(file_path)
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
                "mesh.hdf5/band.hdf5 file {0}: Unknown data format."
                "".format(file_path)
            )

        for idx, q_pos in enumerate(q_pts):
            if np.allclose(q_pos, 0.0, atol=ZERO_TOLERANCE):
                freqs = freqs[idx]
                evecs = evecs[idx].T.reshape(-1, len(freqs) // 3, 3)

                if not np.iscomplex(evecs).any():
                    evecs = evecs.real

                return (freqs, evecs)

    raise RuntimeError(
        "mesh.hdf5/band.hdf5 file {0}: Gamma-point "
        "frequencies/eigenvectors not found.".format(file_path)
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
        cond = f["temperature"][:] == t

        if cond.sum() != 1:
            raise RuntimeError(
                "kappa-m*.hdf5 file {0}: Requested t = {1:.2f} not "
                "found.".format(file_path, t)
            )

        ((t_idx,),) = np.where(cond)

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
            "kappa-m*.hdf5 file {0}: Gamma-point linewidths not found."
            "".format(file_path)
        )


# ---------
# BORN file
# ---------


def hf_dielectric_and_born_from_born(file_path, struct):
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
            "read_hf_dielectric_and_born_from_born() requires the "
            "phonopy.file.IO.parse_born function."
        )

    born_data = parse_BORN(struct.to_phonopy_atoms(), filename=file_path)

    return (
        np.asarray(born_data["dielectric"], dtype=np.float64),
        np.asarray(born_data["born"], dtype=np.float64),
    )
