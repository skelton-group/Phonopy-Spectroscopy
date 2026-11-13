# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Class and routines for storing and working with crystal structures."""

# -------
# Imports
# -------

import warnings

import numpy as np

from functools import lru_cache

from .constants import ZERO_TOLERANCE
from .units import convert_distance_units
from .utility.structure import fractional_to_cartesian_coordinates

try:
    import spglib as spg
except ImportError:
    from pyspglib import spglib as spg

from phonopy.structure.atoms import PhonopyAtoms, atom_data

# ---------
# Constants
# ---------

_ATOMIC_MASSES = {
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226,
    "Pu": 244,
    "Am": 243,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257,
    "Md": 258,
    "No": 259,
    "Lr": 262,
    "Rf": 267,
    "Db": 268,
    "Sg": 269,
    "Bh": 270,
    "Hs": 269,
    "Mt": 278,
    "Ds": 281,
    "Rg": 280,
    "Cn": 285,
    "Uut": 286,  # Uut = Nh
    "Uuq": 289,  # Uuq = Fl
    "Uup": 289,  # Uup = Mc
    "Uuh": 293,  # Uuh = Lv
    "Uus": 294,  # Uus = Ts
    "Uuo": 294,  # Uuo = Og
}
"""Atomic masses of selected elements that are present in the Phonopy
internal database but have their mass set to `None`."""

# ---------
# Functions
# ---------


@lru_cache(maxsize=294)
def lookup_atomic_mass(sym):
    """Look up an atomic mass from an atomic symbol.

    Parameters
    ----------
    sym : str
        Atomic symbol.

    Returns
    -------
    m : float
        Atomic mass (amu).
    """

    sym_check = str(sym).title()

    for _, db_sym, _, db_mass in atom_data:
        if sym_check == db_sym and db_mass is not None:
            return db_mass

    if sym_check in _ATOMIC_MASSES:
        return _ATOMIC_MASSES[sym_check]

    raise ValueError(f'Atomic mass for sym="{sym}" is not available.')


@lru_cache(maxsize=294)
def lookup_atomic_number(sym):
    """Look up an atomic number from an atomic symbol.

    Parameters
    ----------
    sym : str
        Atomic symbol.

    Returns
    -------
    at_num : int
        Atomic number.
    """

    sym_check = str(sym).title()

    for db_at_num, db_sym, _, _ in atom_data:
        if sym_check == db_sym:
            return db_at_num

    raise ValueError(f'Atomic number for sym="{sym}" is not available.')


# ---------------
# Structure class
# ---------------


class Structure:
    def __init__(
        self,
        latt_vecs,
        at_pos,
        at_typ,
        at_m=None,
        conv_trans=None,
    ):
        """Create a new instance of the `Structure` class.

        Parameters
        ----------
        latt_vecs : array_like
            Lattice vectors (shape: `(3, 3)`).
        at_pos : array_like
            Fractional atomic positions (shape: `(N, 3)`).
        at_typ : array_like
            Atom types (shape: `(N,)`).
        at_m : array_like, optional
            Atomic masses (optional, shape: `(N,)`).
        conv_trans : array_like, optional
            Specifies a transformation to the conventional unit cell
            (default: identity matrix).
        """

        latt_vecs = np.array(latt_vecs, dtype=np.float64, copy=True)

        if latt_vecs.shape != (3, 3):
            raise ValueError(
                "latt_vecs must be an array_like with shape (3, 3)."
            )

        at_pos = np.array(at_pos, dtype=np.float64, copy=True)

        # Allow "empty" structures.

        if at_pos.size == 0:
            at_pos = at_pos.reshape((0, 3))

        at_pos = np.atleast_2d(at_pos)

        if at_pos.ndim != 2 or (len(at_pos) > 0 and at_pos.shape[1] != 3):
            raise ValueError("at_pos must be an array_like with shape (N, 3).")

        at_typ = np.array(at_typ, dtype=str, copy=True)

        if at_typ.size == 0:
            at_typ = at_typ.reshape((0,))

        at_typ = np.atleast_1d(at_typ)

        if at_typ.ndim != 1 or len(at_typ) != len(at_pos):
            raise ValueError("at_typ must be an array_like with shape (N,).")

        if at_m is not None:
            at_m = np.array(at_m, dtype=np.float64, copy=True)

            if at_m.size == 0:
                at_m = at_m.reshape((0,))

            at_m = np.atleast_1d(at_m)

            if at_m.ndim != 1 or len(at_m) != len(at_pos):
                raise ValueError("at_m must be an array_like with shape (N,).")
        else:
            at_m = np.array(
                [lookup_atomic_mass(sym) for sym in at_typ], dtype=np.float64
            )

        if (at_m <= 0.0).any():
            raise ValueError("Atomic masses must be non-zero.")

        if conv_trans is not None:
            conv_trans = np.array(conv_trans, dtype=np.float64, copy=True)

            if conv_trans.shape != (3, 3):
                raise ValueError(
                    "If supplied, conv_trans must be an array_like "
                    "with shape (3, 3)."
                )

            # A valid transformation matrix to a conventional cell
            # should have integer elements, although this may not be
            # the case if the matrix has not been specified with
            # sufficient precision.

            max_abs_diff = np.abs(np.rint(conv_trans) - conv_trans).max()

            if (max_abs_diff > ZERO_TOLERANCE).any():
                warnings.warn(
                    "One or more elements in conv_trans deviates from "
                    f"integer values by up to {max_abs_diff:.3e}. This "
                    "could indicate an invalid trandformation matrix "
                    "or insufficient precision.",
                    UserWarning,
                )
        else:
            conv_trans = np.identity(3, dtype=np.float64)

        latt_vecs_conv = conv_trans @ latt_vecs

        latt_vecs.flags.writeable = False
        latt_vecs_conv.flags.writeable = False

        at_pos.flags.writeable = False
        at_typ.flags.writeable = False
        at_m.flags.writeable = False

        conv_trans.flags.writeable = False

        self._v_latt = latt_vecs
        self._v_latt_conv = latt_vecs_conv

        self._at_pos = at_pos
        self._at_typ = at_typ
        self._at_m = at_m

        self._conv_trans = conv_trans

        self._volume = np.abs(np.linalg.det(latt_vecs))

        self._centrosymmetric = None
        self._all_wyckoff = None

    @property
    def lattice_vectors(self):
        """numpy.ndarray : Lattice vectors (shape: `(3, 3)`)."""
        return self._v_latt

    @property
    def conventional_lattice_vectors(self):
        """numpy.ndarray : Lattice vectors of the conventional cell
        (shape: `(3, 3)`."""
        return self._v_latt_conv

    @property
    def atom_positions(self):
        """numpy.ndarray : Atomic positions (shape: `(N, 3)`)."""
        return self._at_pos

    @property
    def atom_types(self):
        """numpy.ndarray : Atom types (shape: `(N,)`)."""
        return self._at_typ

    @property
    def atomic_masses(self):
        """numpy.ndarray : Atomic masses (shape: `(N,)`)."""
        return self._at_m

    @property
    def conventional_transformation_matrix(self):
        """numpy.ndarray : Transformation matrix to convert the
        structure to its conventional cell (shape: `(3, 3)`)."""
        return self._conv_trans

    @property
    def num_atoms(self):
        """int : Number of atoms in the structure."""
        return len(self._at_pos)

    @property
    def volume(self):
        """float : Unit-cell volume."""
        return self._volume

    def reciprocal_lattice_vectors(self, conv=False, two_pi=True):
        r"""Calculate and return the reciprocal lattice vectors.

        Parameters
        ----------
        conv : bool, optional
            If `True`, return the volume of the conventional unit cell
            (default: `False`).
        two_pi : bool, optional
            If `True`, compute the "Physics" definition of the
            reciprocal lattice vectors with a prefactor of 2 \pi
            (default: `True`).

        Returns
        -------
        recip_latt_vec : numpy.ndarray
            Reciprocal lattice vectors (shape: `(3, 3)`).
        """

        v_latt = self._v_latt_conv if conv else self._v_latt
        rec_v_latt = np.linalg.inv(v_latt).T

        if two_pi:
            rec_v_latt *= 2.0 * np.pi

        return rec_v_latt

    def real_space_normal(self, hkl, conv=True):
        """Calculate the real-space normal to the surface with Miller
        index `hkl`.

        Parameters
        ----------
        hkl : array_like
            Integer Miller indices of the surface (shape: `(3,)`).,
        conv : bool, optional
            If `True`, interpret the Miller indices in terms of the
            conventional unit cell (default: `True`).

        Returns
        -------
        norm : numpy.ndarray
            Real-space surface normal in Cartesian coordinates (shape:
            `(3,)`).
        """

        hkl = np.asanyarray(hkl, dtype=int)

        if hkl.shape != (3,):
            raise ValueError("hkl must be an array_like with shape `(3,)`.")

        # Whether we use the "Physics" or "crystallography" definition
        # of the reciprocal lattice doesn't matter because we normalise
        # the result.

        rec_v_latt = self.reciprocal_lattice_vectors(conv=conv, two_pi=False)

        v = np.matmul(hkl, rec_v_latt)
        return v / np.linalg.norm(v)

    def atomic_numbers(self):
        """Return the atomic numbers of the atoms.

        Returns
        -------
        at_nums : numpy.ndarray
            Atomic numbers (shape: `(N,)`).
        """

        return np.array(
            [lookup_atomic_number(sym) for sym in self._at_typ], dtype=int
        )

    def cartesian_positions(self):
        """Return the atomic positions in Cartesian coordinates.

        Returns
        -------
        pos_cart : numpy.ndarray
            Atom positions in Cartesian coordinates (shape: `(N, 3)`).
        """

        return fractional_to_cartesian_coordinates(self._at_pos, self._v_latt)

    def build_supercell(self, dim):
        """Build a supercell of the structure.

        Parameters
        ----------
        dim : array_like
            Integer supercell expansion (shape: `(3,)`).

        Returns
        -------
        sc : Structure
            Supercell.

        Notes
        -----
        The supercell retains the original atom ordering - i.e. for a
        supercell expansion with `m = np.prod(dim)` times the original
        cell volume, atoms `0:m` are the replicas of atom `0` in the
        original cell, atoms `m:2 * m` are replicas of atom `1`, etc.
        """

        dim = np.asanyarray(dim, dtype=int)

        if dim.shape != (3,):
            raise ValueError("dim must be an array_like with shape (3,).")

        if (dim <= 0).any():
            raise ValueError("Supercell dimensions must be non-zero.")

        sc_latt_vecs = dim[:, np.newaxis] * self._v_latt

        # Prepare a set of translations for building the supercell, then
        # generate and scale the translated atom positions.

        grid = np.indices(dim, dtype=int).reshape(3, -1).T

        sc_at_pos = self._at_pos[:, np.newaxis, :] + grid[np.newaxis, :, :]
        sc_at_pos /= dim[np.newaxis, np.newaxis, :]

        mult = np.prod(dim)

        return Structure(
            sc_latt_vecs,
            sc_at_pos.reshape(-1, 3),
            np.repeat(self._at_typ, mult),
            np.repeat(self._at_m, mult),
            conv_trans=self._conv_trans,
        )

    def check_centrosymmetry(self, symprec=1.0e-5):
        """Determine whether the structure is in a centrosymmetric space
        group and all atoms occupy Wyckoff sites with inversion
        symmetry.

        Params
        ------
        symprec : float
            Precision for symmetry analysis (default: 1.0e-5).

        Returns
        -------
        centrosymm : tuple of bool
            Tuple of `(is_centrosymm, all_on_wyckoff)` indicating
            whether the structure is in a centrosymmetric spacegroup
            and whether all the atoms are located on Wyckoff sites.
        """

        # Perform symmetry analysis and get symmetry operations.

        d = spg.get_symmetry_dataset(
            (self._v_latt, self._at_pos, self.atomic_numbers()),
            symprec=symprec,
        )

        op_rs = d.rotations
        op_ts = d.translations

        # Mask inversion operations.

        op_i = -1.0 * np.eye(3, 3)

        mask = np.all(
            np.isclose(op_rs, op_i, atol=ZERO_TOLERANCE), axis=(1, 2)
        )

        if not mask.any():
            # No "global" inversion symmetry.
            return (False, False)

        for r, t in zip(op_rs[mask], op_ts[mask]):
            # (N, 3) @ (3, 3) -> (N, 3) - need to use r.T for right
            # multiplication.

            diff = ((self._at_pos @ r.T) + t) - self._at_pos
            diff = diff - np.round(diff)

            if np.isclose(diff, 0.0, atol=ZERO_TOLERANCE).all(axis=1).all():
                # Inversion operation maps all positions to themselves
                # -> all atoms on sites with inversion symmetry.

                return (True, True)

        return (True, False)

    def to_phonopy_atoms(self, distance_unit="ang"):
        """Return the structure as a `PhonopyAtoms` object.

        Params
        ------
        distance_unit : str, optional
            Specify the distance unit to be used in the returned
            `PhonopyAtoms` object (default: "ang").

        Returns
        -------
        atoms : PhonopyAtoms
            `PhonopyAtoms` object.
        """

        v_latt = self._v_latt

        if distance_unit != "ang":
            v_latt = convert_distance_units(v_latt, "ang", distance_unit)

        return PhonopyAtoms(
            cell=v_latt.tolist(),
            scaled_positions=self.atom_positions.tolist(),
            symbols=self.atom_types.tolist(),
            masses=self.atomic_masses.tolist(),
        )

    @staticmethod
    def from_phonopy_atoms(atoms, distance_unit="ang", **kwargs):
        """Create a `Structure` from a `PhonopyAtoms` object.

        Params
        ------
        atoms : PhonopyAtoms
            `PhonopyAtoms` object.
        distance_unit : str, optional
            Distance unit used in `atoms` (default: "ang").
        kwargs : any, optional
            Optional keyword arguments to the `Structure` class
            constructor.

        Returns
        -------
        struct : Structure
            `Structure` object.
        """

        v_latt = atoms.cell

        if distance_unit != "ang":
            v_latt = convert_distance_units(v_latt, distance_unit, "ang")

        return Structure(
            v_latt,
            atoms.scaled_positions,
            atoms.symbols,
            atoms.masses,
            **kwargs,
        )

    def to_dict(self):
        """Return the internal data as a dictionary of native Python
        types for serialisation.

        Returns
        -------
        d : dict
            Dictionary structure containing internal data as native
            Python types.
        """

        return {
            "lattice_vectors": self._v_latt.tolist(),
            "atom_positions": self._at_pos.tolist(),
            "atom_types": list(self._at_typ),
            "atomic_masses": self._at_m.tolist(),
            "conventional_transformation_matrix": self._conv_trans.tolist(),
        }

    @staticmethod
    def from_dict(d):
        """Create a new `Structure` instance from a dictionary
        generated by `Structure.to_dict()`.

        Parameters
        ----------
        d : dict
            Dictionary generated by `to_dict()`.

        Returns
        -------
        struct : Structure
            `Structure` object constructed from the data in `d`.
        """

        return Structure(
            d["lattice_vectors"],
            d["atom_positions"],
            d["atom_types"],
            at_m=d["atomic_masses"],
            conv_trans=d["conventional_transformation_matrix"],
        )
