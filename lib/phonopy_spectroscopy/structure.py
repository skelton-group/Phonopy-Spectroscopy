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

from .constants import ZERO_TOLERANCE

from .units import convert_distance_units

from .utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
)

from .utility.structure import (
    cartesian_to_fractional_coordinates,
    fractional_to_cartesian_coordinates,
)

_PHONOPY_AVAILABLE = False

try:
    from phonopy.structure.atoms import PhonopyAtoms, atom_data

    _PHONOPY_AVAILABLE = True
except ImportError:
    warnings.warn(
        "Imports from phonopy failed - some functions require phonopy "
        "and will raise exceptions if it is not installed.",
        RuntimeWarning,
    )


# ---------
# Functions
# ---------


def lookup_atomic_mass(symbol):
    """Lookup an atomic mass from an atomic symbol.

    Parameters
    ----------
    symbol : str
        Atomic symbol.

    Returns
    -------
    m : float
        Atomic mass (amu).

    Notes
    -----
    This function requires the `phonopy` package.
    """

    if not _PHONOPY_AVAILABLE:
        raise RuntimeError(
            "lookup_atomic_mass() requires the "
            "phonopy.atoms.atom_data attribute."
        )

    symbol = str(symbol).title()

    for _, db_symbol, _, db_mass in atom_data:
        if symbol == db_symbol and db_mass is not None:
            return db_mass

    raise ValueError(
        'Atomic mass for symbol="{0}" not available in '
        "phonopy.atoms.atom_data.".format(symbol)
    )


def lookup_atomic_number(symbol):
    """Lookup an atomic number from a symbol.

    Parameters
    ----------
    symbol : str
        Atomic symbol.

    Returns
    -------
    num : int
        Atomic number.

    Notes
    -----
    This function requires the `phonopy` package.
    """

    if not _PHONOPY_AVAILABLE:
        raise RuntimeError(
            "lookup_atomic_number() requires the "
            "phonopy.atoms.atom_data attribute."
        )

    for db_at_num, db_symbol, _, _ in atom_data:
        if symbol == db_symbol:
            return db_at_num

    raise ValueError(
        'Atomic number for symbol="{0}" not available in '
        "phonopy.atoms.atom_data.".format(symbol)
    )


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
        cart_to_frac=False,
    ):
        """Create a new instance of the `Structure` class.

        Parameters
        ----------
        latt_vecs : array_like
            Lattice vectors (shape: `(3, 3)`).
        at_pos : array_like
            Atomic positions (shape: `(N, 3)`).
        at_typ : array_like
            Atom types (shape: `(N,)`).
        at_m : array_like, optional
            Atomic masses (optional, shape: `(N,)`).
        conv_trans : array_like, optional
            Specifies a transformation to the conventional unit cell
            (default: identity matrix).
        cart_to_frac : bool, optional
            If `True`, convert `at_pos` from Cartesian to fractional
            coordinates (default: `False`).
        """

        latt_vecs = np_asarray_copy(latt_vecs, dtype=np.float64)

        if not np_check_shape(latt_vecs, (3, 3)):
            raise ValueError(
                "latt_vecs must be an array_like with shape (3, 3)."
            )

        at_pos = np_asarray_copy(at_pos, dtype=np.float64)
        at_typ = np.array([str(typ) for typ in at_typ], dtype=object)

        if len(at_pos) > 0 and not np_check_shape(at_pos, (None, 3)):
            raise ValueError("at_pos must be an array_like with shape (N, 3).")

        n_a = len(at_pos)

        if not np_check_shape(at_typ, (n_a,)):
            raise ValueError("at_typ must be an array_like with shape (N,).")

        if at_m is None:
            at_m = np.array(
                [lookup_atomic_mass(sym) for sym in at_typ], dtype=np.float64
            )

            if (at_m <= 0.0).any():
                warnings.warn(
                    "Atomic mass lookup returned m <= 0 for one or "
                    "more atoms. The atomic masses likely need to be "
                    "specified explicitly with the at_m keyword.",
                    RuntimeWarning,
                )
        else:
            at_m = np_asarray_copy(at_m, dtype=np.float64)

            if not np_check_shape(at_m, (n_a,)):
                raise ValueError(
                    "If supplied, at_m must be an array_like with shape (N,)."
                )

        if (at_m <= 0.0).any():
            raise ValueError("Atomic masses must be larger than zero.")

        if conv_trans is not None:
            conv_trans = np_asarray_copy(conv_trans, dtype=np.float64)

            if not np_check_shape(conv_trans, (3, 3)):
                raise ValueError(
                    "If supplied, prim_trans must be an array_like "
                    "with shape (3, 3)."
                )

            # A valid transformation matrix to a conventional cell
            # should have integer elements, although this may not be
            # the case if the matrix has not been specified with
            # sufficient precision.

            abs_diff = np.abs(np.rint(conv_trans) - conv_trans)

            if (abs_diff > ZERO_TOLERANCE).any():
                warnings.warn(
                    "One or more elements in conv_trans deviates from "
                    "integer values by up to {0:.3e}. This could "
                    "indicate an invalid tranformation matrix or "
                    "insufficient precision.".format(abs_diff.max()),
                    UserWarning,
                )
        else:
            conv_trans = np.identity(3, dtype=np.float64)

        latt_vecs_conv = np.matmul(conv_trans.T, latt_vecs)

        if n_a > 0:
            if cart_to_frac:
                at_pos = cartesian_to_fractional_coordinates(at_pos, latt_vecs)
            else:
                if (np.abs(at_pos) > 1.0).any():
                    warnings.warn(
                        "One or more of at_pos are outside the range "
                        "[-1, 1] expected for fractional coordinates - "
                        "use cart_to_frac=True to convert if needed.",
                        UserWarning,
                    )

        self._v_latt = latt_vecs
        self._v_latt_conv = latt_vecs_conv

        self._at_pos = at_pos
        self._at_typ = at_typ
        self._at_m = at_m

        self._conv_trans = conv_trans

    @property
    def lattice_vectors(self):
        """numpy.ndarray : Lattice vectors (shape: `(3, 3)`)."""
        return np_readonly_view(self._v_latt)

    @property
    def primitive_lattice_vectors(self):
        """numpy.ndarray : Lattive vectors of the primitive cell
        (shape: `(3, 3)`, alias for `lattice_vectors`."""
        return self.lattice_vectors

    @property
    def conventional_lattice_vectors(self):
        """numpy.ndarray : Lattive vectors of the conventional cell
        (shape: `(3, 3)`."""
        return np_readonly_view(self._v_latt_conv)

    @property
    def atom_positions(self):
        """numpy.ndarray : Atomic positions (shape: `(N, 3)`)."""
        return np_readonly_view(self._at_pos)

    @property
    def atom_types(self):
        """numpy.ndarray : Atom types (shape: `(N,)`)."""
        return np_readonly_view(self._at_typ)

    @property
    def atomic_masses(self):
        """numpy.ndarray : Atomic masses (shape: `(N,)`)."""
        return np_readonly_view(self._at_m)

    @property
    def conventional_transformation_matrix(self):
        """numpy.ndarray : Transformation matrix to convert the
        structure to its conventional cell."""
        return np_readonly_view(self._conv_trans)

    @property
    def num_atoms(self):
        """int : Number of atoms in the structure."""
        return self._at_pos.shape[0]

    def volume(self, conv=False):
        """Calculate the unit-cell volume.

        Parameters
        ----------
        conv : bool, optional
            If `True`, return the volume of the conventional unit cell
            (default: `False`).

        Returns
        -------
        v : float
            Unit-cell volume.
        """

        v_latt = self._v_latt_conv if conv else self._v_latt
        return np.abs(np.linalg.det(v_latt))

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
            Recipocal lattice vectors (shape: `(3, 3)`).
        """

        v_latt = self._v_latt_conv if conv else self._v_latt
        rec_v_latt = np.linalg.inv(v_latt).T

        return (2.0 * np.pi * rec_v_latt) if two_pi else rec_v_latt

    def real_space_normal(self, hkl, conv=False):
        """Calculate the real-space normal to the surface with Miller
        index `hkl`.

        Parameters
        ----------
        hkl : array_like
            Integer Miller indices of the surface (shape: `(3,)`).,
        conv : bool, optional
            If `True`, return the volume of the conventional unit cell
            (default: `False`).

        Returns
        -------
        norm : numpy.ndarray
            Real-space sufrace normal in Cartesian coordinates (shape:
            `(3,)`).
        """

        hkl = np.asarray(hkl)

        if not np_check_shape(hkl, (3,)):
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

        # Avoid looking up the same symbol multiple times.

        at_nums_lut = {
            sym: lookup_atomic_number(sym) for sym in np.unique(self._at_typ)
        }

        return np.array([at_nums_lut[sym] for sym in self._at_typ], dtype=int)

    def cartesian_positions(self):
        """Return the atomic positions converted to Cartesian
        coordinates.

        Returns
        -------
        pos_cart : numpy.ndarray
            Atom positions in Cartesian coordinates (shape: `(N, 3)`).
        """

        return fractional_to_cartesian_coordinates(self._at_pos, self._v_latt)

    def to_phonopy_atoms(self, distance_unit="ang"):
        """Return the structure as a `PhonopyAtoms` instance.

        Params
        ------
        distance_unit : str, optional
            Specify the distance unit to be used in the returned
            `PhonopyAtoms` object (default: "ang"").

        Returns
        -------
        atoms : PhonopyAtoms
            `PhonopyAtoms` object containing the structure data.

        Notes
        -----
        This function requires the `phonopy` package.
        """

        if not _PHONOPY_AVAILABLE:
            raise RuntimeError(
                "Structure.to_phonopy_atoms() requires the "
                "phonopy.structure.PhonopyAtoms class."
            )

        v_latt = self._v_latt

        if distance_unit != "ang":
            v_latt = convert_distance_units(v_latt, "ang", distance_unit)

        # The phonopy API uses the idiom "if x" to detect when a
        # parameter x is set, which raises if x is a NumPy array with
        # more than one element.

        return PhonopyAtoms(
            cell=v_latt.tolist(),
            scaled_positions=self.atom_positions.tolist(),
            symbols=self.atom_types.tolist(),
            masses=self.atomic_masses.tolist(),
        )

    @staticmethod
    def from_phonopy_atoms(atoms, distance_unit="ang", conv_trans=None):
        """Create a new `Structure` instance from a `PhonopyAtoms`
        object.

        Params
        ------
        atoms : PhonopyAtoms
            `PhonopyAtoms` object.
        distance_unit : str, optional
            Specify the distance unit used in `atoms` (default: "ang").
        conv_trans : array_like, optional
            Specify the optional `conv_trans` keyword to the `Structure`
            class constructor (default: `None`).

        Returns
        -------
        struct : Structure
            `Structure` object constructed from the data in `atoms`.
        """

        v_latt = atoms.cell

        if distance_unit != "ang":
            v_latt = convert_distance_units(v_latt, distance_unit, "ang")

        return Structure(
            v_latt,
            atoms.scaled_positions,
            atoms.symbols,
            atoms.masses,
            conv_trans=conv_trans,
            cart_to_frac=False,
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
            cart_to_frac=False,
        )
