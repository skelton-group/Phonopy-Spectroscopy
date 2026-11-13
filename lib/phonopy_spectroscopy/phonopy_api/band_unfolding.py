# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Wrappers for using the Phonopy Python API."""

# -------
# Imports
# -------

import numpy as np

from ..interfaces.phonopy_interface import get_distance_unit_for_interface
from ..structure import Structure

from ..utility.numpy_helper import atleast_3d

from ..utility.structure import (
    fractional_to_cartesian_coordinates,
    map_atom_positions,
    invert_atom_map,
    map_qpoints,
)

from phonopy.phonon.band_structure import get_band_qpoints
from phonopy.structure.cells import Primitive
from phonopy.unfolding.core import Unfolding

# ------------------
# BandUnfolder class
# ------------------


class BandUnfolder:
    """Provides an interface to the Phonopy band-unfolding API."""

    def __init__(
        self,
        ph,
        ref_struct,
        sc_mat,
        prim_trans_mat=None,
        atom_map=None,
        ph_calc="vasp",
        symprec=1e-5,
    ):
        """Create a new instance of the `BandUnfolder` class.

        Parameters
        ----------
        ph : phonopy.api_phonopy.Phonopy
            Phonon calculation to unfold.
        ref_struct : Structure
            Reference structure to unfold onto.
        sc_mat : array_like
            Transformation matrix from the reference unit cell to the
            supercell used for the phonon calculation (shape: `(3, 3)`).
        prim_trans_mat : array_like, optional
            Transformation matrix from the reference unit cell to the
            primitive cell (optional, default: `None`).
        atom_map : array_like of int or None, optional
            Mapping of the atom positions in the calculation supercell
            to the positions in the reference structure (`None` attempts
            an automatic mapping).
        ph_calc : str, optional
            Calculator interface used to set up `ph` (default: "vasp").
        symprec : float, optional
            Symmetry tolerance used to set up `ph` (default: 1e-5).

        See Also
        --------
        utility.structure.map_atom_positions
            Create `atom_map` with more control over the algorithm.

        Notes
        -----
        Band unfolding requires at least three components:

        * A phonon calculation.
        * A reference structure specifying the "ideal" atom positions.
        * A map between the atoms in the calculation and reference
          structures.

        In most use cases, the calculation will be performed in a
        supercell with atomic substitutions, vacancies, defects, etc.
        In this case, the reference structure must be an equivalent
        "pristine" supercell, and a transformation matrix to the
        "parent" cell must be specified.

        The mapping can be set using the optional `atom_map` keyword.
        Integers correspond to the indices of atoms in `ref_struct`, and
        `None` specifies that atoms in the calculation that do not match
        any atoms in the reference structure.

        For many situations the map can be generated using the
        `map_atom_positions()` function in the `utility.structure`
        module.

        If `atom_map` is not set, the code will attempt to generate one
        by calling `map_atom_positions()` with the default parameters,
        but this is not guaranteed to work. Calling code will then be
        required to generate and specify the map.

        The supercell is in general related to a unit cell, and possibly
        a primitive cell, by the transformation matrices `sc_mat` and
        `prim_trans`:

        supercell <- `sc_mat` <- unit cell -> `prim_trans` -> primitive cell

        If not specified, `prim_trans` defaults to the identify matrix.

        The composite matrix used to transform the supercell to the
        primitive cell during unfolding is given by
        `sc_mat @ np.linalg.inv(prim_trans)`. This is computed
        automatically and the result is stored in
        `unfolding_supercell_matrix`.

        The unfolding routine performs phonon calculations and requires
        the force constants. Calling code must therefore set the
        `force_constants` property on `ph` before initialisation.

        The calculation structure is taken from `ph`, and the
        calculator used to prepare the phonon calculation must be
        specified to determine whether a distance unit conversion is
        required (the default is "vasp").

        Finally, the primitive cell used for unfolding is obtained from
        `ref_str` and the unfolding supercell matrix described above.
        This involves a symmetry tolerance, which can be set by the
        `symprec` parameter and which should generally be the same as
        used in the phonon calculation (the default is 1e-5).
        """

        if ph.force_constants is None:
            raise RuntimeError(
                "The force_constants property on ph must be set."
            )

        sc_mat = np.array(sc_mat, dtype=np.float64, copy=True)

        if sc_mat.shape != (3, 3):
            raise ValueError("sc_mat must be an array_like with shape (3, 3).")

        if prim_trans_mat is not None:
            prim_trans_mat = np.array(
                prim_trans_mat, dtype=np.float64, copy=True
            )

            if prim_trans_mat.shape != (3, 3):
                raise ValueError(
                    "If supplied, prim_trans must be an array_like "
                    "with shape (3, 3)."
                )
        else:
            # Identity matrix.
            prim_trans_mat = np.identity(3, dtype=np.float64)

        uf_sc_mat = sc_mat @ np.linalg.inv(prim_trans_mat)

        struct = Structure.from_phonopy_atoms(
            ph.supercell,
            distance_unit=get_distance_unit_for_interface(ph_calc),
        )

        # Generate the primitive cell from the reference structure -
        # needed to obtain the reciprocal lattice vectors for generating
        # unfolded band structures.

        prim = Primitive(
            ref_struct.to_phonopy_atoms(),
            np.linalg.inv(uf_sc_mat),
            symprec=symprec,
        )

        prim_struct = Structure.from_phonopy_atoms(
            prim,
            distance_unit=get_distance_unit_for_interface(ph_calc),
        )

        if atom_map is not None:
            # atom_map entries should be either integers or None, so
            # shallow copying should be fine. (The array contents are
            # valudated in the call to invert_atom_map().)
            atom_map = np.array(atom_map, dtype=object, copy=True)
        else:
            atom_map, _ = map_atom_positions(struct, ref_struct)

        # While it is arguably most intuitive to map the calculation
        # structure onto the reference structure, the Phonopy Unfolding
        # class requires the reverse.

        inv_atom_map = invert_atom_map(atom_map, struct, ref_struct)

        # The Unfolding class cannot handle "many -> one" mapping, where
        # multiple atoms in struct map to the same atom in ref_struct.

        for i, inds in enumerate(inv_atom_map):
            if np.ndim(inds) != 0:
                raise RuntimeError(
                    "The BandUnfolder class currently does not support "
                    "many -> one mapping."
                )

        prim_trans_mat.flags.writeable = False
        sc_mat.flags.writeable = False
        uf_sc_mat.flags.writeable = False

        atom_map.flags.writeable = False
        inv_atom_map.flags.writeable = False

        self._phonopy = ph

        self._struct = struct
        self._ref_struct = ref_struct
        self._prim_struct = prim_struct

        self._prim_trans_mat = prim_trans_mat
        self._sc_mat = sc_mat
        self._uf_sc_mat = uf_sc_mat

        self._atom_map = atom_map
        self._inv_atom_map = inv_atom_map

    @property
    def phonopy(self):
        """phonopy.api_phonopy.Phonopy : `Phonopy` object used by this
        class."""
        return self._phonopy

    @property
    def structure(self):
        """Structure : Structure for unfolding."""
        return self._struct

    @property
    def reference_structure(self):
        """Structure : Reference structure for unfolding."""
        return self._ref_struct

    @property
    def primitive_structure(self):
        """Structure : Primitive unit cell for unfolding."""
        return self._prim_struct

    @property
    def supercell_matrix(self):
        """numpy.ndarray : Supercell matrix (shape: `(3, 3)`)."""
        return self._sc_mat

    @property
    def primitive_matrix(self):
        """numpy.ndarray : Primitive transformation matrix (shape:
        `(3, 3)`)."""
        return self._prim_trans_mat

    @property
    def unfolding_supercell_matrix(self):
        """numpy.ndarray : Unfolding supercell matrix (shape: `(3, 3)`)."""
        return self._uf_sc_mat

    @property
    def atom_map(self):
        """numpy.ndarray : Mapping of atoms in `structure` to atoms in
        `reference_structure`."""
        return self._atom_map

    @property
    def inverse_atom_map(self):
        """numpy.ndarray : "Inverse" mapping of atoms in
        `reference_structure` to atoms in `structure`."""
        return self._inv_atom_map

    def _run_unfolding(self, qpts, ref_prim=None):
        """Run the band unfolding for a set of q-points with optional
        Brillouin zone mapping."""

        if ref_prim is not None:
            qpts = map_qpoints(qpts, ref_prim, self._prim_struct)

        unfolding = Unfolding(
            self._phonopy,
            self._uf_sc_mat,
            self._ref_struct.atom_positions,
            self._inv_atom_map,
            qpts,
        )

        unfolding.run()

        return (unfolding.frequencies, unfolding.unfolding_weights)

    def unfold_to_q(self, qpts, ref_prim=None):
        """Unfold to (a) specified q-point(s).

        Parameters
        ----------
        qpts : array_like
            q-point(s) to unfold (shape: `(3,)` or `(M, 3)`).
        ref_prim : Structure, optional
            Reference primitive structure for which the q-point
            coordinates are defined (default: `None`).

        Returns
        -------
        freqs_weights : tuple if numpy.ndarray
            Frequencies and unfolding weights (shape: `(3N,)` or
            `(M, 3N)`, where `N` is the number of atoms in
            `structure`.
        """

        q_pts = np.atleast_1d(np.asanyarray(q_pts, dtype=np.float64))

        ret_single = False

        match q_pts.shape:
            case (3,):
                q_pts = q_pts[np.newaxis, :]
                ret_single = True
            case (_, 3):
                pass
            case _:
                raise ValueError(
                    "q_pts must be an array_like with shape (3,) or (M, 3)."
                )

        freqs, weights = self._run_unfolding(qpts, ref_prim=ref_prim)

        return (
            freqs[0] if ret_single else freqs,
            weights[0] if ret_single else weights,
        )

    def unfold_band_structure(
        self, band_path, num_pts=101, var_seg_len=True, ref_prim=None
    ):
        """Unfold to q-points along a specified band path.

        Parameters
        ----------
        band_path : array_like
            Start/end points of each segment along band path (shape:
            `(M, 2, 3)`).
        num_pts : int, optional
            Number of points between in each segment (default: 101).
        var_seg_len : bool, optional
            If `True`, scale the number of points along segments based
            on the reciprocal lattice vectors (default: `True`).
        ref_prim : Structure, optional
            Reference primitive structure for which the q-point
            coordinates are defined (default: `None`).

        Returns
        -------
        uf_disp : list of tuple of numpy.ndarray
            Tuples of `(q, dists, freqs, weights)` for each segment of
            the unfolded dispersion.
        """

        band_path = atleast_3d(np.asanyarray(band_path, dtype=np.float64))

        if band_path.ndim > 3 or band_path.shape[-2:] != (2, 3):
            raise ValueError(
                "band_path must be an array_like with shape (M, 2, 3)."
            )

        if num_pts <= 1:
            raise ValueError("num_pts must be greater than one.")

        # If var_seg_len is set, pass the reciprocal lattice vectors
        # of the primitive cell to get_band_qpoints() to generate points
        # based on the relative lengths of the segments in the
        # dispersion.

        rec_v_latt = self._prim_struct.reciprocal_lattice_vectors()

        seg_qpts = get_band_qpoints(
            band_path,
            npoints=num_pts,
            rec_lattice=(rec_v_latt.T if var_seg_len else None),
        )

        seg_freqs, seg_weights = [], []

        for qpts in seg_qpts:
            freqs, weights = self._run_unfolding(qpts, ref_prim=ref_prim)

            seg_freqs.append(freqs)
            seg_weights.append(weights)

        # Calculate the (cumulative) distance along the path at each
        # q-point.

        seg_dists = []

        for qpts in seg_qpts:
            qpts = fractional_to_cartesian_coordinates(qpts, rec_v_latt)

            dists = [0.0] + [
                np.linalg.norm(q2 - q1) for q1, q2 in zip(qpts[:-1], qpts[1:])
            ]

            seg_dists.append(np.cumsum(dists, dtype=np.float64))

        for i, dists in enumerate(seg_dists[1:]):
            dists += seg_dists[i][-1]

        return [
            item for item in zip(seg_qpts, seg_dists, seg_freqs, seg_weights)
        ]
