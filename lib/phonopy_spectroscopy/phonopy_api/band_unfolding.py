# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Wrappers for using the Phonopy Python API."""


# -------
# Imports
# -------


import warnings

import numpy as np

from ..structure import (
    cartesian_to_fractional_coordinates,
    fractional_to_cartesian_coordinates,
    Structure,
)

from ..interfaces.phonopy_interface import get_distance_unit_for_interface

from ..utility.numpy_helper import (
    np_asarray_copy,
    np_check_shape,
    np_readonly_view,
    np_expand_dims,
)

from phonopy.structure.cells import Primitive
from phonopy.phonon.band_structure import get_band_qpoints
from phonopy.unfolding.core import Unfolding


# -----------------
# Structure mapping
# -----------------


def _exception_or_user_warning(msg, warn=False):
    """Handle an error message `msg` by raising an `Exception`
    (`warn=False`, default) or a `UserWarning` (`warn=True`).
    """

    if warn:
        warnings.warn(msg, UserWarning)
    else:
        raise Exception(msg)


def map_atom_positions(str_map, str_ref, constraints=None, warn=False):
    """Generate an integer mapping of the closest atomic positions in
    two `Structure` objects.

    Parameters
    ----------
    str_map, str_ref : Structure
        Structures to map.
    constraints : list of (tuple of (list of str or None))
       Specifies pairwise sets of atom types in `str_map` and `str_ref`
       that should be matched (see Notes for examples).
    warn : bool
        Downgrade exceptions raised when potential problems with the
        mapping are detected to warnings (default: `False`).

    Returns
    -------
    map : tuple of numpy.ndarray
        Integer mapping of the atom positions in `str_map` and distances
        to the closest positions in `str_ref`.

    Notes
    -----
    The default behaviour without `constraints` is suitable for a number
    of common matching problems, including e.g. small symmetry-breaking
    distortions, vacancies and atomic substitutions.

    The routine checks for differences of >0.5% in the lengths or >1 deg
    in the angles between lattice vectors and non-unique mapping, and
    raises exceptions if these are found. Setting `warn=True` downgrades
    these to warnings.

    `constraints` can be used to handle more difficult cases. This
    allows the caller to specify pairwise sets of atom types in
    `str_map` and `str_ref` that should be matched. Some examples:

    * `constraints=[(["Ce"], ["Ce"]), (["O"], ["O"])]` specifies that
      atom types should be matched (e.g. for Frenkel defects).
    * `constraints=[(["Sn"], ["Sn"]), (["S", "Se"], ["S"])]` specifies
      that atoms of the same "identity" (cation/anion) should be
      matched (e.g. in alloys).
    * `constraints=[(["H"], [None]), ...]` specifies that some types of
      atom in `str_map` should not be matched (e.g. for impurities).
    * `constraints=[(["H", "C", "N"], ["Pb"]), ...]` specifies that
      multiple atom types in `str_map` should be matched to a single
      atom type in `str_ref` (e.g. for hybrid organic/inorganic
      materials).

    Setting `constraints` allows non-unique mapping, but the parameters
    must cover all atom types in `str_map` including those that do not
    require special handling.

    While the routine should handle a number of typical scenarios, some
    may require further modification. For example:

    * Frenkel defect where the interstitials are far from, and cannot be
      mapped to, the vacancy in the reference structure.
    * Impurities of the same type of the host atoms (e.g. excess
      oxygen).

    Both cases could be approached by setting `constraints` or
    `warn=True` to allow for non-unique mapping, and identifying the
    atoms that require correction based on the returned interatomic
    distances.
    """

    # The approach of mapping the closest atoms assumes the two
    # structures have similar lattice vectors. Check this and issue a
    # warning if this is not the case.

    norms_map = np.linalg.norm(str_map.lattice_vectors, axis=1)
    norms_ref = np.linalg.norm(str_ref.lattice_vectors, axis=1)

    norm_diff = np.abs(norms_map - norms_ref) / norms_ref

    if (norm_diff > 5.0e-3).any():
        _exception_or_user_warning(
            "Maximum difference in lattice vector lengths is {0:.2f}% "
            "> 1%. ".format(100.0 * norm_diff.max()),
            warn=warn,
        )

    thetas = np.zeros((3,), dtype=np.float64)

    for idx in range(3):
        dp = np.dot(str_map.lattice_vectors[idx], str_ref.lattice_vectors[idx])
        thetas[idx] = np.acos(dp / (norms_map[idx] * norms_ref[idx]))

    thetas = np.abs(np.degrees(thetas))

    if (thetas > 1.0).any():
        _exception_or_user_warning(
            "Largest angle between lattice vector is {0:.2f} > 1 deg."
            "".format(thetas.max()),
            warn=warn,
        )

    # Perform the mapping on pairs of groups of indices at a time. This
    # allows constraints on atom types to be applied.

    index_grps = []

    if constraints is not None:
        for map_typs, ref_typs in constraints:
            if len(map_typs) == 0:
                raise ValueError(
                    "Mapping atom types in constraints must be an "
                    "array_like with at least one element."
                )

            inds_map = []

            for sym in map_typs:
                (inds,) = np.where(str_map.atom_types == sym)
                inds_map.extend(inds)

            inds_ref = None

            if ref_typs is not None:
                if len(ref_typs) == 0:
                    raise ValueError(
                        "Reference atom types in constraints must "
                        "either be an array_like with at least one "
                        "element or None."
                    )

                inds_ref = []

                for sym in ref_typs:
                    (inds,) = np.where(str_ref.atom_types == sym)
                    inds_ref.extend(inds)

            index_grps.append((inds_map, inds_ref))
    else:
        index_grps = [
            (
                np.arange(0, str_map.num_atoms, dtype=int),
                np.arange(0, str_ref.num_atoms, dtype=int),
            )
        ]

    # Adjust for differences in unit cells by converting the positions
    # into Cartesian coordinates and then back into fractional
    # coordinates using the lattice vectors of the reference structure.

    str_map_pos_cart = fractional_to_cartesian_coordinates(
        str_map.atom_positions, str_map.lattice_vectors
    )

    str_map_pos_frac_shift = cartesian_to_fractional_coordinates(
        str_map_pos_cart, str_ref.lattice_vectors
    )

    # Perform mapping.

    atom_mapping_dict = {}

    for inds_map, inds_ref in index_grps:
        if inds_ref is not None:
            pos_map = str_map_pos_frac_shift[inds_map]
            pos_ref = str_ref.atom_positions[inds_ref]

            vecs = pos_map[:, np.newaxis, :] - pos_ref[np.newaxis, :, :]

            # Apply periodic boundary conditions.

            vecs[vecs < -0.5] += 1.0
            vecs[vecs >= 0.5] -= 1.0

            # Convert fractional to Cartesian coordinates.

            vecs = np.einsum("ijk,kl", vecs, str_ref.lattice_vectors)

            neighbour_table = np.linalg.norm(vecs, axis=2)

            for i, idx_map in enumerate(inds_map):
                idx = np.argmin(neighbour_table[i])

                atom_mapping_dict[idx_map] = (
                    inds_ref[idx],
                    neighbour_table[i][idx],
                )
        else:
            for idx in inds_map:
                atom_mapping_dict[idx] = (None, None)

    # Check all atoms have been mapped.

    if len(atom_mapping_dict) != str_map.num_atoms:
        raise Exception(
            "Failed to map all atoms in str_map. If constraints were "
            "supplied, check these include all atom types in str_map. "
            "If no constraints were supplied, this is most likely a "
            "bug."
        )

    idx_refs, dists = [], []

    for idx in range(str_map.num_atoms):
        idx_ref, dist = atom_mapping_dict[idx]

        idx_refs.append(idx_ref)
        dists.append(dist)

    # If constraints is not set, check the mapping is unique.

    if constraints is None and len(idx_refs) != len(set(idx_refs)):
        _exception_or_user_warning(
            "Failed to produce a unique mapping - constraints may be "
            "required.",
            warn=warn,
        )

    # Return index array with dtype=object to preserve None if present.

    return (
        np.array(idx_refs, dtype=object),
        np.array(dists, dtype=np.float64),
    )


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
        prim_trans=None,
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
        prim_trans : array_like, optional
            Transformation matrix from the reference unit cell to the
            primitive cell (optional, default: `None`).
        atom_map : array_like of int or None, optional
            Mapping of the atom positions in the calculation supercell
            to the positions in the reference structure (the default
            value of `None` attempts an automatic mapping).
        ph_calc : str, optional
            Calculator interface used to set up `ph` (default: "vasp").
        symprec : float, optional
            Symmetry tolerance used to set up `ph` (default: 1e-5).

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
        `map_atom_positions()` function in this module.

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
        primitive cell during unfolding is given by:
        `np.dot(sc_mat, np.linalg.inv(prim_trans))`

        This is computed automatically and the result is stored in
        `unfolding_supercell_matrix`.

        The unfolding routine performs phonon calculations and requires
        the force constants. Calling code must therefore set the
        `force_constants` property on `ph` before initialisation.

        The calculation structure is taken from `ph`, and the
        "calculator" used to prepare the phonon calculation must be
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

        sc_mat = np_asarray_copy(sc_mat, dtype=np.float64)

        if not np_check_shape(sc_mat, (3, 3)):
            raise ValueError("sc_mat must be an array_like with shape (3, 3).")

        if prim_trans is not None:
            prim_trans = np_asarray_copy(prim_trans, dtype=np.float64)

            if not np_check_shape(prim_trans, (3, 3)):
                raise ValueError(
                    "If supplied, prim_trans must be an array_like "
                    "with shape (3, 3)."
                )
        else:
            # Identity matrix.

            prim_trans = np.eye(3, dtype=np.float64)

        uf_sc_mat = np.dot(sc_mat, np.linalg.inv(prim_trans))

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
            atom_map = np_asarray_copy(atom_map, dtype=object)

            if len(atom_map) != len(struct.num_atoms):
                raise ValueError(
                    "The number of entries in atom_map does not match "
                    "the number of atoms in the supercell used for the "
                    "phonon calculation."
                )

            for i, idx in enumerate(atom_map):
                if idx is not None:
                    idx = int(idx)

                    if idx < 0 or idx >= ref_struct.num_atoms:
                        raise ValueError(
                            "One or more indices in atom_map are "
                            "inconsistent with the number of atoms in "
                            "ref_struct."
                        )

                    atom_map[i] = idx
        else:
            atom_map, _ = map_atom_positions(struct, ref_struct)

        self._phonopy = ph

        self._struct = struct
        self._ref_struct = ref_struct
        self._prim_struct = prim_struct

        self._prim_trans = prim_trans
        self._sc_mat = sc_mat
        self._uf_sc_mat = uf_sc_mat

        self._atom_map = atom_map

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
        return np_readonly_view(self._sc_mat)

    @property
    def primitive_matrix(self):
        """numpy.ndarray : Primitive transformation matrix (shape:
        `(3, 3)`)."""
        return np_readonly_view(self._prim_trans)

    @property
    def unfolding_supercell_matrix(self):
        """numpy.ndarray : Unfolding supercell matrix (shape: `(3, 3)`)."""
        return np_readonly_view(self._uf_sc_mat)

    @property
    def atom_map(self):
        """numpy.ndarray : Mapping of atoms in `structure` to atoms in
        `reference_structure`."""
        return np_readonly_view(self._atom_map)

    def _run_unfolding(self, qpts):
        """Run the band unfolding for a set of q-points."""

        unfolding = Unfolding(
            self._phonopy,
            self._uf_sc_mat,
            self._ref_struct.atom_positions,
            self._atom_map,
            qpts,
        )

        unfolding.run()

        return (unfolding.frequencies, unfolding.unfolding_weights)

    def unfold_to_q(self, qpts):
        """Unfold to (a) specified q-point(s).

        Parameters
        ----------
        qpts : array_like
            q-point(s) to unfold (shape: `(3,)` or `(M, 3)`).

        Returns
        -------
        freqs_weights : tuple if numpy.ndarray
            Frequencies and unfolding weights (shape: `(3N,)` or
            `(M, 3N)`, where `N` is the number of atoms in
            `structure`.
        """

        qpts, n_dim_add = np_expand_dims(
            np.asarray(qpts, dtype=np.float64), (None, 3)
        )

        freqs, weights = self._run_unfolding(qpts)

        return (
            freqs if n_dim_add == 0 else freqs[0],
            weights if n_dim_add == 0 else weights[0],
        )

    def unfold_band_structure(self, band_path, num_pts=101, var_seg_len=True):
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

        Returns
        -------
        uf_disp : list of tuple of numpy.ndarray
            Tuples of `(q, dists, freqs, weights)` for each segment of
            the unfolded dispersion.
        """

        band_path, _ = np_expand_dims(
            np.asarray(band_path, dtype=np.float64), (None, 2, 3)
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
            freqs, weights = self._run_unfolding(qpts)
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
