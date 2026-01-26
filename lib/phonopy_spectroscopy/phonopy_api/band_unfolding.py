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


# -------
# Mapping
# -------


def map_atom_positions(
    map_struct,
    ref_struct,
    max_dist=None,
    type_constraints=None,
    len_tol=0.1,
    ang_tol=1.0,
    allow_non_unique=False,
):
    """Generate an integer mapping of the closest atomic positions in
    two `Structure` objects.

    Parameters
    ----------
    map_struct, ref_struct : Structure
        Structures to map.
    max_dist : float or None, optional
        If set, atoms in `map_struct` that are more than `max_dist` from
        all atoms in `ref_struct` are assumed to have no mapping.
    type_constraints : sequence or None, optional
       If set, specifies pairwise sets of atom types in `map_struct` and
       `ref_struct` to constrain matching (default: `None`).
    allow_non_unique: bool, optional
        Allow multiple atoms in `map_struct` to map to the same atom in
        `ref_struct`.
    tol: float, optional
        Specifies the maximum allowed differences in the metric tensors
        of `map_struct` and `ref_struct` (default: `tol=1.0e-5`).

    Returns
    -------
    map : tuple of list
        Integer mapping of the atom positions in `str_map` and distances
        to the closest positions in `str_ref` (both set to `None` for
        unmapped atoms).

    Notes
    -----
    For a mapping to be physically sound the unit cells of `map_struct`
    and `ref_struct` should be similar. This is checked by comparing the
    (absolute) differences in the lengths of the unit-cell vectors and
    the angles between them to `len_tol` and `ang_tol`.

    The default mapping procedure is designed to give sensible results
    for "simple" cases, including e.g. atomic substitutions and
    vacancies, provided a unique mapping between the atom positions in
    `map_struct` and `ref_struct` can be obtained. This implies that the
    distortions from the "ideal" geometry must be relatively small.

    For more complex cases, the `max_dist`, `type_constraints` and
    `allow_unique` keywords provide more control.

    If `max_dist` is set, atoms in `map_struct` that are more than
    `max_dist` from all atoms in `ref_struct` are mapped to `None`.

    `type_constraints` controls which atom types are mapped to one
    another, and can also be used to specify that some atom types
    in `map_struct` should be mapped to `None`. Constraints are
    specified as pairwise sets of atomic symbols, e.g.
    `type_constraints=[(["Ce"], ["Ce"]), (["O"], ["O"])]`. The second
    "set" in each pair can be subsituted by `None` to specify that the
    atom type(s) in `map_struct` should not be mapped - e.g.
    `type_constraints=[..., (["H"], None)]`. Note that if set
    `type_constraints` must cover all atom types in `map_struct`.

    Finally, non-unique mapping, where multiple atoms in `map_struct`
    map to the same atom in `ref_struct`, can be explicitly allowed by
    setting `allow_non_unique=True`.

    The following are some general comments on how to treat common
    scenarios:

    * Vacancies and atomic substitutions (including e.g. alloys): The
      default behaviour should work in most cases, but
      `type_constraints` may be required if there are significant
      structural distortions.
    * Frenkel defects: If the interstitial is close to the vacancy, it
      may be reasonable to map the interstitial to the (occupied)
      lattice site in `ref_struct`, and the default setup may work. If
      the interstitial is far from the vacancy, it is likely more
      reasonable to map it to `None`, which might be achieved by setting
      an appropriate `max_dist`.
    * Interstitials: These can possibly be handled using `max_dist` to
      map the interstitials to `None`. Alternatively, if the
      interstitial is of a different type to the atoms in `ref_struct`,
      it could be mapped to `None` by setting appropriate
      `type_constraints`.
    * "Many-to-one" mapping: An example of this is mapping the organic
      cation in (CH3NH3)PbI3 to the inorganic cation in CsPbI3. This can
      be achieved with a combination of `type_constraints` and
      `allow_non_unique=True`.

    Finally, while the default behaviour is designed to require user
    intervention if any potential issues are found, we recommend always
    verifying that the mapping is "sane".
    """

    # Check the similarity of the structures by comparing the largest
    # difference in the cell lengths and the largest angles between
    # lattie vectors to the set tolerances.

    if len_tol <= 0.0:
        raise ValueError("len_tol must be > 0.")

    if ang_tol <= 0.0 or ang_tol > 180.0:
        raise ValueError("ang_tol must be > 0 and <= 180.")

    map_norms = np.linalg.norm(map_struct.lattice_vectors, axis=1)
    ref_norms = np.linalg.norm(ref_struct.lattice_vectors, axis=1)

    norm_diff = np.abs(map_norms - ref_norms)

    if (norm_diff > len_tol).any():
        raise Exception(
            "Maximum difference in lattice vector lengths is {0:.3f} "
            "> len_tol = {1:.3f}.".format(norm_diff.max(), len_tol)
        )

    thetas = []

    for idx in range(3):
        dp = np.dot(
            map_struct.lattice_vectors[idx], ref_struct.lattice_vectors[idx]
        )

        cos_theta = np.clip(dp / (map_norms[idx] * ref_norms[idx]), -1.0, 1.0)
        thetas.append(np.arccos(cos_theta))

    thetas = np.abs(np.degrees(thetas))

    if (thetas > ang_tol).any():
        raise Exception(
            "Largest angle between lattice vectors is {0:.2f} > "
            "{1:.2f} deg.".format(thetas.max(), ang_tol)
        )

    # Perform the mapping on pairs of groups of indices at a time - this
    # allows constraints on atom types to be applied.

    index_grps = []

    if type_constraints is not None:
        for map_typs, ref_typs in type_constraints:
            if len(map_typs) == 0:
                raise ValueError(
                    "Mapping atom types in type_constraints must be an "
                    "array_like with at least one element."
                )

            inds_map = []

            for sym in map_typs:
                (inds,) = np.where(map_struct.atom_types == sym)
                inds_map.extend(inds)

            ref_inds = None

            if ref_typs is not None:
                if len(ref_typs) == 0:
                    raise ValueError(
                        "Reference atom types in type_constraints must "
                        "either be an array_like with at least one "
                        "element or None."
                    )

                ref_inds = []

                for sym in ref_typs:
                    (inds,) = np.where(ref_struct.atom_types == sym)
                    ref_inds.extend(inds)

            index_grps.append((inds_map, ref_inds))
    else:
        index_grps = [
            (
                np.arange(0, map_struct.num_atoms, dtype=int),
                np.arange(0, ref_struct.num_atoms, dtype=int),
            )
        ]

    # Adjust for differences in unit cells by converting the positions
    # into Cartesian coordinates and then back into fractional
    # coordinates using the lattice vectors of the reference structure.

    map_pos_cart = fractional_to_cartesian_coordinates(
        map_struct.atom_positions, map_struct.lattice_vectors
    )

    map_pos_frac_shift = cartesian_to_fractional_coordinates(
        map_pos_cart, ref_struct.lattice_vectors
    )

    # Perform mapping.

    atom_mapping_dict = {}

    for inds_map, ref_inds in index_grps:
        if ref_inds is not None:
            pos_map = map_pos_frac_shift[inds_map]
            pos_ref = ref_struct.atom_positions[ref_inds]

            vecs = pos_map[:, np.newaxis, :] - pos_ref[np.newaxis, :, :]

            # Apply periodic boundary conditions.

            vecs[vecs < -0.5] += 1.0
            vecs[vecs >= 0.5] -= 1.0

            # Convert fractional to Cartesian coordinates.

            vecs = np.einsum("ijk,kl", vecs, ref_struct.lattice_vectors)

            neighbour_table = np.linalg.norm(vecs, axis=2)

            for i, idx_map in enumerate(inds_map):
                idx = np.argmin(neighbour_table[i])

                atom_mapping_dict[idx_map] = (
                    ref_inds[idx],
                    neighbour_table[i][idx],
                )
        else:
            for idx in inds_map:
                atom_mapping_dict[idx] = (None, None)

    # If max_dist is set, check distances and set entries where the
    # distance exceeds this to None.

    if max_dist is not None:
        if max_dist <= 0.0:
            raise ValueError("If set, max_dist must be > 0.")

        for idx, (idx_ref, dist) in atom_mapping_dict.items():
            if dist > max_dist:
                atom_mapping_dict[idx] = (None, None)

    # Check all atoms have been mapped.

    if len(atom_mapping_dict) != map_struct.num_atoms:
        raise Exception(
            "Failed to map all atoms in map_struct. If type_constraints"
            "was set, check these include all atom types in "
            "map_struct. If type_constraints was not set, this is most "
            "likely a bug."
        )

    idx_refs, dists = [], []

    for idx in range(map_struct.num_atoms):
        idx_ref, dist = atom_mapping_dict[idx]

        idx_refs.append(idx_ref)
        dists.append(dist)

    # If allow_non_unique is not set, check the mapping is unique.

    if not allow_non_unique is None:
        temp = [idx for idx in idx_refs if idx is not None]

        if len(temp) != len(set(temp)):
            raise Exception(
                "Failed to produce a unique mapping. If a non-unique "
                "mapping is physical this can be explicitly allowed by "
                "setting allow_non_unique=True)."
            )

    # Use dtype=object to ensure None is preserved.

    return (np.array(idx_refs, dtype=object), np.array(dists, dtype=object))


def centred_modulo(a):
    """Perform a "centred modulo" to map values to the range
    [-0.5, 0.5].

    Parameters
    ----------
    a : array_like
        Values to map.

    Returns
    -------
    a_mod : numpy.ndarray
        Mapped values.
    """

    return ((np.asarray(a, dtype=np.float64) + 0.5) % 1.0) - 0.5


def map_qpoints(qpts, ref_struct, map_struct):
    """Map "reduced" q-point(s) defined in the Brillouin zone of a
    reference structure to the Brillouin zone of another structure.

    Parameters
    ----------
    qpts : array_like
        Reduced q-point(s) to map (shape: `(3,)` or `(N, 3)`).
    ref_struct, map_struct : Structure
        Reference structure for which `qpts` are specified and structure
        to map to.

    Returns
    -------
    qpts_map : numpy.ndarray
        Fractional q-point(s) in the Brillouin zone of `map_struct`
        (same shape as `qpts`).
    """

    qpts, n_dim_add = np_expand_dims(
        np.asarray(qpts, dtype=np.float64), (None, 3)
    )

    # No need to include the factor of 2 \pi in the reciprocal lattice
    # vectors as long as we're consistent.

    ref_rec_v_latt = ref_struct.reciprocal_lattice_vectors(two_pi=False)
    map_rec_v_latt = map_struct.reciprocal_lattice_vectors(two_pi=False)

    # Convert the q-points to Cartesian coordinates in the Brillouin
    # zone of the reference structure.

    ref_qpts_cart = fractional_to_cartesian_coordinates(qpts, ref_rec_v_latt)

    # Determine the transformation between the reciprcal lattices of the
    # reference and map structures and rotate the q-points.

    bz_trans_mat = np.dot(np.linalg.inv(ref_rec_v_latt), map_rec_v_latt)

    map_qpts_cart = np.array(
        [np.dot(q, bz_trans_mat) for q in ref_qpts_cart], dtype=np.float64
    )

    # Convert the rotated q-point coordinates back to fractional
    # coordinates in the Brillouin zone of the map structure and apply a
    # centred modulo.

    map_qpts = cartesian_to_fractional_coordinates(
        map_qpts_cart, map_rec_v_latt
    )

    map_qpts = centred_modulo(map_qpts)

    return map_qpts if n_dim_add == 0 else map_qpts[0]


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

            if len(atom_map) != struct.num_atoms:
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

    def _run_unfolding(self, qpts, ref_prim=None, tol=1.0e-5):
        """Run the band unfolding for a set of q-points with optional
        Brillouin zone mapping."""

        if ref_prim is not None:
            qpts = map_qpoints(qpts, ref_prim, self._prim_struct, tol=tol)

        unfolding = Unfolding(
            self._phonopy,
            self._uf_sc_mat,
            self._ref_struct.atom_positions,
            self._atom_map,
            qpts,
        )

        unfolding.run()

        return (unfolding.frequencies, unfolding.unfolding_weights)

    def unfold_to_q(self, qpts, ref_prim=None, tol=1.0e-5):
        """Unfold to (a) specified q-point(s).

        Parameters
        ----------
        qpts : array_like
            q-point(s) to unfold (shape: `(3,)` or `(M, 3)`).
        ref_prim : Structure, optional
            Reference primitive structure for which the q-point
            coordinates are defined (default: `None`).
        tol : float, optional
            Tolerance for checking the calculation and reference
            primitive cells are consistent (default: 1.0e-5).

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

    def unfold_band_structure(
        self,
        band_path,
        num_pts=101,
        var_seg_len=True,
        ref_prim=None,
        tol=1e-5,
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
        tol : float, optional
            Tolerance for checking the calculation and reference
            primitive cells are consistent (default: 1.0e-5).

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
            freqs, weights = self._run_unfolding(
                qpts, ref_prim=ref_prim, tol=tol
            )

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
