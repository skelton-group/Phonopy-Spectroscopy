# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines for working with crystal structures."""


# -------
# Imports
# -------


import glob
import os
import warnings

import numpy as np

from itertools import product

from ..constants import ZERO_TOLERANCE

from ..utility.numpy_helper import np_check_shape, np_expand_dims


# ---------------------
# Coordinate conversion
# ---------------------


def cartesian_to_fractional_coordinates(cart_pos, latt_vecs):
    """Convert positions from Cartesian to fractional coordinates.

    Parameters
    ----------
    cart_pos : array_like
        Atom position or set of positions in Cartesian coordinates
        (shape: `(3,)` or `(N, 3)`).
    latt_vecs : array_like
        Lattice vectors (shape: `(3, 3)`).

    Returns
    -------
    frac_pos : numpy.ndarray
        Atom positions in fractional coordinates (same shape as
        `cart_pos`).
    """

    cart_pos, n_dim_add = np_expand_dims(np.asarray(cart_pos), (None, 3))

    latt_vecs = np.asarray(latt_vecs)

    if not np_check_shape(latt_vecs, (3, 3)):
        raise ValueError("latt_vecs must be an array_like with shape (3, 3).")

    trans_mat = np.linalg.inv(latt_vecs)

    frac_pos = np.zeros_like(cart_pos)

    for i, p in enumerate(cart_pos):
        frac_pos[i] = np.dot(p, trans_mat) % 1.0

    return frac_pos if n_dim_add == 0 else frac_pos[0]


def fractional_to_cartesian_coordinates(frac_pos, latt_vecs):
    """Convert positions from fractional to Cartesian coordinates.

    Parameters
    ----------
    frac_pos : array_like
        Atom position or set of positions in fractional coordinates
        (shape: `(3,)` or `(N, 3)`).
    latt_vecs : array_like
        Lattice vectors (shape: `(3, 3)`).

    Returns
    -------
    cart_pos : numpy.ndarray
        Atom positions in Cartesian coordinates (same shape as
        `frac_pos`).
    """

    frac_pos, n_dim_add = np_expand_dims(np.asarray(frac_pos), (None, 3))

    latt_vecs = np.asarray(latt_vecs)

    if not np_check_shape(latt_vecs, (3, 3)):
        raise ValueError("latt_vecs must be an array_like with shape (3, 3).")

    cart_pos = np.zeros_like(frac_pos)

    for i, p in enumerate(frac_pos):
        cart_pos[i] = np.dot(p, latt_vecs)

    return cart_pos if n_dim_add == 0 else cart_pos[0]


# ---------------------
# Distance calculations
# ---------------------


def calculate_distances_frac(pos, latt_vecs, other_pos=None, ret_vecs=False):
    """Calculate distances between sets of positions in fractional
    coordinates.

    Parameters
    ----------
    pos : array_like
        Positions (shape: `(3,)` or `(N, 3)`).
    latt_vecs : array_like
        Lattice vectors for converting positions to Cartesian
        coordinates (shape: `(3, 3)`).
    other_pos : array_like, optional
        Optional second set of positions (shape: `(3,)` or `(M, 3)`,
        default: `None` implies setting to `pos`).
    ret_vecs : bool, optional
        Return vectors instead of distances (default: `False`).

    Returns
    -------
    dists_or_vecs : numpy.ndarray
        Distances (shape: `(N, N)` or `(N, M)`) or vectors (shape:
        `(N, N, 3)` or `(N, M, 3)`).

    Notes
    -----
    This function can be used in several ways:

    With the default `other_pos=None`, the distances are calculated
    between each pair of positions in `pos`, and the returned distances
    array has shape `(N, N)`.

    `other_pos` can instead be used to specify a second set of
    positions. In this case, the calculation is performed between each
    pair of positions in `pos` and `other_pos`, and the returned
    distances array has shape `(N, M)`.

    In either usage, setting `ret_vecs=True` returns the vectors instead
    of distances between positions (shapes: `(N, N, 3)` or `(N, M, 3)`).
    """

    pos, _ = np_expand_dims(np.asarray(pos, dtype=np.float64), (None, 3))

    latt_vecs = np.asarray(latt_vecs, dtype=np.float64)

    if not np_check_shape(latt_vecs, (3, 3)):
        raise ValueError("latt_vecs must be an array_like with shape (3, 3).")

    if other_pos is not None:
        other_pos, _ = np_expand_dims(
            np.asarray(other_pos, dtype=np.float64), (None, 3)
        )
    else:
        other_pos = pos

    if (np.abs(pos) > 1.0).any() or (
        other_pos is not pos and (np.abs(other_pos) > 1.0).any()
    ):
        warnings.warn(
            "One or more positions are outside the range [-1, 1] "
            "expected for fractional coordinates.",
            UserWarning,
        )

    vecs_frac = pos[:, np.newaxis, :] - other_pos[np.newaxis, :, :]

    # Apply periodic boundary conditions.

    vecs_frac[vecs_frac < -0.5] += 1.0
    vecs_frac[vecs_frac >= 0.5] -= 1.0

    # Convert fractional to Cartesian coordinates.

    vecs = np.einsum("ijk,kl", vecs_frac, latt_vecs)

    return vecs if ret_vecs else np.linalg.norm(vecs, axis=2)


# -----------------
# Structure mapping
# -----------------


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

    Automatic constraints can be applied with `type_constraints="auto"`.
    In this case, atoms in `map_struct` are set to map to the same type
    in `ref_struct`, if present, and `None` otherwise. This covers most
    common usages of constraints.

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
    * Interstitials: These can be handled using `max_dist` to
      map the interstitials to `None`, or, if the
      interstitial is of a different type to the atoms in `ref_struct`,
      using `type_constraints="auto"`.
    * "Many-to-one" mapping: An example of this is mapping the organic
      cation in (CH3NH3)PbI3 to the inorganic cation in CsPbI3. This can
      be achieved with a combination of `type_constraints` and
      `allow_non_unique=True`. Note that `type_constraints="auto"` will
      not work in this case.

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

    # If type_constraints is set to "auto", each atom type in map_struct
    # is mapped to the same type in ref_stuct, if present, or None if
    # the type is not present.

    try:
        if str(type_constraints).lower() == "auto":
            type_constraints = []

            for sym in np.unique(map_struct.atom_types):
                if sym in ref_struct.atom_types:
                    type_constraints.append(([sym], [sym]))
                else:
                    type_constraints.append(([sym], None))
    except ValueError:
        pass

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

            neighbour_table = calculate_distances_frac(
                pos_map, ref_struct.lattice_vectors, pos_ref, ret_vecs=False
            )

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
            if dist is not None and dist > max_dist:
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


# ----------------------
# Brillouin zone mapping
# ----------------------


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

    map_qpts = ((map_qpts + 0.5) % 1.0) - 0.5

    return map_qpts if n_dim_add == 0 else map_qpts[0]
