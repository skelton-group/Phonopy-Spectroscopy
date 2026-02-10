# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines for working with crystal structures."""


# -------
# Imports
# -------


import warnings

import numpy as np

from ..utility.numpy_helper import np_check_shape, np_expand_dims


# ---------------
# Helper routines
# ---------------


def _in_place_centred_modulo(a):
    """Apply an in-place "centred modulo" to fractional coordinates in a
    `numpy.ndarray` object `a`."""

    if not isinstance(a, np.ndarray):
        raise TypeError(
            "a must be a numpy.ndarray (this is most likely a bug)."
        )

    a -= np.round(a)

    return a


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

    frac_pos = np.matmul(cart_pos, np.linalg.inv(latt_vecs))
    frac_pos %= 1.0

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

    cart_pos = np.matmul(frac_pos, latt_vecs)

    return cart_pos if n_dim_add == 0 else cart_pos[0]


# --------------------
# Distance calculation
# --------------------


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

    vecs = pos[:, np.newaxis, :] - other_pos[np.newaxis, :, :]
    vecs = _in_place_centred_modulo(vecs)

    # Apply periodic boundary conditions.

    # vecs[vecs < -0.5] += 1.0
    # vecs[vecs >= 0.5] -= 1.0

    # Convert fractional to Cartesian coordinates.

    vecs = np.matmul(vecs, latt_vecs, out=vecs)

    return vecs if ret_vecs else np.linalg.norm(vecs, axis=2)


# -------------
# Atom grouping
# -------------


def group_atoms(struct, bond_dists=None, default_dist=1.6, atom_inds=None):
    """Group atoms in a structure e.g. into molecules by tracing a
    bonding network.

    Parameters
    ----------
    struct : Structure
        Structure.
    bond_dists : sequence, optional
        Bond distances specified as `(sym, dist)` pairs for all bonds to
        atoms of type `sym`) or `(sym1, sym2, dist)` triples for bonds
        between atoms of type `sym1` and `sym2`).
    default_dist : float, optional
        Default bond distance for bonds not specified in `bond_dists`
        (default: 1.6 Ang).
    atom_inds : array_like or None, optional
        Optionally specify a subset of atoms to "start" groups from -
        other atoms will automatically be included if part of the
        bonding network (default: `None`, group all atoms).

    Returns
    -------
    mol_grp_inds : list of numpy.ndarray
        Indices of atoms in each group.
    """

    if default_dist <= 0.0:
        raise ValueError("default_dist must be > 0.")

    at_typs = struct.atom_types
    at_syms, inv_inds = np.unique(at_typs, return_inverse=True)

    # Given M unique atom types, set up an MxM matrix where each element
    # corresponds to the maximum distance between each pair of types.

    bond_dist_mat = np.zeros((len(at_syms),) * 2, dtype=np.float64)

    if bond_dists is not None:
        # Build a temporary lookup table and fill unassigned distances
        # with default_dist.

        bond_dist_lut = {
            sym: {sym: default_dist for sym in at_syms} for sym in at_syms
        }

        if bond_dists is not None:
            for bond_dist in bond_dists:
                if len(bond_dist) == 2:
                    # (sym, dist) : Set max dist for all bonds to an
                    # atom.

                    sym, dist = bond_dist

                    if sym in at_syms:
                        bond_dist_lut[sym] = {sym: dist for sym in at_syms}

                elif len(bond_dist) == 3:
                    # (sym1, sym2, dist) : Set max dist for a specific
                    # atom pair.

                    sym1, sym2, dist = bond_dist

                    if sym1 in at_syms and sym2 in at_syms:
                        bond_dist_lut[sym1][sym2] = dist
                        bond_dist_lut[sym2][sym1] = dist
                else:
                    raise ValueError(
                        "Entries in bond_dists must specify either "
                        "(sym, dist) pairs or (sym1, sym2, dist) "
                        "triples."
                    )

        # Use lookup table to fill bond-distance matrix.

        for i, sym1 in enumerate(at_syms):
            for j, sym2 in enumerate(at_syms):
                bond_dist_mat[i, j] = bond_dist_lut[sym1][sym2]

    else:
        # If individual bond distances are not specified, set all max
        # dists to default_dist.

        bond_dist_mat[:, :] = default_dist

    if atom_inds is not None:
        for idx in atom_inds:
            if idx < 0 or idx >= struct.num_atoms:
                raise ValueError(
                    "One or more indices in atom_inds are inconsistent "
                    "with the number of atoms in struct."
                )
    else:
        atom_inds = [i for i in range(struct.num_atoms)]

    # Generate a neighbour distance table and use the bond-distance
    # matrix to create a Boolean "bond table".

    dist_table = calculate_distances_frac(
        struct.atom_positions, struct.lattice_vectors, ret_vecs=False
    )

    bond_table = dist_table <= bond_dist_mat[np.ix_(inv_inds, inv_inds)]

    # Group atoms into molecules.

    at_grp_inds = []

    # Keep track of which atoms we've already assigned to molecule
    # groups.

    assigned_inds = set()

    for i in atom_inds:
        if i not in assigned_inds:
            grp_inds = [i]

            while True:
                grp_inds_new = set(grp_inds)

                for idx1 in grp_inds:
                    (inds2,) = np.where(bond_table[idx1])
                    grp_inds_new.update(inds2)

                if len(grp_inds_new) == len(grp_inds):
                    break

                grp_inds = list(grp_inds_new)

            at_grp_inds.append(grp_inds)
            assigned_inds.update(grp_inds)

    return [np.array(grp_inds, dtype=int) for grp_inds in at_grp_inds]


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
    len_tol, ang_tol : float, optional
        Maximum allowed differences in the lengths and angles between
        the lattice vectors of `map_struct` and `ref_struct` (default:
        0.1 Ang, 1 deg).
    allow_non_unique: bool, optional
        Allow multiple atoms in `map_struct` to map to the same atom in
        `ref_struct`.

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

    if not allow_non_unique:
        temp = [idx for idx in idx_refs if idx is not None]

        if len(temp) != len(set(temp)):
            raise Exception(
                "Failed to produce a unique mapping. If a non-unique "
                "mapping is physical this can be explicitly allowed by "
                "setting allow_non_unique=True)."
            )

    # Use dtype=object to ensure None is preserved.

    return (np.array(idx_refs, dtype=object), np.array(dists, dtype=object))


def invert_atom_map(atom_map, map_struct, ref_struct):
    """Invert an integer mapping produced by `map_atom_positions`.

    Parameters
    ----------
    atom_map : array_like
        Atom map.
    map_struct, ref_struct : Structure
        Mapped and reference structures.

    Returns
    -------
    inv_map : numpy.ndarray
        Integer mapping of the atoms in `ref_struct` to atoms in
        `map_struct`.

    Notes
    -----
    Entries can be a single integer (atom in `ref_struct` maps to a
    single atom in `map_struct`), `None` (atom does not map to anything
    in `map_struct`), or a `numpy.ndarray` of integers (atom maps to
    multiple atoms in `map_struct`).
    """

    if len(atom_map) != map_struct.num_atoms:
        raise ValueError(
            "The number of entries in atom_map does not match the "
            "number of atoms in map_struct."
        )

    for i, idx in enumerate(atom_map):
        if idx is not None:
            if idx < 0 or idx >= ref_struct.num_atoms:
                raise ValueError(
                    "One or more indices in atom_map are "
                    "inconsistent with the number of atoms in "
                    "ref_struct."
                )

    inv_atom_map = [[] for _ in range(ref_struct.num_atoms)]

    for idx, ref_idx in enumerate(atom_map):
        if ref_idx is not None:
            inv_atom_map[ref_idx].append(idx)

    for i, inds in enumerate(inv_atom_map):
        if len(inds) == 0:
            inv_atom_map[i] = None
        elif len(inds) == 1:
            inv_atom_map[i] = inds[0]
        else:
            inv_atom_map[i] = np.array(inds, dtype=int)

    return np.array(inv_atom_map, dtype=object)


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
    # reference and map structures and rotate the q-points. Note that we
    # need to work with "normal" column-matrix algebra here.

    trans_mat = np.matmul(map_rec_v_latt.T, np.linalg.inv(ref_rec_v_latt.T))

    # Rotate the Cartesian q-point coordinates, convert to fractional
    # coordinates in the Brillouin zone of map_struct, and apply a
    # centred modulo.

    # The double transpose here is required to converts the q-point
    # coordinates to column vectors and back.

    map_qpts_cart = np.matmul(trans_mat, ref_qpts_cart.T).T

    map_qpts = cartesian_to_fractional_coordinates(
        map_qpts_cart, map_rec_v_latt
    )

    map_qpts = _in_place_centred_modulo(map_qpts)

    diff = map_qpts - qpts
    diff = diff - np.rint(diff)
    assert np.allclose(diff, 0.0)

    return map_qpts if n_dim_add == 0 else map_qpts[0]
