# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for working with crystal structures."""

# -------
# Imports
# -------

import numpy as np

# ---------------------
# Coordinate conversion
# ---------------------


def cartesian_to_fractional_coordinates(cart_pos, latt_vecs):
    """Convert positions from Cartesian to fractional coordinates.

    Parameters
    ----------
    cart_pos : array_like
        Atom position or set of positions in Cartesian coordinates
        (shape: `(..., 3)`).
    latt_vecs : array_like
        Lattice vectors (row-major, shape: `(3, 3)`).

    Returns
    -------
    frac_pos : numpy.ndarray
        Atom positions in fractional coordinates (same shape as
        `cart_pos`).
    """

    cart_pos = np.atleast_1d(np.asanyarray(cart_pos))

    if cart_pos.shape[-1] != 3:
        raise ValueError("cart_pos must be an array_like with shape (..., 3).")

    latt_vecs = np.asanyarray(latt_vecs)

    if latt_vecs.shape != (3, 3):
        raise ValueError("latt_vecs must be an array_like with shape (3, 3).")

    return (cart_pos @ np.linalg.inv(latt_vecs)) % 1.0


def fractional_to_cartesian_coordinates(frac_pos, latt_vecs):
    """Convert positions from fractional to Cartesian coordinates.

    Parameters
    ----------
    frac_pos : array_like
        Atom position or set of positions in fractional coordinates
        (shape: `(..., 3)`).
    latt_vecs : array_like
        Lattice vectors (row-major, shape: `(3, 3)`).

    Returns
    -------
    cart_pos : numpy.ndarray
        Atom positions in Cartesian coordinates (same shape as
        `frac_pos`).
    """

    frac_pos = np.atleast_1d(np.asanyarray(frac_pos))

    if frac_pos.shape[-1] != 3:
        raise ValueError("frac_pos must be an array_like with shape (..., 3).")

    latt_vecs = np.asanyarray(latt_vecs)

    if latt_vecs.shape != (3, 3):
        raise ValueError("latt_vecs must be an array_like with shape (3, 3).")

    return frac_pos @ latt_vecs


# --------------------
# Distance calculation
# --------------------


def calculate_distances_from_fractional_positions(
    pos, latt_vecs, other_pos=None, ret_vecs=False
):
    """Calculate distances between sets of positions in fractional
    coordinates.

    Parameters
    ----------
    pos : array_like
        Positions (shape: `(N, 3)`).
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
        Distances (shape: `(N, M)`) or vectors (shape: `(N, M, 3)`).

    Notes
    -----
    The default `other_pos=None` calculates the distance between each
    pair of positions in `pos` and the return array has shape `(N, N)`.

    If `other_pos` is specified, the distances between each pair of
    positions in `pos` and `other_pos` (shape: `(N, M)`).

    In either case, setting `ret_vecs=True` returns the vectors instead
    of distances (shapes: `(N, N, 3)` or `(N, M, 3)`).
    """

    pos = np.atleast_2d(np.asanyarray(pos))

    if other_pos is not None:
        other_pos = np.atleast_2d(np.asanyarray(other_pos))
    else:
        other_pos = pos

    if (pos.ndim != 2 or pos.shape[-1] != 3) or (
        other_pos.ndim != 2 or other_pos.shape[-1] != 3
    ):
        raise ValueError(
            "pos and other_pos (if supplied) must be array_like with "
            "shape (N, 3) and (M, 3)."
        )

    latt_vecs = np.asanyarray(latt_vecs, dtype=np.float64)

    if latt_vecs.shape != (3, 3):
        raise ValueError("latt_vecs must be an array_like with shape (3, 3).")

    vecs = pos[:, np.newaxis, :] - other_pos[np.newaxis, :, :]
    vecs -= np.round(vecs)

    # Convert fractional to Cartesian coordinates.

    vecs = vecs @ latt_vecs

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
    grp_inds : list of numpy.ndarray
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

    dist_table = calculate_distances_from_fractional_positions(
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
    type_const=None,
    len_tol=0.1,
    ang_tol=1.0,
    allow_multi=False,
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
    type_const : sequence or None, optional
       Specifies pairwise sets of atom types in `map_struct` and
       `ref_struct` to constrain matching (default: `None`).
    len_tol, ang_tol : float, optional
        Maximum allowed differences in the lengths and angles between
        the lattice vectors of `map_struct` and `ref_struct` (default:
        0.1 Ang, 1 deg).
    allow_multi: bool, optional
        Allow multiple atoms in `map_struct` to map to the same atom in
        `ref_struct`.

    Returns
    -------
    map : tuple of list
        Integer mapping of the atom positions in `str_map` and distances
        to the closest positions in `str_ref` (both `None` for unmapped
        atoms).

    Notes
    -----
    A physically sound mapping requires the unit cells of `map_struct`
    and `ref_struct` to be similar. This is checked by comparing the
    differences in the cell lengths and angles against `len_tol` and
    `ang_tol`.

    The default parameters should give sensible results in "simple"
    cases, including e.g. atomic substitutions and vacancies, provided
    the distortions from the "ideal" geometry are relatively small.

    The `max_dist`, `type_const` and `allow_multi` keywords provide more
    control for complex cases.

    If `max_dist` is set, atoms in `map_struct` that are more than
    `max_dist` from all atoms in `ref_struct` are mapped to `None`.

    `type_const` controls which atom types are mapped and can be used to
    specify that some atom types in `map_struct` should be mapped to
    `None`. Constraints are specified as pairwise sets of atomic
    symbols, e.g. `type_const=[(["Ce"], ["Ce"]), (["O"], ["O"])]`.
    Substututing the second set in a pair with `None` indicates that the
    atom type(s) in `map_struct` should not be mapped - e.g.
    `type_const=[..., (["H"], None)]`.

    `type_const` can also be set the the string `"auto"`, which sets
    atoms in `map_struct` to map to the same type in `ref_struct`, if
    present, or `None` otherwise.

    Non-unique mapping, where multiple atoms in `map_struct` map to the
    same atom in `ref_struct`, is allowed by setting `allow_multi=True`.

    The following are suggestions for handling common scenarios:

    * Vacancies and atomic substitutions (including e.g. alloys): The
      default behaviour should work in most cases, but `type_const` may
      be required for significant distortions.
    * Frenkel defects: If the interstitial is close to the vacancy and
      can be mapped to the (occupied) site in `ref_struct`, the defaults
      may work. If the interstitial is far from the vacancy, it is more
      appropriate to map it to `None` using `max_dist`.
    * Interstitials: Use `max_dist` to  map the interstitials to `None`,
      or, if the interstitial is of a different type to the atoms in
      `ref_struct`, set `type_const="auto"`.
    * "Many-to-one" mapping: Manually set `type_const` and
      set `allow_multi=True` (e.g. mapping the organic
      cation in (CH3NH3)PbI3 to the inorganic cation in CsPbI3.)

    The default behaviour is designed to force intervention if potential
    issues are found, but it is strongly recommended to check the
    mapping is "sane".
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

    norm_diffs = np.abs(map_norms - ref_norms)

    if (norm_diffs > len_tol).any():
        d_max = norm_diffs.max()

        raise RuntimeError(
            "Maximum difference in lattice vector lengths is "
            f"{d_max:.3f} > len_tol = {len_tol:.3f}."
        )

    angles = []

    for idx in range(3):
        dp = np.dot(
            map_struct.lattice_vectors[idx], ref_struct.lattice_vectors[idx]
        )

        cos_theta = np.clip(dp / (map_norms[idx] * ref_norms[idx]), -1.0, 1.0)
        angles.append(np.arccos(cos_theta))

    angles = np.abs(np.degrees(angles))

    if (angles > ang_tol).any():
        a_max = angles.max()

        raise RuntimeError(
            f"Largest angle between lattice vectors is {a_max:.2f} > "
            f"{ang_tol:.2f} deg."
        )

    # If type_constraints is set to "auto", each atom type in map_struct
    # is mapped to the same type in ref_stuct, if present, or None if
    # the type is not present.

    if isinstance(type_const, (str, np.str_)):
        type_const_str = str(type_const).lower()

        if type_const_str == "auto":
            type_const = []

            for sym in np.unique(map_struct.atom_types):
                ref_typs = [sym] if sym in ref_struct.atom_types else None
                type_const.append(([sym], ref_typs))
        else:
            raise ValueError(
                f'Unknown string specifier type_const = "{type_const}".'
            )

    # Perform the mapping on pairs of groups of indices at a time - this
    # allows constraints on atom types to be applied.

    index_grps = []

    if type_const is not None:
        for map_typs, ref_typs in type_const:
            if len(map_typs) == 0:
                raise ValueError(
                    "Mapping atom types in type_const must be an "
                    "array_like with at least one element."
                )

            map_inds = np.flatnonzero(np.isin(map_struct.atom_types, map_typs))

            ref_inds = None

            if ref_typs is not None:
                if len(ref_typs) == 0:
                    raise ValueError(
                        "Reference atom types in type_const must be an"
                        "array_like with at least one element or None."
                    )

                ref_inds = np.flatnonzero(
                    np.isin(ref_struct.atom_types, ref_typs)
                )

            index_grps.append((map_inds, ref_inds))
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

    for map_inds, ref_inds in index_grps:
        if ref_inds is not None:
            map_pos = map_pos_frac_shift[map_inds]
            ref_pos = ref_struct.atom_positions[ref_inds]

            neighbour_table = calculate_distances_from_fractional_positions(
                map_pos, ref_struct.lattice_vectors, ref_pos, ret_vecs=False
            )

            inds = np.argmin(neighbour_table, axis=1)
            dists = neighbour_table[np.arange(len(map_inds)), inds]

            for map_idx, ref_idx, dist in zip(map_inds, ref_inds[inds], dists):
                atom_mapping_dict[map_idx] = (ref_idx, dist)
        else:
            for map_idx in map_inds:
                atom_mapping_dict[map_idx] = (None, None)

    # If max_dist is set, check distances and set entries where the
    # distance exceeds this to None.

    if max_dist is not None:
        if max_dist <= 0.0:
            raise ValueError("If set, max_dist must be > 0.")

        for idx, (_, dist) in atom_mapping_dict.items():
            if dist is not None and dist > max_dist:
                atom_mapping_dict[idx] = (None, None)

    # Check all atoms have been mapped.

    if len(atom_mapping_dict) != map_struct.num_atoms:
        raise RuntimeError(
            "Failed to map all atoms in map_struct. If type_const is "
            "set, ensure the constraints include all atom types in "
            "map_struct."
        )

    ref_inds, dists = [], []

    for idx in range(map_struct.num_atoms):
        ref_idx, dist = atom_mapping_dict[idx]

        ref_inds.append(ref_idx)
        dists.append(dist)

    # If allow_multi is not set, check the mapping is unique.

    if not allow_multi:
        temp = [idx for idx in ref_inds if idx is not None]

        if len(temp) != len(set(temp)):
            raise Exception(
                "Failed to obtain a unique mapping. If this is "
                "physical it can be explicitly allowed using "
                "allow_multi=True."
            )

    # Use dtype=object to ensure None is preserved.

    return (np.array(ref_inds, dtype=object), np.array(dists, dtype=object))


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
    Entries in `inv_map` can be single integers (atom in `ref_struct`
    maps to a single atom in `map_struct`), `None` (atom does not map to
    anything) or a sequence of integers (atom maps to multiple atoms`).
    """

    if len(atom_map) != map_struct.num_atoms:
        raise ValueError(
            "Length of atom_map does not match number of atoms in map_struct."
        )

    for idx in atom_map:
        if idx is not None:
            if idx < 0 or idx >= ref_struct.num_atoms:
                raise ValueError(
                    "One or more indices in atom_map are out of range "
                    "for the number of atoms in ref_struct."
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
        Reduced q-point(s) to map (shape: `(3,)` or `(N, 3)`.
    ref_struct, map_struct : Structure
        Reference structure for which `qpts` are specified and structure
        to map to.

    Returns
    -------
    qpts_map : numpy.ndarray
        Fractional q-point(s) in the Brillouin zone of `map_struct`
        (same shape as `qpts`).
    """

    qpts = np.atleast_1d(np.asanyarray(qpts))

    if qpts.shape[-1] != 3:
        raise ValueError(
            "qpts must be an array_like with shape (3,) or (N, 3)."
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

    map_qpts -= np.round(map_qpts)

    return map_qpts
