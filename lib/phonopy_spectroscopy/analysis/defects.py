# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines for analysing phonon calculations on structures with
defects."""


# -------
# Imports
# -------


import warnings

import numpy as np

from ..constants import ZERO_TOLERANCE
from ..structure import Structure

from ..utility.numpy_helper import np_check_shape, np_expand_dims

from ..utility.structure import (
    calculate_centroid,
    calculate_distances_frac,
    group_atoms,
    map_atom_positions,
    invert_atom_map,
)


# ---------------
# Helper routines
# ---------------


def _centroid_or_position(struct, at_inds=None, pos=None, com=False):
    """Return a calculated centroid from a structure and set of atom
    indices, or a specified position."""

    if at_inds is not None:
        at_inds, _ = np_expand_dims(np.asarray(at_inds, dtype=int), (None,))
        return calculate_centroid(struct, at_inds, com=com)

    if pos is None:
        raise Exception("One of atom_idx or pos must be specified.")

    pos = np.asarray(pos, dtype=np.float64)

    if not np_check_shape(pos, (3,)):
        raise ValueError("pos must be an array_like with shape (3,).")

    if (np.abs(pos) > 1.0).any():
        warnings.warn(
            "pos is outside the range [-1, 1] expected for "
            "fractional coordinates.",
            UserWarning,
        )

    return pos


# ----------------------
# Defects and neighbours
# ----------------------


def find_defects(struct, ref_struct, atom_map=None, group=True, **kwargs):
    """Identify the defects in a structure given a reference structure
    and an integer mapping of the atomic positions.

    Parameters
    ----------
    struct, ref_struct : Structure
        Defect and reference (undistorted) structures.
    atom_map : array_like or None, optional
        Integer mapping of the atoms in `struct` to those in
        `ref_struct`, or `None` when no mapping exists (the default
        value of `None` attempts an automatic mapping).
    group : bool, optional
        If `True`, vacancies/interstitials will be grouped using the
        `utility.structure.group_atoms` function.
    **kwargs : any, optional
        Keyword arguments to `utility.structure.group_atoms`.

    Returns
    -------
    subs_ints_vpos : tuple of numpy.ndarray
        Indices or groups of indices (`group=True`) of substituted
        atoms/interstitials, and fractional positions of vacancies in
        `map_struct`.

    See Also
    --------
    utility.structure.map_atom_positions
        Can be used to prepare the `atom_map` required by this function.
    utility.structure.group_atoms
        Used to group atoms when `group=True`.

    Notes
    -----
    This function applies a set of simple heuristics to identify various
    kinds of defects:

    * Atoms in `map_struct` that map to different types in `ref_struct`
      are identified as substitutions.
    * Atoms in `map_struct` that have no mapping to `ref_struct` are
      identified as interstitials.
    * Atoms in `map_struct` that are part of a "group" where multiple
      atom map to the same atom in `ref_struct`.
    * Vacancy centres are identified based on atoms in `ref_struct` that
      are have no mapping from `map_struct`.

    This function requires the `atom_map` to be set up correctly in
    order to identify defects. The default mapping will work for some
    "simple" cases but may fail in more complex ones. For these cases,
    finer control over the mapping can be obtained using the optional
    parameters to the `map_atom_positions` function in the
    `utility.structure` module.

    Atoms can optionally be grouped based on interatomic distances by
    setting `group=True`. This can be used, for example, for
    molecular defects. This uses the `utility.structure.group_atoms`
    function, and the optional parameters to this function can be set
    via keyword arguments.
    """

    if atom_map is not None:
        atom_map = np.asarray(atom_map, dtype=object)
    else:
        atom_map, _ = map_atom_positions(struct, ref_struct)

    # invert_atom_map() handles validation of atom_map.

    inv_atom_map = invert_atom_map(atom_map, struct, ref_struct)

    # Use two heuristics to identify substitutions/interstitials in the
    # defective structure:
    #   1. Atoms of types that are not present in the reference
    #   structure.
    #   2. Atoms that do not map to anything in the reference structure.
    #   3. Atoms that are part of a "group" that map to the same atom
    #   in the reference structure.

    subs_ints_inds = set()

    mask = np.isin(
        struct.atom_types, np.unique(ref_struct.atom_types), invert=True
    )

    (inds,) = np.where(mask)
    subs_ints_inds.update(inds)

    (inds,) = np.where(np.isin(atom_map, [None]))
    subs_ints_inds.update(inds)

    for inds in inv_atom_map:
        if np.ndim(inds) > 0:
            subs_ints_inds.update(inds)

    # Identify vacancies by atoms in the reference structure that do
    # not map to anything in the defective structure.

    (vac_inds,) = np.where(np.isin(inv_atom_map, [None]))

    if group:
        at_grp_inds = group_atoms(
            struct, atom_inds=list(subs_ints_inds), **kwargs
        )

        return (
            np.asarray(at_grp_inds, dtype=object),
            ref_struct.atom_positions[vac_inds],
        )
    else:
        return (
            np.array(list(subs_ints_inds), dtype=int),
            ref_struct.atom_positions[vac_inds],
        )


def find_nearest_neighbours(struct, at_inds=None, pos=None, delta_r=1.0):
    """Identify the nearest neighbours to an atom or arbitrary position
    in a structure using a simple distance search.

    Parameters
    ----------
    struct : Structure
        Structure to analyse.
    at_inds : int, array_like or None, optional
        Index or indices of the atom to analyse (default: `None`).
    pos : array_like, optional
        Position to analyse (shape: `(3,)`) (default: `None`).
    delta_r : float, optional
        Distance from first neighbour to identify additional neighbours
        (default: 1 Ang).

    Returns
    -------
    inds_dists : tuple of numpy.ndarray
        Atom indices and distances of nearest neighbours.

    Notes
    -----
    This function uses a simple distance search to identify nearest
    neighbous.

    The distances of all the atoms in `struct` from a reference
    position, specified implicitly with `atom_idx` or implicitly with
    `pos`, are calculated and the closest non-overlapping atom is
    identified.

    The neighbours are then identified as atoms within a distance range
    of `min_dist` -> `min_dist + delta_r`.
    """

    if at_inds is not None:
        at_inds, _ = np_expand_dims(np.asarray(at_inds, dtype=int), (None,))

    ref_pos = _centroid_or_position(struct, at_inds=at_inds, pos=pos, com=True)

    # ref_pos is a single position, so the distance table returned by
    # calculate_distances_frac() will have shape (1, N).

    (dists,) = calculate_distances_frac(
        ref_pos,
        struct.lattice_vectors,
        other_pos=struct.atom_positions,
        ret_vecs=False,
    )

    all_inds = np.arange(0, struct.num_atoms, dtype=int)

    dist_thr = None

    # If the reference position is defined by an atom position or a
    # centroid, exclude those atoms from the distance threshold.

    if at_inds is not None:
        mask = np.isin(all_inds, at_inds, invert=True)
        dist_thr = dists[mask].min() + delta_r
    else:
        dist_thr = dists.min() + delta_r

    (nn_inds,) = np.where(dists <= dist_thr)

    if at_inds is not None:
        # Exclude atom indices from the neighbour list.

        mask = np.isin(nn_inds, at_inds, invert=True)
        nn_inds = nn_inds[mask]

    return (nn_inds, dists[nn_inds])


# -------------
# Miscellaneous
# -------------


def centre_structure(struct, at_inds=None, pos=None):
    """Centre a structure on an atom or position.

    Parameters
    ----------
    struct : Structure
        Structure to centre.
    at_inds : int, array_like or None, optional
        Index or indices of atom(s) to centre on (default: `None`).
    pos : array_like, optional
        Fractional position to centre on (default: `None`).

    Returns
    -------
    centred_struct : Structure
        Centred structure.
    """

    centre = _centroid_or_position(struct, at_inds=at_inds, pos=pos, com=False)

    # The centre of the cell is (0.5, 0.5, 0.5)

    trans = (0.5 - centre) % 1.0

    return Structure(
        struct.lattice_vectors,
        (struct.atom_positions + trans) % 1.0,
        struct.atom_types,
        at_m=struct.atomic_masses,
        conv_trans=struct.conventional_transformation_matrix,
        cart_to_frac=False,
    )
