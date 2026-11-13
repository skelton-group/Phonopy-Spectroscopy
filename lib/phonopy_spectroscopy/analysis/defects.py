# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for analysing phonon calculations on structures with
defects."""

# -------
# Imports
# -------

import numpy as np

from ..constants import ZERO_TOLERANCE
from ..structure import Structure

from ..utility.structure import (
    calculate_distances_from_fractional_positions,
    group_atoms,
    map_atom_positions,
    invert_atom_map,
)

# ---------------
# Helper routines
# ---------------


def _check_atom_indices(struct, at_inds):
    """Check a set of atom-indices against a structure."""

    if at_inds is not None:
        at_inds = np.atleast_1d(np.asanyarray(at_inds, dtype=int))

        if len(at_inds) == 0 or at_inds.ndim > 1:
            raise ValueError("at_inds must be an array_like with shape (N,).")

        if (at_inds < 0).any() or (at_inds >= struct.num_atoms).any():
            raise ValueError(
                "One or more of at_inds is out of range for the number "
                "of atoms in struct."
            )

        if len(np.unique(at_inds)) != len(at_inds):
            raise ValueError("at_inds contains duplicate indices.")

    return at_inds


def _check_nearest_neighbour_indices(struct, nn_inds):
    """Check a set of nearest-neighbour indices against a structure."""

    nn_inds = np.atleast_1d(np.asanyarray(nn_inds, dtype=int))

    if (nn_inds < 0).any() or (nn_inds > struct.num_atoms).any():
        raise ValueError(
            "One or more of nn_inds is out of range for the number of "
            "atoms in struct."
        )

    if len(np.unique(nn_inds)) != len(nn_inds):
        raise ValueError("nn_inds contains duplicate indices.")

    return nn_inds


def _get_centroid_or_position(struct, at_inds=None, pos=None, com=False):
    """Return a calculated centroid from a structure and set of atom
    indices or a specified position."""

    at_inds = _check_atom_indices(struct, at_inds)

    if at_inds is not None:
        return calculate_centroid_or_centre_of_mass(
            struct.atom_positions[at_inds],
            at_m=(struct.atomic_masses[at_inds] if com else None),
        )

    if pos is None:
        raise Exception("One of atom_idx or pos must be specified.")

    pos = np.asanyarray(pos, dtype=np.float64)

    if pos.shape != (3,):
        raise ValueError("pos must be an array_like with shape (3,).")

    return pos


# ----------------------
# Defects and neighbours
# ----------------------


def calculate_centroid_or_centre_of_mass(at_pos, at_m=None):
    """Calculate the geometric centroid or centre of mass of a set of
    atomic positions.

    Parameters
    ----------
    at_pos : array_like
        Atom positions (shape: `(N, 3)`).
    at_m : array_like or None, optional
        Atomic masses (shape: `(N,)`).

    Returns
    -------
    cent : numpy.ndarray
        Calculated centroid (shape: `(3,)`).
    """

    at_pos = np.atleast_2d(np.asanyarray(at_pos, dtype=np.float64))

    if at_pos.shape[-1] != 3:
        raise ValueError("at_pos must be an array_like with shape (N, 3).")

    if at_m is None:
        return np.mean(at_pos, axis=0)

    at_m = np.atleast_1d(np.asanyarray(at_m, dtype=np.float64))

    if at_m.ndim > 1 or len(at_m) != len(at_pos):
        raise ValueError("at_m must be an array_like with shape (M,).")

    return (at_m[:, np.newaxis] * at_pos).sum(axis=0) / at_m.sum()


def find_defects(struct, ref_struct, atom_map=None, group=True, **kwargs):
    """Identify the defects in a structure given a reference structure
    and an integer mapping of the atomic positions.

    Parameters
    ----------
    struct, ref_struct : Structure
        Defect and reference (undistorted) structures.
    atom_map : array_like or None, optional
        Integer mapping of the atoms in `struct` to those in
        `ref_struct` (`None` attempts an automatic mapping).
    group : bool, optional
        If `True`, substitutions/interstitials will be grouped using the
        `utility.structure.group_atoms` routine.
    **kwargs : any, optional
        Keyword arguments to `utility.structure.group_atoms`.

    Returns
    -------
    subs_ints_vpos : tuple of numpy.ndarray
        Indices or groups of indices (`group=True`) of substituted
        atoms/interstitials and fractional positions of vacancies in
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
      atom map to the same atom in `ref_struct` are identified as
      substitutions/interstitials.
    * Vacancy centres are identified based on atoms in `ref_struct` that
      are have no mapping from `map_struct`.

    This function requires the `atom_map` to be set up correctly. The
    default mapping will work for some "simple" cases but may fail for
    complex defects. For these cases, the optional parameters to the
    `map_atom_positions` routine in the `utility.structure` module
    provides finer control over the algorithm.

    Atoms can optionally be grouped based on interatomic distances by
    setting `group=True`. This can be used, for example, for molecular
    defects. This uses the `utility.structure.group_atoms`
    routine, and the optional parameters to this function can be set via
    keyword arguments.
    """

    if atom_map is not None:
        atom_map = np.asanyarray(atom_map, dtype=object)
    else:
        atom_map, _ = map_atom_positions(struct, ref_struct)

    # invert_atom_map() handles validation of atom_map.

    inv_atom_map = invert_atom_map(atom_map, struct, ref_struct)

    # Use three heuristics to identify substitutions/interstitials in
    # the defective structure:
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

        return (at_grp_inds, ref_struct.atom_positions[vac_inds])
    else:
        return (
            np.array(list(subs_ints_inds), dtype=int),
            ref_struct.atom_positions[vac_inds],
        )


def find_nearest_neighbours(struct, at_inds, pos=None, delta_r=1.0):
    """Identify the nearest neighbours to the atom(s) forming a centre
    or an arbitrary position in a structure using a simple distance
    search.

    Parameters
    ----------
    struct : Structure
        Structure to analyse.
    at_inds : int, array_like or None
        Index or indices of the atom(s) to analyse (shape: `(K,)`), or
        `None` for vacancies.
    pos : array_like, optional
        Position to analyse (shape: `(3,)`) (default: `None`).
    delta_r : float, optional
        Distance from first neighbour to identify additional neighbours
        (default: 1 Ang).

    Returns
    -------
    inds_dists : tuple of numpy.ndarray
        Atom indices and distances of nearest neighbours (shapes:
        `(K',)`).

    Notes
    -----
    This routine uses a simple distance search to identify nearest
    neighbous.

    The distances of all the atoms in `struct` from a reference
    position, specified implicitly with `atom_idx` or explicitly with
    `pos`, are calculated and the closest non-overlapping atom is
    identified.

    The neighbours are then identified as atoms at distances between
    `min_dist` -> `min_dist + delta_r`.
    """

    ref_pos = _get_centroid_or_position(
        struct, at_inds=at_inds, pos=pos, com=True
    )

    # ref_pos is a single position, so the distance table returned by
    # calculate_distances_frac() has shape (1, N).

    (dists,) = calculate_distances_from_fractional_positions(
        ref_pos, struct.lattice_vectors, other_pos=struct.atom_positions
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

    nn_inds = np.flatnonzero(dists <= dist_thr)

    if at_inds is not None:
        # Exclude reference atom indices from the neighbour list.

        mask = np.isin(nn_inds, at_inds, invert=True)
        nn_inds = nn_inds[mask]

    return (nn_inds, dists[nn_inds])


# ----------------
# Characterisation
# ----------------


def off_centring_character(edisps, struct, at_inds, nn_inds, group_corr=False):
    r"""Calculate the "off-centring" character for a set of phonon modes
    and a given defect centre and its neighbouring atoms.

    Parameters
    ----------
    edisps : array_like
        Cartesian displacements (eigendisplacements) (shape: `(N, 3)`
        or `(M, N, 3)`).
    struct : Structure
        Crystal structure.
    at_inds : int, array_like or None
        Index or indices of the atom(s) forming the defect centre
        (shape: `(K,)`), or `None` for a vacancy.
    nn_inds : int or array_like
        Index or indices of nearest neighbour atoms (shape: `(K',)`).
    group_corr : bool, optional
        Reweight the score for "group" defects where the centre contains
        multiple atoms (e.g. a molecule) (default: `False`).

    Returns
    -------
    ocs : numpy.ndarray
        Calculated off-centring scores (shape: `(M,)`).

    Notes
    -----
    The off-centring score is given by:

    .. math::

        \chi^\mathrm{OC}_j =
            \frac{ \left| \langle \boldsymbol{u}_{jk} \rangle_{k \in M} - \langle \boldsymbol{u}_{jk} \rangle_{k \in M^\prime} \right| }{ \sum_{k \in M, M^\prime} \left| \boldsymbol{u}_{jk} \right| }

    When `at_inds` is not set, the defect is treated as a vacancy and
    the terms involving the defect atoms are discarded.
    """

    edisps = np.asanyarray(edisps)

    ret_single = False

    if edisps.ndim == 2:
        edisps = edisps[np.newaxis, :, :]
        ret_single = True

    if edisps.ndim != 3 or edisps.shape[-1] != 3:
        raise ValueError(
            "edisps must be an array_like with shape (N, 3) or (M, N, 3)."
        )

    at_inds = _check_atom_indices(struct, at_inds)
    nn_inds = _check_nearest_neighbour_indices(struct, nn_inds)

    ocs = np.zeros((len(edisps),), dtype=np.float64)

    for i, edisp in enumerate(edisps):
        # Running total of the norms of the defect + neighbour
        # displacements.

        disp_norm_sum = 0.0

        # The defect contribution is only computed if at_inds is set.

        def_w_edisp = 0.0

        if at_inds is not None:
            def_w_edisp = calculate_centroid_or_centre_of_mass(
                edisp[at_inds], at_m=struct.atomic_masses[at_inds]
            )

            disp_norm_sum += np.linalg.norm(edisp[at_inds], axis=-1).sum()

        nn_w_edisp = calculate_centroid_or_centre_of_mass(
            edisp[nn_inds], at_m=struct.atomic_masses[nn_inds]
        )

        disp_norm_sum += np.linalg.norm(edisp[nn_inds], axis=-1).sum()

        if disp_norm_sum < ZERO_TOLERANCE:
            ocs[i] = 0.0
        else:
            ocs[i] = np.linalg.norm(def_w_edisp - nn_w_edisp) / disp_norm_sum

    if group_corr and at_inds is not None:
        # Apply a rewighting to correct the score for "group" defects
        # (e.g. molecules).

        ocs *= len(at_inds)

    return ocs[0] if ret_single else ocs


def radial_breathing_character(edisps, struct, at_inds, nn_inds, pos=None):
    r"""Calculate the "radial-breathing" character for a set of phonon
    modes and a given defect centre and its neighbouring atoms.

    Parameters
    ----------
    edisps : array_like
        Cartesian displacements (eigendisplacements) (shape: `(N, 3)`
        or `(M, N, 3)`).
    struct : Structure
        Crystal structure.
    at_inds : int, array_like or None
        Index or indices of the atom(s) forming the defect centre
        (shape: `(K,)`), or `None` for a vacancy.
    nn_inds : int or array_like
        Index or indices of nearest neighbour atoms (shape: `(K',)`).
    pos : array_like or None, optional
        Position of the defect centre (shape: `(3,)`) (default: `None`).

    Returns
    -------
    rbs : numpy.ndarray
        Calculated radial-breathing scores (shape: `(M,)`).

    Notes
    -----
    The radial-breathing character is given by:

    .. math::

        \chi^\mathrm{B}_j = \frac{ \left| \sum_{k \in M^\prime} \boldsymbol{u}_{jk} \cdot \hat{\boldsymbol{r}}_{k} \right| }{ \sum_{k \in M^\prime} \left| \boldsymbol{u}_{jk} \right| }

    For vacancies, the defect centre must be specified with the `pos`
    keyword.
    """

    edisps = np.asanyarray(edisps)

    ret_single = False

    if edisps.ndim == 2:
        edisps = edisps[np.newaxis, :, :]
        ret_single = True

    if edisps.ndim != 3 or edisps.shape[-1] != 3:
        raise ValueError(
            "edisps must be an array_like with shape (N, 3) or (M, N, 3)."
        )

    centre = _get_centroid_or_position(
        struct, at_inds=at_inds, pos=pos, com=True
    )

    nn_inds = _check_nearest_neighbour_indices(struct, nn_inds)

    # calculate_distances_frac calculates pos - other_pos. To project
    # the eigendisplacements, we need pos = neighbour positions and
    # other_pos = centre. The reshape() drops the extraneous second
    # dimension to make it easier to compute the dot product.

    nn_c_vecs = calculate_distances_from_fractional_positions(
        struct.atom_positions[nn_inds],
        struct.lattice_vectors,
        other_pos=centre,
        ret_vecs=True,
    ).reshape(-1, 3)

    nn_c_vecs /= np.linalg.norm(nn_c_vecs, axis=-1)[:, np.newaxis]

    rbs = np.zeros((len(edisps),), dtype=np.float64)

    for i, edisp in enumerate(edisps):
        nn_edisp = edisp[nn_inds]
        nn_edisp_norms = np.linalg.norm(nn_edisp, axis=-1)

        dot_prods = np.einsum("ij,ij->i", nn_c_vecs, nn_edisp)

        disp_norm_sum = nn_edisp_norms.sum()

        if disp_norm_sum < ZERO_TOLERANCE:
            rbs[i] = 0.0
        else:
            rbs[i] = np.abs(dot_prods.sum()) / disp_norm_sum

    return rbs[0] if ret_single else rbs


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

    centre = _get_centroid_or_position(
        struct, at_inds=at_inds, pos=pos, com=False
    )

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
