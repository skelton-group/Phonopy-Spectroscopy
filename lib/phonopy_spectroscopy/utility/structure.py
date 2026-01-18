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

from ..structure import (
    cartesian_to_fractional_coordinates,
    fractional_to_cartesian_coordinates,
)


# ---------
# Functions
# ---------


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
