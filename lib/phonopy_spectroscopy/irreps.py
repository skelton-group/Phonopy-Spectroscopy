# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Data and routines for handling irreducible representations (irreps)."""


# -------
# Imports
# -------


import warnings

import numpy as np

from .utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
)


# ----
# Data
# ----


_IRREP_ACTIVITIES = {
    # Point group C_1.
    "1": {"ir": ["A"], "raman": ["A"], "all": ["A"]},
    # Point group C_i.
    "-1": {"ir": ["Au"], "raman": ["Ag"], "all": ["Ag", "Au"]},
    # Point group C_2.
    "2": {"ir": ["A", "B"], "raman": ["A", "B"], "all": ["A", "B"]},
    # Point group C_s.
    "m": {"ir": ["A'", "A''"], "raman": ["A'", "A''"], "all": ["A'", "A''"]},
    # Point group C_2h.
    "2/m": {
        "ir": ["Au", "Bu"],
        "raman": ["Ag", "Bg"],
        "all": ["Ag", "Au", "Bg", "Bu"],
    },
    # Point group D_2.
    "222": {
        "ir": ["B1", "B2", "B3"],
        "raman": ["A", "B1", "B2", "B3"],
        "all": ["A", "B1", "B2", "B3"],
    },
    # Point group C_2v.
    "mm2": {
        "ir": ["A1", "B1", "B2"],
        "raman": ["A1", "A2", "B1", "B2"],
        "all": ["A1", "A2", "B1", "B2"],
    },
    # Point group D_2h.
    "mmm": {
        "ir": ["B1u", "B2u", "B3u"],
        "raman": ["Ag", "B1g", "B2g", "B3g"],
        "all": ["Ag", "Au", "B1g", "B1u", "B2g", "B2u", "B3g", "B3u"],
    },
    # Point group C_4.
    "4": {"ir": ["A", "E"], "raman": ["A", "B", "E"], "all": ["A", "B", "E"]},
    # Point group S_4.
    "-4": {"ir": ["B", "E"], "raman": ["A", "B", "E"], "all": ["A", "B", "E"]},
    # Point group C_4h.
    "4/m": {
        "ir": ["Au", "Eu"],
        "raman": ["Ag", "Bg", "Eg"],
        "all": ["Ag", "Au", "Bg", "Bu", "Eg", "Eu"],
    },
    # Point group D_4.
    "422": {
        "ir": ["A2", "E"],
        "raman": ["A1", "B1", "B2", "E"],
        "all": ["A1", "A2", "B1", "B2", "E"],
    },
    # Point group C_4v.
    "4mm": {
        "ir": ["A1", "E"],
        "raman": ["A1", "B1", "B2", "E"],
        "all": ["A1", "A2", "B1", "B2", "E"],
    },
    # Point group D_2d.
    "-42m": {
        "ir": ["B2", "E"],
        "raman": ["A1", "B1", "B2", "E"],
        "all": ["A1", "A2", "B1", "B2", "E"],
    },
    # Point group D_4h.
    "4/mmm": {
        "ir": ["A2u", "Eu"],
        "raman": ["A1g", "B1g", "B2g", "Eg"],
        "all": [
            "A1g",
            "A1u",
            "A2g",
            "A2u",
            "B1g",
            "B1u",
            "B2g",
            "B2u",
            "Eg",
            "Eu",
        ],
    },
    # Point group C_3.
    "3": {"ir": ["A", "E"], "raman": ["A", "E"], "all": ["A", "E"]},
    # Point group C_3i.
    "-3": {
        "ir": ["Au", "Eu"],
        "raman": ["Ag", "Eg"],
        "all": ["Ag", "Au", "Eg", "Eu"],
    },
    # Point group D_3.
    "32": {"ir": ["A2", "E"], "raman": ["A1", "E"], "all": ["A1", "A2", "E"]},
    # Point group C_3v.
    "3m": {"ir": ["A1", "E"], "raman": ["A1", "E"], "all": ["A1", "A2", "E"]},
    # Point group D_3d.
    "-3m": {
        "ir": ["A2u", "Eu"],
        "raman": ["A1g", "Eg"],
        "all": ["A1g", "A1u", "A2g", "A2u", "Eg", "Eu"],
    },
    # Point group C_6.
    "6": {
        "ir": ["A", "E1"],
        "raman": ["A", "E1", "E2"],
        "all": ["A", "B", "E1", "E2"],
    },
    # Point group C_3h.
    "-6": {
        "ir": ["A''", "E'"],
        "raman": ["A'", "E'", "E''"],
        "all": ["A'", "A''", "E'", "E''"],
    },
    # Point group C_6h.
    "6/m": {
        "ir": ["Au", "E1u"],
        "raman": ["Ag", "E1g", "E2g"],
        "all": ["Ag", "Au", "Bg", "Bu", "E1g", "E1u", "E2g", "E2u"],
    },
    # Point group D_6.
    "622": {
        "ir": ["A2", "E1"],
        "raman": ["A1", "E1", "E2"],
        "all": ["A1", "A2", "B1", "B2", "E1", "E2"],
    },
    # Point group C_6v.
    "6mm": {
        "ir": ["A1", "E1"],
        "raman": ["A1", "E1", "E2"],
        "all": ["A1", "A2", "B1", "B2", "E1", "E2"],
    },
    # Point group D_3h.
    "-6m2": {
        "ir": ["A2''", "E'"],
        "raman": ["A1'", "E'", "E''"],
        "all": ["A1'", "A1''", "A2'", "A2''", "E'", "E''"],
    },
    # Point group D_6h.
    "6/mmm": {
        "ir": ["A2u", "E1u"],
        "raman": ["A1g", "E1g", "E2g"],
        "all": [
            "A1g",
            "A1u",
            "A2g",
            "A2u",
            "B1g",
            "B1u",
            "B2g",
            "B2u",
            "E1g",
            "E1u",
            "E2g",
            "E2u",
        ],
    },
    # Point group T.
    "23": {"ir": ["T"], "raman": ["A", "E", "T"], "all": ["A", "E", "T"]},
    # Point group T_h.
    "m-3": {
        "ir": ["Tu"],
        "raman": ["Ag", "Eg", "Tg"],
        "all": ["Ag", "Au", "Eg", "Eu", "Tg", "Tu"],
    },
    # Point group O.
    "432": {
        "ir": ["T1"],
        "raman": ["A1", "E", "T2"],
        "all": ["A1", "A2", "E", "T1", "T2"],
    },
    # Point group T_d.
    "-43m": {
        "ir": ["T2"],
        "raman": ["A1", "E", "T2"],
        "all": ["A1", "A2", "E", "T1", "T2"],
    },
    # Point group O_h.
    "m-3m": {
        "ir": ["T1u"],
        "raman": ["A1g", "Eg", "T2g"],
        "all": [
            "A1g",
            "A1u",
            "A2g",
            "A2u",
            "Eg",
            "Eu",
            "T1g",
            "T2g",
            "T1u",
            "T2u",
        ],
    },
}


# ---------
# Functions
# ---------


def get_irrep_activities(point_group, irrep_type):
    """Return the spectroscopically-active active irreps for a point
    group.

    Parameters
    ----------
    point_group : str
        Point group.
    irrep_type : {"ir", "raman", "all"}
        Type of irrep to return.

    Returns
    -------
    irreps: tuple of (list of str or None)
        Lists of active and inactive irrep symbols for `irrep_type`.
        Returns `(None, None)` if data for `point_group` is not
        available, or if `point_group` foes not have irreps for
        `irrep_type`.
    """

    point_group = str(point_group).lower()
    irrep_type = str(irrep_type).lower()

    # _IRREP_ACTIVITIES is should contain IR/Raman-active and complete
    # irreps for all the point groups that can be associated with
    # the crystallographic spacegroups. To help identify bugs, unknown
    # point groups or missing irreps raise ValueErrors.

    if point_group not in _IRREP_ACTIVITIES:
        raise ValueError(
            'Unknown point_group="{0}" (this may be a bug).'
            "".format(point_group)
        )

    if irrep_type not in _IRREP_ACTIVITIES[point_group]:
        raise ValueError(
            'No irreps of type irrep_type="{0}" for point_group='
            '"{1}" (this may be a bug).'.format(point_group, irrep_type)
        )

    active_irreps = [sym for sym in _IRREP_ACTIVITIES[point_group][irrep_type]]

    inactive_irreps = [
        sym
        for sym in _IRREP_ACTIVITIES[point_group]["all"]
        if sym not in active_irreps
    ]

    return (active_irreps, inactive_irreps)


# ------------
# Irreps class
# ------------


class Irreps:
    """Store and work with a set of mode irreducible representations
    (irreps)."""

    def __init__(self, pt_grp, ir_syms, ir_band_inds):
        """Create a new instance of the `Irreps` class.

        Parameters
        ----------
        pt_grp : str
            Point group symbol.
        ir_syms : array_like
            Symbols of irrep groups (shape: `(N,)`).
        ir_band_inds : array_like
            Indices of bands in irrep groups (shape: `(N, M)`, with `M`
            potentially non-uniform).
        """

        pt_grp = str(pt_grp).lower()

        ir_syms = np.array(
            [str(sym).capitalize() for sym in ir_syms], dtype=object
        )

        # In general, the sets of band indices in ir_band_inds may not
        # be of uniform length, so we cannot store them as a 2D array of
        # dtype=int. If the sets of indices do happen to be uniform,
        # the array needs to be set up "manually" as below to stop
        # NumPy from setting a single data type.

        temp = np.zeros((len(ir_band_inds),), dtype=object)

        for i, inds in enumerate(ir_band_inds):
            temp[i] = np_asarray_copy(inds, dtype=int)

        ir_band_inds = temp

        # Check that we have a set of band indices for each irrep
        # symbol, that each set of indices contains at least one entry,
        # and that there are no duplicate band indices.

        if len(ir_band_inds) != len(ir_syms):
            raise ValueError(
                "ir_band_inds must have the same length as ir_syms."
            )

        for inds in ir_band_inds:
            if len(inds) == 0:
                raise ValueError(
                    "All entries in ir_band_inds must contain at least "
                    "one band index."
                )

        band_inds = set()

        for inds in ir_band_inds:
            for idx in inds:
                if idx in band_inds:
                    raise ValueError(
                        "One or more entries in ir_band_inds has "
                        "duplicate band indices."
                    )

                band_inds.add(idx)

        # If irrep data is available for the point group, check the

        if pt_grp in _IRREP_ACTIVITIES:
            for sym in ir_syms:
                if (
                    sym != "None"
                    and sym not in _IRREP_ACTIVITIES[pt_grp]["all"]
                ):
                    raise ValueError(
                        "{0} is not a valid irrep of point group {1}."
                        "".format(sym, pt_grp)
                    )
        else:
            warnings.warn(
                'Point group "{0}" not recognised. Irrep symbols will '
                'be verified and get_subset() with irreps="ir" or '
                'irreps="raman" will fail.'.format(pt_grp)
            )

        self._pt_grp = pt_grp
        self._ir_syms = ir_syms
        self._ir_band_inds = ir_band_inds

    @property
    def point_group(self):
        """str : Point group symbol."""
        return self._pt_grp

    @property
    def irrep_symbols(self):
        """numpy.ndarray : Irrep group symbols."""
        return np_readonly_view(self._ir_syms)

    @property
    def irrep_band_indices(self):
        """list of numpy.ndarray : Band indices in irrep groups."""
        return [np_readonly_view(inds) for inds in self._ir_band_inds]

    def band_indices_flat(self):
        """Return a flat array of all the band indices covered by the
        irrep groups.

        Returns
        -------
        band_inds : ndarray
            Flat array of band indices.
        """

        band_inds = []

        for inds in self._ir_band_inds:
            band_inds.extend(inds)

        return np.array(band_inds, dtype=int)

    def get_subset(self, band_inds, reset_inds=False):
        """Return a new `Irreps` object for a subset of the bands
        specified by `band_inds`.

        Parameters
        ----------
        band_inds : {"ir", "raman", "all"} or array_like
            Subset of bands to select.
        reset_inds : bool, optional
            If `True`, "reset" the band indices to be continuous from
            zero (default: `False`).

        Returns
        -------
        subset : Irreps
            `Irreps` class containing the selected subset of bands.
        """

        # Determine which band indices to include in the subset.

        subset_inds = None

        band_inds_str = str(band_inds).lower()

        if band_inds_str == "all":
            subset_inds = self.band_indices_flat()

        elif band_inds_str == "ir" or band_inds_str == "raman":
            # This will fail if the point group is not recognised.

            active_irreps, _ = get_irrep_activities(
                self._pt_grp, band_inds_str
            )

            subset_inds = []

            for ir_sym, ir_band_inds in zip(self._ir_syms, self._ir_band_inds):
                if ir_sym in active_irreps:
                    subset_inds.extend(ir_band_inds)

        if subset_inds is None:
            subset_inds = np.asarray(band_inds)

            if not np_check_shape(subset_inds, (None,)):
                raise ValueError(
                    "band_inds must have shape (N,) if specified as an "
                    "array_like."
                )

            if len(subset_inds) != len(set(subset_inds)):
                raise ValueError(
                    "band_inds contains one or more duplicate indices."
                )

        # Determine which irrep groups to include.

        subset_ir_syms, subset_ir_band_inds = [], []

        for ir_sym, ir_band_inds in zip(self._ir_syms, self._ir_band_inds):
            test = np.isin(ir_band_inds, subset_inds)

            if test.all():
                subset_ir_syms.append(ir_sym)
                subset_ir_band_inds.append(ir_band_inds)
            else:
                if test.sum() > 0:
                    raise ValueError(
                        "band_inds cannot cover partial irrep group(s)."
                    )

        # If required, "reset" the band indices in the subset.

        if reset_inds:
            # Lookup table mapping "old" to "new" band indices.

            lut = {
                old: new
                for old, new in zip(
                    sorted(subset_inds), range(len(subset_inds))
                )
            }

            subset_ir_band_inds = [
                [lut[idx] for idx in band_inds]
                for band_inds in subset_ir_band_inds
            ]

        return Irreps(self._pt_grp, subset_ir_syms, subset_ir_band_inds)

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
            "point_group": self._pt_grp,
            "irrep_group_symbols": list(self._ir_syms),
            "irrep_band_indices": list(
                band_inds.tolist() for band_inds in self._ir_band_inds
            ),
        }

    @staticmethod
    def from_dict(d):
        """Create a new `Irreps` instance from a dictionary generated by
        `Irreps.to_dict()`.

        Parameters
        ----------
        d : dict
            Dictionary generated by `to_dict()`.

        Returns
        -------
        irreps : Irreps
            `Irreps` object constructed from the data in `d`.
        """

        return Irreps(
            d["point_group"], d["irrep_group_symbols"], d["irrep_band_indices"]
        )
