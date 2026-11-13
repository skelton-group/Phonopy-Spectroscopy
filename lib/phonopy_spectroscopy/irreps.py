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
"""Spectroscopic activities of irreps of crystallographic pointgroups.
"""

# ---------
# Functions
# ---------


def get_irrep_activities(pt_grp, ir_typ):
    """Return the spectroscopically-active active irreps for a point
    group.

    Parameters
    ----------
    pt_grp : str
        Point group.
    ir_typ : {"ir", "raman", "all"}
        Irrep type to return.

    Returns
    -------
    irreps: tuple of (list of str or None)
        Lists of active and inactive irrep symbols for `irrep_type`.
        Returns `(None, None)` if data for `point_group` is not
        available, or if `point_group` does not have irreps for
        `irrep_type`.
    """

    pt_grp = str(pt_grp).lower()
    ir_typ = str(ir_typ).lower()

    if pt_grp not in _IRREP_ACTIVITIES:
        raise ValueError(
            f'Unknown point_group="{pt_grp}" (this may be a bug).'
        )

    if ir_typ not in _IRREP_ACTIVITIES[pt_grp]:
        raise ValueError(
            f'No irreps of type "{ir_typ}" for point group "{pt_grp}" '
            "(this may be a bug)."
        )

    act_irs = [sym for sym in _IRREP_ACTIVITIES[pt_grp][ir_typ]]

    inact_irs = [
        sym for sym in _IRREP_ACTIVITIES[pt_grp]["all"] if sym not in act_irs
    ]

    return (np.array(act_irs, dtype=str), np.array(inact_irs, dtype=str))


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
            Symbols of irrep groups (shape: `(M,)`). `"None"` or `None`
            indicates unknown irreps.
        ir_band_inds : array_like
            Mapping of bands to irrep groups (shape: `(N,)`).
        """

        pt_grp = str(pt_grp).lower()

        ir_syms = np.atleast_1d(np.array(ir_syms, dtype=str, copy=True))

        for i, sym in enumerate(ir_syms):
            if sym is not None:
                ir_syms[i] = sym.capitalize()
            else:
                ir_syms[i] = "None"

        ir_band_inds = np.atleast_1d(
            np.array(ir_band_inds, dtype=int, copy=True)
        )

        if len(np.unique(ir_band_inds)) != len(ir_syms):
            raise ValueError(
                "The number of unique ir_band_inds is not consistent "
                "with ir_syms."
            )

        if ir_band_inds.min() < 0 or ir_band_inds.max() >= len(ir_syms):
            raise ValueError(
                "One or more ir_band_inds are out of bounds for ir_syms."
            )

        if pt_grp in _IRREP_ACTIVITIES:
            for sym in ir_syms:
                if (
                    sym != "None"
                    and sym not in _IRREP_ACTIVITIES[pt_grp]["all"]
                ):
                    raise ValueError(
                        f'"{sym}" is not a valid irrep of point group {pt_grp}.'
                    )
        else:
            raise ValueError(f'Point group "{pt_grp}" not recognised.')

        # Make arrays readonly.

        ir_syms.flags.writeable = False
        ir_band_inds.flags.writeable = False

        self._pt_grp = pt_grp
        self._ir_syms = ir_syms
        self._ir_band_inds = ir_band_inds

    @property
    def point_group(self):
        """str : Point group symbol."""
        return self._pt_grp

    @property
    def irrep_symbols(self):
        """numpy.ndarray : Irrep group symbols (shape: `(M,)`)."""
        return self._ir_syms

    @property
    def irrep_band_indices(self):
        """numpy.ndarray : Mapping of bands to irrep groups (shape:
        `(N,)`)."""
        return self._ir_band_inds

    def get_subset(self, band_inds):
        """Return a new `Irreps` object for a subset of the bands
        specified by `band_inds`.

        Parameters
        ----------
        band_inds : {"ir", "raman", "all"} or array_like
            Subset of bands to select.

        Returns
        -------
        subset : tuple
            Tuple of `(band_inds, irreps)` containing the band indices
            of the "parent" irreps and a new `Irreps` class for the
            subset of bands.
        """

        ir_syms = self._ir_syms
        ir_band_inds = self._ir_band_inds

        # Determine which band indices to include in the subset.

        sub_inds = None

        if isinstance(band_inds, (str, np.str_)):
            band_inds = band_inds.lower()

            match band_inds:
                case "all":
                    sub_inds = np.arange(0, len(ir_band_inds))
                case "ir" | "raman":
                    act_irs, _ = get_irrep_activities(self._pt_grp, band_inds)

                    # Get the indices of the entries in _ir_band_inds that
                    # point to active irreps.

                    act_inds = np.flatnonzero(np.isin(ir_syms, act_irs))
                    sub_inds = np.flatnonzero(np.isin(ir_band_inds, act_inds))
                case _:
                    raise ValueError(f'Unsupported band_inds="{band_inds}".')
        else:
            # Assume band_inds is an array_like.

            sub_inds = np.atleast_1d(np.asanyarray(band_inds, dtype=int))

            if sub_inds.ndim != 1:
                raise ValueError(
                    "If specified as an array_like, band_inds must "
                    "have shape (N,)."
                )

            if len(np.unique(sub_inds)) != len(sub_inds):
                raise ValueError("band_inds contains duplicate indices.")

            if sub_inds.min() < 0 or sub_inds.max() > len(ir_band_inds):
                raise ValueError("band_inds contains invalid indices.")

        sub_ir_band_inds = ir_band_inds[sub_inds]

        # Make sure sub_inds does not include an incomplete irrep group.

        ir_cnts = np.bincount(ir_band_inds, minlength=len(ir_syms))
        sub_ir_cnts = np.bincount(sub_ir_band_inds, minlength=len(ir_syms))

        if ((sub_ir_cnts > 0) & (sub_ir_cnts != ir_cnts)).any():
            raise ValueError("band_inds cannot cover partial irrep group(s).")

        ir_sym_inds, new_ir_band_inds = np.unique(
            sub_ir_band_inds, return_inverse=True
        )

        return (
            sub_inds,
            Irreps(self._pt_grp, ir_syms[ir_sym_inds], new_ir_band_inds),
        )

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
            "irrep_group_symbols": self._ir_syms.tolist(),
            "irrep_band_indices": self._ir_band_inds.tolist(),
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
