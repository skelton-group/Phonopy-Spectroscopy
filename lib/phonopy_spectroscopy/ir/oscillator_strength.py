# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines implementing single-crystal and powder infrared oscillator
strength simulations."""


# -------
# Imports
# -------


import warnings

import numpy as np

from ..intensity_base import (
    calculate_single_crystal_intensities,
    calculate_powder_intensities,
)


# --------------------------------
# Oscillator strength calculations
# --------------------------------


def _scalar_oscillator_strength(v_i, t, v_s):
    """Calculate the scalar oscillator strength for a mode oscillator
    strength tensor and an incident/scattered polarisation.

    Parameters
    ----------
    t : array_like
        Mode oscillator strength tensor (shape: `(3, 3)`).
    v_i, v_s : array_like
        Incident and scattered light polarisation vectors.

    Returns
    -------
    osc_str : float
        Mode oscillator strength.
    """

    return np.matmul(v_s, np.matmul(t, v_i))


# -----------------------------------
# Single-crystal oscillator strengths
# -----------------------------------


def calculate_single_crystal_oscillator_strengths(
    osc_str, geom, i_pol, s_pol, rot=None
):
    """Calculate the scalar mode oscillator strengths for an infrared
    measurement on a single crystal.

    Parameters
    ----------
    osc_str : array_like
        Mode oscillator strength tensors (shape: `(3, 3)` or
        `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of incident and measured light.
    rot : array_like
        Rotation matrix to apply to oscillator strength tensors.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity (scalar) or set of intensities
        (shape: `(N,)`).

    See Also
    --------
    intensity_base.calculate_single_crystal_intensities : Low-level
        implementation called by this function.
    """

    return calculate_single_crystal_intensities(
        _scalar_oscillator_strength, osc_str, geom, i_pol, s_pol, rot=rot
    )


# ---------------------------
# Powder oscillator strengths
# ---------------------------


def calculate_powder_oscillator_strengths(
    osc_str,
    geom,
    i_pol,
    s_pol,
    po_eta=0.0,
    po_surf_norm=None,
    method="best",
    lc_prec=5,
):
    """Calculate the scalar mode oscillator strengths for an infrared
    measurement on a powder.

    Parameters
    ----------
    osc_str : array_like
        Mode oscillator strength tensors (shape: `(3, 3)` or
        `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisation of incident and measured light.
    po_eta : float, optional
        Fractional excess of crystalites in preferred orientation
        (default: 0.0).
    po_surf_norm : array_like or str, optional
        Surface normal for preferred orientation (default: `None`,
        required if `po_eta` > 0).
    method : {"nquad", "lebedev+circle", "best"}
        Method for calculating intensnties.
    lc_prec : int
        Specifies the precision of the Lebedev/circle quadrature scheme
        for `method="lebedev+circle"`.

    Returns
    -------
    ints : float or numpy.ndarray
        Calcuated intensity (scalar) or set of intensities (shape:
        `(N,)`).

    See Also
    --------
    intensity_base.calculate_powder_intensities : Low-level
        implementation called by this function.
    """

    if method == "best":
        warnings.warn(
            'method="best" is currently not implemented and will '
            'default to method="nquad".',
            RuntimeWarning,
        )

        method = "nquad"

    return calculate_powder_intensities(
        _scalar_oscillator_strength,
        osc_str,
        geom,
        i_pol,
        s_pol,
        po_eta,
        po_surf_norm,
        -1.0 * geom.incident_direction,
        method=method,
        lc_prec=lc_prec,
    )
