# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Low-level routines implementing scalar intensity calculations from
tensor quantities, including powder averaging."""


# -------
# Imports
# -------


import numpy as np

from scipy.integrate import nquad

from .constants import ZERO_TOLERANCE

from .distributions import march_dollase, march_dollase_eta_to_r

from .utility.numpy_helper import np_expand_dims
from .utility.geometry import rotate_tensors, direction_cosine
from .utility.quadrature import lebedev_circle_quad


# -------------------------
# Internal validation/setup
# -------------------------


def _validate_polarisation_and_expand_tensors(t, geom, i_pol, s_pol):
    """Perform common validation and setup for intensity-calculation
    functions.

    Parameters
    ----------
    t : array_like
        Tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of incident and scattered light.

    Returns
    -------
    ret : tuple of (numpy.ndarray, bool)
        `t` with shape `(N, 3, 3)` and a flag set to `True` if `t` has
        multiple tensors, or `False` otherwise.
    """

    if not (
        geom.check_incident_polarisations(i_pol)
        and geom.check_scattered_polarisations(s_pol)
    ):
        raise ValueError(
            "The supplied incident and scattered polarisations are "
            "not possible with the supplied geometry."
        )

    t, n_dim_add = np_expand_dims(np.asarray(t), (None, 3, 3))

    return (t, n_dim_add == 0)


# --------------------------
# Single-crystal intensities
# --------------------------


def calculate_single_crystal_intensities(f, t, geom, i_pol, s_pol, rot=None):
    """Calculate the scalar intensities using a function `f` for
    measurement of a tensor quantity `t` with incident/scattered light
    polarisations `i_pol` and `s_pol`.

    Parameters
    ----------
    f : callable
        Function for calculating intensities.
    t : array_like
        Tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of incident and scattered light.
    rot : array_like, optional
        Optional rotation matrix to apply to tensors (shape: `(3, 3)`).

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).

    Notes
    -----
    `f` is called as `f(v_i, t, v_s)`, where `v_i` and `v_s` are the
    normalised polarisation vectors of the incident and scattered light
    (shape: `(3,)`) and `t` is a tensor (shape: `(3, 3)`).
    """

    # Check polarisations are valid for measurement geometry and
    # "expand" single tensors to a set of tensors.

    ts, ret_multi = _validate_polarisation_and_expand_tensors(
        t, geom, i_pol, s_pol
    )

    # Rotate tensors if required.

    if rot is not None:
        ts = rotate_tensors(ts, rot)

    ints = np.zeros((ts.shape[0],), dtype=np.float64)

    # Weighted sum over polarisation vectors.

    for i, t in enumerate(ts):
        temp = 0.0

        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            temp += w * f(v_i, t, v_s)

        ints[i] = temp

    return ints if ret_multi else ints[0]


# ---------------------------------
# Powder intensities: SciPy nquad()
# ---------------------------------


def calculate_powder_intensities_nquad(f, t, geom, i_pol, s_pol, int_f=None):
    """Calculate the scalar intensities using a function `f` for
    measurement of a tensor quantity `t`, averaged over orientations,
    with incident/scattered light polarisations `i_pol` and `s_pol`. The
    numerical integration is performed using the general SciPy `nquad()`
    routine.

    Parameters
    ----------
    f : callable
        Function for calculating intensities.
    t : array_like
        Tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of incident and scattered light.
    int_f : callable
        Directly specify the integrand function instead of dynamically
        constructing one around `f`.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).

    Notes
    -----
    `f` is called as `f(v_i, t, v_s)`, where `v_i` and `v_s` are the
    normalised polarisation vectors of the incident and scattered light
    (shape: `(3,)`) and `t` is a tensor (shape: `(3, 3)`).

    If `int_f` is specified, `f` is not used and `int_f` is called as
    `inf_f(phi, theta, psi, v_i, t, v_s)`, where `phi`, `theta` and
    `psi` are the Euler angles in radians.
    """

    ts, ret_multi = _validate_polarisation_and_expand_tensors(
        t, geom, i_pol, s_pol
    )

    if int_f is None:
        # "Wrap" f in an integrand function that rotates the tensor and
        # adds the sin(theta) term needed for integrating over spherical
        # polar coordinates and the normalisation factor of
        # 1 / (8 \pi^2).

        def _integrand(phi, theta, psi, v_i, t, v_s):
            rot = direction_cosine(phi, theta, psi)

            return f(v_i, np.matmul(rot, np.matmul(t, rot.T)), v_s) * (
                np.sin(theta) / (8.0 * np.pi**2)
            )

        int_f = _integrand

    ints = np.zeros((ts.shape[0],), dtype=np.float64)

    for i, t in enumerate(ts):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            int, _ = nquad(
                int_f,
                [(0.0, 2.0 * np.pi), (0.0, np.pi), (0.0, 2.0 * np.pi)],
                [v_i, t, v_s],
            )

            ints[i] += w * int

    return ints if ret_multi else ints[0]


def calculate_powder_intensities_with_odf_nquad(
    f, t, geom, i_pol, s_pol, po_eta, po_norm, po_axis, int_f=None
):
    """Calculate the scalar intensities using a function `f` for
    measurement of a tensor quantity `t`, averaged over orientations
    with the March-Dollase distribution function, with
    incident/scattered light polarisations `i_pol` and `s_pol`. The
    numerical integration is performed using the general SciPy `nquad()`
    routine.

    Parameters
    ----------
    f : callable
        Function for calculating intensities.
    t : array_like
        Tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of incident and scattered light.
    po_eta : float
        Crystallite fraction in preferred orientation.
    po_norm, po_axis : array_like
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`).
    int_f : callable
        Directly specify the integrand function instead of dynamically
        constructing one around `f`.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).

    Notes
    -----
    `f` is called as `f(v_i, t, v_s)`, where `v_i` and `v_s` are the
    normalised polarisation vectors of the incident and scattered light
    (shape: `(3,)`) and `t` is a tensor (shape: `(3, 3)`).

    If `int_f` is specified, `f` is not used and `int_f` is called as
    `inf_f(phi, theta, psi, v_i, t, v_s, po_r, po_norm, po_axis)`, where
    `phi`, `theta` and `psi` are the Euler angles in radians and `po_r`,
    `po_norm` and `po_axis` are the parameters for the March-Dollase
    distribution function.
    """

    ts, ret_multi = _validate_polarisation_and_expand_tensors(
        t, geom, i_pol, s_pol
    )

    po_r = march_dollase_eta_to_r(po_eta)

    if int_f is None:
        # Model preferred orientation by including the March-Dollase
        # distribution as a weighting factor in the integrand.

        def _integrand(phi, theta, psi, v_i, t, v_s, po_r, po_norm, po_axis):
            rot = direction_cosine(phi, theta, psi)
            po_alpha = np.arccos(np.dot(np.dot(rot, po_norm), po_axis))

            return (
                march_dollase(po_alpha, po_r)
                * f(v_i, np.matmul(rot, np.matmul(t, rot.T)), v_s)
                * (np.sin(theta) / (8.0 * np.pi**2))
            )

        int_f = _integrand

    ints = np.zeros((ts.shape[0],), dtype=np.float64)

    for i, t in enumerate(ts):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            int, _ = nquad(
                int_f,
                [(0.0, 2.0 * np.pi), (0.0, np.pi), (0.0, 2.0 * np.pi)],
                [v_i, t, v_s, po_r, po_norm, po_axis],
            )

            ints[i] += w * int

    return ints if ret_multi else ints[0]


# -----------------------------------------------
# Powder intensities: Levedev + circle quadrature
# -----------------------------------------------


def calculate_powder_intensities_lebedev_circle(
    f, t, geom, i_pol, s_pol, prec, int_f=None
):
    r"""Calculate the scalar intensities using a function `f` for
    measurement of a tensor quantity `t`, averaged over orientations,
    with incident/scattered light polarisations `i_pol` and `s_pol`. The
    numerical integration is performed using the Lebedev + circle
    quadrature scheme.

    Parameters
    ----------
    f : callable
        Function for calculating intensities.
    t : array_like
        Tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of incident and scattered light.
    prec : int
        Precision of the Lebedev + circle quadrature scheme.
    int_f : callable
        Directly specify the integrand function instead of dynamically
        constructing one around `f`.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).

    Notes
    -----
    `f` is called as `f(v_i, t, v_s)`, where `v_i` and `v_s` are the
    normalised polarisation vectors of the incident and scattered light
    (shape: `(3,)`) and `t` is a tensor (shape: `(3, 3)`).

    If `int_f` is specified, `f` is not used and `int_f` is called as
    `inf_f(phi, theta, psi, v_i, t, v_s)`, where `phi`, `theta` and
    `psi` are the Euler angles in radians.

    For Lebedev + circle quadrature, `int_f` should not include the
    factors of 1 / (8 \pi^2) and sin(\theta) in the integrand, as these
    are accounted for in the quadrature weights.
    """

    ts, ret_multi = _validate_polarisation_and_expand_tensors(
        t, geom, i_pol, s_pol
    )

    if int_f is None:
        # For Lebedev + circle quadrature, the factor of sin(theta) is
        # included in the quadrature weights.

        def _integrand(phi, theta, psi, v_i, t, v_s):
            rot = direction_cosine(phi, theta, psi)
            return f(v_i, np.matmul(rot, np.matmul(t, rot.T)), v_s)

        int_f = _integrand

    ints = np.zeros((ts.shape[0],), dtype=np.float64)

    for i, t in enumerate(ts):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            ints[i] += w * lebedev_circle_quad(
                int_f, prec, n=None, args=[v_i, t, v_s]
            )

    return ints if ret_multi else ints[0]


def calculate_powder_intensities_with_odf_lebedev_circle(
    f, t, geom, i_pol, s_pol, po_eta, po_norm, po_axis, prec, int_f=None
):
    r"""Calculate the scalar intensities using a function `f` for
    measurement of a tensor quantity `t`, averaged over orientations
    with the March-Dollase distribution function, with
    incident/scattered light polarisations `i_pol` and `s_pol`. The
    numerical integration is performed using the Lebedev + circle
    quadrature scheme.

    Parameters
    ----------
    f : callable
        Function for calculating intensities.
    t : array_like
        Tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of incident and scattered light.
    po_eta : float
        Crystallite fraction in preferred orientation.
    po_norm, po_axis : array_like
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`).
    prec : int
        Precision of the Lebedev + circle quadrature scheme.
    int_f : callable
        Directly specify the integrand function instead of dynamically
        constructing one around `f`.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).

    Notes
    -----
    `f` is called as `f(v_i, t, v_s)`, where `v_i` and `v_s` are the
    normalised polarisation vectors of the incident and scattered light
    (shape: `(3,)`) and `t` is a tensor (shape: `(3, 3)`).

    If `int_f` is specified, `f` is not used and `int_f` is called as
    `inf_f(phi, theta, psi, v_i, t, v_s, po_r, po_norm, po_axis)`, where
    `phi`, `theta` and `psi` are the Euler angles in radians and `po_r`,
    `po_norm` and `po_axis` are the parameters for the March-Dollase
    distribution function.

    For Lebedev + circle quadrature, `int_f` should not include the
    factors of 1 / (8 \pi^2) and sin(\theta) in the integrand, as these
    are accounted for in the quadrature weights.
    """

    ts, ret_multi = _validate_polarisation_and_expand_tensors(
        t, geom, i_pol, s_pol
    )

    po_r = march_dollase_eta_to_r(po_eta)

    if int_f is None:

        def _integrand(phi, theta, psi, v_i, t, v_s, po_r, po_norm, po_axis):
            rot = direction_cosine(phi, theta, psi)

            w = march_dollase(
                np.arccos(np.dot(np.dot(rot, po_norm), po_axis)), po_r
            )

            return w * f(v_i, np.matmul(rot, np.matmul(t, rot.T)), v_s)

        int_f = _integrand

    ints = np.zeros((ts.shape[0],), dtype=np.float64)

    for i, t in enumerate(ts):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            ints[i] += w * lebedev_circle_quad(
                int_f, prec, n=None, args=[v_i, t, v_s, po_r, po_norm, po_axis]
            )

    return ints if ret_multi else ints[0]


# -----------------------------------------
# Powder intensities: "dispatcher" function
# -----------------------------------------


def calculate_powder_intensities(
    f,
    t,
    geom,
    i_pol,
    s_pol,
    po_eta=0.0,
    po_norm=None,
    po_axis=None,
    method="best",
    lc_prec=5,
    nquad_int_f=None,
    nquad_odf_int_f=None,
    lc_int_f=None,
    lc_odf_int_f=None,
):
    """Calculate the scalar intensities using a function `f` for
    measurement of a tensor quantity `t`, averaged over orientations
    with the March-Dollase distribution function, with
    incident/scattered light polarisations `i_pol` and `s_pol`. The
    method used for numerical integration is selected using the `method`
    parameter.

    Parameters
    ----------
    f : callable
        Function for calculating intensities.
    t : array_like
        Tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of incident and scattered light.
    po_eta : float
        Crystallite fraction in preferred orientation.
    po_norm, po_axis : array_like
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`).
    method : {"nquad", "lebedev+circle", "best"}
        Method for calculating intensnties.
    lc_prec : int
        Specifies the precision of the Lebedev/circle quadrature scheme
        for `method="lebedev+circle"`.
    nquad_int_f, nquad_odf_int_f, lc_int_f, lc_odf_int_f : callable
        Directly specify the integrand functions to the underlying
        routines for calculating scalar intensities.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).

    Notes
    -----
    `f` is called as `f(v_i, t, v_s)`, where `v_i` and `v_s` are the
    normalised polarisation vectors of the incident and scattered light
    (shape: `(3,)`) and `t` is a tensor (shape: `(3, 3)`).

    For the call signatures of the `*_int_f` parameters, see the
    functions referenced below.

    See Also
    --------
    calculate_powder_intensities_nquad,
    calculate_powder_intensities_with_odf_nquad,
    calculate_powder_intensities_lebedev_circle,
    calculate_powder_intensities_with_odf_lebedev_circle :
        Functions used to perform powder Raman calculations.
    """

    # Determine whether a calculation with a preferred orientation is
    # required.

    pref_orient = np.abs(po_eta) > ZERO_TOLERANCE

    if method == "nquad":
        if pref_orient:
            return calculate_powder_intensities_with_odf_nquad(
                f,
                t,
                geom,
                i_pol,
                s_pol,
                po_eta,
                po_norm,
                po_axis,
                int_f=nquad_odf_int_f,
            )
        else:
            return calculate_powder_intensities_nquad(
                f, t, geom, i_pol, s_pol, int_f=nquad_int_f
            )

    if method == "lebedev+circle":
        if pref_orient:
            return calculate_powder_intensities_with_odf_lebedev_circle(
                f,
                t,
                geom,
                i_pol,
                s_pol,
                po_eta,
                po_norm,
                po_axis,
                lc_prec,
                int_f=lc_odf_int_f,
            )
        else:
            return calculate_powder_intensities_lebedev_circle(
                f, t, geom, i_pol, s_pol, lc_prec, int_f=lc_int_f
            )

    raise ValueError('Unknown method="{0}".'.format(method))
