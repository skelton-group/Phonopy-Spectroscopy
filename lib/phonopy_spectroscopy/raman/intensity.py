# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines implementing single-crystal and powder Raman simulations."""


# -------
# Imports
# -------


import warnings

import numpy as np

from scipy import LowLevelCallable
from scipy.integrate import tplquad

from ..constants import ZERO_TOLERANCE
from ..distributions import march_dollase, march_dollase_eta_to_r
from ..utility.geometry import direction_cosine, rotate_tensors
from ..utility.numpy_helper import np_check_shape, np_expand_dims
from ..utility.quadrature import lebedev_circle_quad

_NUMBA_AVAILABLE = False

try:
    from numba import njit, cfunc, types, carray

    _NUMBA_AVAILABLE = True
except ImportError:
    from ..utility.numba_helper import dummy_njit as njit


# ---------
# Constants
# ---------

_EIGHT_PI_SQUARED = 8.0 * np.pi**2

"""Value of 8 * pi^2. """


# ----------------
# Helper functions
# ----------------


def _validate_polarisation_and_expand_tensors(r_t, geom, i_pol, s_pol):
    """Perform common validation and setup for Raman intensity
    calculations.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.

    Returns
    -------
    ret : tuple of (numpy.ndarray, bool)
        `r_t` with shape `(N, 3, 3)` and a flag set to `True` if `r_t`
        has multiple tensors, or `False` otherwise.
    """

    if not (
        geom.check_incident_polarisations(i_pol)
        and geom.check_scattered_polarisations(s_pol)
    ):
        raise ValueError(
            "The supplied incident and scattered polarisations are "
            "not possible with the supplied geometry."
        )

    r_t, n_dim_add = np_expand_dims(np.asarray(r_t), (None, 3, 3))

    return (r_t, n_dim_add == 0)


def _validate_and_convert_march_dollase_params(po_eta, po_norm, po_axis):
    """Perform common validation and parameter conversion for powder
    Raman intensity calculations using the March-Dollase orientation
    distribution function.

    Parameters
    ----------
    po_eta : float
        Crystallite fraction in the preferred orientation.
    po_norm, po_axis : array_like
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`).

    Returns
    -------
    md_params : tuple of (float, numpy.ndarray, numpy.ndarray)
        Tuple of `(po_r, po_norm, po_axis)`.
    """

    # match_dollase_eta_to_r checks crystallite fraction.

    po_r = march_dollase_eta_to_r(po_eta)

    po_norm = np.asarray(po_norm, dtype=np.float64)
    po_axis = np.asarray(po_axis, dtype=np.float64)

    if not (np_check_shape(po_norm, (3,)) and np_check_shape(po_axis, (3,))):
        raise ValueError(
            "po_norm and po_axis must both be an array_like with shape "
            "(3,)."
        )

    if (
        np.abs(np.linalg.norm(po_norm) - 1.0) > ZERO_TOLERANCE
        or np.abs(np.linalg.norm(po_axis) - 1.0) > ZERO_TOLERANCE
    ):
        raise ValueError(
            "po_norm and po_axis must be non-zero and normalised."
        )

    return (po_r, po_norm, po_axis)


def _match_dtypes(t, v_i, v_s, po_norm=None, po_axis=None):
    """Match the data types of a Raman tensor and a pair of incident/
    scattered polarisation vectors.

    Parameters
    ----------
    t : numpy.ndarray
        Raman tensor (shape: `(3, 3)`).
    v_i, v_s : numpy.ndarray
        Polarisation vectors (shape: `(3,)`).
    po_norm, po_axis : numpy.ndarray or None, optional
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`).

    Returns
    -------
    res : tuple of numpy.ndarray
        `(t, v_i, v_s)` or `(t, v_i, v_s, po_norm, po_axis)` with the
        `numpy.float64` data type if all arguments are real, or the
        `numpy.complex128` type otherwise.

    Notes
    -----
    This function is required to work around a limitation of the Numba.
    """

    if po_norm is not None and po_axis is None:
        raise ValueError(
            "Only one of po_norm/po_axis are set (this is most likely "
            "a bug)."
        )

    dtype = np.float64

    if (
        np.iscomplex(t).any()
        or np.iscomplex(v_i).any()
        or np.iscomplex(v_s).any()
    ):
        dtype = np.complex128

    res = (t.astype(dtype), v_i.astype(dtype), v_s.astype(dtype))

    if po_norm is not None:
        # If supplied, po_norm and po_axis should always be real.

        res = res + (po_norm.astype(dtype), po_axis.astype(dtype))

    return res


# --------------------
# Single-crystal Raman
# --------------------


def calculate_single_crystal_raman_intensities(
    r_t, geom, i_pol, s_pol, rot=None
):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement on a single crystal.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    rot : array_like or None, optional
        Rotation matrix to apply to Raman tensors.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity (scalar) or set of intensities
        (shape: `(N,)`).
    """

    # Check polarisations are valid for measurement geometry and
    # "expand" single tensors to a set of tensors.

    r_t, ret_multi = _validate_polarisation_and_expand_tensors(
        r_t, geom, i_pol, s_pol
    )

    # Rotate tensors if required.

    if rot is not None:
        r_t = rotate_tensors(r_t, rot)

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)

    # Weighted sum over polarisation vectors.

    for i, t in enumerate(r_t):
        temp = 0.0

        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            temp += w * np.abs(np.vdot(v_s, np.dot(t, v_i))).real ** 2

        ints[i] = temp

    return ints if ret_multi else ints[0]


# ------------------------
# Powder Raman: analytical
# ------------------------


def calculate_powder_raman_intensities_analytical(r_t, geom, i_pol, s_pol):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, using the analytical formula.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.

    Returns
    -------
    ints : float or numpy.ndarray
        Calcuated intensity or array of intensities (shape: `(N,)`).

    Notes
    -----
    The formula calculates the intensity from the angle between the
    incident and scattered polarisation, and is only valid if:

    * The Raman tensors are real;
    * The laser polarisation is perpendicular to the collection axis; and
    * The polarisation vectors are real.
    """

    if not (
        geom.check_incident_polarisations(i_pol)
        and geom.check_scattered_polarisations(s_pol)
    ):
        raise ValueError(
            "The supplied incident and scattered polarisation are "
            "not possible with the supplied geometry."
        )

    r_t, n_dim_add = np_expand_dims(np.asarray(r_t), (None, 3, 3))

    # The analytical formula is only valid if the incident polarisation
    # is perpendicular to the collection axis.

    if not i_pol.check_perpendicular(geom.collection_direction):
        raise RuntimeError(
            "The analytical formula can only be used when the "
            "incident polarisation is perpendicular to the collection "
            "axis - use the numerical routines instead."
        )

    # Raman tensors calculated without using the far from resonance
    # (FFR) approximation can in general be complex. Since the
    # analytical formula was derived under the FFR, we only use it
    # for real Raman tensors.

    if np.iscomplex(r_t).any():
        raise RuntimeError(
            "The analytical formula can only be used for real Raman "
            "tensors - use the numerical routines instead."
        )

    # Since the formula was not explicitly derived for complex
    # polarisation vectors, we do not allow this either.

    if i_pol.is_complex or s_pol.is_complex:
        raise RuntimeError(
            "The analytical formula can only be used for real "
            "polarisation vectors - use the numerical routines instead."
        )

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)

    for i, t in enumerate(r_t.real):
        # \alpha^\prime in Porezag and Pederson's notation.

        a_p = np.abs(np.trace(t) / 3.0)

        # (\beta^\prime)^2.

        b_p_2 = (
            (t[0, 0] - t[1, 1]) ** 2
            + (t[0, 0] - t[2, 2]) ** 2
            + (t[1, 1] - t[2, 2]) ** 2
            + 6.0 * (t[0, 1] ** 2 + t[0, 2] ** 2 + t[1, 2] ** 2)
        ) / 2.0

        i_par = a_p**2 + (4.0 / 45.0) * b_p_2
        i_per = (3.0 / 45.0) * b_p_2

        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            cos_chi = np.dot(v_i, v_s) / (
                np.linalg.norm(v_i) * np.linalg.norm(v_s)
            )

            ints[i] += w * (i_per + (i_par - i_per) * cos_chi**2)

    return ints if n_dim_add == 0 else ints[0]


# --------------------------
# Powder Raman: SciPy quad()
# --------------------------


def _setup_int_func_powder_quad(t, v_i, v_s):
    """Generate an integrand function for calculating a powder-average
    Raman intensity using the SciPy `quad()` routine.

    Parameters
    ----------
    t : numpy.ndarray
        Raman tensor (shape: `(3, 3)`).
    v_i, v_s : numpy.ndarray
        Incident and scattered light polarisation vectors.

    Returns
    -------
    int_func : callable or scipy.LowLevelCallable
        Integrand function.

    Notes
    -----
    If Numba is available, the integrand function is JIT-compiled and
    wrapped in a `scipy.LowLevelCallable`.
    """

    if _NUMBA_AVAILABLE:
        # Numba requires dot arguments to have the same dtype.

        t, v_i, v_s = _match_dtypes(t, v_i, v_s)

    @njit(inline="always")
    def _int_func(phi, theta, psi):
        # Ensure r has the same dtype as t.

        r = direction_cosine(phi, theta, psi).astype(t.dtype)

        p = np.vdot(v_s, np.dot(r, np.dot(t, np.dot(r.T, v_i))))

        return (np.abs(p).real ** 2) * (np.sin(theta) / _EIGHT_PI_SQUARED)

    if _NUMBA_AVAILABLE:
        # Wrap integrand function into a SciPy LowLevelCallable.

        @cfunc(types.float64(types.intc, types.CPointer(types.float64)))
        def _int_func_wrapper(n, args):
            a = carray(args, n)
            return _int_func(a[0], a[1], a[2])

        _int_func = LowLevelCallable(_int_func_wrapper.ctypes)

    return _int_func


def calculate_powder_intensities_quad(r_t, geom, i_pol, s_pol):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, using the SciPy `quad()`
    routine.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).
    """

    r_t, ret_multi = _validate_polarisation_and_expand_tensors(
        r_t, geom, i_pol, s_pol
    )

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)

    for i, t in enumerate(r_t):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            f = _setup_int_func_powder_quad(t, v_i, v_s)

            int, _ = tplquad(f, 0.0, 2.0 * np.pi, 0.0, np.pi, 0.0, 2.0 * np.pi)

            ints[i] += w * int

    return ints if ret_multi else ints[0]


def _setup_int_func_powder_md_quad(t, v_i, v_s, po_r, po_norm, po_axis):
    """Generate an integrand function for calculating a powder-average
    Raman intensity, with the March-Dollase orientation distribution
    function, using the SciPy `quad()` routine.

    Parameters
    ----------
    t : array_like
        Raman tensor (shape: `(3, 3)`).
    v_i, v_s : array_like
        Incident and scattered light polarisation vectors.
    po_r : float
        r parameter in the March-Dollase distribution.
    po_norm, po_axis : array_like
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`).

    Returns
    -------
    int_func : callable or scipy.LowLevelCallable
        Integrand function.

    Notes
    -----
    If Numba is available, the integrand function is JIT-compiled and
    wrapped in a `scipy.LowLevelCallable`.
    """

    if _NUMBA_AVAILABLE:
        t, v_i, v_s, po_norm, po_axis = _match_dtypes(
            t, v_i, v_s, po_norm, po_axis
        )

    @njit(inline="always")
    def _int_func(phi, theta, psi):
        r = direction_cosine(phi, theta, psi).astype(t.dtype)
        p = np.vdot(v_s, np.dot(r, np.dot(t, np.dot(r.T, v_i))))
        a = np.arccos(np.dot(np.dot(r, po_norm), po_axis))

        return (
            march_dollase(a, po_r)
            * np.abs(p).real ** 2
            * (np.sin(theta) / _EIGHT_PI_SQUARED)
        )

    if _NUMBA_AVAILABLE:

        @cfunc(types.float64(types.intc, types.CPointer(types.float64)))
        def _int_func_wrapper(n, args):
            a = carray(args, n)
            return _int_func(a[0], a[1], a[2])

        _int_func = LowLevelCallable(_int_func_wrapper.ctypes)

    return _int_func


def calculate_powder_intensities_march_dollase_quad(
    r_t, geom, i_pol, s_pol, po_eta, po_norm, po_axis
):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, and with preferred orientation
    modelled with the March-Dollase orientation distribution function,
    using the SciPy `quad()` routine.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    po_eta : float
        Crystallite fraction in the preferred orientation.
    po_norm, po_axis : array_like
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`).

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).
    """

    r_t, ret_multi = _validate_polarisation_and_expand_tensors(
        r_t, geom, i_pol, s_pol
    )

    po_r, po_norm, po_axis = _validate_and_convert_march_dollase_params(
        po_eta, po_norm, po_axis
    )

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)

    for i, t in enumerate(r_t):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            f = _setup_int_func_powder_md_quad(
                t, v_i, v_s, po_r, po_norm, po_axis
            )

            int, _ = tplquad(
                f,
                0.0,
                2.0 * np.pi,
                0.0,
                np.pi,
                0.0,
                2.0 * np.pi,
            )

        ints[i] += w * int

    return ints if ret_multi else ints[0]


# ------------------------------
# Powder Raman: Lebedev + circle
# ------------------------------


@njit
def _powder_lebedev_int(phi, theta, psi, t, v_i, v_s):
    """Integrand function for calculating powder-averaged Raman
    intensities using Lebedev + circle quadrature.

    Parameters
    ----------
    phi, theta, psi : float
        Euler angles.
    t : array_like
        Raman tensor (shape: `(3, 3)`).
    v_i, v_s : array_like
        Incident and scattered light polarisation vectors.

    Returns
    -------
    int : float
        Scalar Raman intensity.
    """

    r = direction_cosine(phi, theta, psi).astype(t.dtype)

    return (
        np.abs(np.vdot(v_s, np.dot(r, np.dot(t, np.dot(r.T, v_i))))).real ** 2
    )


def calculate_powder_intensities_leb_circ(r_t, geom, i_pol, s_pol, prec):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, using Lebedev + circle
    quadrature.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    prec : int
        Precision of the Lebedev + circle quadrature scheme.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).
    """

    r_t, ret_multi = _validate_polarisation_and_expand_tensors(
        r_t, geom, i_pol, s_pol
    )

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)

    for i, t in enumerate(r_t):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            if _NUMBA_AVAILABLE:
                t, v_i, v_s = _match_dtypes(t, v_i, v_s)

            ints[i] += w * lebedev_circle_quad(
                _powder_lebedev_int, prec, n=None, args=[t, v_i, v_s]
            )

    return ints if ret_multi else ints[0]


@njit
def _powder_lebedev_odf_int(
    phi, theta, psi, v_i, t, v_s, po_r, po_norm, po_axis
):
    """Integrand function for calculating powder-averaged Raman
    intensities, with the March-Dollase orientation distribution
    function, using Lebedev + circle quadrature.

    Parameters
    ----------
    phi, theta, psi : float
        Euler angles.
    t : array_like
        Raman tensor (shape: `(3, 3)`).
    v_i, v_s : array_like
        Incident and scattered light polarisation vectors.
    po_r : float
        r parameter in the March-Dollase distribution.
    po_norm, po_axis : array_like
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`).

    Returns
    -------
    int : float
        Scalar Raman intensity.
    """

    r = direction_cosine(phi, theta, psi).astype(t.dtype)
    po_alpha = np.arccos(np.dot(np.dot(r, po_norm), po_axis))

    return (
        march_dollase(po_alpha, po_r)
        * np.abs(np.vdot(v_s, np.dot(r, np.dot(t, np.dot(r.T, v_i))))).real
        ** 2
    )


def calculate_powder_intensities_with_march_dollase_leb_circ(
    r_t, geom, i_pol, s_pol, po_eta, po_norm, po_axis, prec
):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, and with preferred orientation
    modelled with the March-Dollase orientation distribution function,
    using Lebedev + circle quadrature.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    po_eta : float
        Crystallite fraction in the preferred orientation.
    po_norm, po_axis : array_like
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`).
    prec : int
        Precision of the Lebedev + circle quadrature scheme.

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity or intensities (scalar or shape: `(N,)`).
    """

    r_t, ret_multi = _validate_polarisation_and_expand_tensors(
        r_t, geom, i_pol, s_pol
    )

    po_r, po_norm, po_axis = _validate_and_convert_march_dollase_params(
        po_eta, po_norm, po_axis
    )

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)

    for i, t in enumerate(r_t):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            if _NUMBA_AVAILABLE:
                t, v_i, v_s, po_norm, po_axis = _match_dtypes(
                    t, v_i, v_s, po_norm, po_axis
                )

            ints[i] += w * lebedev_circle_quad(
                _powder_lebedev_odf_int,
                prec,
                n=None,
                args=[v_i, t, v_s, po_r, po_norm, po_axis],
            )

    return ints if ret_multi else ints[0]


def calculate_powder_raman_intensities(
    r_t,
    geom,
    i_pol,
    s_pol,
    po_eta=0.0,
    po_surf_norm=None,
    method="best",
    lc_prec=5,
):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement on a powder.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    po_eta : float, optional
        Fractional excess of crystalites in the preferred orientation
        (default: 0.0).
    po_surf_norm : array_like or str, optional
        Surface normal for the preferred orientation (default: `None`,
        required if `po_eta` > 0).
    method : {"quad", "leb+circ", "best"}
        Method for calculating intensnties.
    lc_prec : int
        Specifies the precision of the Lebedev/circle quadrature scheme
        for `method="leb+circ"`.

    Returns
    -------
    ints : float or numpy.ndarray
        Calcuated intensity (scalar) or set of intensities (shape:
        `(N,)`).
    """

    r_t, _ = np_expand_dims(np.asarray(r_t), (None, 3, 3))

    # Determine whether an energy-dependent Raman calculation and/or
    # a calculation with a preferred orientation are required.

    complex_rt = np.iscomplex(r_t).any()
    pref_orient = po_eta > ZERO_TOLERANCE

    # If the Raman tensors are real, there is no preferred orientation,
    # and the incident polarisation is perpendicular to the collection
    # direction, we can use the analytical formula.

    if (
        method == "best"
        and not complex_rt
        and not pref_orient
        and i_pol.check_perpendicular(geom.collection_direction)
        and not (i_pol.is_complex or s_pol.is_complex)
    ):
        return calculate_powder_raman_intensities_analytical(
            r_t, geom, i_pol, s_pol
        )

    # If a preferred orientation is specified, or if the Lebedev +
    # circle precision is set to the minimum value, the SciPy quad()
    # routine is the "safe" option.

    if method == "best":
        if not pref_orient and lc_prec >= 5:
            method = "leb+circ"
        else:
            method = "quad"

    if method == "quad":
        if not _NUMBA_AVAILABLE:
            warnings.warn(
                'Numerical integration with method="quad" may be '
                "significantly faster if Numba is installed.",
                RuntimeWarning,
            )

        if pref_orient:
            return calculate_powder_intensities_march_dollase_quad(
                r_t,
                geom,
                i_pol,
                s_pol,
                po_eta,
                po_surf_norm,
                -1.0 * geom.incident_direction,
            )
        else:
            return calculate_powder_intensities_quad(r_t, geom, i_pol, s_pol)

    if method == "leb+circ":
        if pref_orient:
            return calculate_powder_intensities_with_march_dollase_leb_circ(
                r_t,
                geom,
                i_pol,
                s_pol,
                po_eta,
                po_surf_norm,
                -1.0 * geom.incident_direction,
                lc_prec,
            )
        else:
            return calculate_powder_intensities_leb_circ(
                r_t, geom, i_pol, s_pol, lc_prec
            )

    raise ValueError('Unknown method: "{0}".'.format(method))
