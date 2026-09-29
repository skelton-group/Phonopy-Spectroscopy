# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines implementing single-crystal and powder Raman simulations."""

# -------
# Imports
# -------

import ctypes
import warnings

import numpy as np

from functools import lru_cache

from scipy import LowLevelCallable
from scipy.integrate import quad, nquad, cubature, qmc_quad
from scipy.stats.qmc import Sobol

from ..constants import ZERO_TOLERANCE
from ..distributions import march_dollase
from ..instrument import Polarisation

from ..utility.geometry import (
    direction_cosine,
    rotate_tensors,
    parse_direction,
    rotation_matrix_from_vectors,
)

from ..utility.numpy_helper import np_check_shape, np_expand_dims

from ..utility.quadrature import (
    circle_circle_euler_angle_quad_rule,
    lebedev_circle_euler_angle_quad_rule,
)

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

"""Value of 8 * pi^2."""

_SCIPY_LLC_CFUNC_SIG = types.float64(
    types.intc,
    types.CPointer(types.float64),
    types.CPointer(types.float64),
)

"""Signature of Numba `@cfuncs` for constructing SciPy
`LowLevelCallable` objects."""

# ----------------
# Helper functions
# ----------------


def _check_raman_tensor_transpose_symmetry(r_t):
    r"""Check whether Raman tensors are transpose symmetric (\alpha =
    \alpha^T).

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).

    Returns
    -------
    trans_sym : bool
        `True` if all tensors in `r_t` are transpose symmetric,
        otherwise `False`.
    """

    r_t, _ = np_expand_dims(r_t, (None, 3, 3))
    return (np.abs(r_t - r_t.swapaxes(2, 1)) < ZERO_TOLERANCE).all()


def _validate_params_and_transform_coords(
    r_t, geom, i_pol, s_pol, po_r=1.0, po_norm=None, po_axis=None
):
    """Perform common validation and setup for Raman intensity calculations.

    Parameters
    ----------
    r_t : array_like
        Raman tensors (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    po_r : float or None, optional
        r parameter for preferred orientation (default: 1.0).
    po_norm, po_axis : array_like or str, optional
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`; must be specified if `po_r != 1.0`).

    Returns
    -------
    params : tuple
        Tuple of `(r_t, i_pol, s_pol, po_r, ret_multi)` with Raman
        tensors (shape: `(N, 3, 3)`), incident/scattered
        `Polarisation`s, r parameter for preferred orientation, and
        a flag for whether to return intensities as a `(N,)` array
        (`ret_multi=True`) or scalar (`ret_multi=False`).

    Notes
    -----
    Setting `po_r=None` is guaranteed to bypasses coordinate
    transformations for preferred orientation calculations.
    """

    # Check polarisations are compatible with the measurement geometry.

    if not (
        geom.check_incident_polarisations(i_pol)
        and geom.check_scattered_polarisations(s_pol)
    ):
        raise ValueError(
            "The supplied incident and scattered polarisations are "
            "not possible with the supplied geometry."
        )

    # "Expand" Raman tensors.

    r_t, n_dim_add = np_expand_dims(np.asarray(r_t), (None, 3, 3))

    # Additional setup for calculations including preferred orientation.

    if po_r is not None and np.abs(po_r - 1.0) > ZERO_TOLERANCE:
        # Check reference axis and normal.

        if po_r <= 0.0:
            raise ValueError("po_r must be > 0.")

        if po_norm is None or po_axis is None:
            raise ValueError(
                "po_norm and po_axis must be specified for "
                "calculations including a preferred orientation."
            )

        po_norm = parse_direction(po_norm)
        po_axis = parse_direction(po_axis)

        # Rotate polarisation vectors to align po_axis with +z.

        r_i = rotation_matrix_from_vectors(po_axis, "+z")

        i_pol = Polarisation([r_i @ v for v in i_pol.vectors], i_pol.weights)
        s_pol = Polarisation([r_i @ v for v in s_pol.vectors], s_pol.weights)

        # Rotate Raman tensors to align po_norm with +z in "initial"
        # geometry.

        r_t = rotate_tensors(r_t, rotation_matrix_from_vectors(po_norm, "+z"))

    return (r_t, i_pol, s_pol, po_r, n_dim_add == 0)


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
        Rotation matrix to apply to Raman tensors (shape: `(3, 3)`).

    Returns
    -------
    ints : numpy.ndarray
        Calcuated intensity/intensities (scalar or shape `(N,)`).
    """

    # Check polarisations are valid for measurement geometry and
    # "expand" single tensors to a set of tensors.

    r_t, i_pol, s_pol, _, ret_multi = _validate_params_and_transform_coords(
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
    measurement with powder averaging using the analytical formula.

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
    incident and scattered polarisation, and is only valid if the Raman
    tensors are symmetric and the incident polarisation is perpendicular
    to the collection direction.
    """

    r_t, i_pol, s_pol, _, ret_multi = _validate_params_and_transform_coords(
        r_t, geom, i_pol, s_pol
    )

    # The analytical formula is only valid if the incident polarisation
    # is perpendicular to the collection axis.

    if not i_pol.check_perpendicular(geom.collection_direction):
        raise RuntimeError(
            "The analytical formula can only be used when the "
            "incident polarisation is perpendicular to the collection "
            "axis - use the numerical routines instead."
        )

    # The analytical formula is derived assuming the Raman tensors are
    # transpose symmetric \alpha = \alpha.T. If this is not the case, we
    # warn rather than raise because this will in most cases probably be
    # due to numerical noise, and small errors should not significantly
    # affect the result.

    if not _check_raman_tensor_transpose_symmetry(r_t):
        warnings.warn(
            "One or more Raman tensors are not transpose symmetric - "
            "the analytical result should be checked against one of "
            "the numerical routines.",
            RuntimeWarning,
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
            cos_chi = np.vdot(v_i, v_s) / (
                np.linalg.norm(v_i) * np.linalg.norm(v_s)
            )

            ints[i] += w * (i_per + (i_par - i_per) * cos_chi**2)

    return ints if ret_multi else ints[0]


# -----------------------
# Powder Raman: numerical
# -----------------------


@njit(fastmath=True, inline="always")
def _powder_raman_int_func(
    phi,
    theta,
    psi,
    t,
    v_i_conj,
    v_s,
    po_r,
):
    """Integrand function for numerical powder Raman intensity
    calculations.

    Parameters
    ----------
    phi, theta, psi : float
        Euler angles.
    t : numpy.ndarray
        Raman tensor (shape: `(3, 3)`).
    v_i_conj, v_s : numpy.ndarray
        Conjugated incident and scattered polarisations (shape: `(3,)`).
    po_r : float or None
        r parameter for preferred orientation.

    Returns
    -------
    i : int
        Raman intensity.

    Notes
    -----
    This is a kernel function designed to be called during numerical
    integration loops, and as such has some design features that calling
    code must be aware of.

    If `po_r != 1` the calculations assumes the instrument geometry and
    Raman tensors have been rotated to align the reference axis and
    preferred direction, respectively, with +z.

    If the incident polarisation is complex, `v_i_conj` must be its
    complex conjugate.

    If the Numba JIT-compiled version is used (recommended), if `v_i`
    and/or `v_s` are complex, `t` must also be complex. Otherwise, the
    polarisation vectors will be cast to real, which may lead to errors.
    """

    r = direction_cosine(phi, theta, psi)

    w = 1.0

    if np.abs(po_r - 1.0) > ZERO_TOLERANCE:
        w *= march_dollase(theta, po_r)

    # Explicit data-type matching is required in for the function to
    # compile with @njit when one of t, v_i_conj or v_s are complex.

    v_i_conj = v_i_conj.astype(t.dtype)
    v_s = v_s.astype(t.dtype)

    r = r.astype(t.dtype)

    return w * np.abs(v_i_conj @ r @ t @ r.T @ v_s) ** 2


@njit(fastmath=True, inline="always")
def _powder_raman_int_func_grid(x, t, v_i_conj, v_s, po_r):
    """Integrand function for numerical powder Raman intensity
    calculations, vectorised over arrays of Euler angles.

    Parameters
    ----------
    x : numpy.ndarray
        Euler angles (shape: `(N, 3)`).
    t : numpy.ndarray
        Raman tensor (shape: `(3, 3)`).
    v_i_conj, v_s : numpy.ndarray
        Conjugated incident and scattered polarisations (shape: `(3,)`).
    po_r : float or None
        r parameter for preferred orientation.

    Returns
    -------
    ints : numpy.ndarray
        Calculated Raman intensities for Euler angles specified in `x`
        (shape: `(N,)`).

    Notes
    -----
    This is a kernel function designed to be called during numerical
    integration loops and has the same requirements as
    `_powder_raman_int_func`.

    See Also
    --------
    _powder_raman_int_func : Scalar Raman intensity calculation called
        from this function.
    """

    n, _ = x.shape

    ints = np.zeros((n,), dtype=np.float64)

    for i in range(n):
        ints[i] = _powder_raman_int_func(
            x[i, 0], x[i, 1], x[i, 2], t, v_i_conj, v_s, po_r
        )

    return ints


# -------------------------------------
# Powder Raman: recursive 1D quadrature
# -------------------------------------

if _NUMBA_AVAILABLE:

    @lru_cache(maxsize=2)
    def _powder_raman_make_quad_int_func_nb(real):
        """Generate a compiled integrand function for numerical
        integration with the SciPy `nquad()` routine using Numba.

        Parameters
        ----------
        real : bool
            `True` if the intensity calculation is to be performed on
            a real Raman tensor and polarisation vectors, otherwise
            `False`.

        Returns
        -------
        int_func : cfunc
            Integrand function.

        Notes
        -----
        The returned function object has the signature:

            `double (double *xx, double *user_data)`

        To be used with `nquad()` the function must be wrapped in a
        `LowLevelCallable` configured with the required `user_data`.

        There are two varints of the compiled function depending on the
        value of `real`, which are generated lazily and cached/reused.

        See Also
        --------
        _powder_raman_int_func_quad : Generate an integrand function
            configured with the required `user_data` for passing
            directly to `nquad()`.
        """

        buf_size = 16 if real else 31

        _cfunc_sig = types.float64(
            types.intc,
            types.CPointer(types.float64),
            types.CPointer(types.float64),
        )

        @cfunc(_SCIPY_LLC_CFUNC_SIG, fastmath=True)
        def _powder_raman_quad_int_func_nb(n, xx, user_data):
            psi, theta, phi = xx[0], xx[1], xx[2]

            buf = carray(user_data, (buf_size,), np.float64)

            if real:
                t = buf[:9].reshape((3, 3))
                v_i_conj = buf[9:12]
                v_s = buf[12:15]
                po_r = buf[15]
            else:
                t = buf[:18].view(np.complex128).reshape((3, 3))
                v_i_conj = buf[18:24].view(np.complex128)
                v_s = buf[24:30].view(np.complex128)
                po_r = buf[30]

            return _powder_raman_int_func(
                phi, theta, psi, t, v_i_conj, v_s, po_r
            ) * (np.sin(theta) / _EIGHT_PI_SQUARED)

        return _powder_raman_quad_int_func_nb


def _powder_raman_make_quad_int_func(t, v_i, v_s, po_r):
    """Generate an integrand function for numerical integration with the
    SciPy `nquad()` routine.

    Parameters
    ----------
    t : numpy.ndarray
        Raman tensor (shape: `(3, 3)`).
    v_i_conj, v_s : numpy.ndarray
        Conjugated incident and scattered polarisations (shape: `(3,)`).
    po_r : float or None
        r parameter for preferred orientation.

    Returns
    -------
    int_func : callable or LowLevelCallable
        Integrand function.

    Notes
    -----
    If Numba is available, a `LowLevelCallable` encapsulating a compiled
    C function will be returned. Otherwise, a Lambda function with the
    same signature will be returned.

    The former shows significantly better performance and should be
    preferred wherever possible.
    """

    if not np_check_shape(t, (3, 3)):
        raise ValueError("t must be an array_like with shape (3, 3).")

    if po_r <= 0.0:
        raise ValueError("po_r must be > 0.")

    v_i_conj = v_i.conj()

    if _NUMBA_AVAILABLE:
        real = not (
            np.iscomplexobj(t) or np.iscomplexobj(v_i) or np.iscomplexobj(v_s)
        )

        int_func = _powder_raman_make_quad_int_func_nb(real=real)

        np_dtype = np.float64 if real else np.complex128

        user_data = [
            t.astype(np_dtype).view(np.float64).ravel(),
            v_i_conj.astype(np_dtype).view(np.float64).ravel(),
            v_s.astype(np_dtype).view(np.float64).ravel(),
            np.array([po_r], dtype=np.float64),
        ]

        user_data = np.concatenate(user_data, dtype=np.float64)

        return LowLevelCallable(
            int_func.ctypes,
            user_data=user_data.ctypes.data_as(ctypes.c_void_p),
            signature="double (int, double *, void *)",
        )
    else:
        return lambda psi, theta, phi: _powder_raman_int_func(
            phi, theta, psi, t, v_i_conj, v_s, po_r
        ) * (np.sin(theta) / _EIGHT_PI_SQUARED)


def calculate_powder_raman_intensities_quad(
    r_t, geom, i_pol, s_pol, po_r=1.0, po_norm=None, po_axis=None
):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, with optional preferred
    orientation modelled using the March-Dollase orientation
    distribution function, using the SciPy `nquad()` routine.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    po_r : float, optional
        r parameter for preferred orientation (default: 1.0).
    po_norm, po_axis : array_like or str, optional
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`; must be specified if `po_r != 1.0`).

    Returns
    -------
    res : tuple of (scalar or numpy.ndarray)
        Tuple of `(ints, errs, func_evals)` with the intensities, error
        estimate and numbers of function evaluations (scalar or shape
        `(N,)`).
    """

    r_t, i_pol, s_pol, po_r, ret_multi = _validate_params_and_transform_coords(
        r_t, geom, i_pol, s_pol, po_r=po_r, po_norm=po_norm, po_axis=po_axis
    )

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)
    errs = np.zeros((r_t.shape[0],), dtype=np.float64)
    func_evals = np.zeros((r_t.shape[0],), dtype=int)

    for i, t in enumerate(r_t):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            f = _powder_raman_make_quad_int_func(t, v_i, v_s, po_r)

            res = nquad(
                f,
                [(0.0, 2.0 * np.pi), (0.0, np.pi), (0.0, 2.0 * np.pi)],
                full_output=True,
                opts={"epsabs": ZERO_TOLERANCE, "epsrel": ZERO_TOLERANCE},
            )

        ints[i] += w * res[0]
        errs[i] += w * res[1]
        func_evals[i] += res[2]["neval"]

    return (
        ints if ret_multi else ints[0],
        errs if ret_multi else errs[0],
        func_evals if ret_multi else func_evals[0],
    )


# ----------------------
# Powder Raman: cubature
# ----------------------


class _PowderRamanCubeIntFunc:
    """Integrand function for numerical powder Raman intensity
    calculations using the SciPy `cubature()` routine."""

    def __init__(self, t, v_i, v_s, po_r=1.0):
        """Create a new instance of the
        `_PowderRamanCubeIntFunc` class.

        Parameters
        ----------
        t : numpy.ndarray
            Raman tensor (shape: `(3, 3)`).
        v_i, v_s : numpy.ndarray
            Incident and scattered polarisations (shape: `(3,)`).
        po_r : float, optional
            r parameter for preferred orientation (default: 1.0).
        """

        if not np_check_shape(t, (3, 3)):
            raise ValueError("t must be an array_like with shape (3, 3).")

        if po_r <= 0.0:
            raise ValueError("po_r must be > 0.")

        self._t = t
        self._v_i_conj = v_i.conj()
        self._v_s = v_s
        self._po_r = po_r

        self._func_evals = 0

    def __call__(self, x):
        """Evaluate the integrand function.

        Parameters
        ----------
        x : numpy.ndarray
            Euler angles (shape: `(N, 3)`).

        Returns
        -------
        ints : numpy.ndarray
            Calculated intensities (shape: `(N,)`).
        """

        self._func_evals += x.shape[0]

        return _powder_raman_int_func_grid(
            x, self._t, self._v_i_conj, self._v_s, self._po_r
        ) * (np.sin(x[:, 1]) / _EIGHT_PI_SQUARED)

    @property
    def function_evaluations(self):
        """int : Number of evaluations of the integrand function since
        instantiation."""

        return self._func_evals


def calculate_powder_raman_intensities_cube(
    r_t,
    geom,
    i_pol,
    s_pol,
    po_r=1.0,
    po_norm=None,
    po_axis=None,
):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, with optional preferred
    orientation modelled with the March-Dollase orientation distribution
    function, using the SciPy `cubature()` routine.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    po_r : float, optional
        r parameter for preferred orientation (default: 1.0).
    po_norm, po_axis : array_like or str, optional
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`; must be specified if `po_r != 1.0`).

    Returns
    -------
    res : tuple of (scalar or numpy.ndarray)
        Tuple of `(ints, errs, func_evals)` with the intensities, error
        estimate and numbers of function evaluations (scalar or shape
        `(N,)`).
    """

    r_t, i_pol, s_pol, po_r, ret_multi = _validate_params_and_transform_coords(
        r_t, geom, i_pol, s_pol, po_r=po_r, po_norm=po_norm, po_axis=po_axis
    )

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)
    errs = np.zeros((r_t.shape[0],), dtype=np.float64)
    func_evals = np.zeros((r_t.shape[0],), dtype=int)

    for i, t in enumerate(r_t):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            f = _PowderRamanCubeIntFunc(t, v_i, v_s, po_r)

            res = cubature(
                f,
                [0.0, 0.0, 0.0],
                [2.0 * np.pi, np.pi, 2.0 * np.pi],
                rule="gk21",
                rtol=ZERO_TOLERANCE,
                atol=ZERO_TOLERANCE,
            )

            ints[i] += w * res.estimate
            errs[i] += w * res.error
            func_evals[i] += f.function_evaluations

    return (
        ints if ret_multi else ints[0],
        errs if ret_multi else errs[0],
        func_evals if ret_multi else func_evals[0],
    )


# ------------------------------
# Powder Raman: Lebedev + circle
# ------------------------------


def calculate_powder_raman_intensities_leb_circ(
    r_t, geom, i_pol, s_pol, p, po_r=1.0, po_norm=None, po_axis=None
):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, with optional preferred
    orientation modelled with the March-Dollase orientation distribution
    function, using Lebedev + circle quadrature.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    p : int
        Precision (order) of the Lebedev + circle quadrature scheme.
    po_r : float, optional
        r parameter for preferred orientation (default: 1.0).
    po_norm, po_axis : array_like or str, optional
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`; must be specified if `po_r != 1.0`).

    Returns
    -------
    res : tuple of (scalar or numpy.ndarray)
        Tuple of `(ints, errs, func_evals)` with the intensities, error
        estimate and numbers of function evaluations (scalar or shape
        `(N,)`).

    Notes
    -----
    The errors are set to `np.nan` because errors are not available for
    Lebedev + circle quadrature.
    """

    r_t, i_pol, s_pol, po_r, ret_multi = _validate_params_and_transform_coords(
        r_t, geom, i_pol, s_pol, po_r=po_r, po_norm=po_norm, po_axis=po_axis
    )

    lc_a, lc_w = lebedev_circle_euler_angle_quad_rule(p)

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)

    for i, t in enumerate(r_t):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            # sin \theta term and normalisation by 8 \pi^2 are included
            # in the Lebedev/circle weights.

            lc_ints = _powder_raman_int_func_grid(
                lc_a, t, v_i.conj(), v_s, po_r
            )

            ints[i] += w * (lc_w * lc_ints).sum()

    # Error estimates are not available for Lebedev/circle quadrature,
    # and the number of function evaluations is fixed.

    errs = np.full((r_t.shape[0],), np.nan, dtype=np.float64)

    func_evals = np.full(
        (r_t.shape[0],),
        lc_a.shape[0] * i_pol.vectors.shape[0] * s_pol.vectors.shape[0],
        dtype=int,
    )

    return (
        ints if ret_multi else ints[0],
        errs if ret_multi else errs[0],
        func_evals if ret_multi else func_evals[0],
    )


# -------------------------------
# Powder Raman: quasi-Monte Carlo
# -------------------------------


def _powder_raman_qmc_int_func(x, t, v_i_conj, v_s, po_r):
    """Integrand function for numerical powder Raman intensity
    calculations using the SciPy `qmc_quad()` routine.

    Parameters
    ----------
    x : numpy.ndarray
        Euler angles (shape: `(3,)` or `(3, N)`).
    t : numpy.ndarray
        Raman tensor (shape: `(3, 3)`).
    v_i_conj, v_s : numpy.ndarray
        Conjugated incident and scattered polarisations (shape: `(3,)`).
    po_r : float
        r parameter for preferred orientation.

    Returns
    -------
    ints : numpy.ndarray
        Calculated intensities (shape: `(N,)`).

    See Also
    --------
    _powder_raman_int_func_grid : Scalar Raman intensity calculation
        called from this function.
    """

    # qmc_quad may pass 1D or 2D arguments (shape: (3,), (3, N)), which
    # both need to be reshaped to (N, 3).

    x = x.reshape((-1, 3)) if x.ndim == 1 else x.T

    return _powder_raman_int_func_grid(x, t, v_i_conj, v_s, po_r) * (
        np.sin(x[:, 1]) / _EIGHT_PI_SQUARED
    )


def calculate_powder_raman_intensities_qmc(
    r_t,
    geom,
    i_pol,
    s_pol,
    n_pts,
    n_est=1,
    po_r=1.0,
    po_norm=None,
    po_axis=None,
):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, with optional preferred
    orientation modelled with the March-Dollase orientation distribution
    function, using the SciPy `qmc_quad()` routine.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    n_pts, n_est : int
        Perform integration with `n_pts` points using `n_est`
        independent repeats to estimate the error.
    po_r : float, optional
        r parameter for preferred orientation (default: 1.0).
    po_norm, po_axis : array_like or str, optional
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`; must be specified if `po_r != 1.0`).

    Returns
    -------
    res : tuple of (scalar or numpy.ndarray)
        Tuple of `(ints, errs, func_evals)` with the intensities, error
        estimate and numbers of function evaluations (scalar or shape:
        `(N,)`).

    Notes
    -----
    With `n_est=1` the standard error cannot be estimated and is usually
    set to `np.nan`.

    The number of function evaluations does not include the initial
    "probing" calls performed by `qmc_quad()` and is therefore
    approximate.
    """

    r_t, i_pol, s_pol, po_r, ret_multi = _validate_params_and_transform_coords(
        r_t, geom, i_pol, s_pol, po_r=po_r, po_norm=po_norm, po_axis=po_axis
    )

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)
    errs = np.zeros((r_t.shape[0],), dtype=np.float64)

    for i, t in enumerate(r_t):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            f = lambda x: _powder_raman_qmc_int_func(
                x, t, v_i.conj(), v_s, po_r
            )

            # Sobol sequences are theoretically faster to generate and
            # sample multidimensional variable spaces more uniformly
            # than the default Halton sequence.

            res = qmc_quad(
                f,
                [0.0, 0.0, 0.0],
                [2.0 * np.pi, np.pi, 2.0 * np.pi],
                n_estimates=n_est,
                n_points=n_pts,
                qrng=Sobol(3, scramble=True),
            )

            ints[i] += w * res.integral
            errs[i] += w * res.standard_error

    func_evals = np.full(
        (r_t.shape[0],),
        n_pts * i_pol.vectors.shape[0] * s_pol.vectors.shape[0],
        dtype=int,
    )

    return (
        ints if ret_multi else ints[0],
        errs if ret_multi else errs[0],
        func_evals if ret_multi else func_evals[0],
    )


# ----------------------------------------
# Powder Raman: mixed 1D quadrature/circle
# ----------------------------------------


@njit(fastmath=True, inline="always")
def _powder_raman_quad_circ_int_func(theta, t, v_i_conj, v_s, po_r, phi_psi):
    r"""Integrand function for numerical powder Raman intensity
    calculations using the mixed quadrature/circle scheme with the SciPy
    `quad()` routine.

    Parameters
    ----------
    theta : float
        Euler angle \theta.
    t : numpy.ndarray
        Raman tensor (shape: `(3, 3)`).
    v_i_conj, v_s : numpy.ndarray
        Conjugated incident and scattered polarisations (shape: `(3,)`).
    po_r : float
        r parameter for preferred orientation.
    phi_psi : numpy.ndarray
        \phi and \psi angles for integration with the "inner" circle
        scheme (shape: `(N, 2)`).

    Returns
    -------
    i : int
        Raman intensity.
    """

    m_m = len(phi_psi)

    x = np.zeros((m_m, 3), dtype=np.float64)

    x[:, 0] = phi_psi[:, 0]
    x[:, 1] = theta
    x[:, 2] = phi_psi[:, 1]

    ints = _powder_raman_int_func_grid(x, t, v_i_conj, v_s, po_r)

    return ints.sum() * (np.sin(theta) / (2.0 * m_m))


if _NUMBA_AVAILABLE:

    @lru_cache
    def _powder_raman_make_quad_circ_int_func_nb(real, m):
        r"""Generate a compiled integrand function for numerical
        integration using the mixed quadrature/circle scheme with the
        SciPy `quad()` routine.

        Parameters
        ----------
        real : bool
            `True` if the intensity calculation is to be performed on
            a real Raman tensor and polarisation vectors, otherwise
            `False`.
        m : int
            Number of points for (order of) the circle rule.

        Returns
        -------
        int_func : cfunc
            Integrand function.

        Notes
        -----
        The returned function object has the signature:

            `double (double *xx, double *user_data)`

        To be used with `quad()` the function must be wrapped in a
        `LowLevelCallable` configured with the required `user_data`.

        There are multiple varints of the compiled function depending on
        the value of `real` and the (hard coded) rule for integrating
        over the \phi and \psi (`m`), which are generated lazily and
        cached/reused.

        See Also
        --------
        _powder_raman_make_quad_circ_int_func : Generate an integrand
            function configured with the required `user_data` for
            passing directly to `quad()`.
        """

        buf_size = 16 + (2 * m * m) if real else 31 + (2 * m * m)

        @cfunc(_SCIPY_LLC_CFUNC_SIG, fastmath=True)
        def _powder_raman_quad_circ_int_func_nb(n, xx, user_data):
            buf = carray(user_data, (buf_size,), np.float64)

            if real:
                t = buf[:9].reshape((3, 3))
                v_i_conj = buf[9:12]
                v_s = buf[12:15]
                po_r = buf[15]
                phi_psi = buf[16:].reshape(m * m, 2)
            else:
                t = buf[:18].view(np.complex128).reshape((3, 3))
                v_i_conj = buf[18:24].view(np.complex128)
                v_s = buf[24:30].view(np.complex128)
                po_r = buf[30]
                phi_psi = buf[31:].reshape(m * m, 2)

            return _powder_raman_quad_circ_int_func(
                xx[0], t, v_i_conj, v_s, po_r, phi_psi
            )

        return _powder_raman_quad_circ_int_func_nb


def _powder_raman_make_quad_circ_int_func(t, v_i, v_s, po_r, m):
    r"""Generate an integrand function for numerical integration using
    the mixed quadrature/circle scheme with the SciPy `quad()` routine.

    Parameters
    ----------
    t : numpy.ndarray
        Raman tensor (shape: `(3, 3)`).
    v_i, v_s : numpy.ndarray
        Incident and scattered polarisations (shape: `(3,)`).
    po_r : float or None
        r parameter for preferred orientation.
    m : int
        Number of points for (order of) the circle rule used to
        integrate over the \phi and \psi angles.

    Returns
    -------
    int_func : callable or LowLevelCallable
        Integrand function.

    Notes
    -----
    If Numba is available, a `LowLevelCallable` encapsulating a compiled
    C function will be returned. Otherwise, a Lambda function with the
    same signature will be returned.

    The former shows significantly better performance and should be
    preferred wherever possible.
    """

    if not np_check_shape(t, (3, 3)):
        raise ValueError("t must be an array_like with shape (3, 3).")

    if po_r <= 0.0:
        raise ValueError("po_r must be > 0.")

    v_i_conj = v_i.conj()

    phi_psi, _ = circle_circle_euler_angle_quad_rule(m)

    if _NUMBA_AVAILABLE:
        real = not (
            np.iscomplexobj(t) or np.iscomplexobj(v_i) or np.iscomplexobj(v_s)
        )

        int_func = _powder_raman_make_quad_circ_int_func_nb(real, m)

        np_dtype = np.float64 if real else np.complex128

        user_data = [
            t.astype(np_dtype).view(np.float64).ravel(),
            v_i_conj.astype(np_dtype).view(np.float64).ravel(),
            v_s.astype(np_dtype).view(np.float64).ravel(),
            np.array([po_r], dtype=np.float64),
            phi_psi.view(np.float64).ravel(),
        ]

        user_data = np.concatenate(user_data, dtype=np.float64)

        return LowLevelCallable(
            int_func.ctypes,
            user_data=user_data.ctypes.data_as(ctypes.c_void_p),
            signature="double (int, double *, void *)",
        )
    else:
        return lambda theta: _powder_raman_quad_circ_int_func(
            theta, t, v_i_conj, v_s, po_r, phi_psi
        )


def calculate_powder_raman_intensities_quad_circ(
    r_t, geom, i_pol, s_pol, m, po_r=1.0, po_norm=None, po_axis=None
):
    r"""Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, with optional preferred
    orientation modelled using the March-Dollase orientation
    distribution function, using the mixed quadrature/circle scheme
    with SciPy `quad()`.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    m : int
        Number of points for (order of) the circle rule used to
        integrate over the \phi and \psi angles.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    po_r : float, optional
        r parameter for preferred orientation (default: 1.0).
    po_norm, po_axis : array_like or str, optional
        Surface normal and reference axis for preferred orientation
        (shape: `(3,)`; must be specified if `po_r != 1.0`).

    Returns
    -------
    res : tuple of (scalar or numpy.ndarray)
        Tuple of `(ints, errs, func_evals)` with the intensities, error
        estimate and numbers of function evaluations (scalar or shape
        `(N,)`).
    """

    r_t, i_pol, s_pol, po_r, ret_multi = _validate_params_and_transform_coords(
        r_t, geom, i_pol, s_pol, po_r=po_r, po_norm=po_norm, po_axis=po_axis
    )

    ints = np.zeros((r_t.shape[0],), dtype=np.float64)
    errs = np.zeros((r_t.shape[0],), dtype=np.float64)
    func_evals = np.zeros((r_t.shape[0],), dtype=int)

    for i, t in enumerate(r_t):
        for v_i, v_s, w in i_pol.combine_with_iter(s_pol):
            f = _powder_raman_make_quad_circ_int_func(t, v_i, v_s, po_r, m)

            res = quad(
                f,
                0.0,
                np.pi,
                full_output=True,
                epsabs=ZERO_TOLERANCE,
                epsrel=ZERO_TOLERANCE,
            )

        ints[i] += w * res[0]
        errs[i] += w * res[1]
        func_evals[i] += m**2 * res[2]["neval"]

    return (
        ints if ret_multi else ints[0],
        errs if ret_multi else errs[0],
        func_evals if ret_multi else func_evals[0],
    )


# ---------------------
# Powder Raman: general
# ---------------------


def calculate_powder_raman_intensities(
    r_t, geom, i_pol, s_pol, po_r=1.0, po_norm=None
):
    """Calculate the scalar Raman intensities for a polarised Raman
    measurement with powder averaging, with optional preferred
    orientation modelled with the March-Dollase orientation
    distribution function, using the best available method.

    Parameters
    ----------
    r_t : array_like
        Raman tensor(s) (shape: `(3, 3)` or `(N, 3, 3)`).
    geom : Geometry
        Measurement geometry.
    i_pol, s_pol : Polarisation
        Polarisations of the incident and scattered light.
    po_r : float, optional
        r parameter for preferred orientation (default: 1.0).
    po_norm : array_like or str, optional
        Surface normal for preferred orientation (shape: `(3,)`; must be
        specified if `po_r != 1.0`).

    Returns
    -------
    ints : float or numpy.ndarray
        Calcuated intensity/intensities (scalar or shape `(N,)`).
    """

    r_t = np.asarray(r_t)

    # Determine whether the Raman tensors are transpose symmetric, and
    # whether the calculation is for an isotropic powder (no preferred
    # orientation).

    trans_sym = _check_raman_tensor_transpose_symmetry(r_t)

    isotropic = np.abs(po_r - 1.0) < ZERO_TOLERANCE

    if isotropic:
        if trans_sym and i_pol.check_perpendicular(geom.collection_direction):
            # Analytical formula.

            return calculate_powder_raman_intensities_analytical(
                r_t, geom, i_pol, s_pol
            )
        else:
            # Lebedev + circle quadrature with p = 5.

            ints, _, _ = calculate_powder_raman_intensities_leb_circ(
                r_t, geom, i_pol, s_pol, p=5
            )

            return ints
    else:
        # Mixed 1D quadrature/circle scheme with m = 8.

        ints, _, _ = calculate_powder_raman_intensities_quad_circ(
            r_t,
            geom,
            i_pol,
            s_pol,
            m=8,
            po_r=po_r,
            po_norm=po_norm,
            po_axis=(-1.0 * geom.collection_direction),
        )

        return ints
