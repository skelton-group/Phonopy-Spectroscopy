# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines for simulating optical spectra."""


# -------
# Imports
# -------


import abc
import warnings

import numpy as np

from ..constants import (
    ZERO_TOLERANCE,
    SPEED_OF_LIGHT,
    VACUUM_PERMITIVITY,
)

from ..units import convert_frequency_units

from ..utility.numpy_helper import (
    np_expand_dims,
    np_check_shape,
)

from ..utility.diagonalisation import reorder_with_branch_tracking


# -------------------
# Optical eigenvalues
# -------------------


def diagonalise_epsilon(eps, branch_tracking=True):
    """Diagonalise a frequency-dependent D-dimensional dielectric
    function and return the eigenvalues and eigenvectors.

    Parameters
    ----------
    eps : array_like
        Dielectric function (shape: `(O,)` or `(O, D, D)`).
    branch_tracking : bool, optional
        If `True`, reorder the optical eigenmodes by performing a
        "branch tracking" over frequencies (default: `True`).

    Returns
    -------
    res : tuple of numpy.ndarray
        `(evals, evecs)` tuple (shapes: `(O, D)`, `(O, D, D)`).

    Notes
    -----
    "0D" dielectric functions with (shape: `(O,)`) will be reshaped to
    `(O, 1, 1)`. In this case, the returned eigenvalues are just `eps`
    (0D) or `eps[:, 0, 0]` (1D), and the eigenvectors are
    `np.ones((len(eps),))`.

    This function applies the following post-processing:

    * If the eigenvectors are all real, the returned `evecs` are
      converted to `np.float64` instead of `np.complex128`.
    """

    eps, _ = np_expand_dims(
        np.asarray(eps, dtype=np.complex128),
        (None, None, None),
        expand_order=(2, 1),
    )

    _, d_1, d_2 = eps.shape

    if d_1 != d_2:
        raise ValueError(
            "The two outermost dimensions of eps must be the same size."
        )

    if d_1 == 1:
        return (eps[:, :, 0], np.ones((d_1, 1, 1), dtype=np.float64))

    evals, evecs = np.linalg.eig(eps)

    # If branch_tracking is set, reorder evals/evecs.

    if branch_tracking:
        evals, evecs = reorder_with_branch_tracking(evals, evecs)

    # If the eigenvectors are real, strip out the complex part for
    # efficiency.

    if not np.iscomplex(evecs).any():
        evecs = np.array(evecs.real, dtype=np.float64)

    return (evals, evecs)


# -------------
# Pellet models
# -------------


def _validate_setup_multiphase_mixture_calcs(eps, fracs):
    """Perform common validation and setup for the mixture modelling
    routines.

    Parameters
    ----------
    eps : list of [float, complex or array_like]
        Dielectric constants/functions.
    fracs : list of float
        Volume fractions.

    Returns
    -------
    res : tuple of (list of float, complex or numpy.ndarray, bool)
        `eps` with dielectric functions converted to `numpy.ndarray`
        objects and a Boolean flag indicating whether to return a
        scalar value.
    """

    if len(fracs) != len(eps):
        raise ValueError("fracs and eps must have the same number of entries.")

    # Determine the shape (O,) for "aligning" the \eps.

    o = None

    for e in eps:
        d = np.ndim(e)

        if d != 0:
            if d != 1:
                raise ValueError(
                    "All eps must be a scalar or an array_like with "
                    "shape (O,)."
                )

            if o is None:
                o = len(e)
            else:
                if len(e) != o:
                    raise ValueError(
                        "If two or more eps are supplied as 1D "
                        "dielectric functions, they must have the same "
                        "shape (O,)."
                    )

    ret_sc = o is None

    if o is None:
        o = 1

    eps_align = []

    for i, e in enumerate(eps):
        if np.ndim(e) == 0:
            eps_align.append(np.full((o,), e, dtype=np.complex128))
        else:
            eps_align.append(np.asarray(e, dtype=np.complex128))

    # Check volume fractions are between 0 and 1, and sum to 1.

    for f in fracs:
        if f < 0.0 or f > 1.0:
            raise ValueError("Volume fractions must be between 0 and 1.")

    if np.abs(np.sum(fracs) - 1.0) > ZERO_TOLERANCE:
        raise ValueError("Volume fractions must sum to one.")

    # Return converted eps and a ret_scalar flat.

    return (eps_align, ret_sc)


def bruggeman_two_phase_mixture(eps_a, eps_b, f_a, f_b):
    """Calculate the effective dielectric response of a two-phase
    mixture using the Bruggeman model.

    Parameters
    ----------
    eps_a, eps_b : float, complex or array_like
        Dielectric constant (scalar) or function (shape: `(O,)`) of the
        two components.
    f_a, f_b : float
        Volume fractions of the two components.

    Returns
    -------
    eps_eff : complex or array_like
        Effective dielectric constant (scalar) or function (shape:
        `(O,)`).
    """

    (eps_a, eps_b), ret_sc = _validate_setup_multiphase_mixture_calcs(
        [eps_a, eps_b], [f_a, f_b]
    )

    # For a two-phase mixture, we need to set up and solve a quadratic
    # equation of the form a x^2 + b x + c = 0.

    coeffs = np.zeros((len(eps_a), 3), dtype=np.complex128)

    coeffs[:, 0] = 2.0

    coeffs[:, 1] = -1.0 * (
        (((3.0 * f_a) - 1.0) * eps_a) + ((3.0 * (1.0 - f_a) - 1.0) * eps_b)
    )

    coeffs[:, 2] = -1.0 * eps_a * eps_b

    eps_eff = np.zeros((len(eps_a),), dtype=np.complex128)

    for i, c in enumerate(coeffs):
        roots = np.roots(c)

        # The root with Im(eps) > 0 corresponds to the material
        # absorbing energy. If there isn't one, we fall back to
        # the root with the largest real part.

        mask = roots.imag > -1.0 * ZERO_TOLERANCE

        if mask.any() > 0:
            roots = roots[mask]

        eps_eff[i] = roots[np.argmax(roots.real)]

    return eps_eff[0] if ret_sc else eps_eff


def bruggeman_three_phase_mixture(eps_a, eps_b, eps_c, f_a, f_b, f_c):
    """Calculate the effective dielectric response of a three-phase
    mixture using the Bruggeman model.

    Parameters
    ----------
    eps_a, eps_b, eps_c : float, complex or array_like
        Dielectric constant (scalar) or function (shape: `(O,)`) of the
        three components.
    f_a, f_b, f_c : float
        Volume fractions of the three components.

    Returns
    -------
    eps_eff : complex or array_like
        Effective dielectric constant (scalar) or function (shape:
        `(O,)`).
    """

    (eps_a, eps_b, eps_c), ret_sc = _validate_setup_multiphase_mixture_calcs(
        [eps_a, eps_b, eps_c], [f_a, f_b, f_c]
    )

    # For a three-phase mixture, we need to set up and solve a cubic
    # equation of the form a x^3 + b x^2 + c x + d = 0.

    coeffs = np.zeros((len(eps_a), 4), dtype=np.complex128)

    coeffs[:, 0] = 4.0

    coeffs[:, 1] = (
        eps_a * (2.0 - 6.0 * f_a)
        + eps_b * (2.0 - 6.0 * f_b)
        + eps_c * (2.0 - 6.0 * f_c)
    )

    coeffs[:, 2] = (
        eps_a * eps_b * (3.0 * f_c - 2.0)
        + eps_a * eps_c * (3.0 * f_b - 2.0)
        + eps_b * eps_c * (3.0 * f_a - 2.0)
    )

    coeffs[:, 3] = -1.0 * eps_a * eps_b * eps_c

    eps_eff = np.zeros((len(eps_a),), dtype=np.complex128)

    for i, c in enumerate(coeffs):
        roots = np.roots(c)

        # Same root-selection logic as for the two-phase case.

        mask = roots.imag > -1.0 * ZERO_TOLERANCE

        if mask.any() > 0:
            roots = roots[mask]

        eps_eff[i] = roots[np.argmax(roots.real)]

    return eps_eff[0] if ret_sc else eps_eff


def bruggeman_multiphase_mixture(
    eps, fracs, tol=1.0e-8, max_iter=1000, damp=1.0e-4
):
    """Calculate the effective dielectric response of a general
    multiphase mixture using the Bruggeman model.

    Parameters
    ----------
    eps : array_like
        Dielectric constants (scalar) or functions (shape: `(O,)`)
        of the components.
    fracs : array_like
        Volume fractions of the components.
    tol : float, optional
        Maximum convergence tolerance (default: 1.0e-8).
    max_iter : int, optional
        Maximum number of iterations for the Newton-Rhapson
        optimisation (default: 1,000)
    damp : float, optional
        Damping constant for determining Newton-Rhapson step (default:
        1.0e-4).

    Returns
    -------
    eps_eff : complex or array_like
        Effective dielectric constant (scalar) or function (shape:
        `(O,)`).
    """

    eps, ret_sc = _validate_setup_multiphase_mixture_calcs(eps, fracs)

    if tol <= 0.0:
        raise ValueError("tol must be > 0.")

    if max_iter <= 0:
        raise ValueError("max_iter must be > 0.")

    if damp < 0.0:
        raise ValueError("damp must be >= 0.")

    # Don't do anything if we don't need to.

    if len(eps) == 1:
        (eps,) = eps
        return eps[0] if ret_sc else eps

    eps = np.array(eps, dtype=np.complex128)
    fracs = np.asarray(fracs, dtype=np.float64)

    eps_eff = (eps * fracs[:, np.newaxis]).sum(axis=0)

    conv = False

    for _ in range(max_iter):
        d = eps + 2.0 * eps_eff[np.newaxis, :]

        g = np.sum(
            (fracs[:, np.newaxis] * (eps - eps_eff[np.newaxis, :])) / d, axis=0
        )

        g_prime = np.sum((-3.0 * fracs[:, np.newaxis] * eps) / d**2, axis=0)

        step = (g * g_prime.conj()) / (np.abs(g_prime) ** 2 + damp)

        eps_eff -= step

        if (np.abs(step) < tol).all():
            conv = True
            break

    if not conv:
        warnings.warn(
            "Newton-Rhapson failed to converge in max_iter = {0:,}."
            "".format(max_iter),
            RuntimeWarning,
        )

    return eps_eff[0] if ret_sc else eps_eff


def bruggeman_mixture(eps, f, eps_mix, f_mix):
    """Calculate the effective dielectric response of a mixture using
    the Bruggeman method.

    Parameters
    ----------
    eps : float, complex or array_like
        Dielectric constant (scalar) or function (shape: `(O,)`) of the
        first component.
    f : float
        Volume fraction of the first component.
    eps_mix : array_like
        Dielectric constants/functions of the other components.
    f_mix : array_like
        Volume fractions of the other components.

    Returns
    -------
    eps_eff : complex or array_like
        Effective dielectric constant (scalar) or function (shape:
        `(O,)`).

    Notes
    -----
    This is a convenience function that automatically selects from three
    implementations of the Bruggeman model depending on the number of
    components in the mixture.

    The Bruggeman model treats all components equally, so `eps`/`f` do
    not need to specify the largest component.

    See Also
    --------
    bruggeman_two_phase_mixture
    bruggeman_three_phase_mixture
    bruggeman_multiphase_mixture
        Implementations of the Bruggeman model for two-, three- and
        multiphase mixtures.
    """

    if len(eps_mix) != len(f_mix):
        raise ValueError(
            "eps_mix and f_mix must have the same number of entries."
        )

    if len(eps_mix) == 0:
        if f < 0.0 or np.abs(f - 1.0) > ZERO_TOLERANCE:
            raise ValueError("For a single component, f must be equal to 1.")

        return eps

    # Each of the bruggeman* methods will perform parameter validation.

    if len(eps_mix) == 1:
        return bruggeman_two_phase_mixture(eps, eps_mix[0], f, f_mix[0])

    if len(eps_mix) == 2:
        return bruggeman_three_phase_mixture(eps, *eps_mix, f, *f_mix)

    else:
        return bruggeman_multiphase_mixture(
            [eps] + [e for e in eps_mix], [f] + [f for f in f_mix]
        )


# ------------------
# Optical properties
# ------------------


def _validate_setup_optical_property_cals(opt_props, x=None, x_units="thz"):
    """Perform common validation and setup for optical-property
    calculation routines.

    Parameters
    ----------
    opt_props : list of array_like
        List of calculated optical properties (shapes: `(O,)` or
        `(O, D)`).
    x : array_like or None, optional
        Frequencies (shape: `(O,)`).
    x_units : str, optional
        Frequency units (default: "thz").

    Returns
    -------
    res : tuple of numpy.ndarray or None
        Tuple of `(opt_props, x)` with a list of optical properties and
        frequencies in rad THz suitable for broadcasting against the
        arrays in `opt_props`.
    """

    # Shape of outer dimension of optical properties.

    if x is not None:
        x = np.asarray(x)

        if not np_check_shape(x, (None,)):
            raise ValueError("x must be an array_like with shape (O,).")

        # Convert x to rad THz.

        x = convert_frequency_units(x, x_units, "rad_thz")

    opt_props = [np.asarray(p) for p in opt_props]

    # Require outer dimension of all optical properties to match the
    # number of frequencies.

    shape_o = len(x) if x is not None else None

    p_1 = opt_props[0]

    if not (
        np_check_shape(p_1, (shape_o,)) or np_check_shape(p_1, (shape_o, None))
    ):
        raise ValueError(
            "Optical properties must be array_likes with shape (O,) or "
            "(O, D)."
        )

    for p_n in opt_props[1:]:
        if not np_check_shape(p_n, p_1.shape):
            raise ValueError(
                "All optical properties must have the same shape."
            )

    # If the optical properties have shape (O, D), "expand" x for
    # broadcasting.

    if x is not None and p_1.ndim == 2:
        x = x[:, np.newaxis]

    return (opt_props, x)


def optical_spectra_from_epsilon(x, eps, x_units):
    r"""Calculate the optical spectra for a frequency-dependent
    dielectric function `eps`: complex refractive index, absorption
    coefficient, reflectivity, optical conductivity and energy loss
    function.

    Parameters
    ----------
    x : array_like
        Frequencies (shape: `(O,)`).
    eps : array_like
        Scalar dielectric function or eigenvalues of a D-dimensional
        dielectric functions in units of relative permittivity (shape:
        `(O,)` or `(O, D)`).
    x_units : str
        Frequency unit.


    Returns
    -------
    res : tuple of numpy.ndarray
        Tuple of `(n_t, a, r, s, l)` (all same shape as `eps`; `a` in
        cm^-1, `s` in S cm^-1).

    Notes
    -----
    The optical spectra are calculated using the formulae below.

    Complex refractive index:

    .. math::

         \tilde{n}(\omega) = \sqrt(\varepsilon(\omega)) = n(\omega) + ik(\omega)

    Absorption coefficient:

    .. math::

        \alpha(\omega) = \frac{2 \omega}{c} k(\omega)

    Reflectivity at normal incidence:

    .. math::

        R(\omega) = \left| \frac{\tilde{n} - 1}{\tilde{n} + 1} \right|^2

    Optical conductivity:

    .. math::

        \sigma(\omega) = -i \omega \varepsilon_0 \left[ \varepsilon(\omega) - 1 \right]

    Loss function:

    .. math::

        L(\omega) = \Im \left[ -\frac{1}{\varepsilon(\omega)} \right]
    """

    (eps,), x = _validate_setup_optical_property_cals(
        [eps], x=x, x_units=x_units
    )

    # Complex refractive index.

    n_t = np.sqrt(eps)

    # Absorption coefficient.

    a = 1.0e-2 * (2.0 * 1.0e12 * x * n_t.imag) / SPEED_OF_LIGHT

    # Reflectivity at normal incidence.

    r = np.abs((n_t - 1) / (n_t + 1)) ** 2

    # Optical conductivity:

    s = 1.0e-2 * -1.0j * 1.0e12 * x * VACUUM_PERMITIVITY * (eps - 1.0)

    # Energy loss function.

    l = (-1.0 / eps).imag

    return (n_t, a, r, s, l)


def complex_phase_factor(x, n_t, t, x_units="thz"):
    r"""Calculate the complex phase factor from a complex refractive
    index and sample thickness.

    Params
    ------
    x : array_like
        Frequencies in `x_units` (shape: `(O,)`).
    n_t : array_like
        Complex refractive index (shape: `(O,)`).
    t : float
        Sample thickness in mm.
    x_units : st, optional
        Units of `x` (default: "thz").

    Returns
    -------
    phi : array_like
        Complex phase factor (shape: `(O,)`).

    Notes
    -----
    The complex phase factor is calculated as:

    .. math::

        \phi(\omega) = \exp \left[ \frac{i \omega t}{c} \tilde{n}(\omega) \right]
    """

    (n_t,), x = _validate_setup_optical_property_cals(
        [n_t], x=x, x_units=x_units
    )

    return np.exp(-1.0j * ((1.0e12 * x * t * 1.0e-3) / SPEED_OF_LIGHT) * n_t)


def intrinsic_transmission_absorbance(a, t):
    r"""Calculate the "intrinsic" (Beer-Lambert) transmission and
    absorbance from an absorption coefficient and sample thickness.

    Params
    ------
    a : array_like
        Frequency-dependent absorption coefficient (shape: `(O,)` or
        `(O, D)`).
    t : float
        Sample thickness (same distance units as `a`).

    Returns
    -------
    res : tuple of numpy.ndarray
        Tuple of `(trans_t, abs_t)` (shape: `(O,)`).

    Notes
    -----
    The intrinsic absorbance and transmission are calculated as follows:

    .. math::

         A(\omega) = \alpha(\omega) t

    .. math::

         T(\omega) = \exp \left[ - A(\omega) \right ]
    """

    (a,), _ = _validate_setup_optical_property_cals([a])

    if t <= 0.0:
        raise ValueError("t must be >= 0.")

    abs_t = a * t
    trans_t = np.exp(-1.0 * abs_t)

    return (trans_t, abs_t)


def normal_transmission_absorbance(a, r, t):
    r"""Calculate the "normal" (single reflection) transmission and
    absorbance from an absorption coefficient, reflectivity and sample
    thickness.

    Params
    ------
    a, r : array_like
        Frequency-dependent absorption coefficient and reflectivity
        (shape: `(O,)` or `(O, D)`).
    t : float
        Sample thickness (same distance units as `a`).

    Returns
    -------
    res : tuple of numpy.ndarray
        Tuple of `(abs_t, trans_t)` (shape: `(O,)`).

    Notes
    -----
    The normal transmission and absorbance are calculated as follows:

    .. math::

         T(\omega) = \left[ 1 - R(\omega) \right]^2 \exp \left[ -\alpha(\omega) t \right ]

    .. math::

         A(\omega) = -  \log_{10} \left[T(\omega) \right]
    """

    (a, r), _ = _validate_setup_optical_property_cals([a, r])

    if t <= 0.0:
        raise ValueError("t must be >= 0.")

    trans_t = (1.0 - r) ** 2 * np.exp(-1.0 * a * t)

    return (trans_t, -1.0 * np.log10(trans_t))


def incoherent_transmission_absorbance(a, r, t):
    r"""Calculate the "incoherent" (multiple reflection) transmission
    and absorbance from an absorption coefficient, reflectivity and
    sample thickness.

    Params
    ------
    a, r : array_like
        Frequency-dependent absorption coefficient and reflectivity
        (shape: `(O,)` or `(O, D)`).
    t : float
        Sample thickness (same distance units as `a`).

    Returns
    -------
    res : tuple of numpy.ndarray
        Tuple of `(abs_t, trans_t)` (shape: `(O,)`).

    Notes
    -----
    The transmission and absorbance are calculated as follows:

    .. math::

         T(\omega) = \frac{\left[ 1 - R(\omega) \right]^2 \exp \left[ -\alpha(\omega) t \right ]}{1 - R^2(\omega) \exp \left[ -2 \alpha(\omega) t \right]}

    .. math::

         A(\omega) = -  \log_{10} \left[T(\omega) \right]
    """

    (a, r), _ = _validate_setup_optical_property_cals([a, r])

    if t <= 0.0:
        raise ValueError("t must be >= 0.")

    trans_t = ((1.0 - r) ** 2 * np.exp(-1.0 * a * t)) / (
        1.0 - r**2 * np.exp(-2.0 * a * t)
    )

    return (trans_t, -1.0 * np.log10(trans_t))
