# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for modelling the dielectric functions of multiphase
mixtures."""

# -------
# Imports
# -------

import warnings

import numpy as np

from ..constants import ZERO_TOLERANCE

# ---------
# Functions
# ---------


def _validate_setup_multiphase_mixture(eps, fracs):
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


def bruggeman_two_phase_scalar(eps_a, eps_b, f_a, f_b):
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

    (eps_a, eps_b), ret_sc = _validate_setup_multiphase_mixture(
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


def bruggeman_three_phase_scalar(eps_a, eps_b, eps_c, f_a, f_b, f_c):
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

    (eps_a, eps_b, eps_c), ret_sc = _validate_setup_multiphase_mixture(
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


def bruggeman_multiphase_scalar(
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

    eps, ret_sc = _validate_setup_multiphase_mixture(eps, fracs)

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


def bruggeman_scalar(eps, f, eps_mix, f_mix):
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
        return bruggeman_two_phase_scalar(eps, eps_mix[0], f, f_mix[0])

    if len(eps_mix) == 2:
        return bruggeman_three_phase_scalar(eps, *eps_mix, f, *f_mix)

    else:
        return bruggeman_multiphase_scalar(
            [eps] + [e for e in eps_mix], [f] + [f for f in f_mix]
        )
