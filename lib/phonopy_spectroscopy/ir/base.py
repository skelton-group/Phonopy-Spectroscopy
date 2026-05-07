# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Low-level functions and base classes for simulating optical spectra."""


# -------
# Imports
# -------


import abc
import warnings

import numpy as np
import pandas as pd

from ..constants import (
    ZERO_TOLERANCE,
    SPEED_OF_LIGHT,
    VACUUM_PERMITIVITY,
)

from ..spectrum_base import SpectrumBase

from ..units import convert_frequency_units

from ..utility.numpy_helper import (
    np_expand_dims,
    np_check_shape,
    np_readonly_view,
    np_discard_imag_if_real,
)

from ..utility.diagonalisation import reorder_with_branch_tracking


# ------------------
# Optical eigenmodes
# ------------------


def optical_eigenmodes_from_epsilon(eps, branch_tracking=True):
    """Diagonalise a frequency-dependent D-dimensional dielectric
    function and return the eigenvalues and eigenvectors (optical
    eigenmodes) with optional branch tracking.

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

    evecs = np_discard_imag_if_real(evecs)

    return (evals, evecs)


# -------------------
# Multiphase mixtures
# -------------------


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


def optical_properties_from_epsilon(x, eps, x_units):
    r"""Calculate the optical properties for a frequency-dependent
    dielectric function `eps`: complex refractive index, absorption
    coefficient, optical conductivity and electron energy loss
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
        Tuple of `(n_t, a, s, l)` (all same shape as `eps`; `a` in
        cm^-1, `s` in S cm^-1).

    Notes
    -----
    The optical properties are calculated using the formulae below.

    Complex refractive index:

    .. math::

         \tilde{n}(\omega) = \sqrt(\varepsilon(\omega)) = n(\omega) + ik(\omega)

    Absorption coefficient:

    .. math::

        \alpha(\omega) = \frac{2 \omega}{c} k(\omega)

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

    # Optical conductivity:

    s = 1.0e-2 * -1.0j * 1.0e12 * x * VACUUM_PERMITIVITY * (eps - 1.0)

    # Energy loss function.

    l = (-1.0 / eps).imag

    return (n_t, a, s, l)


def optical_spectra_from_optical_properties(n_t, a, t, n_f=1.0, n_b=1.0):
    r"""Calculate the "standard" optical properties from a complex
    refractive index, absorption coefficient, sample thickness and
    refractive indices of the front/back (entrance/exit)
    media: intrinsic (Beer-Lambert) absorbance, incoherent "single"
    (infinite bulk) and total reflectivity at normal incidence, and
    incoherent transmission.

    Parameters
    ----------
    n_t : array_like
        Complex refractive index.
    a : array_like
        Absorption coefficient (cm^-1).
    t : float
        Sample thickness (mm).
    n_f, n_b : float, optional
        Refractive indices of the front and back media (default: 1.0 =
        vacuum ~ air).

    Returns
    -------
    res : tuple of numpy.ndarray
        Tuple of `(a_int, r_s, r_t, t)` (all same shape as `a`/`n_t`).

    Notes
    -----
    The optical spectra are calculated using the formulae below:

    Intrinsic absorbance:

    .. math::

        A_\mathrm{int} = \exp \left[- \alpha t \right]

    Relectivity from the front and back surfaces at normal incidence:

    .. math::

        R_\mathrm{f} = \left| \frac{n_f - \tilde{n}}{n_f + \tilde{n}} \right|^2

    .. math::

        R_\mathrm{b} = \left| \frac{\tilde{n} - n_b}{\tilde{n} + n_b} \right|^2

    Incoherent "single" reflectivity:

    .. math::

        R_\mathrm{s} = R_\mathrm{f}

    Incoherent total reflectivity:

    .. math::

        R_\mathrm{t} = R_\mathrm{f} + \frac{(1 - R_\mathrm{f})^2 R_\mathrm{b} A_\mathrm{int}^2}{1 - R_\mathrm{f} R_\mathrm{b} A_\mathrm{int}^2}

    Incoherent total transmission:

    .. math::

        T = \frac{(1 - R_\mathrm{f}) (1 - R_\mathrm{b}) A_\mathrm{int}}{1 - R_\mathrm{f} R_\mathrm{b} A_\mathrm{int}^2}
    """

    (n_t, a), _ = _validate_setup_optical_property_cals([n_t, a])

    # Reflectivity of the front and back surfaces at normal incidence.

    r_f = np.abs((n_f - n_t) / (n_f + n_t)) ** 2
    r_b = np.abs((n_t - n_b) / (n_t + n_b)) ** 2

    # Incoherent reflectivity, transmission and absorptance.

    a_int = np.exp(-1.0 * a * 1.0e-1 * t)

    r_incoh = r_f + (
        ((1.0 - r_f) ** 2 * r_b * a_int**2) / (1.0 - r_f * r_b * a_int**2)
    )

    t_incoh = ((1.0 - r_f) * (1.0 - r_b) * a_int) / (
        1.0 - r_f * r_b * a_int**2
    )

    return (a_int, r_f, r_incoh, t_incoh)


# -------------------------
# OpticalSpectrumBase class
# -------------------------


# -------------------------
# OpticalSpectrumBase class
# -------------------------


class OpticalSpectrumBase(SpectrumBase, abc.ABC):
    """Abstract base class for optical spectra."""

    def __init__(self, x, x_units="thz", t=1.0, n_f=1.0, n_b=1.0):
        """Create a new instance of the `OpticalSpectrumBase` class.

        Parameters
        ----------
        t : float, optonal
            Thickness in mm (default: 1 mm).
        n_f, n_b : float, optional
            Refractive indices of the front (indicent) and back (exit)
            media (default: 1.0 = vacuum ~ air).
        """

        super(OpticalSpectrumBase, self).__init__(x=x, x_units=x_units)

        if t <= 0.0:
            raise ValueError("t must be > 0.")

        self._t = t
        self._n_f = n_f
        self._n_b = n_b

        self._ref_s = None
        self._ref_t = None
        self._trans = None
        self._aps = None
        self._abs = None

    @abc.abstractmethod
    def _init_single_reflectivity(self):
        """Set the `_r_s` field."""

        raise NotImplementedError(
            "_init_single_reflectivity() must be implemented in "
            "derived classes."
        )

    @abc.abstractmethod
    def _init_total_reflectivity_and_transmission(self):
        """Set the `_ref_t` and `_trans` fields."""

        raise NotImplementedError(
            "_init_total_reflectivity_and_transmission() must be "
            "implemented in derived classes."
        )

    def _lazy_init_single_reflectivity(self):
        """Lazy initialisation of single reflectivity."""

        if self._ref_s is not None:
            self._init_single_reflectivity()

    def _lazy_init_total_reflectivity_and_transmission(self):
        """Lazy initialisation of total reflectivity, transission, and
        associated quantities."""

        if self._ref_t is not None:
            self._init_total_reflectivity_and_transmission()

            self._aps = 1.0 - (self._ref_t + self._trans)
            self._abs = -1.0 * np.log10(self._trans)

    @property
    def sample_thickness(self):
        """float : Sample thickness in mm."""
        return self._t

    @property
    def front_medium_refractive_index(self):
        """float : Refractive index of the incident (front) medium."""
        return self._n_f

    @property
    def back_medium_refractive_index(self):
        """float : Refractive index of the exit (back) medium."""
        return self._n_b

    @property
    def single_reflectivity(self):
        """numpy.ndarray : Single (infinite bulk) reflectivity (shape:
        `(O,)`)."""

        self._lazy_init_single_reflectivity()
        return np_readonly_view(self._ref_s)

    @property
    def total_reflectivity(self):
        """numpy.ndarray : Total reflectivity (shape: `(O,)`)."""

        self._lazy_init_total_reflectivity_and_transmission()
        return np_readonly_view(self._ref_t)

    @property
    def transmission(self):
        """numpy.ndarray : Transmission (shape: `(O,)`)."""

        self._lazy_init_total_reflectivity_and_transmission()
        return np_readonly_view(self._trans)

    @property
    def absorptance(self):
        """numpy.ndarray : Absorptance (shape: `(O,)`)."""

        self._lazy_init_total_reflectivity_and_transmission()
        return np_readonly_view(self._aps)

    @property
    def absorbance(self):
        """numpy.ndarray : Decadic absorbance (shape: `(O,)`)."""

        self._lazy_init_total_reflectivity_and_transmission()
        return np_readonly_view(self._abs)

    def spectrum(self):
        """Return the optical spectra as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the optical spectra.
        """

        d = {
            "freq_energy": self._x,
            "single_reflectivity": self.single_reflectivity,
            "total_reflectivity": self.total_reflectivity,
            "transmission": self.transmission,
            "absorptance": self.absorptance,
            "absorbance": self.absorbance,
        }

        return pd.DataFrame(d)
