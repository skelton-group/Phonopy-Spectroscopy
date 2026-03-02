# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines implementing the calculation of optical spectra from
frequency-dependent dielectric functions."""


# -------
# Imports
# -------


from abc import ABC

import numpy as np
import pandas as pd

from ..constants import (
    ZERO_TOLERANCE,
    SPEED_OF_LIGHT,
    VACUUM_PERMITIVITY,
)

from .dielectric_function import (
    _EPSILON_UNIT_TEXT_LABEL,
    _EPSILON_UNIT_PLOT_LABEL,
)

from ..spectrum_base import SpectrumBase
from ..units import convert_frequency_units

from ..utility.numpy_helper import (
    np_readonly_view,
    np_asarray_copy,
    np_expand_dims,
    np_check_shape,
)


# ---------
# Constants
# ---------


_ALPHA_UNIT_TEXT_LABEL = r"\alpha / cm^-1"
_ALPHA_UNIT_PLOT_LABEL = r"$\alpha$ / cm$^{-1}$"

_SIGMA_UNIT_TEXT_LABEL = r"\sigma / S cm^-1"
_SIGMA_UNIT_PLOT_LABEL = r"$\sigma$ / S cm$^{-1}$"


# -------------------
# Optical eigenvalues
# -------------------


def diagonalise_epsilon(eps):
    """Diagonalise a frequency-dependent D-dimensional dielectric
    function and return the eigenvalues and eigenvectors.

    Parameters
    ----------
    eps : array_like
        Dielectric function (shape: `(O,)` or `(O, D, D)`).

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
        np_asarray_copy(eps, dtype=np.complex128),
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

    # If the eigenvectors are real, strip out the complex part for
    # efficiency.

    if not np.iscomplex(evecs).any():
        evals = np.array(evals.real, dtype=np.float64)

    return (evals, evecs)


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

    s = -1.0j * x * VACUUM_PERMITIVITY * (eps - 1.0)

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

         A(\omega) = - \ln \left[T(\omega) \right]
    """

    (a, r), _ = _validate_setup_optical_property_cals([a, r])

    if t <= 0.0:
        raise ValueError("t must be >= 0.")

    trans_t = (1.0 - r) ** 2 * np.exp(-1.0 * a * t)

    return (trans_t, -1.0 * np.log(trans_t))


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

         T(\omega) = \frac{\left[ 1 - R(\omega) \right]^2 \exp \left[ -\alpha(\omega) t \right ]}{1 - R(\omega)^2 \exp \left[ -2 \alpha(\omega) t \right]}

    .. math::

         A(\omega) = - \ln \left[T(\omega) \right]
    """

    (a, r), _ = _validate_setup_optical_property_cals([a, r])

    if t <= 0.0:
        raise ValueError("t must be >= 0.")

    trans_t = ((1.0 - r) ** 2 * np.exp(-1.0 * a * t)) / (
        1.0 - r**2 * np.exp(-2.0 * a * t)
    )

    return (trans_t, -1.0 * np.log(trans_t))


# -------------------------
# OpticalSpectrumBase class
# -------------------------


class OpticalSpectrumBase(ABC):
    """Abstract base class for optical spectra."""

    # These properties are shared by three different classes, so an
    # abstract base class is used to avoid repetition.

    @property
    def epsilon_unit_text_label(self):
        """str : Dielectric constant unit label suitable for plain-text
        output."""
        return _EPSILON_UNIT_TEXT_LABEL

    @property
    def epsilon_unit_plot_label(self):
        """str : Dielectric constant unit label suitable for plotting
        (contains TeX strings)."""
        return _EPSILON_UNIT_PLOT_LABEL

    @property
    def absorption_coefficient_unit_text_label(self):
        """str : Absorption coefficient unit label suitable for
        plain-text output."""
        return _ALPHA_UNIT_TEXT_LABEL

    @property
    def absorption_coefficient_unit_plot_label(self):
        """str : Absorption coefficient unit label suitable for plotting
        (contains TeX strings)."""
        return _ALPHA_UNIT_PLOT_LABEL

    @property
    def optical_conductivity_unit_text_label(self):
        """str : Optical conductivity unit label suitable for plain-text
        output."""
        return _SIGMA_UNIT_TEXT_LABEL

    @property
    def optical_conductivity_unit_plot_label(self):
        """str : Optical conductivity unit label suitable for plotting
        (contains TeX strings)."""
        return _SIGMA_UNIT_PLOT_LABEL


# ------------------------------
# OpticalEigenmodeSpectrum class
# ------------------------------


class OpticalEigenmodeSpectrum(SpectrumBase, OpticalSpectrumBase):
    """Class for finding the optical eigenmodes of a D-dimensional
    dielectric function and simulating the optical spectra."""

    def __init__(self, x, eps, x_units="thz", t=1.0):
        """Create a new instance of the `OpticalSpectrumBase` class.

        Parameters
        ----------
        x : array_like
            Frequencies in `x_units` (shape: `(O,)`).
        eps : array_like
            D-dimensional dielectric function in units of relative
            permittivity (shape: `(O,) or `(O, D, D)`).
        x_units : str, optional
            Frequency unit (default: "thz").
        t : float or None, optional
            Sample thickness for calculating absorbance and transmission
            in mm (default: 1 mm).
        """

        # Checks and validates x and x_units, and sets the _x and
        # _x_units fields.

        super(OpticalEigenmodeSpectrum, self).__init__(x=x, x_units=x_units)

        if t <= 0.0:
            raise ValueError("t must be >= 0.")

        # Diagonalise dielectric function to find the optical
        # eigenmodes. diagonalise_epsilon() checks the shape of eps and
        # handles "0D" or 1D dielectric functions.

        evals, evecs = diagonalise_epsilon(eps)

        if len(evals) != len(self._x):
            raise ValueError(
                "oe_evals must be an array_like with shape (O,) or "
                "(O, D, D)."
            )

        self._eps = eps
        self._mode_evals = evals
        self._mode_evecs = evecs

        self._t = t

        # Initialise fields for eigenmode optical properties.

        self._mode_n_t = None
        self._mode_a = None
        self._mode_r = None
        self._mode_s = None
        self._mode_l = None

        self._mode_phi = None

        self._mode_trans_int = None
        self._mode_trans_norm = None
        self._mode_trans_incoh = None

        self._mode_abs_int = None
        self._mode_abs_norm = None
        self._mode_abs_incoh = None

    def _lazy_init_optical_properties(self):
        """Lazy initialisation of optical properties."""

        if self._mode_n_t is None:
            n_t, a, r, s, l = optical_spectra_from_epsilon(
                self._x, self._mode_evals, self._x_units
            )

            self._mode_n_t = n_t
            self._mode_a = a
            self._mode_r = r
            self._mode_s = s
            self._mode_l = l

    def _lazy_init_transmission_absorbance(self):
        """Lazy initialisation of absorbance and transmission."""

        self._lazy_init_optical_properties()

        if self._mode_abs_int is None:
            # Thickness in cm for consistency with absorption
            # coefficient.

            t = self._t * 1.0e-1

            self._mode_phi = complex_phase_factor(
                self._x, self._mode_n_t, t, x_units=self._x_units
            )

            self._mode_trans_int, self._mode_abs_int = (
                intrinsic_transmission_absorbance(self._mode_a, t)
            )

            self._mode_trans_norm, self._mode_abs_norm = (
                normal_transmission_absorbance(self._mode_a, self._mode_r, t)
            )

            trans, abs = incoherent_transmission_absorbance(
                self._mode_a, self._mode_r, t
            )

            self._mode_trans_incoh = trans
            self._mode_abs_incoh = abs

    @property
    def epsilon(self):
        """numpy.ndarray : Dielectric function (shape: `(O, D, D)`)."""
        return np_readonly_view(self._eps)

    @property
    def num_dims(self):
        """int : Dimensionality (number of eigenmodes)."""
        return self._mode_evals.shape[1]

    @property
    def mode_eigenvalues(self):
        r"""numpy.ndarray : Eigenvalues of the optical eigenmodes in
        \eps_0 (shape: `(O, D)`).
        """
        return np_readonly_view(self._mode_evals)

    @property
    def mode_eigenvectors(self):
        """numpy.ndarray : Eigenvectors of the optical eigenmodes
        (shape: `(O, D, D)`) (column-major format - the `j`th
        eigenvector at the `i`th frequency is
        `epsilon_eigenvectors[i, :, j]`)."""

        return np_readonly_view(self._mode_evecs)

    @property
    def mode_eigenvectors_row(self):
        """numpy.ndarray : Eigenvectors of the optical eigenmodes
        (shape: `(O, D, D)`) (row-major format - the `j`th eigenvector
        at the `i`th frequency is `epsilon_eigenvectors[i, j]`)
        (this is more natural for many operations)."""

        evecs = np.swap_axes(self._mode_evecs, 1, 2)
        return np_readonly_view(evecs)

    @property
    def mode_refractive_index(self):
        """numpy.ndarray : Complex refractive indices n + ik of the
        optical eigenmodes (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return np_readonly_view(self._mode_n_t)

    @property
    def mode_absorption_coefficient(self):
        """numpy.ndarray : Absorption coefficients of the optical
        eigenmodes in cm^-1 (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return np_readonly_view(self._mode_a)

    @property
    def mode_reflectivity(self):
        """numpy.ndarray : Reflectivity of the optical eigenmodes
        (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return np_readonly_view(self._mode_r)

    @property
    def mode_optical_conductivity(self):
        """numpy.ndarray : Optical conductivity of the optical
        eigenmodes in S cm^-1 (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return np_readonly_view(self._mode_s)

    @property
    def mode_energy_loss_function(self):
        """numpy.ndarray : Energy loss functions of the optical
        eigenmodes (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return np_readonly_view(self._mode_l)

    @property
    def sample_thickness(self):
        """float : Sample thickness in mm."""
        return self._t

    @property
    def mode_complex_phase_factor(self):
        """numpy.ndarray : Complex phase factor of the optical
        eigenmodes at `sample_thickness` (shape: `(O, D))."""

        self._lazy_init_transmission_absorbance()
        return np_readonly_view(self._mode_phi)

    @property
    def mode_intrinsic_transmission(self):
        """numpy.ndarray : Intrinsic (Beer-Lambert) transmission of the
        optical eigenmodes at `sample_thickness` (shape: `(O, D)`)."""

        self._lazy_init_transmission_absorbance()
        return np_readonly_view(self._mode_trans_int)

    @property
    def mode_normal_transmission(self):
        """numpy.ndarray : Normal (single-reflection) transmission of
        the optical eigenmodes at `sample_thickness` (shape: `(O, D)`).
        """

        self._lazy_init_transmission_absorbance()
        return np_readonly_view(self._mode_trans_norm)

    @property
    def mode_incoherent_transmission(self):
        """numpy.ndarray : Incoherent (multiple-reflection) transmission
        of the optical eigenmodes at `sample_thickness` (shape:
        `(O, D)`)."""

        self._lazy_init_transmission_absorbance()
        return np_readonly_view(self._mode_trans_incoh)

    @property
    def mode_intrinsic_absorbance(self):
        """numpy.ndarray : Intrinsic (Beer-Lambert) absorbance of the
        optical eigenmodes at `sample_thickness` (shape: `(O, D)`)."""

        self._lazy_init_transmission_absorbance()
        return np_readonly_view(self._mode_abs_int)

    @property
    def mode_normal_absorbance(self):
        """numpy.ndarray : Normal (single-reflection) absorbance of the
        optical eigenmodes at `sample_thickness` (shape: `(O, D)`)."""

        self._lazy_init_transmission_absorbance()
        return np_readonly_view(self._mode_abs_norm)

    @property
    def mode_incoherent_absorbance(self):
        """numpy.ndarray : Incoherent (multiple-reflection) absorbance
        of the optical eigenmodes at `sample_thickness` (shape:
        `(O, D)`)."""

        self._lazy_init_transmission_absorbance()
        return np_readonly_view(self._mode_abs_incoh)

    def spectrum(self):
        """Return the optical eigenmode eigenvalues and optical spectra
        as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the dielectric function and derived
            quantities.
        """

        self._lazy_init_optical_properties()
        self._lazy_init_transmission_absorbance()

        d = {"freq_energy": self._x}

        for i in range(self.num_dims):
            n = i + 1

            for k, prop in [
                ("epsilon_re", self._mode_evals.real),
                ("epsilon_im", self._mode_evals.imag),
                ("refractive_index", self._mode_n_t.real),
                ("extinction_coefficient", self._mode_n_t.imag),
                ("absorption_coefficient", self._mode_a),
                ("reflectivity", self._mode_r),
                ("optical_conductivity_re", self._mode_s.real),
                ("optical_conductivity_im", self._mode_s.imag),
                ("energy_loss_function", self._mode_l),
                ("phase_factor_re", self._mode_phi.real),
                ("phase_factor_im", self._mode_phi.imag),
                ("intrinsic_transmission", self._mode_trans_int),
                ("intrinsic_absorbance", self._mode_trans_int),
                ("normal_transmission", self._mode_trans_norm),
                ("normal_absorbance", self._mode_abs_norm),
                ("incoherent_transmission", self._mode_trans_incoh),
                ("incoherent_absorbance", self._mode_abs_incoh),
            ]:
                d["{0}_{1}".format(k, n)] = prop[:, i]

        return pd.DataFrame(d)

    @staticmethod
    def from_infrared_dielectric_function(eps_ir, t=1.0):
        """Create an `OpticalEigenmodeSpectrum` from an
        `InfraredDielectricFunction object.

        Parameters
        ----------
        eps_ir : InfraredDielectricFunction
            Infrared dielectric function.
        t : float, optonal
            Thickness in mm (default: 1 mm).

        Returns
        -------
        oe_sp : OpticalEigenmodeSpectrum
            Optical eigenmode spectrum.
        """

        return OpticalEigenmodeSpectrum(
            eps_ir.x, eps_ir.epsilon, x_units=eps_ir.x_units, t=t
        )


# ---------------------
# OpticalSpectrum class
# ---------------------


class OpticalSpectrum(SpectrumBase, OpticalSpectrumBase):
    """Class for simulating unpolarised optical spectra from optical
    eigenmode spectra."""

    def __init__(self, oe_sp, eps_eels=None, ave_w=None):
        """Create a new instance of the `UnpolarisedOpticalSpectrum`
        class.

        Parameters
        ----------
        oe_sp : OpticalEigenmodeSpectrum
            Optical eigenmode spectrum.
        eps_eels : array_like or None, optional
            For 2D optical eigenmode spectra, specify the "inaccessible"
            diagonal element of the dielectric function for calculating
            the electron energy loss spectrum (EELS) (shape: `(O,)`,
            default: `None`).
        ave_w : array_like or None, optional
            Averaging weights for computing mode-average properties
            (shape: `(O, D)`, default: equal weights).
        """

        # This effectively "copies" the x-axis related quantities from
        # the initialising OpticalEigenmodeSpectrum.

        super(OpticalSpectrum, self).__init__(x=oe_sp.x, x_units=oe_sp.x_units)

        self._oe_sp = oe_sp

        if ave_w is not None:
            # Converts (O,) -> (O, 1).

            ave_w, _ = np_expand_dims(
                np_asarray_copy(ave_w, dtype=np.float64), (None, None)
            )

            if not np_check_shape(ave_w, (len(oe_sp.x), oe_sp.num_dims)):
                raise ValueError(
                    "If supplied, ave_w must be an array_like with "
                    "shape (O, D)."
                )

            if not np.allclose(ave_w.sum(axis=-1), 1.0, atol=ZERO_TOLERANCE):
                raise ValueError(
                    "If supplied, the weights in ave_w for each "
                    "frequency must sum to 1."
                )
        else:
            ave_w = (
                np.ones((len(oe_sp.x), oe_sp.num_dims), dtype=np.float64)
            ) / oe_sp.num_dims

        self._eps_eff = (oe_sp.mode_eigenvalues * ave_w).sum(axis=-1)

        self._n_t = (oe_sp.mode_refractive_index * ave_w).sum(axis=-1)
        self._a = (oe_sp.mode_absorption_coefficient * ave_w).sum(axis=-1)
        self._r = (oe_sp.mode_reflectivity * ave_w).sum(axis=-1)
        self._s = (oe_sp.mode_optical_conductivity * ave_w).sum(axis=-1)
        self._l = (oe_sp.mode_energy_loss_function * ave_w).sum(axis=-1)

        # Recalculate absorbance from averaged transmission.

        trans_int = (oe_sp.mode_intrinsic_transmission * ave_w).sum(axis=-1)
        trans_norm = (oe_sp.mode_normal_transmission * ave_w).sum(axis=-1)
        trans_incoh = (oe_sp.mode_incoherent_transmission * ave_w).sum(axis=-1)

        self._trans_int = trans_int
        self._trans_norm = trans_norm
        self._trans_incoh = trans_incoh

        self._abs_int = -1.0 * np.log(trans_int)
        self._abs_norm = -1.0 * np.log10(trans_norm)
        self._abs_incoh = -1.0 * np.log10(trans_incoh)

        self._eels = None

        if eps_eels is not None:
            if self._oe_sp.num_dims != 2:
                raise ValueError(
                    "eps_eels should only be specified when "
                    "oe_sp.num_dims is 2."
                )

            eps_eels = np.asarray(eps_eels, dtype=np.complex128)

            if not np_check_shape(eps_eels, (len(oe_sp.x),)):
                raise ValueError(
                    "If supplied, eps_eels must be an array_like with "
                    "shape (O,)."
                )

            self._eels = (-1.0 / eps_eels).imag

    @property
    def optical_eigenmode_spectrum(self):
        """OpticalEigenmodeSpectrum : Optical eigenmodes used to
        generate the spectra."""
        return self._oe_sp

    @property
    def epsilon(self):
        r"""numpy.ndarray : Dielectric function in \eps_0 (shape:
        (`O,`))."""
        return np_readonly_view(self._eps_eff)

    @property
    def refractive_index(self):
        """numpy.ndarray : Complex refractive index n + ik (shape:
        `(O,)`)."""
        return np_readonly_view(self._n_t)

    @property
    def absorption_coefficient(self):
        """numpy.ndarray : Absorption coefficient in cm^-1 (shape:
        `(O,)`)."""
        return np_readonly_view(self._a)

    @property
    def reflectivity(self):
        """numpy.ndarray : Reflectivity at normal incidence (shape:
        `(O,)`)."""
        return np_readonly_view(self._r)

    @property
    def optical_conductivity(self):
        """numpy.ndarray : Optical conductivity in S m^-1 (shape:
        `(O,)`)."""
        return np_readonly_view(self._s)

    @property
    def energy_loss_function(self):
        """numpy.ndarray : Energy loss function (shape: `(O,)`)."""
        return np_readonly_view(self._l)

    @property
    def sample_thickness(self):
        """float : Sample thickness in mm."""
        return self._oe_sp.sample_thickness

    @property
    def intrinsic_transmission(self):
        """numpy.ndarray : Intrinsic (Beer-Lambert) transmission at
        `sample_thickness` (shape: `(O, D)`)."""
        return np_readonly_view(self._trans_int)

    @property
    def normal_transmission(self):
        """numpy.ndarray : Normal (single-reflection) transmission at
        `sample_thickness` (shape: `(O, D)`).
        """
        return np_readonly_view(self._trans_norm)

    @property
    def incoherent_transmission(self):
        """numpy.ndarray : Incoherent (multiple-reflection) transmission
        at `sample_thickness` (shape: `(O, D)`)."""
        return np_readonly_view(self._trans_incoh)

    @property
    def intrinsic_absorbance(self):
        """numpy.ndarray : Intrinsic (Beer-Lambert) absorbance at
        `sample_thickness` (shape: `(O,)`)."""
        return np_readonly_view(self._abs_int)

    @property
    def normal_absorbance(self):
        """numpy.ndarray : Normal (single-reflection) absorbance at
        `sample_thickness` (shape: `(O, D)`)."""
        return np_readonly_view(self._abs_norm)

    @property
    def incoherent_absorbance(self):
        """numpy.ndarray : Incoherent (multiple-reflection) absorbance
        at `sample_thickness` (shape: `(O, D)`)."""
        return np_readonly_view(self._abs_incoh)

    @property
    def electron_energy_loss_function(self):
        """numpy.ndarray or None : Electron energy loss function (EELS)
        if specified."""

        if self._eels is not None:
            return np_readonly_view(self._eels)
        else:
            return None

    def spectrum(self):
        """Return the simulated dielectric function and derived
        quantities as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the dielectric function and derived
            quantities.
        """

        d = {
            "freq_energy": self._oe_sp.x,
            "epsilon_re": self._eps_eff.real,
            "epsilon_im": self._eps_eff.imag,
            "refractive_index": self._n_t.real,
            "extinction_coefficient": self._n_t.imag,
            "absorption_coefficient": self._a,
            "reflectivity": self._r,
            "optical_conductivity_re": self._s.real,
            "optical_conductivity_im": self._s.imag,
            "energy_loss_function": self._l,
            "intrinsic_transmission": self._trans_int,
            "normal_transmission": self._trans_norm,
            "incoherent_transmission": self._trans_incoh,
            "intrinsic_absorbance": self._abs_int,
            "normal_absorbance": self._abs_norm,
            "incoherent_absorbance": self._abs_incoh,
        }

        return pd.DataFrame(d)


# -----------------------------------
# InputPolarisedOpticalSpectrum class
# -----------------------------------


class InputPolarisedOpticalSpectrum(OpticalSpectrum):
    """Class for simulating optical spectra measured in a collinear
    geometry with polarised incident light, using the optical eigenmodes
    of the "accessible" block of polarisations."""

    def __init__(self, oe_sp, i_pol, eps_eels=None):
        """Create a new instance of the `InputPolarisedOpticalSpectrum`
        class.

        Parameters
        ----------
        oe_sp : OpticalEigenmodeSpectrum
            Optical eigenmode spectrum if the "accessible" 2x2 block of
            the dielectric tensor (must be 2D, i.e. `num_dims == 2`).
        i_pol : Polarisation
            Polarisation of the incident light (must be defined in the
            x/y plane).
        eps_eels : array_like or None, optional
            Specify the "inaccessible" diagonal element of the
            dielectric function for calculating the electron energy loss
            spectrum (EELS) (shape: `(O,)`, default: `None`).

        Notes
        -----
        This workflow implemented in this class assumes a collinear
        measurement along the z-axis and performs calculations using the
        optical eigenmodes of the accessible 2x2 (x/y) block of the
        dielectric function. Therefore, `oe_sp` must be 2D
        (`num_dims == 2`), and vectors in `i_pol` must be in the x/y
        plane.
        """

        if oe_sp.num_dims != 2:
            raise ValueError("oe_sp must be 2D (num_dims = 2).")

        if not i_pol.check_perpendicular("z"):
            raise ValueError(
                "All vectors in i_pol must be defined in the x/y plane "
                "(perpendicular to z)."
            )

        # For polarised incident light we use "power coupling" where
        # the optical spectra are weighted by the absolute squared
        # projections of the polarisation vector onto the inverse of the
        # mode eignvector matrix.

        # (The inverse is used because the dielectric tensors may be
        # non-Hermitian and the eigenvectors are not guaranteed to be
        # orthogonal.)

        sq_coeffs = np.zeros((len(oe_sp.x), oe_sp.num_dims), dtype=np.float64)

        inv_evecs = np.linalg.inv(oe_sp.mode_eigenvectors)

        for v, w in i_pol.iter_v_w():
            sq_coeffs += w * np.abs(np.matmul(inv_evecs, v[:2])) ** 2

        super(InputPolarisedOpticalSpectrum, self).__init__(
            oe_sp, eps_eels=eps_eels, ave_w=sq_coeffs
        )

        self._i_pol = i_pol

    @property
    def incident_polarisation(self):
        """Polarisation : Polarisation of the incident light."""
        return self._i_pol
