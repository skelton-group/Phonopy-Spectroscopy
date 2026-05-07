# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines for finding the optical eigenmodes of complex dielectric
functions and calculating the associated optical properties and spectra.
"""


# -------
# Imports
# -------


import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

from .dielectric_function import (
    _EPSILON_UNIT_TEXT_LABEL,
    _EPSILON_UNIT_PLOT_LABEL,
)

from .base import (
    optical_eigenmodes_from_epsilon,
    bruggeman_mixture,
    optical_properties_from_epsilon,
    optical_spectra_from_optical_properties,
    OpticalSpectrumBase,
)

from ..constants import ZERO_TOLERANCE

from ..instrument import Polarisation

from ..spectrum_base import SpectrumBase

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


# -----------------------
# OpticalEigenmodes class
# -----------------------


class OpticalEigenmodes(SpectrumBase):
    """Class for finding the optical eigenmodes of a D-dimensional
    dielectric function and simulating the optical spectra."""

    def __init__(
        self,
        x,
        eps,
        x_units="thz",
        branch_tracking=True,
        p_vol_frac=1.0,
        p_binder_eps=1.0,
        p_den=1.0,
    ):
        """Create a new instance of the `OpticalEigenmodes` class.

        Parameters
        ----------
        x : array_like
            Frequencies in `x_units` (shape: `(O,)`).
        eps : array_like
            D-dimensional dielectric function in units of relative
            permittivity (shape: `(O,) or `(O, D, D)`).
        x_units : str, optional
            Frequency unit (default: "thz").
        branch_tracking : bool, optional
            If `True`, reorder the optical eigenmodes by performing a
            "branch tracking" over frequencies (default: `True`).
        p_vol_frac : float, optional
            Volume fraction of material in a pellet (default: 1.0).
        p_binder_eps : float or tuple of numpy.ndarray, optional
            Dielectric constant or tuple of `(x, eps_x)` specifying the
            frequency-dependent dielectric function of the pellet
            "binder" material (default: 1.0 = vacuum ~ air).
        p_den : float, optional
            Density of the pellet (default: 1.0).

        Notes
        -----
        Specifying `p_vol_frac` < 1.0 or `p_den` < 1.0 simulates a
        pellet by modifying the dielectric response along the principal
        axes (i.e. optical eigenmode eigenvalues) according to the
        Bruggeman model.

        With `p_den=1.0`, the pellet is treated as a two-phase system
        with the dielectric constant or scalar (1D) frequency-dependent
        dielectric function of the binder set with `p_binder_eps`.

        With `p_den` < 1.0, the pellet is treated as a three-phase
        system with vacuum ~ air as the third medium.

        Different parameter combinations can therefore be used to model
        low-density pure powders, fully dense pellets, and low-density
        pellets.

        If `p_binder_eps` specifies a frequency-dependent dielectric
        function, it will be interpolated to the same frequency axis as
        `eps`. It must therefore be in the same `x_units` as `eps`, and
        must span the same or wider frequency range.
        """

        # Checks and validates x and x_units, and sets the _x and
        # _x_units fields.

        super(OpticalEigenmodes, self).__init__(x=x, x_units=x_units)

        # Diagonalise dielectric function to find the optical
        # eigenmodes. diagonalise_epsilon() checks the shape of eps and
        # handles "0D" or 1D dielectric functions.

        evals, evecs = optical_eigenmodes_from_epsilon(
            eps, branch_tracking=branch_tracking
        )

        if len(evals) != len(self._x):
            raise ValueError(
                "oe_evals must be an array_like with shape (O,) or "
                "(O, D, D)."
            )

        if p_vol_frac < 0.0 or p_vol_frac > 1.0:
            raise ValueError("p_vol_frac must be between 0 and 1.")

        if p_den < 0.0 or p_den > 1.0:
            raise ValueError("p_den must be between 0 and 1.")

        if p_vol_frac < 1.0 or p_den < 1.0:
            # Modify the eigenvalues for a mixed pellet.

            eps_mix, f_mix = [], []

            if p_vol_frac < 1.0:
                # Include the "host" specified by p_binder_eps in the
                # mixture.

                if np.ndim(p_binder_eps) != 0:
                    x, eps = p_binder_eps

                    if x.min() > self._x.min() or x.max() < self._x.max():
                        raise ValueError(
                            "If eps_binder specifies a frequency-dependent "
                            "dielectric function, it must be calculated "
                            "over the same or a wider energy range than "
                            "the material dielectric function."
                        )

                    interp_func = interp1d(x, eps, kind="linear")
                    eps_mix.append(interp_func(self._x))
                else:
                    eps_mix.append(p_binder_eps)

                f_mix.append(p_den * (1.0 - p_vol_frac))

            if p_den < 1.0:
                p_vol_frac = p_den * p_vol_frac

                # Include air in the mixture.

                eps_mix.append(1.0)
                f_mix.append(1.0 - p_den)

            evals_new = np.zeros_like(evals, dtype=np.complex128)

            for i in range(evals.shape[1]):
                # bruggeman_mixture() automatically selects the two-
                # or three-phase analytical model depending on the
                # number of additional components in eps_mix/f_mix.

                evals_new[:, i] = bruggeman_mixture(
                    evals[:, i], p_vol_frac, eps_mix, f_mix
                )

            evals = evals_new

        self._mode_evals = evals
        self._mode_evecs = evecs

        # Initialise fields for eigenmode optical properties.

        self._mode_n_t = None
        self._mode_a = None
        self._mode_s = None
        self._mode_l = None

        self._dual_basis = None

    def _lazy_init_optical_properties(self):
        """Lazy initialisation of optical properties."""

        if self._mode_n_t is None:
            n_t, a, s, l = optical_properties_from_epsilon(
                self._x, self._mode_evals, self._x_units
            )

            self._mode_n_t = n_t
            self._mode_a = a
            self._mode_s = s
            self._mode_l = l

    def _lazy_init_dual_basis(self):
        """Lazy initialisation of dual basis."""

        if self._dual_basis is None:
            self._dual_basis = np.linalg.inv(self._mode_evecs)

    @property
    def num_dims(self):
        """int : Dimensionality (number of eigenmodes)."""
        return self._mode_evals.shape[1]

    @property
    def eigenvalues(self):
        r"""numpy.ndarray : Eigenvalues of the optical eigenmodes in
        \eps_0 (shape: `(O, D)`).
        """
        return np_readonly_view(self._mode_evals)

    @property
    def eigenvectors(self):
        """numpy.ndarray : Eigenvectors of the optical eigenmodes
        (shape: `(O, D, D)`) (column-major format - the `j`th
        eigenvector at the `i`th frequency is
        `epsilon_eigenvectors[i, :, j]`)."""

        return np_readonly_view(self._mode_evecs)

    @property
    def eigenvectors_row(self):
        """numpy.ndarray : Eigenvectors of the optical eigenmodes
        (shape: `(O, D, D)`) (row-major format - the `j`th eigenvector
        at the `i`th frequency is `epsilon_eigenvectors[i, j]`)
        (this is more natural for many operations)."""

        evecs = np.swapaxes(self._mode_evecs, 1, 2)
        return np_readonly_view(evecs)

    @property
    def dual_basis(self):
        """numpy.ndarray : Inverse of the mode eigenvectors (dual basis)
        for eigenmode projections (shape: `(O, D, D)`)."""

        self._lazy_init_dual_basis()
        return np_readonly_view(self._dual_basis)

    @property
    def refractive_index(self):
        """numpy.ndarray : Complex refractive indices n + ik of the
        optical eigenmodes (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return np_readonly_view(self._mode_n_t)

    @property
    def absorption_coefficient(self):
        """numpy.ndarray : Absorption coefficients of the optical
        eigenmodes in cm^-1 (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return np_readonly_view(self._mode_a)

    @property
    def optical_conductivity(self):
        """numpy.ndarray : Optical conductivity of the optical
        eigenmodes in S cm^-1 (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return np_readonly_view(self._mode_s)

    @property
    def energy_loss_function(self):
        """numpy.ndarray : Energy loss functions of the optical
        eigenmodes (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return np_readonly_view(self._mode_l)

    @property
    def eigenvalue_unit_text_label(self):
        """str : Mode polarisability (eigenvalue) unit label suitable
        for plain-text output.
        """
        return _EPSILON_UNIT_TEXT_LABEL

    @property
    def eigenvalue_unit_plot_label(self):
        """str : Mode polarisability (eigenvalue) unit label suitable
        for plotting (contains TeX strings)."""
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

    def spectrum(self):
        """Return the optical eigenmode eigenvalues and optical spectra
        as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the dielectric function and derived
            quantities.
        """

        d = {"freq_energy": self._x}

        for i in range(self.num_dims):
            n = i + 1

            for k, prop in [
                ("epsilon_re", self._mode_evals.real),
                ("epsilon_im", self._mode_evals.imag),
                ("refractive_index", self._mode_n_t.real),
                ("extinction_coefficient", self._mode_n_t.imag),
                ("absorption_coefficient", self._mode_a),
                ("optical_conductivity_re", self._mode_s.real),
                ("optical_conductivity_im", self._mode_s.imag),
                ("energy_loss_function", self._mode_l),
            ]:
                d["{0}_{1}".format(k, n)] = prop[:, i]

        return pd.DataFrame(d)

    def unpolarised_eigenmode_average_optical_spectrum(
        self, t=1.0, n_f=1.0, n_b=1.0
    ):
        """Calculate an eigenmode-average optical spectrum for
        unpolarised incident and detected light.

        Parameters
        ----------
        t : float, optonal
            Thickness in mm (default: 1 mm).
        n_f, n_b : float, optional
            Refractive indices of the front (indicent) and back (exit)
            media (default: 1.0 = vacuum ~ air).

        Returns
        -------
        sp : EigenmodeAverageOpticalSpectrum
            Eigenmode-average optical spectrum.
        """

        return EigenmodeAverageOpticalSpectrum(self, t=t, n_f=n_f, n_b=n_b)

    def standard_polarised_eigenmode_average_optical_spectrum(
        self, i_pol=None, d_pol=None, t=1.0, n_f=1.0, n_b=1.0
    ):
        """Calculate an eigenmode-average optical spectrum for polarised
        incident and/or detected light in a standard "on axis"
        transmission geometry.

        Parameters
        ----------
        i_pol, d_pol : Polarisation or None, optional
            Polarisations of incident and detected light (default:
            `None` = unpolarised incident/detected light).
        t : float, optonal
            Thickness in mm (default: 1 mm).
        n_f, n_b : float, optional
            Refractive indices of the front (indicent) and back (exit)
            media (default: 1.0 = vacuum ~ air).

        Returns
        -------
        sp : EigenmodeAverageOpticalSpectrum
            Eigenmode-average optical spectrum.

        Notes
        -----
        This function is only valid when `num_dims` is 2 or 3. For the
        2D case, the third component of the polarisation vectors is
        discarded.

        This means that calling this function with the default `i_pol`
        and `s_pol`
        """

        if self.num_dims != 2 and self.num_dims != 3:
            raise RuntimeError(
                "Polarised eigenmode-average optical spectra can only "
                "be calculated from the eigenmodes of 2D or 3D "
                "dielectric functions."
            )

        # Default polarisations.

        if i_pol is None:
            i_pol = Polarisation.unpolarised_incident("z")

        if d_pol is None:
            d_pol = Polarisation.unpolarised_collected("z")

        ave_w = np.zeros((len(self._x), self.num_dims), dtype=np.float64)

        for v_i, v_d, w in i_pol.combine_with_iter(d_pol):
            # "Truncate" polarisation vectors for 2D case.

            v_i, v_d = v_i[: self.num_dims], v_d[: self.num_dims]

            w_i = np.einsum("nij,j->ni", self.dual_basis, v_i)
            w_d = np.einsum("j,nij->ni", v_d.conj(), self._mode_evecs)

            ave_w += w * np.abs(w_i) ** 2 * np.abs(w_d) ** 2

        return EigenmodeAverageOpticalSpectrum(
            self, t=t, n_f=n_f, n_b=n_b, ave_w=ave_w
        )

    @staticmethod
    def from_infrared_dielectric_function(eps_ir, **kwargs):
        """Create an `OpticalEigenmodeSpectrum` from an
        `InfraredDielectricFunction object.

        Parameters
        ----------
        eps_ir : InfraredDielectricFunction
            Infrared dielectric function.
        **kwargs : any
            Other optional keywords to the `OpticalEigenmodeSpectrum`
            constructor.

        Returns
        -------
        oe_sp : OpticalEigenmodeSpectrum
            Optical eigenmode spectrum.
        """

        return OpticalEigenmodes(
            eps_ir.x, eps_ir.epsilon, x_units=eps_ir.x_units, **kwargs
        )


# -------------------------------------
# EigenmodeAverageOpticalSpectrum class
# -------------------------------------


class EigenmodeAverageOpticalSpectrum(OpticalSpectrumBase):
    """Class for simulating unpolarised optical spectra from optical
    eigenmode spectra."""

    def __init__(self, oe_sp, t=1.0, n_f=1.0, n_b=1.0, ave_w=None):
        """Create a new instance of the
        `EigenmodeAverageOpticalSpectrum` class.

        Parameters
        ----------
        oe_sp : OpticalEigenmodeSpectrum
            Optical eigenmode spectrum.
        t : float, optonal
            Thickness in mm (default: 1 mm).
        n_f, n_b : float, optional
            Refractive indices of the front (indicent) and back (exit)
            media (default: 1.0 = vacuum ~ air).
        ave_w : array_like or None, optional
            Averaging weights for computing mode-average properties
            (shape: `(O, D)`, default: equal weights).
        """

        # This effectively "copies" the x-axis related quantities from
        # the initialising OpticalEigenmodeSpectrum.

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

        super(EigenmodeAverageOpticalSpectrum, self).__init__(
            x=oe_sp.x, x_units=oe_sp.x_units, t=t, n_f=n_f, n_b=n_b
        )

        self._oe_sp = oe_sp
        self._ave_w = ave_w

        self._abs_int = None

    def _init_optical_spectra(self):
        """Set the `a_int` and base class `_r_s`, `_r_t` and `_t`
        fields."""

        a_int, r_s, r_t, t = optical_spectra_from_optical_properties(
            self._oe_sp.refractive_index,
            self._oe_sp.absorption_coefficient,
            self._t,
            self._n_f,
            self._n_b,
        )

        self._abs_int = (a_int * self._ave_w).sum(axis=1)
        self._ref_s = (r_s * self._ave_w).sum(axis=1)
        self._ref_t = (r_t * self._ave_w).sum(axis=1)
        self._trans = (t * self._ave_w).sum(axis=1)

    def _init_single_reflectivity(self):
        """Implements the base class `_init_single_reflectivity()`
        abstract method."""

        self._init_optical_spectra()

    def _init_total_reflectivity_and_transmission(self):
        """Implements the base class
        `_init_total_reflectivity_and_transmission()` abstract method.
        """

        self._init_optical_spectra()

    @property
    def optical_eigenmodes(self):
        """OpticalEigenmodes : Optical eigenmodes used to generate the
        spectra."""
        return self._oe_sp

    @property
    def averaging_weights(self):
        """numpy.ndarray : Averaging weights for optical eigenmodes
        (shape: `(O, D)`)."""
        return np_readonly_view(self._ave_w)

    @property
    def intrinsic_absorbance(self):
        """numpy.ndarray : Intrinsic (Beer-Lambert) absorbance (shape:
        `(O,)`.)."""

        self._init_total_reflectivity_and_transmission()
        return np_readonly_view(self._abs_int)

    def spectrum(self):
        """Return the optical spectra as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the optical spectra.
        """

        # Insert additional intrinsic_absorbance column into DataFrame
        # generated by base class method.

        df = super(EigenmodeAverageOpticalSpectrum, self).spectrum()
        df.insert(0, "intrinsic_absorbance", self.intrinsic_absorbance)

        return df
