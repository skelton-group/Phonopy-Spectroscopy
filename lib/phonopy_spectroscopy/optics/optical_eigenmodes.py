# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Objects and routines for finding the optical eigenmodes of complex
dielectric functions and calculating the associated optical properties
and spectra."""

# -------
# Imports
# -------

import numpy as np
import pandas as pd

from .spectrum import OpticalSpectrumBase

from ..constants import ZERO_TOLERANCE, SPEED_OF_LIGHT, VACUUM_PERMITIVITY
from ..instrument import SPPolarisation
from ..units import convert_frequency_units

# --------------------------
# Optical properties/spectra
# --------------------------


def optical_properties_from_epsilon(e, eps, e_units):
    r"""Calculate the "standard" optical properties from an
    energy-dependent dielectric function: complex refractive index,
    absorption coefficient, optical conductivity and electron energy
    loss function.

    Parameters
    ----------
    e : array_like
        Frequencies (shape: `(O,)`).
    eps : array_like
        Scalar dielectric function or array_like of dielectric functions
        in units of relative permittivity (shape: `(O, ...)`).
    e_units : str
        Energy unit.

    Returns
    -------
    res : tuple of numpy.ndarray
        Tuple of `(n_t, a, s, l)` (shapes: `(O, ...)`; `a` in cm^-1 and
        `s` in S cm^-1).

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

    e = np.atleast_1d(np.asanyarray(e))

    if len(e) == 0 or e.ndim > 1:
        raise ValueError("e must be an array_like with shape (O,).")

    eps = np.atleast_1d(np.asanyarray(eps))

    if eps.shape[0] != len(e):
        raise ValueError("eps must be an array_like with shape (O, ...).")

    # Convert e to rad.  THz and reshape to broadcast against eps.

    e = convert_frequency_units(e, e_units, "rad_thz")
    e = e.reshape((len(e),) + (1,) * (eps.ndim - 1))

    # Complex refractive index.

    n_t = np.sqrt(eps)

    # Absorption coefficient.

    a = 1.0e-2 * (2.0 * 1.0e12 * e * n_t.imag) / SPEED_OF_LIGHT

    # Optical conductivity:

    s = 1.0e-2 * -1.0j * 1.0e12 * e * VACUUM_PERMITIVITY * (eps - 1.0)

    # Energy loss function.

    l = (-1.0 / eps).imag

    return (n_t, a, s, l)


def optical_spectra_from_optical_properties(n_t, a, t, n_f=1.0, n_b=1.0):
    r"""Calculate the "standard" optical properties from a complex
    refractive index, absorption coefficient, sample thickness and
    refractive indices of the front/back (entrance/exit)
    media: incoherent intrinsic (Beer-Lambert) transmission, "single"
    (infinite bulk) and total reflectivity at normal incidence, and
    (total) transmission.

    Parameters
    ----------
    n_t : array_like
        Complex refractive index (shape: `(O, ...)`).
    a : array_like
        Absorption coefficient in cm^-1 (shape: `(O, ...)`).
    t : float
        Sample thickness/path length (mm).
    n_f, n_b : float, optional
        Refractive indices of the front and back media (default: 1.0 =
        vacuum ~ air).

    Returns
    -------
    res : tuple of numpy.ndarray
        Tuple of `(t_int, r_s, r_t, t)` (shapes: `(O, ...)`).

    Notes
    -----
    The optical spectra are calculated using the formulae below:

    Intrinsic transmission:

    .. math::

        T_\mathrm{int} = \exp \left[- \alpha t \right]

    Relectivity at the front and back surfaces at normal incidence:

    .. math::

        R_\mathrm{f} = \left| \frac{n_f - \tilde{n}}{n_f + \tilde{n}} \right|^2

    .. math::

        R_\mathrm{b} = \left| \frac{\tilde{n} - n_b}{\tilde{n} + n_b} \right|^2

    Incoherent "single" reflectivity:

    .. math::

        R_\mathrm{s} = R_\mathrm{f}

    Incoherent total reflectivity:

    .. math::

        R_\mathrm{t} = R_\mathrm{f} + \frac{(1 - R_\mathrm{f})^2 R_\mathrm{b} T_\mathrm{int}^2}{1 - R_\mathrm{f} R_\mathrm{b} T_\mathrm{int}^2}

    Incoherent total transmission:

    .. math::

        T = \frac{(1 - R_\mathrm{f}) (1 - R_\mathrm{b}) T_\mathrm{int}}{1 - R_\mathrm{f} R_\mathrm{b} T_\mathrm{int}^2}
    """

    n_t = np.atleast_1d(np.asanyarray(n_t))
    a = np.atleast_1d(np.asanyarray(a))

    if a.shape != n_t.shape:
        raise ValueError("n_t and a must have the same shape.")

    # Reflectivity at the front and back surfaces at normal incidence.

    r_f = np.abs((n_f - n_t) / (n_f + n_t)) ** 2
    r_b = np.abs((n_t - n_b) / (n_t + n_b)) ** 2

    # Incoherent intrinsic transmission, total reflectivity and
    # transmission.

    t_i = np.exp(-1.0 * a * 1.0e-1 * t)

    r_t = r_f + (
        ((1.0 - r_f) ** 2 * r_b * t_i**2) / (1.0 - r_f * r_b * t_i**2)
    )

    t_t = ((1.0 - r_f) * (1.0 - r_b) * t_i) / (1.0 - r_f * r_b * t_i**2)

    return (t_i, r_f, r_t, t_t)


# -----------------------
# OpticalEigenmodes class
# -----------------------


class OpticalEigenmodes:
    """Class for finding the optical eigenmodes of a D-dimensional
    dielectric function."""

    absorption_coefficient_unit_text_label = r"\alpha / cm^-1"
    """str : Absorption coefficient unit label suitable for plain-text
    output."""

    absorption_coefficient_unit_plot_label = r"$\alpha$ / cm$^{-1}$"
    """str : Absorption coefficient unit label suitable for plotting
    (contains TeX strings)."""

    optical_conductivity_unit_text_label = r"\sigma / S cm^-1"
    """str : Optical conductivity unit label suitable for plain-text
    output."""

    optical_conductivity_unit_plot_label = r"$\sigma$ / S cm$^{-1}$"
    """str : Optical conductivity unit label suitable for plotting
    (contains TeX strings)."""

    def __init__(self, eps, branch_tracking=True):
        """Create a new instance of the `OpticalEigenmodes` class.

        Parameters
        ----------
        eps : DielectricFunction
            D-dimensional dielectric function.
        branch_tracking : bool, optional
            If `True`, reorder the optical eigenmodes by performing a
            "branch tracking" over frequencies (default: `True`).
        """

        # Diagonalise dielectric function and get optical eigenmodes.

        evals, evecs = eps.eigenmodes(branch_tracking=branch_tracking)

        evals.flags.writeable = False
        evecs.flags.writeable = False

        self._eps = eps

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
                self._eps.energy, self._mode_evals, self._eps.energy_units
            )

            n_t.flags.writeable = False
            a.flags.writeable = False
            s.flags.writeable = False
            l.flags.writeable = False

            self._mode_n_t = n_t
            self._mode_a = a
            self._mode_s = s
            self._mode_l = l

    def _lazy_init_dual_basis(self):
        """Lazy initialisation of dual basis."""

        if self._dual_basis is None:
            if self._eps.num_dims == 1:
                self._dual_basis = self._mode_evecs
            else:
                dual_basis = np.linalg.inv(self._mode_evecs)
                dual_basis.flags.writeable = False

                self._dual_basis = dual_basis

    @property
    def dielectric_function(self):
        """DielectricFunction : Underlying dielectric function."""
        return self._eps

    @property
    def num_dims(self):
        """int : Dimensionality (number of eigenmodes)."""
        return self._eps.num_dims

    @property
    def energy(self):
        """numpy.ndarray : Energies in `energy_units` (shape: `(O,)`)."""
        return self._eps.energy

    @property
    def energy_units(self):
        """str : Units of `energy`."""
        return self._eps.energy_units

    @property
    def eigenvalues(self):
        r"""numpy.ndarray : Optical eigenmode eigenvalues in \eps_0
        (shape: `(O, D)`)."""
        return self._mode_evals

    @property
    def eigenvectors(self):
        """numpy.ndarray : Eigenvectors of the optical eigenmodes
        (shape: `(O, D, D)`) (column-major format - the `j`th
        eigenvector at the `i`th frequency is
        `epsilon_eigenvectors[i, :, j]`)."""
        return self._mode_evecs

    @property
    def dual_basis(self):
        """numpy.ndarray : Inverse of the mode eigenvectors (dual basis)
        for eigenmode projections (shape: `(O, D, D)`)."""

        self._lazy_init_dual_basis()
        return self._dual_basis

    @property
    def refractive_index(self):
        """numpy.ndarray : Complex refractive indices n + ik of the
        optical eigenmodes (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return self._mode_n_t

    @property
    def absorption_coefficient(self):
        """numpy.ndarray : Absorption coefficients of the optical
        eigenmodes in cm^-1 (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return self._mode_a

    @property
    def optical_conductivity(self):
        """numpy.ndarray : Optical conductivity of the optical
        eigenmodes in S cm^-1 (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return self._mode_s

    @property
    def energy_loss_function(self):
        """numpy.ndarray : Energy loss functions of the optical
        eigenmodes (shape: `(O, D)`)."""

        self._lazy_init_optical_properties()
        return self._mode_l

    @property
    def energy_unit_text_label(self):
        return self._eps.energy_unit_text_label

    @property
    def energy_unit_plot_label(self):
        return self._eps.energy_unit_plot_label

    @property
    def eigenvalue_unit_text_label(self):
        """str : Mode polarisability (eigenvalue) unit label suitable
        for plain-text output.
        """
        return self._eps.epsilon_unit_text_label

    @property
    def eigenvalue_unit_plot_label(self):
        """str : Mode polarisability (eigenvalue) unit label suitable
        for plotting (contains TeX strings)."""
        return self._eps.epsilon_unit_plot_label

    def eigenmode_optical_properties(self):
        """Return the eigenmode eigenvalues and optical properties as
        a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the eigenvalues and optical
            properties.
        """

        self._lazy_init_optical_properties()

        d = {"freq_energy": self.energy}

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
                d[f"{k}_{n}"] = prop[:, i]

        return pd.DataFrame(d)


# ------------------------------------
# SingleCrystalOpticalEigenmodes class
# ------------------------------------


class SingleCrystalOpticalEigenmodes(OpticalEigenmodes):
    """Class for finding the optical eigenmodes of the accessible
    block of the dielectric function of an oriented single crystal and
    simulating on-axis optical spectra."""

    def __init__(self, eps, branch_tracking=True):
        """Create a new instance of the `SingleCrystalOpticalEigenmodes`
        class.

        Parameters
        ----------
        eps : DielectricFunction
            2D dielectric function.
        branch_tracking : bool, optional
            If `True`, reorder the optical eigenmodes by performing a
            "branch tracking" over frequencies (default: `True`).
        """

        if eps.num_dims != 2:
            raise ValueError(
                "SingleCrystalOpticalEigenmodes objects must be "
                "initialised with a 2D dielectric function - use "
                "OpticalEigenmodes instead."
            )

        super().__init__(eps, branch_tracking=branch_tracking)

    def incoherent_on_axis_optical_spectrum(
        self, i_pol=None, d_pol=None, t=1.0, n_f=1.0, n_b=1.0
    ):
        """Calculate an incoherent optical spectrum for a standard
        "on axis" geometry.

        Parameters
        ----------
        i_pol, d_pol : SPPolarisation or None, optional
            Polarisations of incident and detected light (default:
            unpolarised incident/detected light).
        t : float, optonal
            Sample thickness (path length) in mm (default: 1 mm).
        n_f, n_b : float, optional
            Refractive indices of the front (indicent) and back (exit)
            media (default: 1.0 = vacuum ~ air).

        Returns
        -------
        sp : EigenmodeAverageOpticalSpectrum
            Eigenmode-average optical spectrum.
        """

        # Default polarisations.

        if i_pol is None:
            i_pol = SPPolarisation.unpolarised_incident()

        if d_pol is None:
            d_pol = SPPolarisation.unpolarised_collected()

        # Convert to 3D polarisation vectors using the standard s = y
        # and p = x convention for on-axis geometries.

        i_pol = i_pol.to_polarisation(q="+z", n="+z")
        d_pol = d_pol.to_polarisation(q="+z", n="+z")

        ave_w = np.zeros((len(self.energy), 2), dtype=np.float64)

        for v_i, v_d, w in i_pol.combine_with_iter(d_pol):
            # "Truncate" polarisation vectors for 2D projection.

            v_i, v_d = v_i[:2], v_d[:2]

            w_i = np.einsum("nij,j->ni", self.dual_basis, v_i)
            w_d = np.einsum("j,nij->ni", v_d.conj(), self._mode_evecs)

            # Incoherent summation.

            ave_w += w * np.abs(w_i) ** 2 * np.abs(w_d) ** 2

        return EigenmodeAverageOpticalSpectrum(
            self, t=t, n_f=n_f, n_b=n_b, ave_w=ave_w
        )


# -------------------------------------
# EigenmodeAverageOpticalSpectrum class
# -------------------------------------


class EigenmodeAverageOpticalSpectrum(OpticalSpectrumBase):
    """Class for simulating incoherent optical spectra from optical
    eigenmodes."""

    def __init__(self, oes, t, n_f, n_b, ave_w=None):
        """Create a new instance of the
        `EigenmodeAverageOpticalSpectrum` class.

        Parameters
        ----------
        oe_sp : OpticalEigenmodes
            Optical eigenmodes of a D-dimensional dielectric function.
        t : float
            Sample thickness (path length) in mm (default: 1 mm).
        n_f, n_b : float, optional
            Refractive indices of the front (indicent) and back (exit)
            media (default: 1.0 = vacuum ~ air).
        ave_w : array_like or None, optional
            Averaging weights for computing mode-average properties
            (shape: `(O, D)`, default: equal weights).
        """

        super().__init__(
            oes.energy, t=t, n_f=n_f, n_b=n_b, e_units=oes.energy_units
        )

        if ave_w is not None:
            ave_w = np.atleast_1d(np.array(ave_w, dtype=np.float64, copy=True))

            if ave_w.ndim == 1:
                # Assume (O,) for optical eigenmodes of a 1D dielectric
                # function.

                ave_w = ave_w[:, np.newaxis]

            if ave_w.shape != (len(oes.energy), oes.num_dims):
                raise ValueError(
                    "ave_w must be an array_like with shpe (O, D)."
                )

            if not np.allclose(ave_w.sum(axis=1), 1.0, atol=ZERO_TOLERANCE):
                raise ValueError("ave_w at each frequency must sum to 1.")
        else:
            # Equal weight.

            ave_w = (
                np.ones((len(oes.energy), oes.num_dims), dtype=np.float64)
            ) / oes.num_dims

        ave_w.flags.writeable = False

        self._oes = oes
        self._ave_w = ave_w

        self._trans_int = None
        self._aps_int = None
        self._abs_int = None

    def _init_optical_spectra(self):
        """Set the `_trans_int`, `abs_int` and `aps_int` fields, and the
        base class `_r_s`, `_r_t` and `_trans` fields."""

        if self._trans_int is None:
            t_i, r_s, r_t, t_t = optical_spectra_from_optical_properties(
                self._oes.refractive_index,
                self._oes.absorption_coefficient,
                self._t,
                self._n_f,
                self._n_b,
            )

            self._trans_int = (t_i * self._ave_w).sum(axis=1)
            self._ref_s = (r_s * self._ave_w).sum(axis=1)
            self._ref_t = (r_t * self._ave_w).sum(axis=1)
            self._trans = (t_t * self._ave_w).sum(axis=1)

    def _init_single_reflectivity(self):
        """Implements the base class `_init_single_reflectivity()`
        abstract method."""

        self._init_optical_spectra()

    def _init_total_reflectivity_and_transmission(self):
        """Implements the base class
        `_init_total_reflectivity_and_transmission()` abstract method.
        """

        self._init_optical_spectra()

        if self._aps_int is None:
            aps_int = 1.0 - self._trans_int
            abs_int = -1.0 * np.log10(self._trans_int)

            aps_int.flags.writeable = False
            abs_int.flags.writeable = False

            self._aps_int = aps_int
            self._abs_int = abs_int

    @property
    def optical_eigenmodes(self):
        """OpticalEigenmodes : Optical eigenmodes used to generate the
        spectra."""
        return self._oes

    @property
    def averaging_weights(self):
        """numpy.ndarray : Averaging weights for optical eigenmodes
        (shape: `(O, D)`)."""
        return self._ave_w

    @property
    def intrinsic_transmission(self):
        """numpy.ndarray : Intrinsic (Beer-Lambert) transmission (shape:
        `(O,)`.)."""

        self._init_total_reflectivity_and_transmission()
        return self._trans_int

    @property
    def intrinsic_absorptance(self):
        """numpy.ndarray : Intrinsic (Beer-Lambert) absorptance (shape:
        `(O,)`)."""

        self._init_total_reflectivity_and_transmission()
        return self._aps_int

    @property
    def intrinsic_absorbance(self):
        """numpy.ndarray : Intrinsic (Beer-Lambert) absorbance (shape:
        `(O,)`."""

        self._init_total_reflectivity_and_transmission()
        return self._abs_int

    def spectrum(self):
        """Return the optical spectra as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the optical spectra.
        """

        # Insert additional intrinsic_transmission column into DataFrame
        # generated by base class method.

        df = super().spectrum()

        df.insert(1, "intrinsic_transmission", self.intrinsic_transmission)
        df.insert(2, "intrinsic_absorptance", self.intrinsic_absorptance)
        df.insert(3, "intrinsic_absorbance", self.intrinsic_absorbance)

        return df
