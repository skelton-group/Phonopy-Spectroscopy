# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Classes and routines for simulating Raman spectra."""


# -------
# Imports
# -------


import warnings

import pandas as pd
import numpy as np

from ..constants import (
    ZERO_TOLERANCE,
    AMU_TO_KG,
    PLANCK_CONSTANT_J,
    SPEED_OF_LIGHT,
)

from ..distributions import lorentzian, phonon_occupation_number
from ..spectrum_base import GammaPhononSpectrumBase
from ..units import nm_to_thz

from ..utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
    np_expand_dims,
)


# ---------
# Constants
# ---------


_RAMAN_INTENSITY_PREFACTOR = (
    (1.0e10**2)
    * (2.0 * (np.pi**2))
    * (PLANCK_CONSTANT_J / (SPEED_OF_LIGHT**4))
    * ((1.0e-10**4) / AMU_TO_KG)
)


# ---------
# Functions
# ---------


def adjust_gamma_phonons_for_spectrum_type(
    freqs, ints, lws, irrep_syms=None, spectrum_type="stokes"
):
    """Adjust a set of band frequencies, intensities, linewidths and,
    optionally, irreps for different types of Raman spectrum.

    Parameters
    ----------
    freqs : array_like
        Frequencies (shape: `(N,)`).
    ints : array_like
        Intensities for `N` bands (shape: `(N,)` or `N` modes and `M`
        calculations (shape: `(N, M)`).
    lws : array_like or None, optional
        Linewidths (shape: `(N,)`).
    irrep_syms : array_like or None, optional
        Irreducible representation (irrep) symbols (shape: `(N,)`).
    spectrum_type : {"stokes", "anti-stokes", "both"}
        Type of Raman spectrum (default "stokes").

    Returns
    -------
    phonons : tuple of (numpy.ndarray or None)
        `(freqs, ints, lws, irreps)` tuple with input data adjusted for
        `spectrum_type`.
    """

    freqs = np.asarray(freqs)

    if not np_check_shape(freqs, (None,)):
        raise ValueError("freqs must be an array_like with shape (N,).")

    ints, n_dim_add = np_expand_dims(np.asarray(ints), (len(freqs), None))

    if lws is not None:
        lws = np.asarray(lws)

        if not np_check_shape(lws, (len(freqs),)):
            raise ValueError("lws must be an array_like with shape (N,).")

    if irrep_syms is not None:
        irrep_syms = np.asarray(irrep_syms)

        if not np_check_shape(irrep_syms, (len(freqs),)):
            raise ValueError(
                "irrep_syms must be an array_like with shape (N,)."
            )

    spectrum_type = str(spectrum_type).lower()

    if spectrum_type == "stokes":
        # Nothing to do - return input data.
        return (
            freqs,
            ints if n_dim_add == 0 else ints.reshape((-1,)),
            lws,
            irrep_syms,
        )

    if spectrum_type == "anti-stokes" or spectrum_type == "both":
        # Reverse data for anti-stokes spectrum.

        freqs_as = (-1.0 * freqs)[::-1]
        ints_as = ints[::-1, :]
        lws_as = lws[::-1] if lws is not None else None
        irrep_syms_as = irrep_syms[::-1] if irrep_syms is not None else None

        if spectrum_type == "anti-stokes":
            return (
                freqs_as,
                ints_as if n_dim_add == 0 else ints_as.reshape((-1,)),
                lws_as,
                irrep_syms_as,
            )

        else:
            ints_both = np.concatenate((ints_as, ints), axis=0)

            irrep_syms_both = None

            if irrep_syms is not None:
                irrep_syms_both = np.concatenate(
                    (irrep_syms_as, irrep_syms), axis=0
                )

            return (
                np.concatenate((freqs_as, freqs), axis=0),
                ints_both if n_dim_add == 0 else ints_both.reshape((-1,)),
                np.concatenate((lws_as, lws), axis=0),
                irrep_syms_both,
            )

    raise ValueError('Unsupported spectrum_type="{0}".'.format(spectrum_type))


def modulate_intensities(freqs, ints, w=None, t=None):
    """Modulate Raman band intensities for a measurement wavelength
    and/or temperature.

    Parameters
    ----------
    freqs : array_like
        Frequencies (shape: `(N,)`).
    ints : array_like
        Intensities for `N` bands (shape: `(N,)` or `N` modes and `M`
        calculations (shape: `(N, M)`).
    w, t : float, optional
        Measurement wavelength in nm and temperature in K.

    Returns
    -------
    mod_ints : numpy.ndarray
        Modulated intensities (same shape as `ints`).
    """

    freqs = np.asarray(freqs)

    if not np_check_shape(freqs, (None,)):
        raise ValueError("freqs must be an array_like with shape (N,).")

    ints, n_dim_add = np_expand_dims(np.asarray(ints), (len(freqs), None))

    if w is not None:
        if w <= 0.0:
            raise ValueError("If supplied, w must be non-zero and positive.")

        f_i = nm_to_thz(w)

        scale = _RAMAN_INTENSITY_PREFACTOR * (
            ((1.0e12 * (f_i - freqs)) ** 4) / (1.0e12 * np.abs(freqs))
        )

        ints *= scale[:, np.newaxis]

    if t is not None:
        if t <= 0.0:
            raise ValueError("If supplied, t must be non-zero and positive.")

        occ_nums = phonon_occupation_number(np.abs(freqs), t)

        mask = freqs >= 0.0
        inv_mask = np.logical_not(mask)

        # The weighting factors are different for Stokes and
        # anti-Stokes branches.

        if mask.any():
            ints[mask, :] *= occ_nums[mask, np.newaxis] + 1.0

        if inv_mask.any():
            ints[inv_mask, :] *= occ_nums[inv_mask, np.newaxis]

    return ints if n_dim_add == 0 else ints.reshape((-1,))


# -----------------------
# RamanSpectrumBase class
# -----------------------


class RamanSpectrumBase(GammaPhononSpectrumBase):
    """Base class for generating simulated Raman spectra from sets of
    frequencies, band intensities and linewidths."""

    def __init__(
        self,
        freqs,
        ints,
        lws,
        irreps=None,
        w=None,
        t=None,
        spectrum_type="stokes",
        **kwargs
    ):
        r"""Create a new instance of the `RamanSpectrum` class.

        Parameters
        ----------
        freqs : array_like
            Frequencies in THz (shape: `(N,)`).
        ints : array_like
            Band intensities in Ang^4 / sqrt(amu) (shape: `(N,)` for 1D
            spectra, or `(N, M)` for 2D spectra).
        lws : array_like
            Linewidths in THz (shape: `(N,)`).
        irreps : Irreps or None, optional
            `Irreps` object assigning bands to irrep groups.
        w : float or None, optional
            Measurement wavelength in nm to calculate intensity
            modulation envelope (default: None).
        t : float or None, optional
            Temperature in K to calculate phonon occupation numbers for
            intensity modulation envelope (default: None)
        spectrum_type : {"stokes", "anti-stokes", "both"}
            Type of spectrum (default: "stokes").
        **kwargs : any
            Keyword arguments to the `GammaPhononSpectrumBase`
            constructor.

        See Also
        --------
        spectrum.GammaPhononSpectrumBase
            Base class for `RamanSpectrum`.

        Notes
        -----
        During initialisation, if irreps are supplied the
        frequencies/linewidths and intensities of degenerate bands are
        averaged and summed, respectively, so the outer dimensions of
        properties (`N`) may differ from those supplied to the
        constructor.
        """

        freqs = np_asarray_copy(freqs, dtype=np.float64)

        if not np_check_shape(freqs, (None,)):
            raise ValueError("freqs must be an array_like with shape `(N,)`.")

        ints, n_dim_add = np_expand_dims(
            np_asarray_copy(ints, dtype=np.float64), (len(freqs), None)
        )

        lws = np_asarray_copy(lws, dtype=np.float64)

        if not np_check_shape(lws, (len(freqs),)):
            raise ValueError("lws must be an array_like with shape `(N,)`.")

        # If we have any imaginary modes with non-zero intensity, it
        # not cause issues in the calculations but it may give
        # unphysical results in some cases.

        mask = np.logical_and(freqs < 0.0, np.abs(freqs) > ZERO_TOLERANCE)

        if (np.abs(ints[mask, :]) > ZERO_TOLERANCE).any():
            warnings.warn(
                "Imaginary modes with non-zero intensity may lead to "
                "unphysical results, particularly if applying "
                "temperature modulation to the intensities and/or "
                "simulating anti-Stokes spectra.",
                UserWarning,
            )

        # If irreps are supplied, average frequencies/linewidths and sum
        # intensities.

        irrep_syms = None

        if irreps is not None:
            ir_band_inds = irreps.irrep_band_indices

            freqs = np.array(
                [np.mean(freqs[inds]) for inds in ir_band_inds],
                dtype=np.float64,
            )

            ints = np.array(
                [np.sum(ints[inds, :], axis=0) for inds in ir_band_inds],
                dtype=np.float64,
            )

            lws = np.array(
                [np.mean(lws[inds]) for inds in ir_band_inds], dtype=np.float64
            )

            irrep_syms = np_asarray_copy(irreps.irrep_symbols, dtype=object)

        # Adjust freqs, ints, lws and irrep_syms for spectrum type.

        freqs, ints, lws, irrep_syms = adjust_gamma_phonons_for_spectrum_type(
            freqs,
            ints,
            lws,
            irrep_syms=irrep_syms,
            spectrum_type=spectrum_type,
        )

        # Apply intensity modulation.

        ints = modulate_intensities(freqs, ints, w, t)

        # Call the GammaPhononSpectrumBase constructor to handle
        # "x-axis"-related intialisation.

        super(RamanSpectrumBase, self).__init__(
            freqs, lws, irrep_syms=irrep_syms, **kwargs
        )

        # Store additional parameters.

        self._ints = ints

        self._w = w
        self._t = t

        self._spectrum_type = spectrum_type

        # Set _sp to default value.

        self._sp = None

    def _lazy_init_spectrum(self):
        """Simulate spectrum on first call to `y`, `z` or `spectrum()`."""

        if self._sp is None:
            x = self.x

            _, num_sp = self._ints.shape

            sp = np.zeros((len(x), num_sp), dtype=np.float64)

            for sp_idx in range(self._ints.shape[1]):
                for f, i, lw in zip(
                    self.raman_shifts, self._ints[:, sp_idx], self.linewidths
                ):
                    sp[:, sp_idx] += lorentzian(x, i, f, lw)

            self._sp = sp

    @property
    def frequencies(self):
        """numpy.ndarray : Phonon frequencies (shape: `(N,)`)."""
        # If the spectrum includes anti-Stokes branches, the frequencies
        # are actually Raman shifts and may include negative values. To
        # avoid confusion, we hide the base class property and intercept
        # it with a warning.

        if (
            self._spectrum_type == "anti-stokes"
            or self._spectrum_type == "both"
        ):
            warnings.warn(
                "The values returned by frequencies are Raman shifts "
                "and not phonon frequencies (to avoid this warning, "
                "use the raman_shifts property alias instead).",
                RuntimeWarning,
            )

        return super(RamanSpectrumBase, self).frequencies

    @property
    def raman_shifts(self):
        """numpy.ndarray : Raman shifts (shape: `(N,)`)."""
        return super(RamanSpectrumBase, self).frequencies

    @property
    def spectrum_type(self):
        """{"stokes", "anti-stokes", "both"} : Type of spectrum."""
        return self._spectrum_type

    @property
    def measurement_wavelength(self):
        """float or None : Measurement wavelength."""
        return self._w

    @property
    def measurement_temperature(self):
        """float or None : Measurement temperature."""
        return self._t

    @property
    def _intensity_unit_text_label(self):
        """str : Intensity unit label suitable for plain-text output."""

        # If a measurement wavelength was set, the intensities are
        # converted to differential cross sections. If not, the
        # intensities are "raw" values.

        return (
            "d(sigma)/d(omega) / (Ang^2 sr^-1)"
            if self._w is not None
            else "I^Raman / (Ang^4 amu^-1)"
        )

    @property
    def _intensity_unit_plot_label(self):
        """str : Intensity unit label suitable for plotting (may contain
        TeX strings)."""

        # (See comment on _intensity_unit_text_label.)

        return (
            r"$d \sigma / d \Omega$ / ($\mathrm{\AA}^2$ sr$^{-1}$)"
            if self._w is not None
            else r"$I^\mathrm{Raman}$ / ($\mathrm{AA}^4$ amu$^{-1}$)"
        )

    def _get_data_frame_column_headers(self):
        """Return a set of column headers for the `pandas.DataFrame`
        objects created by `peak_table()` and `spectrum()`.

        Returns
        -------
        col_hdrs : list of str
            Column headers.
        """

        hdr_base = "cross_sect" if self._w is not None else "int"

        if self._ints.shape[1] == 1:
            return [hdr_base]
        else:
            return [
                "{0}_{1}".format(hdr_base, i + 1)
                for i in range(self._ints.shape[1])
            ]

    def peak_table(self):
        """Return the peak table as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the peak table.
        """

        d = {
            "freq_energy": self._freqs,
            "linewidth": self._lws,
            "irrep": self._irrep_syms,
        }

        for i, h in enumerate(self._get_data_frame_column_headers()):
            d[h] = self._ints[:, i]

        return pd.DataFrame(d)

    def spectrum(self):
        """Return the simulated spectrum as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the spectrum.
        """

        self._lazy_init_spectrum()

        d = {"freq_energy": self.x}

        for i, h in enumerate(self._get_data_frame_column_headers()):
            d[h] = self._sp[:, i]

        return pd.DataFrame(d)


# ---------------------
# RamanSpectrum1D class
# ---------------------


class RamanSpectrum1D(RamanSpectrumBase):
    """Class for generating simulated 1D Raman spectra from sets of
    frequencies, band intensities and linewidths."""

    def __init__(
        self,
        freqs,
        ints,
        lws,
        irreps=None,
        w=None,
        t=None,
        spectrum_type="stokes",
        x_range=None,
        x_res=None,
        x_units="thz",
        **kwargs
    ):
        """Create a new instance of the `RamanSpectrum1D` class.

        Parameters
        ----------
        freqs : array_like
            Frequencies in THz (shape: `(N,)`).
        ints : array_like
            Band intensities in Ang^4 / sqrt(amu) (shape: `(N,)`).
        lws : array_like
            Linewidths in THz (shape: `(N,)`).
        irreps : Irreps or None, optional
            `Irreps` object assigning bands to irrep groups.
        w : float or None, optional
            Measurement wavelength in nm to calculate intensity
            modulation envelope (default: None).
        t : float or None, optional
            Temperature in K to calculate phonon occupation numbers for
            intensity modulation envelope (default: None)
        spectrum_type : {"stokes", "anti-stokes", "both"}
            Type of spectrum (default: "stokes").
        x_range : tuple of float or None, optional
            `(min, max)` range of x-axis of simulated spectrum in
            `x_units` (default: automatically determined).
        x_res : float or None, optional
            Resolution of x-axis in `x_units` (default: automatically
            determined).
        x_units : str or None, optional
            x-axis units of simulated spectrum (default: `"thz"`).
        **kwargs : any
            Keyword arguments to parent class constructor(s).
        """

        # Perform base class initialisation.

        super(RamanSpectrum1D, self).__init__(
            freqs,
            ints,
            lws,
            irreps=irreps,
            w=w,
            t=t,
            spectrum_type=spectrum_type,
            x_range=x_range,
            x_res=x_res,
            x_units=x_units,
            **kwargs
        )

        # Check supplied intensities are 1D.

        if not np_check_shape(self._ints, (len(self._freqs), 1)):
            raise ValueError(
                "For 1D spectra, ints must be a an array_like with "
                "shape (N,)."
            )

    @property
    def intensities(self):
        """numpy.ndarray : Band intensities (shape: `(N,)`)."""
        return np_readonly_view(self._ints.reshape((-1,)))

    @property
    def y(self):
        """numpy.ndarray : Spectrum (shape: `(O,)`.)"""
        self._lazy_init_spectrum()
        return np_readonly_view(self._sp.reshape((-1,)))

    @property
    def y_unit_text_label(self):
        """str : y-axis label suitable for plain-text output."""
        return self._intensity_unit_text_label

    @property
    def y_unit_plot_label(self):
        """str: y-axis label suitable for plotting (may contain TeX
        strings)."""
        return self._intensity_unit_plot_label


# ---------------------
# RamanSpectrum2D class
# ---------------------


class RamanSpectrum2D(RamanSpectrumBase):
    """Class for generating simulated 2D Raman spectra from sets of
    frequencies, band intensities and linewidths."""

    def __init__(
        self,
        freqs,
        ints,
        lws,
        d2_axis_vals,
        d2_unit_text_label,
        irreps=None,
        w=None,
        t=None,
        spectrum_type="stokes",
        x_range=None,
        x_res=None,
        x_units="thz",
        d2_unit_plot_label=None,
        d2_col_hdrs=None,
        **kwargs
    ):
        """Create a new instance of the `RamanSpectrum2D` class.

        Parameters
        ----------
        freqs : array_like
            Frequencies in THz (shape: `(N,)`).
        ints : array_like
            Band intensities in Ang^4 / sqrt(amu) (shape: `(N,)`).
        lws : array_like
            Linewidths in THz (shape: `(N,)`).
        d2_axis_vals : array_like
            Values of secondary axis (shape: `(M,)`).
        d2_unit_text_label : str
            Text label for secondary axis units.
        irreps : Irreps or None, optional
            `Irreps` object assigning bands to irrep groups.
        w : float or None, optional
            Measurement wavelength in nm to calculate intensity
            modulation envelope (default: None).
        t : float or None, optional
            Temperature in K to calculate phonon occupation numbers for
            intensity modulation envelope (default: None)
        spectrum_type : {"stokes", "anti-stokes", "both"}
            Type of spectrum (default: "stokes").
        x_range : tuple of float or None, optional
            `(min, max)` range of x-axis of simulated spectrum in
            `x_units` (default: automatically determined).
        x_res : float or None, optional
            Resolution of x-axis in `x_units` (default: automatically
            determined).
        x_units : str or None, optional
            x-axis units of simulated spectrum (default: `"thz"`).
        d2_unit_plot_label : str or None, optional
            Plot label for secondary axis units (default:
            `d2_unit_text_label`).
        d2_col_hdrs : array_like or None, optional
            Column headings for `pandas.DataFrame` objects (shape:
            `(M,)`; default: automatically determined).
        **kwargs : any
            Keyword arguments to parent class constructor(s).
        """

        # Pass args/kwargs through to base class.

        super(RamanSpectrum2D, self).__init__(
            freqs,
            ints,
            lws,
            irreps=irreps,
            w=w,
            t=t,
            spectrum_type=spectrum_type,
            x_range=x_range,
            x_res=x_res,
            x_units=x_units,
            **kwargs
        )

        d2_axis_vals = np_asarray_copy(d2_axis_vals)

        if not np_check_shape(d2_axis_vals, (self._ints.shape[1],)):
            raise ValueError(
                "d2_axis_vals must be an array_like with shape `(M,)`."
            )

        d2_unit_text_label = str(d2_unit_text_label)

        if d2_unit_plot_label is not None:
            d2_unit_plot_label = str(d2_unit_plot_label)
        else:
            d2_unit_plot_label = d2_unit_text_label

        if d2_col_hdrs is not None:
            d2_col_hdrs = np.array([str(h) for h in d2_col_hdrs], dtype=object)

            if not np_check_shape(d2_col_hdrs, (self._ints.shape[1],)):
                raise ValueError(
                    "If supplied, d2_col_hdrs must be an array_like "
                    "with shape `(M,)`."
                )

        self._d2_axis_vals = d2_axis_vals

        self._d2_unit_text_label = d2_unit_text_label
        self._d2_unit_plot_label = d2_unit_plot_label

        self._d2_col_hdrs = d2_col_hdrs

    def _get_data_frame_column_headers(self):
        """Overrides base class method to return column headers for
        `pandas.DataFrame` objects supplied by during construction.

        Returns
        ------
        col_hdrs : list of str
            Column headers.
        """

        if self._d2_col_hdrs is not None:
            return self._d2_col_hdrs

        return super()._get_data_frame_column_headers()

    @property
    def intensities(self):
        """numpy.ndarray : Band intensities (shape: `(N, M)`)."""
        return np_readonly_view(self._ints)

    @property
    def y(self):
        """numpy.ndarray : y values (shape: `(M,)`)."""
        return np_readonly_view(self._d2_axis_vals)

    @property
    def z(self):
        """numpy.ndarray : z values (shape: `(O, M)`)."""

        self._lazy_init_spectrum()
        return np_readonly_view(self._sp)

    @property
    def y_unit_text_label(self):
        """str : y-axis unit label suitable for plain-text output."""
        return self._d2_unit_text_label

    @property
    def y_unit_plot_label(self):
        """str: y-axis unit label suitable for plotting (may contain TeX
        strings.)"""
        return self._d2_unit_plot_label

    @property
    def z_unit_text_label(self):
        """str : z-axis unit label suitable for plain-text output."""
        return self._intensity_unit_text_label

    @property
    def z_unit_plot_label(self):
        """str : z-axis unit label suitable for plotting (may contain
        TeX strings)."""
        return self._intensity_unit_plot_label
