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

from itertools import product

from ..constants import (
    AMU_TO_KG,
    PLANCK_CONSTANT_J,
    SPEED_OF_LIGHT,
    ZERO_TOLERANCE,
)

from ..distributions import lorentzian, phonon_occupation_number

from ..irreps import Irreps

from ..units import (
    convert_frequency_energy_units,
    nm_to_thz,
    get_frequency_energy_unit_text_label,
    get_frequency_energy_unit_plot_label,
)

from ..utility.spectrum import auto_range

# ---------
# Constants
# ---------

_RAMAN_INTENSITY_PREFACTOR = (
    (1.0e10**2)
    * (2.0 * (np.pi**2))
    * (PLANCK_CONSTANT_J / (SPEED_OF_LIGHT**4))
    * ((1.0e-10**4) / AMU_TO_KG)
)
"""float : Prefactor for converting Raman intensities in Ang^4 / amu to 
cross sections in Ang^2 / sr."""

# ---------
# Functions
# ---------


def adjust_band_params_for_spectrum_type(
    shifts, ints, lws, ir_syms=None, sp_type="stokes"
):
    """Adjust band frequencies, intensities, linewidths and, optionally,
    irreps for different types of Raman spectra.

    Parameters
    ----------
    freqs : scalar or array_like
        Frequencies (scalar or shape: `(..., N)`).
    ints : scalar or array_like
        Intensities (scalar or shape: `(..., N)`).
    lws : scalar, array_like or None, optional
        Linewidths (scalar or shape: `(..., N)`).
    ir_syms : scalar, array_like or None, optional
        Irreducible representation (irrep) symbols (scalar or shape:
        `(..., N)`).
    sp_type : {"stokes", "anti-stokes", "both"}
        Type of Raman spectrum (default "stokes").

    Returns
    -------
    params : tuple of (numpy.ndarray or None)
        `(shifts, ints, lws, irreps)` tuple with input data adjusted for
        `sp_type` and broadcast to a common shape.

    Notes
    -----
    The array shapes of `freqs` and `lws`/`ir_syms` if provided are
    matched to the shape of `ints` according to standard NumPy
    brodcasting rules.
    """

    ints = np.asanyarray(ints)
    shifts = np.asanyarray(shifts)

    if lws is not None:
        lws = np.asanyarray(lws)

    if ir_syms is not None:
        ir_syms = np.asanyarray(ir_syms)

    sp_type = str(sp_type).lower()

    match str(sp_type).lower():
        case "stokes":
            pass

        case "anti-stokes":
            # Negate frequencies and reverse band ordering.

            shifts = -1.0 * shifts[..., ::-1]
            ints = ints[..., ::-1]

            if lws is not None:
                lws = lws[..., ::-1]

            if ir_syms is not None:
                ir_syms = ir_syms[..., ::-1]

        case "both":
            # Concatenate stokes/anti-Stokes data. For this, we need at
            # least 1D arrays.

            shifts, ints = np.atleast_1d(shifts, ints)

            if (shifts < 0.0).any():
                warnings.warn(
                    "Imaginary frequencies may give unexpected results for "
                    'sp_type="both".',
                    RuntimeWarning,
                )

            shifts = np.concatenate(
                (-1.0 * shifts[..., ::-1], shifts), axis=-1
            )
            ints = np.concatenate((ints[..., ::-1], ints), axis=-1)

            if lws is not None:
                lws = np.atleast_1d(lws)
                lws = np.concatenate((lws[..., ::-1], lws), axis=-1)

            if ir_syms is not None:
                ir_syms = np.atleast_1d(ir_syms)

                ir_syms = np.concatenate(
                    (ir_syms[..., ::-1], ir_syms), axis=-1
                )

        case _:
            raise ValueError(f'Unknown sp_type="{sp_type}".')

    # Broadcast shapes against ints.

    ints_shape = ints.shape

    return (
        np.broadcast_to(shifts, ints_shape),
        ints,
        np.broadcast_to(lws, ints_shape) if lws is not None else None,
        np.broadcast_to(ir_syms, ints_shape) if ir_syms is not None else None,
    )


def intensities_to_cross_sections(s, i, w, t, s_units="thz"):
    """Convert "raw" Raman intensities in Ang^4 / amu to differential
    cross sections in Ang^2 / sr.

    Parameters
    ----------
    s : scalar or array_like
        Raman shifts in `s_units` (scalar or shape: `(..., O)`).
    i : scalar or array_like
        Intensities in Ang^4 / amu (scalar or shape: `(..., O)`).
    w, t : float
        Measurement wavelength in nm and temperature in K.
    s_units : str, optional
        Units of `s` (default: "thz").

    Returns
    -------
    xs : numpy.ndarray
        Cross sections in Ang^2 / sr (same shape as `i`).

    Notes
    -----
    `s` must be broadcastable to the same shape as `i`.
    """

    if w <= 0.0 or t <= 0.0:
        raise ValueError("w and t must be > 0.")

    i = np.asanyarray(i)

    # Convert x to THz and broadcast against y.

    s = np.broadcast_to(
        convert_frequency_energy_units(np.asanyarray(s), s_units, "thz"),
        i.shape,
    )

    # Power modulation.

    f_i = nm_to_thz(w)

    scale = _RAMAN_INTENSITY_PREFACTOR * (
        ((1.0e12 * (f_i - s)) ** 4) / (1.0e12 * np.abs(s))
    )

    # Temperature modulation.

    n = phonon_occupation_number(np.abs(s), t)

    t_w = np.where(s >= 0.0, n + 1.0, n)

    return i * scale * t_w


# ---------------------
# NDRamanSpectrum class
# ---------------------


class NDRamanSpectrum:
    """Class for simulating N-dimensional Raman spectra with optional
    averaging over crystallite orientations."""

    intensity_unit_text_label = "I / (Ang^4 amu^-1)"
    """str : Intensity unit label suitable for plain-text output."""

    cross_section_unit_text_label = "d(sigma)/d(omega) / (Ang^2 sr^-1)"
    """str : Cross section unit label suitable for plain-text output."""

    cross_section_unit_plot_label = (
        r"$d \sigma / d \Omega$ / ($\mathrm{\AA}^2$ sr$^{-1}$)"
    )
    """str : Cross section unit label suitable for plotting (contains
    TeX strings)."""

    def __init__(
        self,
        freqs,
        ints,
        lws,
        w,
        t,
        irreps=None,
        ir_syms="?",
        sp_type="stokes",
        sp_units="thz",
        sp_range=None,
        sp_res=None,
        or_a=None,
        or_w=None,
    ):
        r"""Create a new instance of the `NDRamanSpectrum` class.

        Parameters
        ----------
        freqs : scalar or array_like
            Frequencies (shape: `(..., N)` or `(..., L, N)`).
        ints : scalar or array_like
            Intensities (shape: `(..., N)` or `(..., L, N)`).
        lws : scalar, array_like or None, optional
            Linewidths (scalar or shape: `(..., N)` or `(..., L, N)`).
        w, t : float
            Measurement wavelength (nm) and temperature (K).
        irreps : Irreps or None, optional
            Irreducible representations (irreps) of the bands included
            in the spectrum (assumed to apply to all dimensions for
            > 1D spectra, overrides `ir_syms`).
        ir_syms : scalar or array_like, optional
            Irrep symbols of the bands in the spectrum (scalar or shape:
            `(..., N)`, default: "?").
        sp_type : {"stokes", "anti-stokes", "both"}, optional
            Spectrum type (default: "stokes").
        sp_units : str, optional
            Raman shift units of simulated spectrum (default: "thz").
        sp_range : tuple of float, optional
            `(min, max)` for simulated spectrum (default: automatic).
        sp_res : float, optional
            Resolution of simulated spectrum (default: automatic).
        or_a : array_like or None, optional
            Specifies the three Euler angles \phi, \theta and \psi for
            each orientation (L) (shape: `(L, 3)`).
        or_w : array_like or None, optional
            Optional averaging weights for orientations (L) (shape:
            `(L,)`, defaults to 1 / L).

        Notes
        -----
        If `irreps` are supplied they are used to average the
        frequencies/linewidths and sum the intensities of degenerate
        irreps.
        """

        # Enforce minimum ndim and dtype, and copy data to ensure
        # "ownership".

        freqs = np.atleast_1d(np.array(freqs, dtype=np.float64, copy=True))
        ints = np.atleast_1d(np.array(ints, dtype=np.float64, copy=True))
        lws = np.atleast_1d(np.array(lws, dtype=np.float64, copy=True))

        # Handle irrep summing/averaging if needed.

        if isinstance(irreps, Irreps):
            # "Reduce" N -> M bands by averaging
            # frequencies/linewidths and summing intensities.

            inds = irreps.irrep_band_indices
            ir_syms = np.copy(irreps.irrep_symbols)

            n = len(inds)
            m = len(ir_syms)

            unique_inds = np.unique(inds)

            # Broadcast the outer dimensions of freqs, ints and lws
            # against ir_band_inds.

            freqs = np.broadcast_to(freqs, freqs.shape[:-1] + (n,))
            ints = np.broadcast_to(ints, ints.shape[:-1] + (n,))
            lws = np.broadcast_to(lws, lws.shape[:-1] + (n,))

            # Average freqs/linewidths and sum intensities.

            freqs_new = np.zeros(freqs.shape[:-1] + (m,), dtype=np.float64)
            ints_new = np.zeros(ints.shape[:-1] + (m,), dtype=np.float64)
            lws_new = np.zeros(lws.shape[:-1] + (m,), dtype=np.float64)

            for idx in unique_inds:
                freqs_new[..., idx] = np.mean(freqs[..., inds == idx], axis=-1)
                ints_new[..., idx] = ints[..., inds == idx].sum(axis=-1)
                lws_new[..., idx] = np.mean(lws[..., inds == idx], axis=-1)

            freqs, ints, lws = freqs_new, ints_new, lws_new
        else:
            ir_syms = np.atleast_1d(np.array(ir_syms, dtype=str, copy=True))

        # Convert frequnecy units. This will check that sp_units are
        # supported.

        sp_units = str(sp_units).lower()

        freqs = convert_frequency_energy_units(freqs, "thz", sp_units)
        lws = convert_frequency_energy_units(lws, "thz", sp_units)

        # adjust_band_params_for_spectrum_type() checks sp_type and
        # handles shape matching via broadcasting.

        sp_type = str(sp_type).lower()

        shifts, ints, lws, ir_syms = adjust_band_params_for_spectrum_type(
            freqs, ints, lws, ir_syms, sp_type=sp_type
        )

        # modulate_intensities_for_spectrum() checks w and t.

        w, t = float(w), float(t)

        if w <= 0.0 or t <= 0.0:
            raise ValueError("w and t must be > 0.")

        # Handle parameters for optional averaging over crystallite
        # orientations.

        if or_a is not None:
            if ints.ndim < 2:
                raise ValueError(
                    "Averaging over orientations requires at least two "
                    "dimensional input."
                )

            or_a = np.atleast1d(np.array(or_a, dtype=np.float64, copy=True))

            if or_a.shape != (ints.shape[-2], 3):
                raise ValueError(
                    "or_a must be an array_like with shape (L, 3)."
                )

            if or_w is not None:
                or_w = np.atleast_1d(
                    np.array(or_w, dtype=np.float64, copy=True)
                )
                or_w = np.broadcast_to(or_w, len(or_a))

                if (or_w.sum() - 1.0) < ZERO_TOLERANCE:
                    raise ValueError("or_w must sum to 1.")
            else:
                or_w = np.full(len(or_a), 1.0 / len(or_a), dtype=np.float64)
        else:
            # Set or_w to None if or_a is not specified.

            or_w = None

        # Set frequency grid for simulating spectra. auto_range() will
        # check and use sp_range and sp_res if set.

        (s_min, s_max), s_res = auto_range(
            shifts,
            lws,
            f_range=sp_range,
            f_res=sp_res,
        )

        s = np.arange(s_min, s_max + s_res / 10.0, s_res)

        # Mark arrays/views as readonly.

        freqs.flags.writeable = False
        ints.flags.writeable = False
        lws.flags.writeable = False
        ir_syms.flags.writeable = False

        if or_a is not None:
            or_a.flags.writeable = False
            or_w.flags.writeable = False

        s.flags.writeable = False

        self._band_shifts = freqs
        self._band_ints = ints
        self._band_lws = lws
        self._band_ir_syms = ir_syms

        self._band_xs = None

        self._w = w
        self._t = t
        self._sp_type = sp_type

        self._or_a = or_a
        self._or_w = or_w

        self._s = s
        self._xs = None

        self._s_units = sp_units

    def _lazy_set_band_x_sects(self):
        """Compute band cross sections and set the _band_xs field on
        first call to `band_cross_sections`."""

        if self._band_xs is None:
            band_xs = intensities_to_cross_sections(
                self._band_shifts,
                self._band_ints,
                self._w,
                self._t,
                s_units=self._s_units,
            )

            band_xs.flags.writeable = False

            self._band_xs = band_xs

    def _lazy_set_xs(self):
        """Compute cross sections and set the _xs field on first call to
        `cross_sections`."""

        if self._xs is None:
            # lorentzian() accepts variable batch dimensions and returns
            # (O,) + (...) -> (..., O). We then sum over the
            # second-to-last axis to sum the spectra over bands.

            xs = lorentzian(
                self._s,
                self._band_ints,
                self._band_shifts,
                self._band_lws,
            ).sum(axis=-2)

            # If we're averaging over orientations, we perform a
            # weighted sum over the last axis.

            if self._or_w is not None:
                xs = np.average(xs, axis=-2, weights=self._or_w)

            xs = intensities_to_cross_sections(
                self._s, xs, self._w, self._t, s_units=self._s_units
            )

            xs.flags.writeable = False

            self._xs = xs

    def _get_column_hdrs_for_spectral_dims(self, hdrs=None):
        """Get column headers for the product over spectral dimensions
        for outputting peak tables and spectra.

        Params
        ------
        hdrs : array_like or None
            Headers for each column in the product over spectral
            dimensions.

        Notes
        -----
        The number of spectral dimensions is given by the `num_dims`
        property.
        """

        n_d = self.num_spectral_dims

        if hdrs is not None:
            hdrs = np.atleast_1d(np.asanyarray(hdrs, dtype=str))

            if (n_d == 0 and len(hdrs) > 0) or (
                (len(hdrs) != np.prod(self._band_ints.shape[:n_d]))
            ):
                raise ValueError(
                    "hdrs must only be specify headers for each column "
                    "in the product over outer dimensions."
                )
        else:
            if n_d > 0:
                # Default to labelling with indices.

                hdrs = [
                    "_".join(str(i) for i in idx + 1)
                    for idx in np.ndindex(*self._band_ints.shape[:n_d])
                ]
            else:
                # Single spectral dimension -> empty array.
                hdrs = []

        return hdrs

    @property
    def band_raman_shifts(self):
        """numpy.ndarray : Band Raman shifts (shape: `(..., N)`)."""
        return self._band_shifts

    @property
    def band_intensities(self):
        """numpy.ndarray : Band intensities (shape: `(..., N)`)."""
        return self._band_ints

    @property
    def band_cross_sections(self):
        """numpy.ndarray : Band differential cross sections (shape:
        `(..., N)`)."""

        self._lazy_set_band_x_sects()
        return self._band_xs

    @property
    def band_linewidths(self):
        """numpy.ndarray : Band linewidths (shape: `(..., N)`)."""
        return self._band_lws

    @property
    def band_irrep_symbols(self):
        """numpy.ndarray : Band irreducible representation (irrep)
        symbols (shape: `(..., N)`)."""
        return self._band_ir_syms

    @property
    def measurement_wavelength(self):
        """float : Measurement wavelength (nm)."""
        return self._w

    @property
    def measurement_temperature(self):
        """float : Measurement temperature (K)."""
        return self._t

    @property
    def spectrum_type(self):
        """str : Spectrum type."""
        return self._sp_type

    @property
    def orientation_angles(self):
        """numpy.ndarray or None : Euler angles of crystallite orientations."""
        return self._or_a

    @property
    def orientation_averaging_weights(self):
        """numpy.ndarray or None : Averaging weights for crystallite
        orientations."""
        return self._or_w

    @property
    def raman_shift(self):
        """numpy.ndarray : Raman shifts in `frequency_units` (shape:
        `(O,)`)."""
        return self._s

    @property
    def cross_section(self):
        """numpy.ndarray : Cross sections in Ang^2 / sr (shape:
        `(..., O)`)."""

        self._lazy_set_xs()
        return self._xs

    @property
    def raman_shift_units(self):
        """str : Raman shift units."""
        return self._s_units

    @property
    def orientation_average(self):
        """bool : `True` if the spectrum is an average over crystallite
        orientations, otherwise `False`."""
        return self._or_w is not None

    @property
    def num_spectral_dims(self):
        """int : Number of spectral dimensions (...)."""

        d = self._band_ints.ndim
        return d - 2 if self.orientation_average else d - 1

    @property
    def raman_shift_unit_text_label(self):
        """str : Raman shift unit label suitable for plain-text output."""
        return get_frequency_energy_unit_text_label(self._s_units)

    @property
    def raman_shift_unit_plot_label(self):
        """str : Raman shift unit label suitable for plotting (may
        contain TeX strings)."""
        return get_frequency_energy_unit_plot_label(self._s_units)

    def peak_table(self, or_stack=True):
        """Generate a peak table with the band Raman shifts,
        intensities/cross sections, linewidths and irrep symbols as a
        Pandas `DataFrame`.

        Params
        ------
        or_stack : bool
            If `True`, stack orientations into a single axis to match
            the data shape of the spectra (default: `True`).

        Returns
        -------
        df : DataFrame
            Peak table.
        """

        shifts = self.band_raman_shifts
        ints = self.band_intensities
        x_sects = self.band_cross_sections
        lws = self.band_linewidths
        ir_syms = self._band_ir_syms

        if self.orientation_average and or_stack:
            shape = shifts.shape[:-2] + (-1,)

            shifts = shifts.reshape(shape)
            ints = ints.reshape(shape)
            x_sects = x_sects.reshape(shape)
            ir_syms = ir_syms.reshape(shape)

        data = {}

        if ints.ndim == 1:
            # Single set of peak parameters.

            data["raman_shift"] = shifts
            data["intensity"] = ints
            data["cross_section"] = x_sects
            data["linewidth"] = lws
            data["irrep"] = ir_syms
        else:
            # Multiple sets of peak parameters.

            hdrs = None

            if self.orientation_average and not or_stack:
                oa_hdrs = [
                    "or{0}".format(i + 1) for i in range(len(self._or_a))
                ]

                if self.num_spectral_dims > 0:
                    prod = product(
                        self._get_column_hdrs_for_spectral_dims(), oa_hdrs
                    )

                    hdrs = ["{0}_{1}".format(h, oa_h) for h, oa_h in prod]
                else:
                    hdrs = oa_hdrs
            else:
                hdrs = self._get_column_hdrs_for_spectral_dims()

            for idx, h in zip(np.ndindex(ints.shape), hdrs):
                data[f"raman_shift_{h}"] = shifts[idx]
                data[f"intensity_{h}"] = ints[idx]
                data[f"cross_section_{h}"] = x_sects[idx]
                data[f"linewidth_{h}"] = lws[idx]
                data[f"irrep_{h}"] = ir_syms[idx]

        return pd.DataFrame(data)

    def spectrum(self):
        """Return the simulated spectra as a Pandas `DataFrame`.

        Returns
        -------
        df : DataFrame
            Spectrum.
        """

        data = {"raman_shift": self.raman_shift}

        xs = self.cross_section

        if xs.ndim == 1:
            # 1D spectrum.
            data["cross_section"] = xs
        else:
            # ND spectrum.

            for inds, s in zip(
                np.ndindex(*xs.shape),
                self._get_column_hdrs_for_spectral_dims(),
            ):
                data[f"cross_section_{s}"] = xs[inds]

        return pd.DataFrame(data)

    def stack_with(self, *spectra, sp_type=None, sp_units=None):
        """Stack the calling spectrum with one or more other spectra
        into a new multi-dimensional spectrum with `num_dims + 1`
        spectral dimensions.

        Params
        ------
        *spectra : NDRamanSpectrum
            Spectra to concatenate with.
        sp_type : {"stokes", "anti-stokes", "both"} or None, optional
            Spectrum type (default: same as calling spectrum).
        sp_units : str or None, optional
            Frequency units of simulated spectrum (default: same as
            calling spectrum).

        Returns
        -------
        sp : NDRamanSpectrum
            Concatenated spectrum.

        Notes
        -----
        For spectra to be "stackable" the following conditions must
        hold:

        * The band parameters must have the same shape.
        * The measurement wavelength and temperature must be the same.
        * The spectra must be of the same type (i.e., Stokes,
          anti-Stokes or both).
        * If the spectra are orientational averages, they must be
          averaged over the same set of angles with the same averaging
          weights.

        An "inclusive" grid of frequencies (Raman shifts) is found from
        the minimum/maximum frequency and smallest step across all the
        spectra.

        To correctly initialise the returned NDRamanSpectrum object, the
        frequencies/linewidths need to be converted to THz and the the
        transformations for `spectrum_type="anti-stokes"` and
        `spectrum_type="both"` need to be reversed. This can lead to
        unnecessary work if the calling object and spectra are only
        required to generate the stacked spectrum. To avoid this, the
        spectrum type and frequency units of the calling spectrum can be
        overridden with the `sp_type` and `sp_units` keywords.
        """

        for sp in spectra:
            # The frequencies, linewidths and intensities all have the
            # same shape, so we only need to check one.

            if sp.band_intensities.shape != self._band_ints.shape:
                raise ValueError(
                    "All spectra must have the same data shape as the "
                    "calling spectrum."
                )

            if (
                sp.measurement_wavelength != self._w
                or sp.measurement_temperature != self._t
            ):
                raise ValueError(
                    "All spectra must have the same measurement "
                    "wavelength and temperature as the calling spectrum."
                )

            if sp.spectrum_type != self._sp_type:
                raise ValueError(
                    "All spectra must be of the same spectrum_type as "
                    "the calling spectrum."
                )

            if self.orientation_average:
                if not sp.orientation_average or not (
                    np.allclose(
                        self._or_a, sp.orientation_angles, atol=ZERO_TOLERANCE
                    )
                    and np.allclose(
                        self._or_w,
                        sp.orientation_averaging_weights,
                        atol=ZERO_TOLERANCE,
                    )
                ):
                    raise ValueError(
                        "If the calling spectrum is (not) an "
                        "orientational average, all spectra must also "
                        "(not) be an orientational average."
                    )

        f_min, f_max = self._s.min(), self._s.max()
        f_res = np.mean(self._s[1:] - self._s[:-1])

        f_min = f_res * np.floor(f_min / f_res)
        f_max = f_res * np.ceil(f_max / f_res)

        s = (len(spectra) + 1,) + self._band_shifts.shape

        freqs = np.zeros(s, dtype=np.float64)
        ints = np.zeros(s, dtype=np.float64)
        lws = np.zeros(s, dtype=np.float64)
        ir_syms = np.zeros(s, dtype=str)

        for i, sp in enumerate((self,) + spectra):
            freqs[i] = convert_frequency_energy_units(
                sp.band_raman_shifts, sp.raman_shift_units, "thz"
            )

            ints[i] = sp.band_intensities

            lws[i] = convert_frequency_energy_units(
                sp.band_linewidths, sp.raman_shift_units, "thz"
            )

            ir_syms[i] = sp.band_irrep_symbols

        # If sp_type is not Stokes, "undo" the data transformation to
        # ensure correct handling in the new spectrum.

        match self._sp_type:
            case "stokes":
                pass
            case "anti-stokes":
                freqs = -1.0 * freqs[..., ::-1]
                ints = ints[..., ::-1]
                lws = lws[..., ::-1]
                ir_syms = ir_syms[..., ::-1]
            case "both":
                idx = self._band_shifts.shape[-1] // 2

                freqs = freqs[..., idx:]
                ints = ints[..., idx:]
                lws = lws[..., idx:]
                ir_syms = ir_syms[..., idx:]
            case _:
                raise RuntimeError(
                    f'Unknown spectrum_type "{self._sp_type}". (This '
                    "may be a bug.)"
                )

        if sp_type is None:
            sp_type = self._sp_type

        if sp_units is None:
            sp_units = self._s_units

        return NDRamanSpectrum(
            freqs,
            ints,
            lws,
            self._w,
            self._t,
            irreps=None,
            ir_syms=ir_syms,
            sp_type=sp_type,
            sp_units=sp_units,
            sp_range=(f_min, f_max),
            sp_res=f_res,
            or_a=self._or_a,
            or_w=self._or_w,
        )


# -------------------
# RamanSpectrum class
# -------------------


class RamanSpectrum(NDRamanSpectrum):
    """Class for simulating standard 1D Raman spectra."""

    def __init__(
        self,
        freqs,
        ints,
        lws,
        w,
        t,
        irreps=None,
        ir_syms="?",
        sp_type="stokes",
        sp_units="thz",
        sp_range=None,
        sp_res=None,
        or_a=None,
        or_w=None,
    ):
        r"""Create a new instance of the `NDRamanSpectrum` class.

        Parameters
        ----------
        freqs : scalar or array_like
            Frequencies (shape: `(N,)` or `(L, N)`).
        ints : scalar or array_like
            Intensities (shape: `(N,)` or `(L, N)`).
        lws : scalar, array_like or None, optional
            Linewidths (scalar or shape: `(N,)` or `(L, N)`).
        w, t : float
            Measurement wavelength (nm) and temperature (K).
        irreps : Irreps or None, optional
            Irreducible representations (irreps) of the bands included
            in the spectrum (assumed to apply to all dimensions for
            > 1D spectra, overrides ir_syms).
        ir_syms : scalar or array_like, optional
            Irrep symbols of the bands in the spectrum (scalar or shape:
            `(N,)`, default: "?").
        sp_type : {"stokes", "anti-stokes", "both"}. optional
            Spectrum type (default: "stokes").
        sp_units : str, optional
            Frequency units of simulated spectrum (default: "thz").
        sp_range : tuple of float, optional
            `(min, max)` for simulated spectrum (default: automatic).
        sp_res : float, optional
            Resolution of simulated spectrum (default: automatic).
        or_a : array_like or None, optional
            Specifies the three Euler angles \phi, \theta and \psi for
            each orientation (L) (shape: `(L, 3)`).
        or_w : array_like or None, optional
            Optional averaging weights for orientations (L) (shape:
            `(L,)`, defaults to 1 / L).
        """

        super().__init__(
            freqs,
            ints,
            lws,
            w,
            t,
            irreps=irreps,
            ir_syms=ir_syms,
            sp_type=sp_type,
            sp_units=sp_units,
            sp_range=sp_range,
            sp_res=sp_res,
            or_a=or_a,
            or_w=or_w,
        )

        if self.num_spectral_dims != 0:
            raise ValueError(
                "Standard Raman spectra must not be multidimensional "
                "- use NDRamanSpectrum of another derived class "
                "instead."
            )


# --------------------------------
# AngleResolvedRamanSpectrum class
# --------------------------------


class AngleResolvedRamanSpectrum(NDRamanSpectrum):
    """Class for simulating angle-resolved 2D Raman spectra."""

    def __init__(
        self,
        freqs,
        ints,
        lws,
        w,
        t,
        rot_angles,
        irreps=None,
        ir_syms="None",
        sp_type="stokes",
        sp_units="thz",
        sp_range=None,
        sp_res=None,
        or_a=None,
        or_w=None,
        angle_unit_text_label="angle / deg",
        angle_unit_plot_label=r"Angle / $^\circ$",
    ):
        r"""Create a new instance of the `AngleResolvedRamanSpectrum`
        class.

        Parameters
        ----------
        freqs : array_like
            Frequencies (shape: `(N,)`, `(M, N)` or `(M, L, N)`).
        ints : array_like
            Intensities (shape: `(M, L, N)`).
        lws : scalar or array_like
            Linewidths (scalar or shape: `(N,)`, `(M, N)` or
            `(M, L, N)`).
        w, t : float
            Measurement wavelength (nm) and temperature (K).
        rot_angles : array_like
            Rotation angles in degrees (shape: `(M,)`).
        irreps : Irreps or None, optional
            Irreducible representations (irreps) of the bands included
            in the spectrum (assumed to apply to all dimensions for
            > 1D spectra, overrides ir_syms).
        ir_syms : scalar, array_like or None, optional
            Irrep symbols of the bands included in the spectrum
            (scalar or shape: `(..., N)`, default: "None").
        sp_type : {"stokes", "anti-stokes", "both"}. optional
            Spectrum type (default: "stokes").
        sp_units : str, optional
            Frequency units of simulated spectrum (default: "thz").
        sp_range : tuple of float, optional
            `(min, max)` for simulated spectrum (default: automatic).
        sp_res : float, optional
            Resolution of simulated spectrum (default: automatic).
        or_a : array_like or None, optional
            Specifies the three Euler angles \phi, \theta and \psi for
            each orientation (L) (shape: `(L, 3)`).
        or_w : array_like or None, optional
            Optional averaging weights for orientations (L) (shape:
            `(L,)`, defaults to 1 / L).
        angle_unit_text_label, angle_unit_plot_label : str, optional
            Angle unit labels for plain-text output and plotting.
        """

        super().__init__(
            freqs,
            ints,
            lws,
            w,
            t,
            irreps=irreps,
            ir_syms=ir_syms,
            sp_type=sp_type,
            sp_units=sp_units,
            sp_range=sp_range,
            sp_res=sp_res,
            or_a=or_a,
            or_w=or_w,
        )

        if self.num_spectral_dims != 1:
            raise ValueError(
                "Angle resolved Raman spectra must have a single "
                "spectral dimension - use NDRamanSpectrum of another "
                "derived class instead."
            )

        rot_angles = np.atleast_1d(
            np.array(rot_angles, dtype=np.float64, copy=True)
        )

        if rot_angles.ndim != 1 or len(rot_angles) != len(self._band_ints):
            raise ValueError("angles must be an array_like with shape (M,).")

        rot_angles.flags.writeable = False

        self._rot_angles = rot_angles

        self._angle_unit_text_label = str(angle_unit_text_label)
        self._angle_unit_plot_label = str(angle_unit_plot_label)

    def _get_column_hdrs_for_spectral_dims(self, hdrs=None):
        if hdrs is None:
            hdrs = ["{0:.2f}".format(a) for a in self._rot_angles]

        return super()._get_column_hdrs_for_spectral_dims(hdrs)

    @property
    def rotation_angles(self):
        """numpy.ndarray : Angles in degrees (shape: `(M,)`)."""
        return self._rot_angles

    @property
    def rotation_angle_unit_text_label(self):
        """str : Angle unit label suitable for plain-text output."""
        return self._angle_unit_text_label

    @property
    def rotation_angle_unit_plot_label(self):
        """str : Angle unit label suitable for plotting (contains TeX
        strings)."""
        return self._angle_unit_plot_label
