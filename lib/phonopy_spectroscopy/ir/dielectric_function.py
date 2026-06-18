# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Class implementing the calculation of the infrared dielectric
function."""

# -------
# Imports
# -------

import numpy as np
import pandas as pd

from ..constants import VASP_INFRARED_DIELECTIC_TO_RELATIVE_PERMITTIVITY
from ..distributions import lorentz_oscillator
from ..spectrum_base import GammaPhononSpectrumBase
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

_EPSILON_UNIT_TEXT_LABEL = r"\eps / \eps_0"
_EPSILON_UNIT_PLOT_LABEL = r"$\epsilon$ / $\epsilon_0$"

# --------------------------------
# InfraredDielectricFunction class
# --------------------------------


class InfraredDielectricFunction(GammaPhononSpectrumBase):
    """Base class for generating simulated infrared dielectric functions
    from sets of frequencies, mode oscillator strengths and linewidths.
    """

    def __init__(
        self, freqs, osc_strs, lws, cell_volume, eps_inf, irreps=None, **kwargs
    ):
        """Create a new instance of the `InfraredDielectricFunction`
        class.

        Parameters
        ----------
        freqs : array_like
            Frequencies in THz (shape: `(N,)`).
        osc_strs : array_like
            Dipole oscillator strengths in e^2 amu^-1 (shape:
            `(N, 3, 3)`) .
        lws : array_like
            Linewidths in THz (shape: `(N,)`).
        cell_volume : float
            Unit-cell volume in Ang^3.
        eps_inf : array_like
            High-frequency dielectric constant in units of relative
            permittivity (shape: `(3, 3)`).
        irreps : Irreps or None, optional
            `Irreps` object assigning bands to irrep groups.
        **kwargs
            Keyword arguments to the `GammaPhononSpectrumBase`
            constructor.
        """

        osc_strs, _ = np_expand_dims(
            np_asarray_copy(osc_strs, dtype=np.complex128), (None, 3, 3)
        )

        if not np.iscomplex(osc_strs).all():
            osc_strs = osc_strs.real

        # If irreps are supplied, average frequencies/linewidths and sum
        # mode oscillator strengths.

        irrep_syms = None

        if irreps is not None:
            ir_band_inds = irreps.irrep_band_indices

            freqs = np.array(
                [np.mean(freqs[inds]) for inds in ir_band_inds],
                dtype=np.float64,
            )

            osc_strs = np.array(
                [
                    np.sum(osc_strs[inds, :, :], axis=0)
                    for inds in ir_band_inds
                ],
                dtype=osc_strs.dtype,
            )

            lws = np.array(
                [np.mean(lws[inds]) for inds in ir_band_inds], dtype=np.float64
            )

            irrep_syms = np_asarray_copy(irreps.irrep_symbols, dtype=object)

        if cell_volume <= 0.0:
            raise ValueError("cell_volume must be positive and non-zero.")

        eps_inf = np_asarray_copy(eps_inf, dtype=np.float64)

        if not np_check_shape(eps_inf, (3, 3)):
            raise ValueError(
                "If supplied, eps_inf must be an array_like with shape "
                "(3, 3)."
            )

        # The high-frequency dielectric constant should be real.

        if np.iscomplex(eps_inf).any():
            raise ValueError("If supplied, eps_inf must be real.")

        # Call the GammaPhononSpectrumBase constructor to handle
        # "x-axis"-related intialisation.

        super(InfraredDielectricFunction, self).__init__(
            freqs, lws, irrep_syms=irrep_syms, **kwargs
        )

        # Store parameters.

        self._osc_strs = osc_strs
        self._cell_volume = cell_volume

        self._eps_inf = eps_inf

        # Initialise _eps field for "lazy" evaluation.

        self._eps = None

    def _lazy_init_epsilon(self):
        """Lazy initialisation of dielectric function."""

        if self._eps is None:
            # The conversion factor to relative permittivity assumes
            # oscillator strengths in e^2 / amu, volumes in Ang^2 and
            # frequencies in THz.

            x, freqs, lws = self.x, self.frequencies, self.linewidths

            if self._x_units != "thz":
                x = convert_frequency_units(x, self._x_units, "thz")
                freqs = convert_frequency_units(freqs, self._x_units, "thz")
                lws = convert_frequency_units(lws, self._x_units, "thz")

            eps = lorentz_oscillator(x, self._osc_strs[0], freqs[0], lws[0])

            for osc_str, freq, lw in zip(
                self._osc_strs[1:], freqs[1:], lws[1:]
            ):
                eps += lorentz_oscillator(x, osc_str, freq, lw)

            # Convert to relative permittivity.

            eps *= (
                VASP_INFRARED_DIELECTIC_TO_RELATIVE_PERMITTIVITY
                / self._cell_volume
            )

            # If a high-frequency dielectric constant was supplied
            # during initialisation, add it to the dielectric function.

            if self._eps_inf is not None:
                # Note that _eps_inf may be a ("regular") float.
                eps += self._eps_inf[np.newaxis, :, :]

            self._eps = eps

    @property
    def oscillator_strengths(self):
        """numpy.ndarray : Mode oscillator strengths in e^2 amu^-1
        (shape: '(N, 3, 3)')."""
        return np_readonly_view(self._osc_strs)

    @property
    def epsilon_inf(self):
        r"""numpy.ndarray : High-frequency dielectric constant in \eps_0
        (shape: `(3, 3)`)."""
        return np_readonly_view(self._eps_inf)

    @property
    def epsilon(self):
        r"""numpy.ndarray : Complex dielectric function in \eps_0
        (shape: `(O, 3, 3)`)."""

        self._lazy_init_epsilon()
        return np_readonly_view(self._eps)

    @property
    def mode_oscillator_strength_unit_text_label(self):
        """str : Mode oscillator strength unit label suitable for
        plain-text output."""

        return "e^2 amu^-1"

    @property
    def epsilon_unit_text_label(self):
        """str : Dielectric function unit label suitable for plain-text
        output."""

        return _EPSILON_UNIT_TEXT_LABEL

    @property
    def epsilon_unit_plot_label(self):
        """str : Dielectric function unit label suitable for plotting
        (contains TeX strings)."""

        return _EPSILON_UNIT_PLOT_LABEL

    def peak_table(self):
        """Return the peak table as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the peak table.
        """

        d = {"freq_energy": self._freqs, "linewidth": self._lws}

        if self._irrep_syms is not None:
            d["irrep"] = self._irrep_syms
        else:
            d["irrep"] = ["None"] * len(self._freqs)

        for i, j, label in InfraredDielectricFunction._DF_COL_INDS_HDRS:
            d["osc_str_{0}".format(label)] = self._osc_strs[:, i, j]

        return pd.DataFrame(d)

    def spectrum(self):
        """Return the simulated dielectric function as a Pandas
        `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the dielectric function.
        """

        self._lazy_init_epsilon()

        d = {"freq_energy": self.x}

        for i, j, label in InfraredDielectricFunction._DF_COL_INDS_HDRS:
            d["epsilon_re_{0}".format(label)] = self._eps.real[:, i, j]

        for i, j, label in InfraredDielectricFunction._DF_COL_INDS_HDRS:
            d["epsilon_im_{0}".format(label)] = self._eps.imag[:, i, j]

        return pd.DataFrame(d)

    _DF_COL_INDS_HDRS = [
        (0, 0, "xx"),
        (1, 1, "yy"),
        (2, 2, "zz"),
        (0, 1, "xy"),
        (0, 2, "xz"),
        (1, 2, "yz"),
    ]

    """list of tuples of (int, int, str) : Indices and column labels
    for constructing Pandas `DataFrame` objects from the unique
    components of the mode oscillator strengths/tensor dielectric
    function."""
