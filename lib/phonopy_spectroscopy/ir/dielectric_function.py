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

from ..constants import (
    VASP_INFRARED_DIELECTIC_TO_RELATIVE_PERMITTIVITY,
    ZERO_TOLERANCE,
)

from ..distributions import lorentz_oscillator
from ..optics.response_function import DielectricFunction
from ..units import convert_frequency_units

from ..utility.numpy_helper import atleast_3d
from ..utility.spectrum import auto_range

# --------------------------------
# InfraredDielectricFunction class
# --------------------------------


class InfraredDielectricFunction(DielectricFunction):
    """Class for generating simulated infrared dielectric functions
    from sets of frequencies, mode oscillator strengths and linewidths.
    """

    mode_oscillator_strength_unit_text_label = "e^2 amu^-1"
    """str : Mode oscillator strength unit label suitable for plain-text
    output."""

    def __init__(
        self,
        freqs,
        osc_strs,
        lws,
        cell_volume,
        eps_inf,
        irreps=None,
        ir_syms="?",
        f_units="thz",
        f_range=None,
        f_res=None,
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
        irreps : Irreps or None, optional
            Irreducible representations (irreps) of the bands included
            in the spectrum (overrides `ir_syms`).
        ir_syms : scalar, array_like or None, optional
            Irrep symbols of the bands included in the spectrum (scalar
            or shape: `(N,)`, default: "?").
        f_units : str
            Frequency units (default: "thz").
        f_range : tuple of float, optional
            `(min, max)` for simulated dielectric function (default:
            automatic).
        f_res : float, optional
            Resolution of simulated dielectric function (default:
            automatic).
        """

        freqs = np.atleast_1d(np.array(freqs, dtype=np.float64, copy=True))

        osc_strs = atleast_3d(
            np.array(osc_strs, dtype=np.complex128, copy=True)
        )

        lws = np.atleast_1d(np.array(lws, dtype=np.float64, copy=True))

        # Broadcast freqs, osc_strs and lws to a common N.

        n = max(len(freqs), osc_strs.shape[0], len(lws))

        freqs = np.broadcast_to(freqs, (n,))
        osc_strs = np.broadcast_to(osc_strs, (n, 3, 3))
        lws = np.broadcast_to(lws, (n,))

        # Phonon mode dipole oscillator strengths will usually be real,
        # but complex input is possible (and physical) in other use
        # cases.

        osc_strs = np.real_if_close(osc_strs, tol=ZERO_TOLERANCE)

        # If irreps are supplied, average frequencies/linewidths and sum
        # mode oscillator strengths.

        if irreps is not None:
            # "Reduce" N -> M bands by averaging
            # frequencies/linewidths and summing intensities.

            inds = irreps.irrep_band_indices
            ir_syms = np.copy(irreps.irrep_symbols)

            if len(inds) != n:
                raise ValueError(
                    "irreps must assign all N bands to M irrep groups."
                )

            # Average freqs/linewidths and sum intensities.

            m = len(ir_syms)

            freqs_new = np.zeros((m,), dtype=np.float64)
            osc_strs_new = np.zeros((m, 3, 3), dtype=np.float64)
            lws_new = np.zeros((m,), dtype=np.float64)

            for idx in np.unique(inds):
                freqs_new[idx] = np.mean(freqs[inds == idx])
                osc_strs_new[idx] = osc_strs[inds == idx].sum(axis=0)
                lws_new[idx] = np.mean(lws[inds == idx])

            freqs, osc_strs, lws = freqs_new, osc_strs_new, lws_new
        else:
            # Broadcast ir_syms against N.

            ir_syms = np.broadcast_to(
                np.atleast_1d(np.array(ir_syms, dtype=str)), (n,)
            )

        if cell_volume <= 0.0:
            raise ValueError("cell_volume must be > 0.")

        eps_inf = np.atleast_2d(np.array(eps_inf, dtype=np.float64, copy=True))

        if eps_inf.shape != (3, 3):
            raise ValueError(
                "eps_inf must be an array_like with shape (3, 3)."
            )

        # Determine energy axis range.

        (e_min, e_max), e_res = auto_range(
            freqs,
            lws,
            f_range=f_range,
            f_res=f_res,
            f_pad_mult=5.0,
            f_res_min_pts=1000,
        )

        e = np.arange(e_min, e_max + e_res / 10.0, e_res)

        eps = lorentz_oscillator(e, osc_strs, freqs, lws).sum(axis=0)
        eps *= VASP_INFRARED_DIELECTIC_TO_RELATIVE_PERMITTIVITY / cell_volume
        eps += eps_inf

        # Convert units.

        freqs = convert_frequency_units(freqs, "thz", f_units)
        lws = convert_frequency_units(lws, "thz", f_units)

        e = convert_frequency_units(e, "thz", f_units)

        super().__init__(e, eps, e_units=f_units)

        # Store additional data.

        freqs.flags.writeable = False
        osc_strs.flags.writeable = False
        lws.flags.writeable = False
        ir_syms.flags.writeable = False
        eps_inf.flags.writeable = False

        self._freqs = freqs
        self._osc_strs = osc_strs
        self._lws = lws
        self._ir_syms = ir_syms
        self._eps_inf = eps_inf

    @property
    def mode_frequencies(self):
        """numpy.ndarray : Mode frequencies in `energy_units`
        (shape: `(N,)`)."""
        return self._freqs

    @property
    def mode_oscillator_strengths(self):
        """numpy.ndarray : Mode oscillator strengths in e^2 amu^-1
        (shape: '(N, 3, 3)')."""
        return self._osc_strs

    @property
    def mode_linewidths(self):
        """numpy.ndarray : Mode linewidths in `energy_units` (shape:
        `(N,)`)."""
        return self._lws

    @property
    def mode_irrep_symbols(self):
        """numpy.ndarray : Mode irrep symbols (shape: `(N,)`)."""
        return self._ir_syms

    @property
    def epsilon_inf(self):
        r"""numpy.ndarray : High-frequency dielectric constant in \eps_0
        (shape: `(3, 3)`)."""
        return self._eps_inf

    @property
    def frequency(self):
        """numpy.ndarray : Frequencies in `frequency_units` (alias for
        base class `energy`)."""
        return self.energy

    @property
    def frequency_units(self):
        """str : Frequency units (alias for base class `energy_units`)."""
        return self.energy_units

    @property
    def frequency_unit_text_label(self):
        """str : Frequency unit label suitable for plain-text output
        (alias for base class `energy_unit_text_label`)."""
        return self.energy_unit_text_label

    @property
    def frequency_unit_plot_label(self):
        """str : Frequency unit label suitable for plotting (may contain
        TeX strings, alias for base class `frequency_unit_plot_label`).
        """
        return self.energy_unit_plot_label

    def peak_table(self):
        """Return the peak table as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the peak table.
        """

        d = {"freq_energy": self._freqs}

        osc_str = self._osc_strs

        for idx, x in enumerate("xyz"):
            for idy, y in enumerate("xyz"):
                d[f"osc_str_{x}{y}"] = osc_str[:, idx, idy]

        d["linewidth"] = self._lws
        d["irrep"] = self._ir_syms

        return pd.DataFrame(d)
