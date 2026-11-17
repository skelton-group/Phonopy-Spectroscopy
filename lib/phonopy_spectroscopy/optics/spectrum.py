# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Base classes for simulated optical spectra."""

# -------
# Imports
# -------

import abc

import numpy as np
import pandas as pd

from ..units import (
    get_supported_frequency_energy_units,
    get_frequency_energy_unit_text_label,
    get_frequency_energy_unit_plot_label,
)

# -------------------------
# OpticalSpectrumBase class
# -------------------------


class OpticalSpectrumBase(abc.ABC):
    """Abstract base class for optical spectra."""

    def __init__(self, e, t, n_f, n_b, e_units="ev"):
        """Create a new instance of the `OpticalSpectrumBase` class.

        Parameters
        ----------
        e : array_like
            Frequencies in `e_units` (shape: `(O,)`).
        t : float
            Thickness (path length) in mm (default: 1 mm).
        n_f, n_b : float, optional
            Refractive indices of the front (indicent) and back (exit)
            media.
        e_units : str, optional
            Energy units (default: "ev").
        """

        e = np.atleast_1d(np.array(e, dtype=np.float64, copy=True))

        if e.ndim > 1 or len(e) == 0:
            raise ValueError("f must be an array_like with shape (O,).")

        if t <= 0.0:
            raise ValueError("t must be > 0.")

        if e_units not in get_supported_frequency_energy_units():
            raise ValueError(f'f_units="{e_units}" are not supported.')

        e.flags.writeable = False

        self._e = e

        self._t = t
        self._n_f = n_f
        self._n_b = n_b

        self._e_units = e_units

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

        if self._ref_s is None:
            self._init_single_reflectivity()

        self._ref_s.flags.writeable = False

    def _lazy_init_total_reflectivity_and_transmission(self):
        """Lazy initialisation of total reflectivity, transission, and
        associated quantities."""

        if self._ref_t is None or self._trans is None:
            self._init_total_reflectivity_and_transmission()

        self._ref_t.flags.writeable = False
        self._trans.flags.writeable = False

        if self._aps is None or self._abs is None:
            aps = 1.0 - (self._ref_t + self._trans)
            abs = -1.0 * np.log10(self._trans)

            aps.flags.writeable = False
            abs.flags.writeable = False

            self._aps = aps
            self._abs = abs

    @property
    def energy(self):
        """numpy.ndarray : Energies in `e_units` (shape: `(O,)`)."""
        return self._e

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
    def energy_units(self):
        """str : Energy units."""
        return self._e_units

    @property
    def energy_unit_text_label(self):
        """str : Energy unit label suitable for plain-text output."""
        return get_frequency_energy_unit_text_label(self._e_units)

    @property
    def energy_unit_plot_label(self):
        """str : Energy unit label suitable for plotting (may contain
        TeX strings)."""
        return get_frequency_energy_unit_plot_label(self._e_units)

    @property
    def single_reflectivity(self):
        """numpy.ndarray : Single (infinite bulk) reflectivity (shape:
        `(O,)`)."""

        self._lazy_init_single_reflectivity()
        return self._ref_s

    @property
    def total_reflectivity(self):
        """numpy.ndarray : Total reflectivity (shape: `(O,)`)."""

        self._lazy_init_total_reflectivity_and_transmission()
        return self._ref_t

    @property
    def transmission(self):
        """numpy.ndarray : Transmission (shape: `(O,)`)."""

        self._lazy_init_total_reflectivity_and_transmission()
        return self._trans

    @property
    def absorptance(self):
        """numpy.ndarray : Absorptance (shape: `(O,)`)."""

        self._lazy_init_total_reflectivity_and_transmission()
        return self._aps

    @property
    def absorbance(self):
        """numpy.ndarray : Decadic absorbance (shape: `(O,)`)."""

        self._lazy_init_total_reflectivity_and_transmission()
        return self._abs

    def spectrum(self):
        """Return the optical spectra as a Pandas `DataFrame`.

        Returns
        -------
        df : pandas.DataFrame
            `DataFrame` containing the optical spectra.
        """

        d = {
            "energy": self._e,
            "single_reflectivity": self.single_reflectivity,
            "total_reflectivity": self.total_reflectivity,
            "transmission": self.transmission,
            "absorptance": self.absorptance,
            "absorbance": self.absorbance,
        }

        return pd.DataFrame(d)
