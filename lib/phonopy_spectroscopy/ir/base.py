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

from ..spectrum_base import SpectrumBase

from ..utility.numpy_helper import np_readonly_view

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

        if self._ref_s is None:
            self._init_single_reflectivity()

    def _lazy_init_total_reflectivity_and_transmission(self):
        """Lazy initialisation of total reflectivity, transission, and
        associated quantities."""

        if self._ref_t is None:
            self._init_total_reflectivity_and_transmission()

        if self._aps is None:
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
