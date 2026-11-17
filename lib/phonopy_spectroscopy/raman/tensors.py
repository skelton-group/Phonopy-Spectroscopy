# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Class for storing and interpolating energy-dependent Raman tensors."""

# -------
# Imports
# -------

import warnings

import numpy as np

from scipy.interpolate import interp1d

from ..constants import ZERO_TOLERANCE
from ..gamma_phonons import gamma_phonons_from_dict, FullPolarGammaPhonons
from ..units import nm_to_ev
from ..utility.geometry import parse_direction
from ..utility.io_helper import complex_array_to_dict, dict_to_complex_array

# ------------------
# RamanTensors class
# ------------------


class RamanTensors:
    """Class to store a phonon calculation with Raman tensors and get
    Raman tensors for experimental measurement setups."""

    def __init__(self, ph_calc, r_t, e=0.0, band_inds=None):
        """Create a new instance of the `RamanTensors` class.

        Parameters
        ----------
        ph_calc : GammaPhonons or PolarGammaPhonons
            Phonon calculation.
        r_t : array_like
            Single Raman tensor, `N` tensors or `N` tensors at `M`
            photon energies (shapes: `(3, 3)`, `(N, 3, 3)` or
            `(N, M, 3, 3)`).
        e : scalar or array_like, optional
            Photon energies in eV (scalar or shape: `(M,)`).
        band_inds : scalar, array_like or None, optional
            Band indices for which the Raman tensors were calculated
            (default: assume all bands).
        acc_mode_tol : str or float, optional
            Tolerance for identifying the acoustic modes to check
            whether the calculated bands form a complete basis
            (default: "auto").

        Notes
        -----
        `e` is required when `r_t` have shape `(N, M, 3, 3)` and must be
        monotonically increasing. When `e` is not specified, Raman
        tensors with shape `(3, 3)` or `(N, 3, 3)` are assumed to be
        calculated at E = 0.

        Note that the default energy for `get_tensors_for_measurement()`
        is zero (i.e., the far-from-resonance approximation). This will
        fail if `e` does not include E = 0 and therefore produces a
        warning.
        """

        r_t = np.array(r_t, dtype=np.complex128, copy=True)

        match r_t.shape:
            case (3, 3):
                r_t = r_t[np.newaxis, np.newaxis, :, :]
            case (_, 3, 3):
                r_t = r_t[:, np.newaxis, :, :]
            case (_, _, 3, 3):
                pass
            case _:
                raise ValueError(
                    "r_t must be an array_like with shape (3, 3), "
                    "(N, 3, 3) or (N, M, 3, 3)."
                )

        e = np.atleast_1d(np.array(e, dtype=np.float64, copy=True))

        if len(e) != r_t.shape[1]:
            raise ValueError("e must have shape (M,).")

        if (e < 0.0).any():
            raise ValueError("e cannot be negative.")

        if len(e) > 1 and (e[1:] < e[:-1]).any():
            raise ValueError("e must be monotonically increasing.")

        if np.abs(e[0]) > ZERO_TOLERANCE:
            warnings.warn(
                "r_t does not include tensors for the "
                "far-from-resonance approximation (E = 0). Calls to "
                "get_tensors_for_measurement() with the default e=0.0 "
                "will fail.",
                RuntimeWarning,
            )

        if band_inds is not None:
            band_inds = np.atleast_1d(
                np.array(band_inds, dtype=int, copy=True)
            )

            if (band_inds < 0).any() or (band_inds >= ph_calc.num_modes).any():
                raise ValueError(
                    "One or more indices in band_inds are invalid."
                )
        else:
            band_inds = np.arange(0, ph_calc.num_modes, dtype=int)

        if r_t.shape[0] != len(band_inds):
            raise ValueError(
                "band_inds with shape (N,) must be supplied, or r_t "
                "must provide tensors for all bands in ph_calc."
            )

        # If ph_calc has irreps, generate a new Irreps object for the
        # subset of bands in the calculation. If band_inds spans
        # multiple irrep groups, this will raise.

        irreps = None

        if ph_calc.irreps is not None:
            _, irreps = ph_calc.irreps.get_subset(band_inds)

        # Determine whether generating Raman tensors with NAC is
        # possible - this requires that ph_calc is a
        # FullPolarGammaPhonons and that the tensors are calculated for
        # all Raman-active bands.

        supports_nac = (
            isinstance(ph_calc, FullPolarGammaPhonons)
            and np.isin(ph_calc.raman_band_indices, band_inds).all()
        )

        # Mark arrays as readonly.

        e.flags.writeable = False
        r_t.flags.writeable = False
        band_inds.flags.writeable = False

        self._ph_calc = ph_calc

        self._e = e
        self._r_t = r_t
        self._band_inds = band_inds

        self._band_irreps = irreps

        self._supports_nac = supports_nac

        # Fields for maintaining a cache key/data.

        self._cache_k = None
        self._cache_d = None

    def _get_tensor_energy(self, e, w):
        """Determine the energy at which to calculate Raman tensors.

        Params
        ------
        e : float
            Energy.
        w : float or None
            Wavelength (overrides `e`).

        Returns
        -------
        e_rt : float
            Energy at which to return Raman tensors.
        """

        if e < 0.0:
            raise ValueError("e must be >= 0.")

        if w is not None:
            if w <= 0.0:
                raise ValueError("w must be > 0.")

            e = nm_to_ev(w)

        return e

    @property
    def structure(self):
        """Structure : Crystal structure."""
        return self._ph_calc.structure

    @property
    def phonon_calculation(self):
        """GammaPhonons : Phonon calculation."""
        return self._ph_calc

    @property
    def raman_tensors(self):
        """numpy.ndarray : Raman tensors in Ang^2 / sqrt(amu) (shape:
        `(N, M, 3, 3)`)."""
        return self._r_t

    @property
    def energies(self):
        """numpy.ndarray : Photon energies in eV (shape: `(M,)`)."""
        return self._e

    @property
    def band_indices(self):
        """numpy.ndarray : Indices of the bands for which Raman tensors
        are available."""
        return self._band_inds

    @property
    def irreps(self):
        """Irreps or None : Irreducible representations (irreps) of the
        bands for which Raman tensors are available."""
        return self._band_irreps

    @property
    def num_modes(self):
        """int : Number of Raman tensors, N."""
        return self._r_t.shape[0]

    @property
    def num_energies(self):
        """int : Number of photon energies, M."""
        return self._r_t.shape[1]

    @property
    def energy_dependent(self):
        """bool : `True` if the Raman tensors are energy-dependent,
        otherwise `False`."""
        return len(self._e) > 1

    @property
    def supports_nac_tensors(self):
        """bool : `True` if generating tensors with non-analytical
        corrections (NAC) is supported, otherwise `False`."""
        return self._supports_nac

    def get_tensors(self, e=0.0, w=None):
        """Get Raman tensors at the supplied photon energy/wavelength.

        Parameters
        ----------
        e : float, optional
            Measurement photon energy (default: 0.0).
        w : float or None, optional
            Measurement wavelength (default: `None`, overrides `e` if
            specified).

        Returns
        -------
        r_t : numpy.ndarray
            Raman tensors (shape: `(N, 3, 3)`).

        Notes
        -----
        The default behaviour of this function is to return Raman
        tensors at E = 0. This is the typical starting point for
        simulating spectra using the far-from-resonance approximation.
        """

        e_rt = self._get_tensor_energy(e, w)

        # Return cached tensors if possible.

        if self._cache_k is not None and self._cache_d is not None:
            q_c, e_rt_c = self._cache_k

            if q_c is None and np.isclose(e_rt, e_rt_c, atol=ZERO_TOLERANCE):
                _, r_t = self._cache_d
                return r_t.copy()

        r_t = None

        # If e matches a stored energy exactly, return the tensors
        # calculated for that energy. Otherwise, attempt a linear
        # interpolation.

        mask = np.isclose(self._e, e_rt, atol=ZERO_TOLERANCE)

        if mask.any():
            # np.nonzero() should always find one index because we
            # require e to be monotonically increasing.

            idx = np.flatnonzero(mask)[0]
            r_t = self._r_t[:, idx, :, :].copy()
        else:
            if len(self._e) == 1:
                raise RuntimeError(
                    "Cannot interpolate tensors calculated at a single energy."
                )

            if e_rt < self._e.min() or e_rt > self._e.max():
                raise ValueError(
                    "e/w is outside the bounds of the tensor energies."
                )

            f = interp1d(self._e, self._r_t, kind="linear", axis=1)
            r_t = f(e_rt)

        # At some energies, and in particular at e = 0, the Raman
        # tensors are likely to be real. In this case, discarding the
        # real part allows further calculations to use simplified
        # arithmetic.

        r_t = np.real_if_close(r_t, tol=ZERO_TOLERANCE)

        self._cache_k = (None, e_rt)
        self._cache_d = (self._ph_calc, r_t.copy())

        return r_t

    def get_tensors_with_nac(self, q, e=0.0, w=None):
        """Get Raman tensors at the supplied photon energy/wavelength
        including a non-analytical correction (NAC).

        Parameters
        ----------
        q : array_like or str
            Approach direction for calculating NAC.
        e : float, optional
            Measurement photon energy (default: 0.0).
        w : float or None, optional
            Measurement wavelength (default: `None`, overrides `e` if
            specified).

        Returns
        -------
        res : tuple of (PolarGammaPhononsBase, numpy.ndarray)
            Tuple of `(ph_calc, r_t)` containing a NAC-modified phonon
            calculation and the calculated Raman tensors.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `q`.

        Notes
        -----
        This function function applies the NAC to the underlying phonon
        calculation, then uses the modified eigenvectors to map the
        Raman tensors of the unmodified calcualation to the new basis.

        This assumes that the Raman tensors have been calculated for all
        Raman-active bands, and the other modes are assumed to be
        inactive.
        """

        if not self._supports_nac:
            raise RuntimeError(
                "Generating Raman tensors with a NAC is not supported "
                "because the underlying phonon_calculation does not "
                "support this and/or because band_indices do not "
                "include all the Raman-active modes."
            )

        q = parse_direction(q)
        e_rt = self._get_tensor_energy(e, w)

        # Return cached tensors if possible.

        if self._cache_k is not None and self._cache_d is not None:
            q_c, e_rt_c = self._cache_k

            if q_c is not None:
                if np.allclose(q, q_c, atol=ZERO_TOLERANCE) and (
                    np.isclose(e_rt, e_rt_c, atol=ZERO_TOLERANCE)
                ):
                    ph_calc, r_t = self._cache_d
                    return (ph_calc, r_t.copy())

        r_t = None

        # Generate a new phonon calculation with the NAC correction.

        ph_calc_nac = self._ph_calc.gamma_phonons_with_nac(parse_direction(q))

        # Get a complete set of tensors, assuming the tensors for modes
        # not included in the calculation have zero activity.

        r_t = np.zeros((self._ph_calc.num_modes, 3, 3), dtype=np.complex128)
        r_t[self._band_inds] = self.get_tensors(e=e_rt, w=None)

        # Project eigenvectors of split modes onto original ones.

        proj = ph_calc_nac.eigenvector_projection(self._ph_calc)

        # Generate Raman tensors for corrected modes from a sum over
        # the original modes weighted by the eigenvector projection
        # coefficients.

        r_t_new = np.real_if_close(
            np.einsum("ij,jkl->ikl", proj, r_t), tol=ZERO_TOLERANCE
        )

        # Cache phonon calculation and projected Raman tensors.

        self._cache_k = (q, e_rt)
        self._cache_d = (ph_calc_nac, r_t_new.copy())

        return (ph_calc_nac, r_t_new)

    def to_dict(self):
        """Return the internal data as a dictionary of native Python
        types for serialisation.

        Returns
        -------
        d : dict
            Dictionary structure containing internal data as native
            Python types.
        """

        return {
            "phonon_calculation": self._ph_calc.to_dict(),
            "raman_tensors": complex_array_to_dict(self._r_t),
            "energies": self._e.tolist(),
            "band_indices": self._band_inds.tolist(),
        }

    @staticmethod
    def from_dict(d):
        """Create a new `RamanTensors` instance from a dictionary
        generated by `RamanTensors.to_dict()`.

        Parameters
        ----------
        d : dict
            Dictionary generated by `to_dict()`.

        Returns
        -------
        r_t : RamanTensors
            `RamanTensors` object constructed from the data in `d`.
        """

        return RamanTensors(
            gamma_phonons_from_dict(d["phonon_calculation"]),
            dict_to_complex_array(d["raman_tensors"]),
            e=d["energies"],
            band_inds=d["band_indices"],
        )
