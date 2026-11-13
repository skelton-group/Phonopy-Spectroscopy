# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Class for storing and working with Gamma-point phonon calculations."""

# -------
# Imports
# -------

import warnings

import numpy as np

from .constants import (
    VASP_TO_THZ,
    VASP_NAC_PREFACTOR,
    VASP_DIELECTRIC_TO_RELATIVE_PERMITTIVITY,
    ZERO_TOLERANCE,
)

from .irreps import Irreps
from .structure import Structure
from .utility.geometry import parse_direction
from .utility.io_helper import complex_array_to_dict, dict_to_complex_array
from .utility.numpy_helper import atleast_3d
from .utility.quadrature import unit_sphere_lebedev_quad_rule

from phonopy import Phonopy
from phonopy.harmonic.dynmat_to_fc import DynmatToForceConstants

# ---------
# Functions
# ---------


def mode_effective_charges(edisps, born_charges):
    """Calculate the mode effective charges from "eigendisplacements"
    (eigenvectors after division by sqrt(mass)) and Born
    effective-charge tensors.

    Parameters
    ----------
    edisps : array_like
        Eigendisplacements (shape: `(..., N, 3)`).
    born_charges : array_like
        Born effective-charge tensors (shape: `(N, 3, 3)`).

    Returns
    -------
    mode_eff_chg : numpy.ndarray
        Mode effective charges (shape: `(..., 3,)`).
    """

    edisps, born_charges = np.asanyarray(edisps), np.asanyarray(born_charges)

    if edisps.ndim < 2 or edisps.shape[-1] != 3:
        raise ValueError(
            "edisps must be an array_like with shape (..., N, 3)."
        )

    if (born_charges.ndim < 3) or (
        born_charges.shape[-3:] != (edisps.shape[-2], 3, 3)
    ):
        raise ValueError(
            "born_charges must be an array_like with shape (..., N, 3, 3)."
        )

    return np.einsum("...nj, nij -> ...i", edisps, born_charges)


def mode_oscillator_strengths(mode_eff_chg):
    """Calculate the mode dipole oscillator strengths from a set of
    mode effective charges.

    Parameters
    ----------
    mode_eff_chg : array_like
        Mode effective charges (shape: `(..., 3)`).

    Returns
    -------
    mode_osc_str : numpy.ndarray
        Mode oscillator strengths (shape: `(..., 3, 3)`).
    """

    mode_eff_chg = np.asanyarray(mode_eff_chg)

    if mode_eff_chg.ndim == 0 or mode_eff_chg.shape[-1] != 3:
        raise ValueError(
            "mode_eff_chg must be an array_like with shape (..., 3)."
        )

    return np.einsum("...i, ...j -> ...ij", mode_eff_chg, mode_eff_chg.conj())


def get_gamma_phonons(
    struct,
    freqs,
    evecs,
    eps_inf=None,
    born_charges=None,
    lws=None,
    irreps=None,
    t=None,
):
    r"""Construct and return a `GammaPhonons`, `PolarGammaPhonons` or
    `FullPolarGammaPhonons` for a given set of calculation data.

    Parameters
    ----------
    struct : Structure
        Crystal structure.
    freqs : array_like
        Phonon frequencies in THz (shape: `(M,)`).
    evecs : array_like
        Phonon eigenvectors in mass-weighted units of sqrt(amu)
        (shape: `(M, N, 3)`).
    eps_inf : array_like or None, optional
        High-frequency dielectric constant \eps_inf (shape:
        `(3, 3)`, must be supplied of `born_charges` is set).
    born_charges : array_like or None, optional
        Born effective charge tensors (shape: `(N, 3, 3)`).
    lws : array_like or None, optional
        Phonon linewidths in THz (shape: `(M,)`) or `None`.
    irreps : Irreps or None, optional
        `Irreps` object specifying the point group and assigning
        bands to irrep groups.
    t : float or None, optional
        Calculation temperature (default: `None`).

    Returns
    -------
    gamma_ph : GammaPhonons, PolarGammaPhonons or FullPolarGammaPhonons
        Phonon calculation.

    Notes
    -----
    If `freqs` and `evecs` contain a "complete" phonon calculation with
    M = 3N frequencies/eigenvectors and `born_charges` are available,
    the function returns a `FullPolarGammaPhonons` object. For "partial"
    calculations with M < 3N, it returns a `PolarGammaPhonons` or
    `GammaPhonons` depending on whether `born_charges` are available.
    """

    args = (struct, freqs, evecs)
    kwargs = {"lws": lws, "irreps": irreps, "t": t}

    if born_charges is not None:
        if eps_inf is None:
            raise ValueError(
                "eps_inf must be set if born_charges are supplied."
            )

        if len(freqs) == 3 * struct.num_atoms:
            # All *GammaPhonons constructors check the lengths of freqs
            # and evecs are consistent.

            return FullPolarGammaPhonons(
                *args, eps_inf, born_charges, **kwargs
            )
        else:
            return PolarGammaPhonons(*args, eps_inf, born_charges, **kwargs)

    return GammaPhonons(*args, **kwargs)


def gamma_phonons_from_dict(d):
    """Create a new `GammaPhonons`, `PolarGammaPhonons` or
    `FullPolarGammaPhonons` object from a dictionary generated by the
    `to_dict()` method of any of these objects.

    Parameters
    ----------
    d : dict
        Dictionary generated by `to_dict()`.

    Returns
    -------
    gamma_ph : GammaPhonons, PolarGammaPhonons or FullPolarGammaPhonons
        Phonon calculation.

    See Also
    --------
    get_gamma_phonons :
        Construct and return a `GammaPhonon`, `PolarGammaPhonons` or
        `FullPolarGammaPhonons` object depending on the available data.
    """

    struct = Structure.from_dict(d["structure"])

    freqs = d["frequencies"]
    evecs = dict_to_complex_array(d["eigenvectors"])
    lws = d["linewidths"]

    irreps = None

    if "irreps" in d:
        irreps = Irreps.from_dict(d["irreps"])

    t = d["temperature"]

    eps_inf = d.pop("eps_inf", None)
    born_charges = d.pop("born_charges", None)

    return get_gamma_phonons(
        struct,
        freqs,
        evecs,
        eps_inf=eps_inf,
        born_charges=born_charges,
        lws=lws,
        irreps=irreps,
        t=t,
    )


# ------------------
# GammaPhonons class
# ------------------


class GammaPhonons:
    """Class for storing and working with Gamma-point phonon
    calculations."""

    def __init__(self, struct, freqs, evecs, lws=None, irreps=None, t=None):
        r"""Create a new instance of the `GammaPhonons` class.

        Parameters
        ----------
        struct : Structure
            Crystal structure.
        freqs : array_like
            Phonon frequencies in THz (shape: `(M,)`).
        evecs : array_like
            Phonon eigenvectors in mass-weighted units of sqrt(amu)
            (shape: `(M, N, 3)`).
        lws : array_like or None, optional
            Phonon linewidths in THz (shape: `(M,)`) or `None`.
        irreps : Irreps or None, optional
            `Irreps` object specifying the point group and assigning
            bands to irrep groups.
        t : float or None, optional
            Optional temperature at which the calculation was performed
            (default: `None`).
        acc_mode_tol : str or float, optional
            Tolerance for identifying the acoustic modes (default:
            "auto").
        """

        n_a = struct.num_atoms

        if n_a == 0:
            raise ValueError("struct must contain at least one atom.")

        freqs = np.atleast_1d(np.array(freqs, dtype=np.float64, copy=True))

        evecs = atleast_3d(np.array(evecs, dtype=np.complex128, copy=True))

        if evecs.shape != (len(freqs), n_a, 3):
            raise ValueError(
                "evecs must be an array_like with shape (M, N, 3)."
            )

        if len(evecs) > 3 * n_a:
            raise ValueError("A phonon calculation cannot have M > 3N modes.")

        evecs = np.real_if_close(evecs, tol=ZERO_TOLERANCE)

        if lws is not None:
            # Avoids spurious warning about small complex parts.

            lws = np.atleast_1d(np.array(lws, dtype=np.float64, copy=True))

            if lws.shape != (len(freqs),):
                raise ValueError("lws must be an array_like with shape (M,).")

            # Negative linewidths are unphysical, but small negative
            # values are possible due to numerical noise.

            if (lws < 0.0).any():
                min_lw = lws.min()

                if np.abs(min_lw) > ZERO_TOLERANCE:
                    warnings.warn(
                        "One or more linewidths are negative and will "
                        "be converted to absoute values (min linewidth "
                        f"is {min_lw:.3e}).",
                        UserWarning,
                    )

                lws = np.abs(lws)

        if irreps is not None:
            if len(irreps.irrep_band_indices) != len(evecs):
                raise ValueError(
                    "irreps must assign all M bands to irrep groups."
                )

        if t is not None and t <= 0.0:
            raise ValueError("t must be > 0.")

        freqs.flags.writeable = False
        evecs.flags.writeable = False

        if lws is not None:
            lws.flags.writeable = False

        self._struct = struct
        self._freqs = freqs
        self._evecs = evecs
        self._lws = lws
        self._irreps = irreps
        self._t = t

        self._edisps = None

        self._acc_band_inds = None
        self._ir_band_inds = None
        self._raman_band_inds = None
        self._imag_band_inds = None

    def _lazy_get_band_activities(self):
        """Determine band spectroscopic activities and set the
        _acc_band_inds, _ir_band_inds and _raman_band_inds fields."""

        if self._acc_band_inds is None:
            acc_band_inds = None

            if self.num_modes == 3 * self._struct.num_atoms:
                self._lazy_calc_eigendisplacements()

                # Per-mode mean displacements (M, 1, 3).

                mean_disps = np.mean(self._edisps, axis=1, keepdims=True)

                # Sums of square residual displacements relative to mean
                # (M,).

                sq_res = np.sum((self._edisps - mean_disps) ** 2, axis=(1, 2))

                # Combine the sum of square residuals with the squared
                # frequencies into a consolidated score where a lower
                # value indicates a more "acoustic-like" mode where all
                # the atoms have similar displacements and the absolute
                # frequency is close to zero.

                score = sq_res + self._freqs**2

                # Take the three modes with the lowest value.

                acc_band_inds = np.argsort(score)[:3]
            else:
                # Empty array.

                acc_band_inds = np.array([], dtype=int)

            # If irreps are available, identify the IR- and Raman-active
            # modes. Otherwise, assume all non-acoustic modes are
            # active.

            ir_band_inds, raman_band_inds = None, None

            if self._irreps is not None:
                ir_band_inds, _ = self._irreps.get_subset("ir")

                mask = np.isin(ir_band_inds, acc_band_inds, invert=True)
                ir_band_inds = ir_band_inds[mask]

                raman_band_inds, _ = self._irreps.get_subset("raman")

                mask = np.isin(raman_band_inds, acc_band_inds, invert=True)
                raman_band_inds = raman_band_inds[mask]
            else:
                mask = np.isin(
                    np.arange(0, self.num_modes, dtype=int),
                    acc_band_inds,
                    invert=True,
                )

                ir_band_inds = raman_band_inds = np.flatnonzero(mask)

            # Find the indices of any imaginary modes excluding the
            # acoustic modes if possible.

            imag_band_inds = np.flatnonzero(self._freqs < 0.0)

            mask = np.isin(imag_band_inds, acc_band_inds, invert=True)
            imag_band_inds = imag_band_inds[mask]

            acc_band_inds.flags.writeable = False
            ir_band_inds.flags.writeable = False
            raman_band_inds.flags.writeable = False
            imag_band_inds.flags.writeable = False

            self._acc_band_inds = acc_band_inds
            self._ir_band_inds = ir_band_inds
            self._raman_band_inds = raman_band_inds
            self._imag_band_inds = imag_band_inds

    def _lazy_calc_eigendisplacements(self):
        """Calculate the eigendisplacements and set the _edisps field."""

        if self._edisps is None:
            sqrt_m = np.sqrt(self._struct.atomic_masses)
            edisps = self._evecs / sqrt_m[np.newaxis, :, np.newaxis]

            edisps.flags.writeable = False
            self._edisps = edisps

    @property
    def structure(self):
        """Structure : Crystal structure."""
        return self._struct

    @property
    def frequencies(self):
        """numpy.ndarray : Phonon frequencies (shape: `(M,)`)."""
        return self._freqs

    @property
    def eigenvectors(self):
        """numpy.ndarray : Phonon eigenvectors (shape: `(M, N, 3)`)."""
        return self._evecs

    @property
    def linewidths(self):
        """numpy.ndarray or None : Phonon linewidths (shape: `(M,)`)."""
        return self._lws

    @property
    def irreps(self):
        """Irreps or None : Irreducible representations (irreps)."""
        return self._irreps

    @property
    def temperature(self):
        """float or None : Calculation temperature."""
        return self._t

    @property
    def num_modes(self):
        """int : Number of modes."""
        return len(self._freqs)

    @property
    def eigendisplacements(self):
        """numpy.ndarray : Phonon eigendisplacements (eigenvectors
        divided by sqrt(mass)) (shape: `(M, N, 3)`)."""

        self._lazy_calc_eigendisplacements()
        return self._edisps

    @property
    def acoustic_band_indices(self):
        """numpy.ndarray : Indices of the acoustic bands (shape: `(3,)`,
        empty for "partial" calculations with M < 3N)."""

        self._lazy_get_band_activities()
        return self._acc_band_inds

    @property
    def ir_band_indices(self):
        """numpy.ndarray : Indices of the M' <= M infrared (IR)-active
        bands exclusding the acoustic modes (shape: `(M',)`)."""

        self._lazy_get_band_activities()
        return self._ir_band_inds

    @property
    def raman_band_indices(self):
        """numpy.ndarray : Indices of the M'' <= M Raman-active bands
        excluding the acoustic modes (shape: `(M',)`)."""

        self._lazy_get_band_activities()
        return self._raman_band_inds

    @property
    def imaginary_band_indices(self):
        """numpy.ndarray : Indices of the M''' bands with imaginary
        frequencies (if present) excluding the acoustic modes (shape:
        `(M''',)`)."""

        self._lazy_get_band_activities()
        return self._imag_band_inds

    def eigenvector_projection(self, gamma_ph):
        """Project the phonon eigenevectors onto those of another
        `GammaPhonon` object.

        Parameters
        ----------
        other : GammaPhonons
            Second `GammaPhonons` object.

        Returns
        -------
        proj : numpy.ndarray
            Eigenvector projections (shape: `(M, M'')`).

        Notes
        -----
        `proj` is arranged such that `proj[i, j]` is the projection of
        the eigenvector of the `i`th mode of the calling `GammaPhonons`
        onto the `j`th mode of `gamma_ph`.
        """

        n_a = self._struct.num_atoms

        if gamma_ph.structure.num_atoms != n_a:
            raise ValueError(
                "The structure in gamma_ph must have the same number "
                "of atoms as the calling object."
            )

        evecs_ref = self._evecs.reshape((-1, 3 * n_a))
        evecs_cmp = gamma_ph.eigenvectors.reshape((-1, 3 * n_a))

        return np.inner(evecs_ref.conj(), evecs_cmp)

    def to_dict(self):
        """Return the internal data as a dictionary of native Python
        types for serialisation.

        Returns
        -------
        d : dict
            Dictionary structure containing internal data as native
            Python types.
        """

        lws = self._lws.tolist() if self._lws is not None else None
        irreps = self._irreps.to_dict() if self._irreps is not None else None

        return {
            "structure": self._struct.to_dict(),
            "frequencies": self._freqs.tolist(),
            "eigenvectors": complex_array_to_dict(self._evecs),
            "linewidths": lws,
            "irreps": irreps,
            "temperature": self._t,
        }

    @staticmethod
    def from_dict(d):
        """Create a new `GammaPhonons` instance from a dictionary
        generated by `GammaPhonons.to_dict()`.

        Parameters
        ----------
        d : dict
            Dictionary generated by `to_dict()`.

        Returns
        -------
        gamma_phonons : GammaPhonons
            `GammaPhonons` object constructed from the data in `d`.
        """

        irreps = None

        if d["irreps"] is not None:
            irreps = Irreps.from_dict(d["irreps"])

        return GammaPhonons(
            Structure.from_dict(d["structure"]),
            d["frequencies"],
            dict_to_complex_array(d["eigenvectors"]),
            lws=d["linewidths"],
            irreps=irreps,
            t=d["temperature"],
        )


# ---------------------------
# PolarGammaPhononsBase class
# ---------------------------


class PolarGammaPhonons(GammaPhonons):
    """Class for storing and working with a Gamma-point phonon
    calculation on a polar crystal."""

    def __init__(
        self,
        struct,
        freqs,
        evecs,
        eps_inf,
        born_charges,
        lws=None,
        irreps=None,
        t=None,
    ):
        r"""Create a new instance of the `PolarGammaPhononsBase` class.

        Parameters
        ----------
        struct : Structure
            Crystal structure.
        freqs : array_like
            Phonon frequencies in THz (shape: `(M,)`).
        evecs : array_like
            Phonon eigenvectors in mass-weighted units of sqrt(amu)
            (shape: `(M, N, 3)`).
        eps_inf : array_like
            High-frequency dielectric constant \eps_inf (shape:
            `(3, 3)`).
        born_charges : array_like
            Born effective charge tensors (shape: `(N, 3, 3)`).
        lws : array_like or None, optional
            Phonon linewidths in THz (shape: `(M,)`) or `None`.
        irreps : Irreps or None, optional
            `Irreps` object specifying the point group and assigning
            bands to irrep groups.
        t : float or None, optional
            Calculation temperature (default: `None`).
        """

        super(PolarGammaPhonons, self).__init__(
            struct, freqs, evecs, lws=lws, irreps=irreps, t=t
        )

        eps_inf = np.array(eps_inf, dtype=np.float64, copy=True)

        if eps_inf.shape != (3, 3):
            raise ValueError(
                "eps_inf must be an array_like with shape (3, 3)."
            )

        born_charges = np.array(born_charges, dtype=np.float64, copy=True)

        if born_charges.shape != (self._struct.num_atoms, 3, 3):
            raise ValueError(
                "born_charges must be an array_like with shape (N, 3, 3)."
            )

        eps_inf.flags.writeable = False
        born_charges.flags.writeable = False

        self._eps_inf = eps_inf
        self._born_charges = born_charges

        self._mode_eff_chg = None
        self._mode_osc_str = None

    def _lazy_calc_mode_effective_charges(self):
        """Calculate the mode effective charges on first call to
        `mode_effective_charges`, `pop_frequency` or
        `_lazy_calc_mode_oscillator_strengths`."""

        if self._mode_eff_chg is None:
            self._lazy_calc_eigendisplacements()

            mode_eff_chg = mode_effective_charges(
                self._edisps, self._born_charges
            )

            mode_eff_chg.flags.writeable = False
            self._mode_eff_chg = mode_eff_chg

    def _lazy_calc_mode_oscillator_strengths(self):
        """Calculate the mode oscillator strengths on first call to
        `mode_oscillator_strengths` or `dielectric_function`."""

        if self._mode_osc_str is None:
            self._lazy_calc_mode_effective_charges()

            mode_osc_str = mode_oscillator_strengths(self._mode_eff_chg)

            mode_osc_str.flags.writeable = False
            self._mode_osc_str = mode_osc_str

    @property
    def epsilon_inf(self):
        r"""numpy.ndarray : High-frequency dielectric constant \eps_inf
        in units of relative permittivity (shape: `(3, 3)`)."""
        return self._eps_inf

    @property
    def born_effective_charges(self):
        """numpy.ndarray : Born effective-charge tensors in e (shape:
        `(N, 3, 3)`)."""
        return self._born_charges

    @property
    def mode_effective_charges(self):
        """numpy.ndarray : Mode effective charges in e / sqrt(amu)
        (shape: `(M, 3)`)."""

        self._lazy_calc_mode_effective_charges()
        return self._mode_eff_chg

    @property
    def mode_oscillator_strengths(self):
        """numpy.ndarray : Mode oscillator strengths in e^2 / amu
        (shape: `(M, 3, 3)`)."""

        self._lazy_calc_mode_oscillator_strengths()
        return self._mode_osc_str

    def to_dict(self):
        """Return the internal data as a dictionary of native Python
        types for serialisation.

        Returns
        -------
        d : dict
            Dictionary structure containing internal data as native
            Python types.
        """

        d = super().to_dict()

        d["eps_inf"] = self._eps_inf.tolist()
        d["born_charges"] = self._born_charges.tolist()

        return d

    @staticmethod
    def from_dict(d):
        """Create a new `PolarGammaPhonons` instance from a dictionary
        generated by `PolarGammaPhonons.to_dict()`.

        Parameters
        ----------
        d : dict
            Dictionary generated by `to_dict()`.

        Returns
        -------
        gamma_ph : PolarGammaPhonons
            `PolarGammaPhonons` object constructed from the data in `d`.
        """

        irreps = None

        if d["irreps"] is not None:
            irreps = Irreps.from_dict(d["irreps"])

        return PolarGammaPhonons(
            Structure.from_dict(d["structure"]),
            d["frequencies"],
            dict_to_complex_array(d["eigenvectors"]),
            d["eps_inf"],
            d["born_charges"],
            lws=d["linewidths"],
            irreps=irreps,
            t=d["temperature"],
        )


# -----------------------
# PolarGammaPhonons class
# -----------------------


class FullPolarGammaPhonons(PolarGammaPhonons):
    """Class for storing and working with Gamma-point phonon
    calculations on a polar crystals with a complete set of 3N modes."""

    def __init__(
        self,
        struct,
        freqs,
        evecs,
        eps_inf,
        born_charges,
        lws=None,
        irreps=None,
        t=None,
    ):
        r"""Create a new instance of the `PolarGammaPhonons` class.

        Parameters
        ----------
        struct : Structure
            Crystal structure.
        freqs : array_like
            Phonon frequencies in THz (shape: `(3N,)`).
        evecs : array_like
            Phonon eigenvectors in mass-weighted units of sqrt(amu)
            (shape: `(3N, N, 3)`).
        eps_inf : array_like
            High-frequency dielectric constant \eps_inf (shape:
            `(3, 3)`).
        born_charges : array_like
            Born effective charge tensors (shape: `(N, 3, 3)`).
        lws : array_like or None, optional
            Phonon linewidths in THz (shape: `(3N,)`) or `None`.
        irreps : Irreps or None, optional
            `Irreps` object specifying the point group and assigning
            bands to irrep groups.
        t : float or None, optional
            Optional temperature at which the calculation was performed
            (default: `None`).
        """

        super(FullPolarGammaPhonons, self).__init__(
            struct,
            freqs,
            evecs,
            eps_inf,
            born_charges,
            lws=lws,
            irreps=irreps,
            t=t,
        )

        # Many features of this class rely on a complete set of
        # Gamma-point frequencies and eigenvectors.

        if len(self._freqs) != 3 * self._struct.num_atoms:
            raise ValueError(
                "PolarGammaPhonons objects must have a complete set of "
                "3N phonon frequencies and eigenvectors."
            )

        self._eps_ionic = None

    def _get_dynmat_evecs_evals(self):
        r"""Reconstruct the eigenvalues and eigengectors of the
        dynamical matrix D(q = \Gamma) from the phonon frequencies and
        eigenvectors.

        Returns
        -------
        res : tuple of numpy.ndarray
            Tuple of `(evals, evecs)` (shapes: `(3N,)` and `(3N, 3N)`,
            "VASP" units).
        """

        dim = self.num_modes

        evals = np.copysign((self._freqs / VASP_TO_THZ) ** 2, self._freqs)

        # Reshape eigenvectors to "flat" column format.

        evecs = np.reshape(self._evecs, (dim, dim)).T

        return (evals, evecs)

    def _lazy_calc_epsilon_ionic(self):
        r"""Calculate the ionic contribution to the static dielectric
        constant \eps_ionic on first call to `epsilon_ionic` or
        `epsilon_static`."""

        evals, evecs = self._get_dynmat_evecs_evals()

        # Construct a Phonopy object to obtain the primitive cell and
        # "supercell".

        phonopy = Phonopy(self._struct.to_phonopy_atoms(), np.eye(3, 3))

        # Construct a DynmatToForceConstants object to reverse transform
        # the dynamical matrix to the corresponding force constants.

        d2f = DynmatToForceConstants(phonopy.primitive, phonopy.supercell)

        d2f.create_dynamical_matrices(
            eigenvalues=[evals], eigenvectors=[evecs]
        )

        d2f.run()

        # "Flatten" fc2 to a (3 n_a) x (3 n_a) matrix; the dimensions
        # need to be transposed to obtain the correct block structure
        # after reshaping.

        fc2 = d2f.force_constants

        n_a = self._struct.num_atoms
        fc2 = fc2.transpose(0, 2, 1, 3).reshape((3 * n_a, 3 * n_a))

        # Invert Hessian. Testing suggests fc2 is generally badly
        # conditioned, and np.linalg.pinv() handles this much better
        # than np.linalg.inv().

        inv_h = np.linalg.pinv(fc2)

        # Reshaping Born charges array allows eps_ionic to be written as
        # a matrix triple product z_mat.T @ inv_h @ z_mat.

        z_mat = self._born_charges.reshape(3 * n_a, 3)

        eps_ionic = (
            VASP_DIELECTRIC_TO_RELATIVE_PERMITTIVITY
            * (z_mat.T @ inv_h @ z_mat)
            / self._struct.volume
        )

        eps_ionic.flags.writeable = False
        self._eps_ionic = eps_ionic

    @property
    def epsilon_ionic(self):
        r"""numpy.ndarray : Ionic contribution to dielectric constant
        \eps_ionic in units of relative permittivity (shape: `(3, 3)`).
        """

        self._lazy_calc_epsilon_ionic()
        return self._eps_ionic

    @property
    def epsilon_static(self):
        r"""numpy.ndarray : Static dielectric constant
        \eps_static = \eps_inf + \eps_ionic in units of relative
        permittivity (shape: `(3, 3)`)."""

        self._lazy_calc_epsilon_ionic()
        return self._eps_inf + self._eps_ionic

    def pop_frequency(self, p="high"):
        """Calculate the so-called polar-optic phonon (POP) frequency
        using the vectors from a Lebedev quadrature rule to average
        over the unit sphere.

        Parameters
        ----------
        p : int or str
            Integer precision of the Lebedev quadrature rule or one of
            {"min", "low", "medium", "high", "ultra", "max"} (default:
            "high").

        Returns
        -------
        pop_freq : float
            POP frequency in THz.

        Notes
        -----
        The method implemented here follows that in the AMSET code,[1]_
        and the default `p` is chosen based on this.

        References
        ----------
        .. [1] A. M. Ganose, J. Park, A. Faghaninia, R. Woods-Robinson,
           K. A. Persson and A. Jain, Nature Comm. 12, 2222 (2021), DOI:
           10.1038/s41467-021-22440-5
        """

        self._lazy_calc_mode_effective_charges()

        # Calculate the POP frequency over IR-active bands excluding
        # imaginary modes, if present.

        mask = np.isin(
            self.ir_band_indices, self.imaginary_band_indices, invert=True
        )

        band_inds = self.ir_band_indices[mask]

        if len(band_inds) < len(self.ir_band_indices):
            warnings.warn(
                "Imaginary modes are ignored when calculating the "
                "POP frequency.",
                RuntimeWarning,
            )

        # Vectors and weights for numerical integration.

        q_v, q_w = unit_sphere_lebedev_quad_rule(p, angles=False)

        mode_w = np.zeros(len(band_inds), dtype=np.float64)

        for i, (freq, eff_chg) in enumerate(
            zip(self._freqs[band_inds], self._mode_eff_chg[band_inds])
        ):
            mode_w[i] = sum(
                w * np.abs(np.matmul(v, eff_chg)) for v, w in zip(q_v, q_w)
            ) / np.sqrt(np.abs(freq))

        return ((mode_w * self._freqs[band_inds]) / mode_w.sum()).sum()

    def gamma_phonons_with_nac(self, q):
        r"""Recalculate the phonon frequencies and eigenvectors with a
        non-analytical correction applied to the dynamical matrix.

        Parameters
        ----------
        q : array_like or str
            Incident wavevector for applying NAC.

        Returns
        -------
        gamma_ph : GammaPhonons
            `GammaPhonons` object updated with the corrected frequencies
            and, if applicable, linewidths.

        See Also
        --------
        utility.geometry.parse_direction : Accepted inputs for `q`.

        Notes
        -----
        While the eigenvectors of the original calculation must be real,
        the eigenvectors of the the dynamical matrix with the
        non-analytical correction may be complex. This can be determined
        from the `dtype` of the `eigenvectors` property of the returned
        `GammaPhonons` object.

        If the calling object has linewidths, the linewidths of the
        corrected modes are determined by modifying the eigenvalues of
        the uncorrected dynamical matrix to
        :math:` \lambda_j = \omega_j^2 - i \omega_j \Gamma_j`. The
        modified frequencies and linewidths are then determined from the
        real and imaginary parts of the eigenvalues of the corrected
        matrix.

        The returned `GammaPhonons` does not have irreps, regardless of
        whether the calling object does - this is because LO/TO
        splitting can mix modes with different character, which makes
        the assignment of irreps somewhat questionable. In the case of
        minimal mixing, it may be possible to match modes from the
        original and corrected calculations and assign irreps using
        eigenvector projection.
        """

        q = parse_direction(q)

        # Reconstruct dynamical matrix from frequencies and
        # eigenvectors.

        evals, evecs = self._get_dynmat_evecs_evals()

        if self._lws is not None:
            evals = np.array(evals, dtype=np.complex128)

            # Make sure we copy _freqs and _lws as they are marked
            # readonly.

            freqs = np.array(self._freqs, dtype=np.complex128, copy=True)
            lws = np.array(self._lws, dtype=np.complex128, copy=True)

            freqs[freqs < 0] *= 1.0j
            lws[lws < 0] *= 1.0j

            evals += (1.0j * freqs * lws) / (VASP_TO_THZ**2)

        dynmat = np.matmul(evecs, np.matmul(np.diag(evals), evecs.T))

        # Compute and add NAC correction.

        q_dot_z = (q @ self._born_charges).ravel()
        dynmat_corr = np.outer(q_dot_z, q_dot_z)

        inv_sqrt_m = np.repeat(1.0 / np.sqrt(self._struct.atomic_masses), 3)

        dynmat += dynmat_corr * (
            VASP_NAC_PREFACTOR
            / (self._struct.volume * np.dot(q, np.matmul(self._eps_inf, q)))
            * (inv_sqrt_m[:, np.newaxis] * inv_sqrt_m[np.newaxis, :])
        )

        # Diagonalise corrected D(q) to find frequencies and
        # eigenvectors.

        evals_new, evecs_new = np.linalg.eig(dynmat)

        # With complex eigenvalues the output of numpy.linalg.eig is
        # unlikely to be sorted by the frequencies (= real part of the
        # eigenvalues), and so need resorting to match "typical" phonon
        # calculations.

        inds = np.argsort(evals_new.real)

        evals_new = evals_new[inds]
        evecs_new = evecs_new[:, inds]

        # Prepare a new GammaPhonons object with corrected modes.

        freqs_new = (
            np.copysign(np.sqrt(np.abs(evals_new.real)), evals_new.real)
            * VASP_TO_THZ
        )

        lws_new = None

        if self._lws is not None:
            lws_new = (
                evals_new.imag / np.sqrt(evals_new.real.astype(np.complex128))
            ) * VASP_TO_THZ

            # If we have imaginary modes the linewidths could acquire
            # a significant complex component, which is unphysical and
            # will cause problems in downstream calculations. This needs
            # a looser tolerance than the default ZERO_TOLERANCE to
            # avoid spurious warnings from acoustic modes.

            lws_new = np.real_if_close(lws_new, tol=(100.0 * ZERO_TOLERANCE))

            if np.iscomplexobj(lws_new):
                warnings.warn(
                    "One or more NAC-corrected linewidths have complex "
                    "components that will be discarded (max. abs. "
                    f"imag. = {np.abs(lws_new.imag).max():.3e}.",
                    RuntimeWarning,
                )

                lws_new = lws_new.real

        # evecs_new is in column-major format -> reorder columns,
        # transpose to row-major format, then reshape to (3N, N, 3).

        evecs_new = evecs_new.T.reshape(self._evecs.shape)

        gamma_ph = PolarGammaPhonons(
            self._struct,
            freqs_new,
            evecs_new,
            self._eps_inf,
            self._born_charges,
            lws=lws_new,
            irreps=None,
        )

        return gamma_ph
