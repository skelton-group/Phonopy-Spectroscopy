# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Class for storing and working with phonon calculations."""


# -------
# Imports
# -------


import warnings

import numpy as np

from scipy.optimize import linear_sum_assignment

from .constants import (
    VASP_TO_THZ,
    VASP_NAC_PREFACTOR,
    VASP_DIELECTRIC_TO_RELATIVE_PERMITTIVITY,
    ZERO_TOLERANCE,
)

from .irreps import Irreps
from .structure import Structure

from .utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
)

_PHONOPY_AVAILABLE = False

try:
    from phonopy import Phonopy
    from phonopy.harmonic.dynmat_to_fc import DynmatToForceConstants

    _PHONOPY_AVAILABLE = True
except ImportError:
    warnings.warn(
        "Imports from Phonopy failed - some functions require Phonopy "
        "and will raise exceptions if it is not installed.",
        RuntimeWarning,
    )


# ------------------
# GammaPhonons class
# ------------------


class GammaPhonons:
    """Class for storing and working with a Gamma-point phonon
    calculation."""

    def __init__(self, struct, freqs, evecs, lws=None, irreps=None):
        r"""Create a new instance of the `GammaPhonons` class.

        Parameters
        ----------
        struct : Structure
            Crystal structure.
        freqs : array_like
            Phonon frequencies in THz (shape: `(3N,)`).
        evecs : array_like
            Phonon eigenvectors in mass-weighted units of sqrt(amu)
            (shape: `(3N, N, 3)`).
        lws : array_like or None, optional
            Phonon linewidths in THz (shape: `(3N,)`) or `None`.
        irreps : Irreps or None, optional
            `Irreps` object specifying the point group and assigning
            bands to irrep groups.
        """

        n_a = struct.num_atoms

        if n_a == 0:
            raise ValueError(
                'struct cannot be "empty" and must contain at least one atom.'
            )

        freqs = np_asarray_copy(freqs, dtype=np.float64)

        if not np_check_shape(freqs, (3 * n_a,)):
            raise ValueError("freqs must be an array_like with shape (3N,).")

        evecs = np_asarray_copy(evecs, dtype=np.float64)

        if not np_check_shape(evecs, (3 * n_a, n_a, 3)):
            raise ValueError(
                "evecs must be an array_like with shape (3N, N, 3)."
            )

        # Eigenvectors are in general complex, but Gamma-point
        # eigenvectors must be real.

        if np.iscomplex(evecs).any():
            if (np.abs(evecs.imag) > ZERO_TOLERANCE).any():
                raise ValueError("Gamma-point eigenvectors should be real.")

            evecs = evecs.real

        if lws is not None:
            lws = np_asarray_copy(lws, dtype=np.float64)

            if not np_check_shape(lws, (3 * n_a,)):
                raise ValueError(
                    "If supplied, lws must be an array_like with shape (3N,)."
                )

            if (lws < 0.0).any():
                raise ValueError("Linewidths cannot be negative.")

        if irreps is not None:
            ir_band_inds_flat = irreps.band_indices_flat()

            if len(ir_band_inds_flat) != 3 * n_a:
                raise ValueError(
                    "If supplied, irreps must assign all bands to "
                    "irrep groups."
                )

            if ir_band_inds_flat.max() >= (3 * n_a):
                raise RuntimeError(
                    "One or more band indices in irreps are not "
                    "compatible with the number of modes in the "
                    "phonon calculation."
                )

        self._struct = struct
        self._freqs = freqs
        self._evecs = evecs
        self._lws = lws
        self._irreps = irreps

    @property
    def structure(self):
        """Structure : Crystal structure."""
        return self._struct

    @property
    def frequencies(self):
        """numpy.ndarray : Phonon frequencies (shape: `(3N,)`)."""
        return np_readonly_view(self._freqs)

    @property
    def eigenvectors(self):
        """numpy.ndarray : Phonon eigenvectors (shape: `(3N, N, 3)`)."""
        return np_readonly_view(self._evecs)

    @property
    def linewidths(self):
        """numpy.ndarray or None : Phonon linewidths (shape: `(3N,)`)."""
        return np_readonly_view(self._lws) if self._lws is not None else None

    @property
    def num_modes(self):
        """int : Number of modes."""
        return len(self._freqs)

    @property
    def has_linewidths(self):
        """bool : `True` if linewidths are available, otherwise `False`."""
        return self._lws is not None

    @property
    def has_irreps(self):
        """bool : `True` if irredicuble representations (irreps) are
        available, otherwise `False`.
        """
        return self._irreps is not None

    @property
    def irreps(self):
        """Irreps or None : `Irreps` object with the point group and
        irrep symbols and indices of band groups."""
        return self._irreps

    def get_acoustic_mode_indices(self):
        """Return the band indices of the acoustic modes.

        Returns
        -------
        band_inds : numpy.ndarray
            Band indices of the acoustic modes.
        """

        # If we have irreps, select irrep group(s) for which the
        # average frequency is closest to f = 0 until we have chosen
        # sufficient groups to cover three modes.

        if self._irreps is not None:
            ir_ave_freqs = [
                np.mean(self._freqs[band_inds])
                for band_inds in self._irreps.irrep_band_indices
            ]

            subset_band_inds = []

            for _, band_inds in sorted(
                zip(ir_ave_freqs, self._irreps.irrep_band_indices)
            ):
                subset_band_inds.extend(band_inds)

                if len(subset_band_inds) == 3:
                    return np.array(subset_band_inds, dtype=int)

            raise RuntimeError(
                "Unable to select a set of acoustic modes spanning "
                "complete irrep groups. This may indicate an issue "
                "with the phonon calculation."
            )

        # If not, find the three modes with frequencies closest to
        # f = 0.

        return np.argsort(np.abs(self._freqs))[:3]

    def eigendisplacements(self):
        r"""Return the phonon eigendisplacements (eigenvectors divided
        by sqrt(mass)).

        Returns
        -------
        edisps : numpy.ndarray
            Eigendisplacements (shape: `(3N, N, 3)`).
        """

        sqrt_m = np.sqrt(self._struct.atomic_masses)
        return self._evecs / sqrt_m[np.newaxis, :, np.newaxis]

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
            "eigenvectors": self._evecs.tolist(),
            "linewidths": lws,
            "irreps": irreps,
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
            d["eigenvectors"],
            lws=d["linewidths"],
            irreps=irreps,
        )


# -----------------------
# PolarGammaPhonons class
# -----------------------


class PolarGammaPhonons(GammaPhonons):
    """Class for storing and working with a Gamma-point phonon
    calculation, including non-analytical corrections in polar
    compounds."""

    def __init__(
        self,
        struct,
        freqs,
        evecs,
        eps_inf,
        born_charges,
        lws=None,
        irreps=None,
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
        """

        super(PolarGammaPhonons, self).__init__(
            struct, freqs, evecs, lws=lws, irreps=irreps
        )

        eps_inf = np_asarray_copy(eps_inf, dtype=np.float64)

        if not np_check_shape(eps_inf, (3, 3)):
            raise ValueError(
                "eps_inf must be an array_like with shape (3, 3)."
            )

        born_charges = np_asarray_copy(born_charges, dtype=np.float64)

        if not np_check_shape(born_charges, (self._struct.num_atoms, 3, 3)):
            raise ValueError(
                "born_charges must be an array_like with shape (N, 3, 3)."
            )

        self._eps_inf = eps_inf
        self._born_charges = born_charges

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

        if not _PHONOPY_AVAILABLE:
            raise RuntimeError(
                "Calculating epsilon_ionic and epsilon_static requires "
                "the phonopy.Phonopy and "
                "phonopy.harmonic.dynmat_to_fc.DynmatToForceConstants "
                "classes."
            )

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

        eps_ionic = np.zeros((3, 3), dtype=np.float64)

        n_dof, _ = inv_h.shape

        for i in range(n_dof):
            i_at, i_dir = i // 3, i % 3

            for j in range(n_dof):
                j_at, j_dir = j // 3, j % 3

                for a in range(3):
                    for b in range(3):
                        eps_ionic[a, b] += (
                            self._born_charges[i_at][i_dir, a]
                            * inv_h[i, j]
                            * self._born_charges[j_at][j_dir, b]
                        )

        self._eps_ionic = (
            VASP_DIELECTRIC_TO_RELATIVE_PERMITTIVITY
            * eps_ionic
            / self._struct.volume()
        )

    @property
    def epsilon_inf(self):
        r"""numpy.ndarray : High-frequency dielectric constant \eps_inf
        (shape: `(3, 3)`)."""
        return np_readonly_view(self._eps_inf)

    @property
    def epsilon_ionic(self):
        r"""numpy.ndarray : Ionic contribution to dielectric constant
        \eps_ionic (shape: `(3, 3)`)."""

        self._lazy_calc_epsilon_ionic()
        return np_readonly_view(self._eps_ionic)

    @property
    def epsilon_static(self):
        r"""numpy.ndarray : Static dielectric constant
        \eps_static = \eps_inf + \eps_ionic (shape: `(3, 3)`)."""

        self._lazy_calc_epsilon_ionic()
        return self._eps_inf + self._eps_ionic

    @property
    def born_effective_charges(self):
        """numpy.ndarray : Born effective-charge tensors (shape:
        `(N, 3, 3)`)."""
        return np_readonly_view(self._born_charges)

    def gamma_phonons_with_nac(self, q):
        """Recalculate the phonon frequencies and eigenvectors with a
        non-analytical correction applied to the dynamical matrix.

        Parameters
        ----------
        q : array_like or str
            Incident wavevector for applying NAC.

        Returns
        -------
        res : tuple of (GammaPhonons, numpy.ndarray)
            Tuple of `(gamma_ph, evec_proj)` containing a `GammaPhonons`
            object, updated with the corrected frequencies and
            eigenvectors, and a matrix of the absolute projections of
            the corrected eigenvectors onto the uncorrected ones (shape:
            `(3N, 3N)`).

        Notes
        -----
        The phonon modes from the corrected dynamical matrix are
        reordered to match the original modes using eigenvector
        projection. This is necessary in order to match the corrected
        modes to the linewidths and irreps of the original calculation,
        but it does mean the frequencies may not be in order.

        If the calling object has linewidths, these are passed through
        to the `GammaPhonons` object in the return tuple.

        If the calling object has irreps, these are "flattened" so that
        degenerate modes, which may be split by the correction, retain
        their irrep symbol but are no longer grouped. (This is to avoid
        them being incorrectly averaged by degeneracy-handling routines
        elsewhere.)
        """

        # Reconstruct dynamical matrix from frequencies and
        # eigenvectors.

        evals, evecs = self._get_dynmat_evecs_evals()
        dynmat = np.matmul(evecs, np.matmul(np.diag(evals), evecs.T))

        # Compute and add NAC correction.

        dynmat_corr = np.zeros_like(dynmat, dtype=np.float64)

        inv_sqrt_m = np.repeat(1.0 / np.sqrt(self._struct.atomic_masses), 3)

        q_dot_z = np.matmul(q, self._born_charges)

        dim = self.num_modes

        for i in range(dim):
            i_at = i // 3
            i_dir = i % 3

            for j in range(dim):
                j_at = j // 3
                j_dir = j % 3

                dynmat_corr[i, j] = q_dot_z[i_at, i_dir] * q_dot_z[j_at][j_dir]

        dynmat += dynmat_corr * (
            VASP_NAC_PREFACTOR
            / (self._struct.volume() * np.dot(q, np.matmul(self._eps_inf, q)))
            * (inv_sqrt_m[:, np.newaxis] * inv_sqrt_m[np.newaxis, :])
        )

        # Diagonalise corrected D(q) to find frequencies and
        # eigenvectors.

        evals_new, evecs_new = np.linalg.eigh(dynmat)

        # Compute eigenvector projections.

        evec_proj = np.abs(np.dot(evecs.T, evecs_new))

        # Match original and corrected modes using projections.

        _, col_inds = linear_sum_assignment(1.0 - evec_proj)

        # Prepare a new GammaPhonons object with corrected modes.

        evals_new = evals_new[col_inds]

        freqs_new = np.copysign(
            np.sqrt(np.abs(evals_new)) * VASP_TO_THZ, evals_new
        )

        # evecs_new is in column-major format -> reorder columns,
        # transpose to row-major format, then reshape to (3N, N, 3).

        evecs_new = evecs_new[:, col_inds].T.reshape(self._evecs.shape)

        irreps_new = None

        if self._irreps is not None:
            irreps_new = Irreps(
                self._irreps.point_group,
                self._irreps.symbols_flat(),
                [[i] for i in range(dim)],
            )

        gamma_ph = GammaPhonons(
            self._struct,
            freqs_new,
            evecs_new,
            lws=self._lws,
            irreps=irreps_new,
        )

        # Reformat eigenvector projections so that rows -> corrected
        # modes after reordering and columns -> original modes.

        evec_proj = evec_proj[:, col_inds].T

        return (gamma_ph, evec_proj)

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
            "eigenvectors": self._evecs.tolist(),
            "eps_inf": self._eps_inf.tolist(),
            "born_charges": self._born_charges.tolist(),
            "linewidths": lws,
            "irreps": irreps,
        }

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
            d["eigenvectors"],
            d["eps_inf"],
            d["born_charges"],
            lws=d["linewidths"],
            irreps=irreps,
        )
