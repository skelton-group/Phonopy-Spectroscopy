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

from .utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
    np_expand_dims,
)

from .utility.quadrature import unit_sphere_lebedev_quad_rule

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


# ---------
# Functions
# ---------


def mode_effective_charges(edisps, born_charges):
    """Calculate the mode effective charges from a set of
    "eigendisplacements" (eigenvectors after division by sqrt(mass)) and
    Born effective-charge tensors.

    Parameters
    ----------
    edisps : array_like
        Eigendisplacements (shape: `(N, 3)` or `(M, N, 3)`).
    born_charges : array_like
        Born effective-charge tensors (shape: `(N, 3, 3)`).

    Returns
    -------
    mode_eff_chg : numpy.ndarray
        Mode effective charges (shape: `(3,)` or `(M, 3)`).
    """

    born_charges = np.asarray(born_charges)

    if not np_check_shape(born_charges, (None, 3, 3)):
        raise ValueError(
            "born_charges must be an array_like with shape (N, 3, 3)."
        )

    n, _, _ = born_charges.shape

    edisps, n_dim_add = np_expand_dims(edisps, (None, n, 3))

    mode_eff_chg = np.einsum("mnj, nij -> mi", edisps, born_charges)

    return mode_eff_chg if n_dim_add == 0 else mode_eff_chg[0]


def mode_oscillator_strengths(mode_eff_chg):
    """Calculate the mode dipole oscillator strengths from a set of
    mode effective charges.

    Parameters
    ----------
    mode_eff_chg : array_like
        Mode effective charges (shape: `(3,)` or `(M, 3)`).

    Returns
    -------
    mode_osc_str : numpy.ndarray
        Mode oscillator strengths (shape: `(3, 3)` or `(M, 3, 3)`).
    """

    mode_eff_chg, n_dim_add = np_expand_dims(mode_eff_chg, (None, 3))

    mode_osc_str = np.einsum(
        "mi, mj -> mij", mode_eff_chg, mode_eff_chg.conj()
    )

    return mode_osc_str if n_dim_add == 0 else mode_osc_str


# ------------------
# GammaPhonons class
# ------------------


class GammaPhonons:
    """Class for storing and working with a Gamma-point phonon
    calculation."""

    def __init__(self, struct, freqs, evecs, lws=None, irreps=None, t=None):
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
        t : float or None, optional
            Optional temperature at which the calculation was performed
            (default: `None`).
        """

        n_a = struct.num_atoms

        if n_a == 0:
            raise ValueError(
                'struct cannot be "empty" and must contain at least one atom.'
            )

        freqs = np_asarray_copy(freqs, dtype=np.float64)

        if not np_check_shape(freqs, (3 * n_a,)):
            raise ValueError("freqs must be an array_like with shape (3N,).")

        evecs = np_asarray_copy(evecs, dtype=np.complex128)

        if not np_check_shape(evecs, (3 * n_a, n_a, 3)):
            raise ValueError(
                "evecs must be an array_like with shape (3N, N, 3)."
            )

        # If the eigenvectors are real, discard the complex part.

        if not np.iscomplex(evecs).any():
            evecs = evecs.real

        if lws is not None:
            lws = np_asarray_copy(lws, dtype=np.float64)

            if not np_check_shape(lws, (3 * n_a,)):
                raise ValueError(
                    "If supplied, lws must be an array_like with shape (3N,)."
                )

            # Negative linewidths are unphysical, but small negative
            # values are possible due to numerical noise.

            if (lws < 0.0).any():
                min_lw = lws.min()

                if min_lw < -1.0 * ZERO_TOLERANCE:
                    warnings.warn(
                        "One or more linewidths are negative and will "
                        "be converted to absoute values (min linewidth "
                        "is {0:.3e}).".format(min_lw),
                        UserWarning,
                    )

                lws = np.abs(lws)

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

        if t is not None:
            if t <= 0.0:
                raise ValueError("If specified, t must be greater than zero.")

        self._struct = struct
        self._freqs = freqs
        self._evecs = evecs
        self._lws = lws
        self._irreps = irreps
        self._t = t

        self._edisps = None
        self._acc_mode_inds = None

    def _lazy_calc_eigendisplacements(self):
        """Calculate the eigendisplacements (eigenvectors divided by
        sqrt(mass)) and set the _edisps field."""

        if self._edisps is None:
            sqrt_m = np.sqrt(self._struct.atomic_masses)
            self._edisps = self._evecs / sqrt_m[np.newaxis, :, np.newaxis]

    def _lazy_find_acoustic_mode_inds(self):
        """Find the indices of the acoustic modes and set the
        _acc_mode_inds field."""

        if self._acc_mode_inds is None:
            self._lazy_calc_eigendisplacements()

            # Per-mode mean displacements (3N, 1, 3).

            mean_disps = np.mean(self._edisps, axis=1, keepdims=True)

            # Sums of square residual displacements relative to mean
            # (3N,).

            sq_res = np.sum((self._edisps - mean_disps) ** 2, axis=(1, 2))

            # Norms of the residuals.

            res_norms = np.linalg.norm(mean_disps.reshape(-1, 3), axis=1)

            # Combine all three into a consolidated "score".

            score = (
                sq_res + self._freqs**2 + 1.0 / (res_norms + ZERO_TOLERANCE)
            )

            acc_mode_inds = np.argsort(score)[:3]

            if self._irreps is not None:
                # If we have irreps, the acoustic modes must span
                # "complete" irrep groups.
                try:
                    subset = self._irreps.get_subset(acc_mode_inds)
                except ValueError:
                    raise RuntimeError(
                        "The acoustic modes do not span complete irrep "
                        "groups. This likely indicates an issue with "
                        "the irreps, or may be a bug."
                    )

            self._acc_mode_inds = acc_mode_inds

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
    def irreps(self):
        """Irreps or None : `Irreps` object with the point group and
        irrep symbols and indices of band groups."""
        return self._irreps

    @property
    def temperature(self):
        """float or None : Temperature at which the calculation was
        performed."""
        return self._t

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
    def eigendisplacements(self):
        """numpy.ndarray : phonon eigendisplacements (eigenvectors
        divided by sqrt(mass)) (shape: `(3N, N, 3)`)."""

        self._lazy_calc_eigendisplacements()
        return self._edisps

    @property
    def acoustic_mode_indices(self):
        """numpy.ndarray : Band indices of the acoustic modes."""

        self._lazy_find_acoustic_mode_inds()
        return np_readonly_view(self._acc_mode_inds)

    @property
    def imaginary_mode_indices(self):
        """numpy.ndarray : Band indices of imaginary modes, if present."""

        (inds,) = np.where(self._freqs < 0.0)

        # Acoustic modes often have small imaginary frequencies.

        self._lazy_find_acoustic_mode_inds()
        mask = np.isin(inds, self._acc_mode_inds, invert=True)

        return inds[mask]

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
            Eigenvector projections (shape: `(3N, 3N)`).

        Notes
        -----
        `proj` is arranged such that `proj[i, j]` is the projection of
        the eigenvector of the `i`th mode of the calling `GammaPhonons`
        onto the `j`th mode of `gamma_ph`.
        """

        n = self.num_modes

        if gamma_ph.num_modes != n:
            raise ValueError(
                "gamma_ph must have the same number of modes as the "
                "calling object."
            )

        evecs_ref = self._evecs.reshape((n, n))
        evecs_cmp = gamma_ph.eigenvectors.reshape(gamma_ph, (n, n))

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

        evecs = {"real": self._evecs.real.tolist()}

        if np.iscomplexobj(self._evecs):
            evecs["imag"] = self._evecs.imag.tolist()

        return {
            "structure": self._struct.to_dict(),
            "frequencies": self._freqs.tolist(),
            "eigenvectors": evecs,
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

        evecs = np.asarray(d["eigenvectors"]["real"], dtype=np.float64)

        if "imag" in d["eigenvectors"]:
            evecs_imag = np.asarray(
                d["eigenvectors"]["imag"], dtype=np.float64
            )

            if not np.equal(np.shape(evecs), np.shape(evecs_imag)).all():
                raise ValueError(
                    '"eigenvectors" key contains "real" and "imag" '
                    "keys with different array shapes."
                )

            evecs = evecs + 1.0j * evecs_imag

        return GammaPhonons(
            Structure.from_dict(d["structure"]),
            d["frequencies"],
            evecs,
            lws=d["linewidths"],
            irreps=irreps,
            t=d["temperature"],
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

        super(PolarGammaPhonons, self).__init__(
            struct, freqs, evecs, lws=lws, irreps=irreps, t=t
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

        self._mode_eff_chg = None
        self._mode_osc_str = None

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

    def _lazy_calc_mode_effective_charges(self):
        """Calculate the mode effective charges on first call to
        `mode_effective_charges`, `pop_frequency` or
        `_lazy_calc_mode_oscillator_strengths`."""

        self._lazy_calc_eigendisplacements()

        if self._mode_eff_chg is None:
            self._mode_eff_chg = mode_effective_charges(
                self._edisps, self._born_charges
            )

    def _lazy_calc_mode_oscillator_strengths(self):
        """Calculate the mode oscillator strengths on first call to
        `mode_oscillator_strengths` or `dielectric_function`."""

        if self._mode_osc_str is None:
            self._lazy_calc_mode_effective_charges()
            self._mode_osc_str = mode_oscillator_strengths(self._mode_eff_chg)

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
        in units of relative permittivity (shape: `(3, 3)`)."""
        return np_readonly_view(self._eps_inf)

    @property
    def epsilon_ionic(self):
        r"""numpy.ndarray : Ionic contribution to dielectric constant
        \eps_ionic in units of relative permittivity (shape: `(3, 3)`).
        """

        self._lazy_calc_epsilon_ionic()
        return np_readonly_view(self._eps_ionic)

    @property
    def epsilon_static(self):
        r"""numpy.ndarray : Static dielectric constant
        \eps_static = \eps_inf + \eps_ionic in units of relative
        permittivity (shape: `(3, 3)`)."""

        self._lazy_calc_epsilon_ionic()
        return self._eps_inf + self._eps_ionic

    @property
    def born_effective_charges(self):
        """numpy.ndarray : Born effective-charge tensors in e (shape:
        `(N, 3, 3)`)."""
        return np_readonly_view(self._born_charges)

    @property
    def mode_effective_charges(self):
        """numpy.ndarray : Mode effective charges in e / sqrt(amu)
        (shape: `(3N, 3)`)."""

        self._lazy_calc_mode_effective_charges()
        return np_readonly_view(self._mode_eff_chg)

    @property
    def mode_oscillator_strengths(self):
        """numpy.ndarray : Mode oscillator strengths in e^2 / amu
        (shape: `(3N, 3, 3)`)."""

        self._lazy_calc_mode_oscillator_strengths()
        return np_readonly_view(self._mode_osc_str)

    def pop_frequency(self, lebedev_prec=53, active_only=True):
        """Calculate the so-called polar-optic phonon (POP) frequency
        using the vectors from a Lebedev quadrature rule to average
        over the unit sphere.

        Parameters
        ----------
        lebedev_prec : int
            Precision of the Lebedev quadrature rule (default: 53).
        active_only : bool
            If `True`, and if the underlying Gamma-point phonon
            calculation has irreps, calculate the POP frequency using
            only the infrared-active modes (default: `True`).

        Returns
        -------
        pop_freq : float
            POP frequency in THz.

        Notes
        -----
        The method implemented here follows that in the AMSET code,[1]_
        and the default `lebedev_prec` is chosen based on this.

        References
        ----------
        .. [1] A. M. Ganose, J. Park, A. Faghaninia, R. Woods-Robinson,
           K. A. Persson and A. Jain, Nature Comm. 12, 2222 (2021), DOI:
           10.1038/s41467-021-22440-5
        """

        self._lazy_calc_mode_effective_charges()

        band_inds = None

        if active_only:
            if self._irreps is not None:
                band_inds = self._irreps.get_subset("ir").band_indices_flat()
        else:
            band_inds = np.arange(0, len(self._freqs), dtype=int)

        # Exclude acoustic and imaginary modes.

        self._lazy_find_acoustic_mode_inds()

        mask = np.isin(band_inds, self._acc_mode_inds, invert=True)

        imag_mode_inds = self.imaginary_mode_indices

        if len(imag_mode_inds) > 0:
            warnings.warn(
                "Imaginary modes are ignored when calculating the "
                "POP frequency.",
                RuntimeWarning,
            )

            mask = np.logical_or(
                mask, np.isin(band_inds, imag_mode_inds, invert=True)
            )

        band_inds = band_inds[mask]

        # Vectors and weights for numerical integration.

        q_v, q_w = unit_sphere_lebedev_quad_rule(lebedev_prec, ret="vectors")

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

        # Reconstruct dynamical matrix from frequencies and
        # eigenvectors.

        evals, evecs = self._get_dynmat_evecs_evals()

        if self._lws is not None:
            evals = np.asarray(evals, dtype=np.complex128)

            freqs = np.asarray(self._freqs, dtype=np.complex128)
            lws = np.asarray(self._lws, dtype=np.complex128)

            freqs[freqs < 0] *= 1.0j
            lws[lws < 0] *= 1.0j

            evals += (1.0j * freqs * lws) / (VASP_TO_THZ**2)

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
                evals_new.imag
                / np.sqrt(evals_new.real.astype(np.complex128))
                * VASP_TO_THZ
            )

        # evecs_new is in column-major format -> reorder columns,
        # transpose to row-major format, then reshape to (3N, N, 3).

        evecs_new = evecs_new.T.reshape(self._evecs.shape)

        gamma_ph = GammaPhonons(
            self._struct,
            freqs_new,
            evecs_new,
            lws=lws_new.real,
            irreps=None,
        )

        # Check none of the linewidths of the none-acoustic modes had a
        # significant imaginary part.

        mask = np.isin(
            np.arange(0, gamma_ph.num_modes, dtype=int),
            gamma_ph.acoustic_mode_indices,
            invert=True,
        )

        max_imag = np.abs(lws_new.imag[mask]).max()

        if max_imag > ZERO_TOLERANCE:
            warnings.warn(
                "Maximum imaginary part of linewidths is {0:.3e} "
                "> ZERO_TOLERANCE = {1:.3e}."
                "".format(max_imag, ZERO_TOLERANCE),
                RuntimeWarning,
            )

        return gamma_ph

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
            "temperature": self._t,
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
            t=d["temperature"],
        )
