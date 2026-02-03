# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""High-level `InfraredCalculation` object providing an API for using
Gamma-point phonon and Born effective-charge calculations to generate
simulated infrared (IR) spectra and related quantities."""


# -------
# Imports
# -------


import warnings

import numpy as np

from .dielectric_func import (
    TensorInfraredDielectricFunction,
    ScalarInfraredDielectricFunction,
)

from .oscillator_strength import (
    calculate_single_crystal_oscillator_strengths,
    calculate_powder_oscillator_strengths,
)

from ..constants import (
    ZERO_TOLERANCE,
    DIELECTRIC_TO_RELATIVE_PERMITTIVITY,
)

from ..phonon import GammaPhonons

from ..utility.geometry import rotation_matrix_from_vectors

from ..utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
)

from ..utility.quadrature import unit_sphere_lebedev_quad_rule


# -------------------------
# InfraredCalculation class
# -------------------------


class InfraredCalculation:
    r"""Combine a `GammaPhonons` object, a set of Born effective
    charges, and, optionally, a high-frequency dielectric constant
    \eps_inf, to generate simulated infrared (IR) spectra and related
    quantities."""

    def __init__(self, gamma_ph, born_charges, eps_inf=None):
        r"""Create a new instance of the `InfraredCalculation` class.

        Parameters
        ----------
        gamma_ph : GammaPhonons
            Gamma-point phonon calculation.
        born_charges : array_like
            Born effective charge tensors (shape: `(N, 3, 3)`).
        eps_inf : array_like or None, optional
            High-frequency dielectric constant \eps_inf (shape:
            `(3, 3)`, default: `None`).
        """

        born_charges = np_asarray_copy(born_charges, dtype=np.float64)

        if not np_check_shape(
            born_charges, (gamma_ph.structure.num_atoms, 3, 3)
        ):
            raise ValueError(
                "born_charges must be an array_like with shape (N, 3, 3)."
            )

        if eps_inf is not None:
            eps_inf = np_asarray_copy(eps_inf, dtype=np.float64)

            if not np_check_shape(eps_inf, (3, 3)):
                raise ValueError(
                    "If supplied, eps_inf must be an array_like with "
                    "shape (3, 3)."
                )

        self._gamma_ph = gamma_ph
        self._born_charges = born_charges

        self._eps_inf = eps_inf
        self._eps_ionic = None

        self._mode_eff_chg = None
        self._mode_osc_str = None

    def _lazy_calc_epsilon_ionic(self):
        r"""Calculate the ionic contribution to the static dielectric
        constant \eps_ionic on first call to `epsilon_ionic`,
        `epsilon_static` or `dielectric_function`."""

        # Invert Hessian. Testing suggests h is generally badly
        # conditioned, and np.linalg.pinv() handles this much better
        # than np.linalg.inv().

        inv_h = np.linalg.pinv(self._gamma_ph.hessian())

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
            DIELECTRIC_TO_RELATIVE_PERMITTIVITY
            * eps_ionic
            / self._gamma_ph.structure.volume()
        )

    def _lazy_calc_mode_effective_charges(self):
        """Calculate the mode effective charges on first call to
        `mode_effective_charges`, `pop_frequency` or
        `_lazy_calc_mode_oscillator_strengths`."""

        if self._mode_eff_chg is None:
            mode_eff_chg = np.zeros(
                (self._gamma_ph.num_modes, 3), dtype=np.float64
            )

            for i, edisp in enumerate(self._gamma_ph.eigendisplacements()):
                temp = np.zeros_like(edisp)

                for j in range(self._gamma_ph.structure.num_atoms):
                    temp[j] = np.matmul(self._born_charges[j], edisp[j])

                mode_eff_chg[i] = temp.sum(axis=0)

            self._mode_eff_chg = mode_eff_chg

    def _lazy_calc_mode_oscillator_strengths(self):
        """Calculate the mode oscillator strengths on first call to
        `mode_oscillator_strengths` or `dielectric_function`."""

        if self._mode_osc_str is None:
            self._lazy_calc_mode_effective_charges()

            mode_osc_str = np.zeros(
                (self._gamma_ph.num_modes, 3, 3), dtype=np.float64
            )

            for i in range(self._gamma_ph.num_modes):
                mode_osc_str[i, :, :] = np.outer(
                    self._mode_eff_chg[i],
                    self._mode_eff_chg[i],
                )

            self._mode_osc_str = mode_osc_str

    @property
    def structure(self):
        """Structure : Underlying `Structure` object."""
        return self._gamma_ph.structure

    @property
    def gamma_phonons(self):
        """GammaPhonons : Underlying `GammaPhonons` object."""
        return self._gamma_ph

    @property
    def born_effective_charges(self):
        """numpy.ndarray : Born effective-charge tensors (shape:
        `(N, 3, 3)`)."""
        return np_readonly_view(self._born_charges)

    @property
    def epsilon_inf(self):
        r"""numpy.ndarray or None : High-frequency dielectric constant
        \eps_inf (shape: `(3, 3)`)."""

        if self._eps_inf is not None:
            return np_readonly_view(self._eps_inf)

        return None

    @property
    def epsilon_ionic(self):
        r"""numpy.ndarray : Ionic contribution to dielectric constant
        \eps_ionic (shape: `(3, 3)`)."""

        self._lazy_calc_epsilon_ionic()
        return np_readonly_view(self._eps_ionic)

    @property
    def epsilon_static(self):
        r"""numpy.ndarray or None : Static dielectric constant
        \eps_static = \eps_inf + \eps_ionic (shape: `(3, 3)`)."""

        if self._eps_inf is not None:
            self._lazy_calc_epsilon_ionic()
            return self._eps_inf + self._eps_ionic

        return None

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

    def _get_calc_params(self, lw, add_eps_inf, active_only):
        r"""Determine parameters for infrared dielectric function
        calculations.

        Parameters
        ----------
        lw : float or None
            Uniform linewidth or scale factor for calculated linewidths,
            depending on whether calculation has linewidths.
        add_eps_inf : bool, optional
            If `True`, add the high-frequency dielectric constant
            \eps_inf to the dielectric function.
        active_only : bool
            If `True`, and if the underlying Gamma-point phonon
            calculations has irreps, return parameters for the subset of
            IR-active modes.

        Returns
        -------
        params : dict
            Dictionary of calculation parameters with the following: `{
                "linewidths": numpy.ndarray,
                "epsilon_inf" : numpy.ndarray
                }`
        """

        params = {}

        # If requested, and if the Gamma-point phonon calculation has
        # irreps, identify the IR-active modes.

        band_inds = None

        if active_only:
            if self._gamma_ph.has_irreps:
                band_inds = self._gamma_ph.irreps.get_subset(
                    "ir", reset_inds=False
                ).band_indices_flat()

                # The acoustic modes share irreps with the IR-active
                # modes and need to be explicitly excluded.

                mask = np.isin(
                    band_inds,
                    self._gamma_ph.get_acoustic_mode_indices(),
                    invert=True,
                )

                band_inds = band_inds[mask]

        if band_inds is None:
            band_inds = list(range(self._gamma_ph.num_modes))

        # Frequencies.

        params["frequencies"] = self._gamma_ph.frequencies[band_inds]

        # Oscillator strengths.

        self._lazy_calc_mode_oscillator_strengths()

        params["oscillator_strengths"] = self._mode_osc_str[band_inds]

        # Linewidths.

        if lw is None:
            lw = 1.0 if self._gamma_ph.linewidths is not None else 0.5
        else:
            if lw < ZERO_TOLERANCE:
                raise ValueError("lw cannot be zero or negative.")

        params["linewidths"] = (
            lw * self._gamma_ph.linewidths[band_inds]
            if self._gamma_ph.has_linewidths
            else lw * np.ones((len(band_inds),), dtype=np.float64)
        )

        # Irreps.

        irreps = None

        if self._gamma_ph.has_irreps:
            irreps = self._gamma_ph.irreps.get_subset(
                band_inds, reset_inds=True
            )

        params["irreps"] = irreps

        # High-frequency dielectric constant.

        eps_inf = None

        if add_eps_inf:
            eps_inf = self._eps_inf

        params["epsilon_inf"] = eps_inf

        return params

    def tensor_dielectric_function(
        self, lw=None, add_eps_inf=True, active_only=True, **kwargs
    ):
        r"""Simulate the tensor infrared dielectric function.

        Parameters
        ----------
        lw : float or None, optional
            Uniform linewidth or scale factor for calculated linewidths
            (defaults: 0.5 THz uniform linewidth or scale factor of
            1.0, depending on whether calculation has linewidths).
        add_eps_inf : bool, optional
            If `True`, add the high-frequency dielectric constant
            \eps_inf, if available, to the dielectric function (default:
            `True`).
        active_only : bool, optional
            If `True`, and if the underlying Gamma-point phonon
            calculation has irreps, simulate the dielectric function
            using only the infrared-active modes (default: `True`).

        **kwargs : any
            Keyword arguments to the `TensorInfraredDielectricFunction`
            constructor.

        Returns
        -------
        dielectric_func : TensorInfraredDielectricFunction
            Simulated dielectric function.

        See Also
        --------
        ir.dielectric_func.TensorInfraredDielectricFunction
            Object returned by this function.
        """

        # Determine linewidths.

        params = self._get_calc_params(lw, add_eps_inf, active_only)

        # Calculate mode oscillator strengths.

        self._lazy_calc_mode_oscillator_strengths()

        return TensorInfraredDielectricFunction(
            params["frequencies"],
            params["oscillator_strengths"],
            params["linewidths"],
            self._gamma_ph.structure.volume(),
            irreps=params["irreps"],
            eps_inf=params["epsilon_inf"],
            **kwargs,
        )

    def scalar_single_crystal_dielectric_function(
        self,
        hkl,
        geom,
        i_pol,
        s_pol,
        rot=None,
        lw=None,
        add_eps_inf=True,
        active_only=True,
        **kwargs
    ):
        r"""Simulate the scalar infrared dielectric function of a single
        crystal.

        Parameters
        ----------
        hkl : tuple of int
            Miller index of crystal surface to orient along incident
            direction.
        geom : Geometry
            Measurement geometry.
        i_pol, s_pol : Polarisation
            Polarisation of incident and scattered light.
        rot : array_like or None, optional
            Rotation matrix to realign crystal after `hkl` rotation
            (default: `None`).
        lw : float or None, optional
            Uniform linewidth or scale factor for calculated linewidths
            (defaults: 0.5 THz uniform linewidth or scale factor of
            1.0, depending on whether calculation has linewidths).
        add_eps_inf : bool, optional
            If `True`, add the high-frequency dielectric constant
            \eps_inf, if available, to the dielectric function (default:
            `True`).
        active_only : bool, optional
            If `True`, and if the underlying Gamma-point phonon
            calculation has irreps, simulate the dielectric function
            using only the infrared-active modes (default: `True`).
        **kwargs : any
            Keyword arguments to the `ScalarInfraredDielectricFunction`
            constructor.

        Returns
        -------
        dielectric_func : ScalarInfraredDielectricFunction
            Simulated dielectric function.

        See Also
        --------
        ir.dielectric_func.ScalarInfraredDielectricFunction
            Object returned by this function.
        """

        params = self._get_calc_params(lw, add_eps_inf, active_only)

        r = rotation_matrix_from_vectors(
            self._gamma_ph.structure.real_space_normal(hkl, conv=True),
            -1.0 * geom.incident_direction,
        )

        if rot is not None:
            if not np_check_shape(rot, (3, 3)):
                raise ValueError(
                    "If supplied, rot must be an array_like with shape "
                    "(3, 3)."
                )

            r = np.matmul(rot, r)

        self._lazy_calc_mode_oscillator_strengths()

        osc_str = calculate_single_crystal_oscillator_strengths(
            params["oscillator_strengths"], geom, i_pol, s_pol, rot=r
        )

        # If adding \eps_inf to the simulated dielectric function, we
        # need to compute the scalar value in the same way as the scalar
        # oscillator strengths.

        eps_inf = params["epsilon_inf"]

        if eps_inf is not None:
            eps_inf = calculate_single_crystal_oscillator_strengths(
                eps_inf, geom, i_pol, s_pol, rot=rot
            )

        return ScalarInfraredDielectricFunction(
            params["frequencies"],
            osc_str,
            params["linewidths"],
            self._gamma_ph.structure.volume(),
            irreps=params["irreps"],
            eps_inf=eps_inf,
            **kwargs,
        )

    def scalar_powder_dielectric_function(
        self,
        geom,
        i_pol,
        s_pol,
        po_hkl=None,
        po_eta=0.0,
        lw=None,
        add_eps_inf=True,
        active_only=True,
        method="best",
        lc_prec=5,
        **kwargs
    ):
        r"""Simulate the scalar dielectric function of a powder.

        Parameters
        ----------
        geom : Geometry
            Measurement geometry.
        i_pol, s_pol : str, Polarisation or list of Polarisation
            Polarisation(s) of incident and scattered light. The
            scattered polarisation may also be specified by one of
            {"parallel", "cross", "sum"}.
        po_hkl : tuple of int or None, optional
            Miller index of preferred orientation (default: None).
        po_eta : float, optional
            Fraction of crystallites with the preferred orientation
            (default: 0.0).
        lw : float or None, optional
            Uniform linewidth or scale factor for calculated linewidths
            (defaults: 0.5 THz uniform linewidth or scale factor of
            1.0, depending on whether calculation has linewidths).
        add_eps_inf : bool, optional
            If `True`, add the high-frequency dielectric constant
            \eps_inf, if available, to the dielectric function (default:
            `True`).
        active_only : bool, optional
            If `True`, and if the underlying Gamma-point phonon
            calculation has irreps, simulate the dielectric function
            using only the infrared-active modes (default: `True`).
        method : {"nquad", "lebedev+circle", "best"}, optional
            Method for powder averaging (default: `"best"`).
        lc_prec : int, optional
            Precision of the Lebedev + circle numerical quadrature
            scheme (default: 5).
        **kwargs : any
            Keyword arguments to the `ScalarInfraredDielectricFunction`
            constructor.

        Returns
        -------
        dielectric_func : ScalarInfraredDielectricFunction
            Simulated dielectric function.

        See Also
        --------
        ir.dielectric_func.ScalarInfraredDielectricFunction
            Object returned by this function.
        """

        params = self._get_calc_params(lw, add_eps_inf, active_only)

        po_surf_norm = None

        if po_hkl is not None:
            po_surf_norm = self._gamma_ph.structure.real_space_normal(
                po_hkl, conv=True
            )

        self._lazy_calc_mode_oscillator_strengths()

        osc_str = calculate_powder_oscillator_strengths(
            params["oscillator_strengths"],
            geom,
            i_pol,
            s_pol,
            po_eta=po_eta,
            po_surf_norm=po_surf_norm,
            method=method,
            lc_prec=lc_prec,
        )

        eps_inf = params["epsilon_inf"]

        if eps_inf is not None:
            eps_inf = calculate_powder_oscillator_strengths(
                eps_inf,
                geom,
                i_pol,
                s_pol,
                po_eta=po_eta,
                po_surf_norm=po_surf_norm,
                method=method,
                lc_prec=lc_prec,
            )

        return ScalarInfraredDielectricFunction(
            params["frequencies"],
            osc_str,
            params["linewidths"],
            self._gamma_ph.structure.volume(),
            irreps=params["irreps"],
            eps_inf=eps_inf,
            **kwargs,
        )

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

        if active_only and self._gamma_ph.has_irreps:
            band_inds = self._gamma_ph.irreps.get_subset(
                "ir"
            ).band_indices_flat()
        else:
            band_inds = list(range(self._gamma_ph.num_modes))

        mode_w = np.zeros((len(band_inds),), dtype=np.float64)

        q_v, q_w = unit_sphere_lebedev_quad_rule(lebedev_prec, ret="vectors")

        freqs = self._gamma_ph.frequencies

        (inds,) = np.where(freqs < -1.0 * ZERO_TOLERANCE)
        acc_inds = self._gamma_ph.get_acoustic_mode_indices()

        for idx in inds:
            if idx in band_inds and idx not in acc_inds:
                warnings.warn(
                    "Imaginary modes are ignored when calculating the "
                    "POP frequency.",
                    RuntimeWarning,
                )

            break

        for i, (freq, eff_chg) in enumerate(
            zip(freqs[band_inds], self._mode_eff_chg[band_inds])
        ):
            # Ignore modes with zero frequency and imaginary modes.

            if freq > ZERO_TOLERANCE:
                mode_w[i] = sum(
                    w * np.abs(np.matmul(v, eff_chg)) for v, w in zip(q_v, q_w)
                ) / np.sqrt(np.abs(freq))

        return ((mode_w * freqs[band_inds]) / mode_w.sum()).sum()

    def to_dict(self):
        """Return the internal data as a dictionary of native Python
        types for serialisation.

        Returns
        -------
        d : dict
            Dictionary structure containing internal data as native
            Python types.
        """

        eps_inf = self._eps_inf.tolist() if self._eps_inf is not None else None

        return {
            "gamma_phonons": self._gamma_ph.to_dict(),
            "born_charges": self._born_charges.tolist(),
            "epsilon_inf": eps_inf,
        }

    @staticmethod
    def from_dict(d):
        """Create a new `InfraredCalculation` instance from a dictionary
        generated by `InfraredCalculation.to_dict()`.

        Parameters
        ----------
        d : dict
            Dictionary generated by `to_dict()`.

        Returns
        -------
        calc : InfraredCalculation
            `InfraredCalculation` object constructed from the data in
            `d`.
        """

        return InfraredCalculation(
            GammaPhonons.from_dict(d["gamma_phonons"]),
            d["born_charges"],
            d["epsilon_inf"],
        )
