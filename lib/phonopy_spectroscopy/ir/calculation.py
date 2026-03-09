# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""High-level `InfraredCalculation` object providing an API for
generating simulated infrared (IR) dielectric functions and optical
spectra."""


# -------
# Imports
# -------


import warnings

import numpy as np

from .dielectric_function import InfraredDielectricFunction

from .spectrum import (
    OpticalEigenmodeSpectrum,
    OpticalSpectrum,
    InputPolarisedOpticalSpectrum,
)

from ..constants import ZERO_TOLERANCE
from ..phonon import PolarGammaPhonons
from ..utility.geometry import rotation_matrix_from_vectors, rotate_tensors
from ..utility.numpy_helper import np_readonly_view
from ..utility.quadrature import unit_sphere_lebedev_quad_rule


# -------------------------
# InfraredCalculation class
# -------------------------


class InfraredCalculation:
    r"""Use a `PolarGammaPhonons` object to generate simulated infrared
    (IR) spectra."""

    def __init__(self, ph_calc):
        r"""Create a new instance of the `InfraredCalculation` class.

        Parameters
        ----------
        gamma_ph : PolarGammaPhonons
            Gamma-point phonon calculation including a high-frequency
            dielectric constant \eps_inf and Born effective charges.
        """

        self._ph_calc = ph_calc

        self._mode_eff_chg = None
        self._mode_osc_str = None

    def _lazy_calc_mode_effective_charges(self):
        """Calculate the mode effective charges on first call to
        `mode_effective_charges`, `pop_frequency` or
        `_lazy_calc_mode_oscillator_strengths`."""

        if self._mode_eff_chg is None:
            mode_eff_chg = np.zeros(
                (self._ph_calc.num_modes, 3), dtype=np.float64
            )

            for i, edisp in enumerate(self._ph_calc.eigendisplacements()):
                temp = np.zeros_like(edisp)

                for j in range(self._ph_calc.structure.num_atoms):
                    temp[j] = np.matmul(
                        self._ph_calc.born_effective_charges[j], edisp[j]
                    )

                mode_eff_chg[i] = temp.sum(axis=0)

            self._mode_eff_chg = mode_eff_chg

    def _lazy_calc_mode_oscillator_strengths(self):
        """Calculate the mode oscillator strengths on first call to
        `mode_oscillator_strengths` or `dielectric_function`."""

        if self._mode_osc_str is None:
            self._lazy_calc_mode_effective_charges()

            mode_osc_str = np.zeros(
                (self._ph_calc.num_modes, 3, 3), dtype=np.float64
            )

            for i in range(self._ph_calc.num_modes):
                mode_osc_str[i, :, :] = np.outer(
                    self._mode_eff_chg[i],
                    self._mode_eff_chg[i],
                )

            self._mode_osc_str = mode_osc_str

    @property
    def structure(self):
        """Structure : Underlying `Structure` object."""
        return self._ph_calc.structure

    @property
    def phonon_calculation(self):
        """PolarGammaPhonons : Underlying `PolarGammaPhonons` object."""
        return self._ph_calc

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

    def dielectric_function(
        self, lw=None, active_only=True, hkl=None, rot=None, **kwargs
    ):
        r"""Simulate the tensor infrared dielectric function.

        Parameters
        ----------
        lw : float or None, optional
            Uniform linewidth or scale factor for calculated linewidths
            (defaults: 0.5 THz uniform linewidth or scale factor of
            1.0, depending on whether calculation has linewidths).
        active_only : bool, optional
            If `True`, and if the underlying Gamma-point phonon
            calculation has irreps, simulate the dielectric function
            using only the infrared-active modes (default: `True`).
        hkl : array_like of int or None, optional
            Reorient crystal so the normal of the surface with the
            specified Miller index is oriented antiparallel to the
            incident direction (default: `None`).
        rot : array_like or None, optional
            Optional rotation to reorient the crystal - if `hkl` is also
            set, the rotation is applied after the `hkl` reorientation.
        **kwargs : any
            Keyword arguments to the `InfraredDielectricFunction`
            constructor.

        Returns
        -------
        eps_ir : InfraredDielectricFunction
            Simulated dielectric function.

        See Also
        --------
        ir.dielectric_function.InfraredDielectricFunction
            Object returned by this function.
        """

        # If requested, and if the Gamma-point phonon calculation has
        # irreps, identify the IR-active modes.

        band_inds = None

        if active_only:
            if self._ph_calc.has_irreps:
                band_inds = self._ph_calc.irreps.get_subset(
                    "ir", reset_inds=False
                ).band_indices_flat()

                # The acoustic modes share irreps with the IR-active
                # modes and need to be explicitly excluded.

                mask = np.isin(
                    band_inds,
                    self._ph_calc.get_acoustic_mode_indices(),
                    invert=True,
                )

                band_inds = band_inds[mask]

        if band_inds is None:
            band_inds = list(range(self._ph_calc.num_modes))

        # Linewidths.

        if lw is not None and lw < ZERO_TOLERANCE:
            raise ValueError("lw cannot be zero or negative.")

        lws = None

        if self._ph_calc.has_linewidths:
            lws = self._ph_calc.linewidths[band_inds]

            if lw is not None:
                lws = lw * lws
        else:
            if lw is None:
                lw = 0.5

            lws = lw * np.ones((len(band_inds),), dtype=np.float64)

        # Oscillator strengths.

        self._lazy_calc_mode_oscillator_strengths()
        osc_strs = self._mode_osc_str

        # \eps_inf.

        eps_inf = self._ph_calc.epsilon_inf

        # Rotate oscillator strengths and \eps_inf if required.

        if hkl is not None or rot is not None:
            r = rotation_matrix_from_vectors(
                self._ph_calc.structure.real_space_normal(hkl, conv=True), "z"
            )

            rot = r if rot is None else np.matmul(rot, r)

        if rot is not None:
            osc_strs = rotate_tensors(osc_strs, rot)
            eps_inf = rotate_tensors(eps_inf, rot)

        # Irreps.

        irreps = None

        if self._ph_calc.has_irreps:
            irreps = self._ph_calc.irreps.get_subset(
                band_inds, reset_inds=True
            )

        return InfraredDielectricFunction(
            self._ph_calc.frequencies[band_inds],
            osc_strs[band_inds],
            lws,
            self._ph_calc.structure.volume(),
            self._ph_calc.epsilon_inf,
            irreps=irreps,
            **kwargs,
        )

    def powder_optical_spectrum_ema(self, t=1.0, diag_eps=True, **kwargs):
        """Simulate the optical spectrum of a powder using the
        effective-medium approximation

        Parameters
        ----------
        t : float, optional
            Sample thickness in mm (default: 1 mm).
        diag_eps : bool, optional
            Diagonalise the dielectric function before taking the scalar
            average (default: `True`).
        **kwargs : any
            Optional arguments to `dielectric_function`.

        Returns
        -------
        sp : OpticalSpectrum
            Simulated spectrum.

        See Also
        --------
        ir.spectrum.OpticalSpectrum
            Object returned by this function.
        """

        eps_ir = self.dielectric_function(**kwargs)

        eps_eff = None

        if diag_eps:
            # Average of the optical eigenmode eigenvalues.

            oe_sp = OpticalEigenmodeSpectrum(
                eps_ir.x, eps_ir.epsilon, x_units=eps_ir.x_units
            )

            eps_eff = np.mean(oe_sp.mode_eigenvalues, axis=-1)
        else:
            # Average of the diagonal elements.

            eps_eff = np.trace(eps_ir.epsilon, axis1=1, axis2=2) / 3.0

        oe_sp = OpticalEigenmodeSpectrum(
            eps_ir.x, eps_eff, x_units=eps_ir.x_units, t=t
        )

        return OpticalSpectrum(oe_sp)

    def powder_optical_spectrum_average_eigenmodes(self, t=1.0, **kwargs):
        """Simulate the optical spectrum of a powder by averaging the
        properties of the optical eigenmodes.

        Parameters
        ----------
        t : float, optional
            Sample thickness in mm (default: 1 mm).
        **kwargs : any
            Optional arguments to `dielectric_function`.

        Returns
        -------
        sp : OpticalSpectrum
            Simulated spectrum.

        See Also
        --------
        ir.spectrum.OpticalSpectrum
            Object returned by this function.
        """

        oe_sp = OpticalEigenmodeSpectrum.from_infrared_dielectric_function(
            self.dielectric_function(**kwargs), t=t
        )

        return OpticalSpectrum(oe_sp)

    def single_crystal_unpolarised_optical_spectrum(
        self, hkl, t=1.0, rot=None, **kwargs
    ):
        """Simulate the optical spectrum of a single-crystal in a
        collinear geometry along the z-axis with unpolarised incident
        light.

        Parameters
        ----------
        hkl : array_like of int
            Miller index of the surface to orient antiparallel to the
            incident direction.
        t : float, optional
            Sample thickness in mm (default: 1 mm).
        rot : array_like or None, optional
            Optional rotation to reorient the crystal after the `hkl`
            rotation.
        **kwargs : any
            Optional arguments to `dielectric_function`.

        Returns
        -------
        sp : OpticalSpectrum
            Simulated spectrum.

        See Also
        --------
        ir.spectrum.OpticalSpectrum
            Object returned by this function.
        """

        eps_ir = self.dielectric_function(hkl=hkl, rot=rot, **kwargs)

        oe_sp = OpticalEigenmodeSpectrum(
            eps_ir.x, eps_ir.epsilon[:, :2, :2], x_units=eps_ir.x_units, t=t
        )

        return OpticalSpectrum(oe_sp, eps_eels=eps_ir.epsilon[:, 2, 2])

    def single_crystal_input_polarised_optical_spectrum(
        self, hkl, i_pol, t=1.0, rot=None, **kwargs
    ):
        """Simulate the optical spectrum of a single-crystal in a
        collinear geometry along the z-axis with polarised incident
        light.

        Parameters
        ----------
        hkl : array_like of int
            Miller index of the surface to orient antiparallel to the
            incident direction.
        i_pol : Polarisation
            Polarisation of incident light.
        t : float, optional
            Sample thickness in mm (default: 1 mm).
        rot : array_like or None, optional
            Optional rotation to reorient the crystal after the `hkl`
            rotation.
        **kwargs : any
            Optional arguments to `dielectric_function`.

        Returns
        -------
        sp : InputPolarisedOpticalSpectrum
            Simulated spectrum.

        Notes
        -----
        This routine assumes a collinear measurement geometry with the
        incident and collected light along the z-axis. The polarisation
        must therefore be defined in the x/y plane.

        See Also
        --------
        ir.spectrum.InputPolarisedOpticalSpectrum
            Object returned by this function.
        """

        eps_ir = self.dielectric_function(hkl=hkl, rot=rot, **kwargs)

        oe_sp = OpticalEigenmodeSpectrum(
            eps_ir.x, eps_ir.epsilon[:, :2, :2], x_units=eps_ir.x_units, t=t
        )

        return InputPolarisedOpticalSpectrum(
            oe_sp, i_pol, eps_eels=eps_ir.epsilon[:, 2, 2]
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

        if active_only and self._ph_calc.has_irreps:
            band_inds = self._ph_calc.irreps.get_subset(
                "ir"
            ).band_indices_flat()
        else:
            band_inds = list(range(self._ph_calc.num_modes))

        mode_w = np.zeros((len(band_inds),), dtype=np.float64)

        q_v, q_w = unit_sphere_lebedev_quad_rule(lebedev_prec, ret="vectors")

        freqs = self._ph_calc.frequencies

        (inds,) = np.where(freqs < -1.0 * ZERO_TOLERANCE)
        acc_inds = self._ph_calc.get_acoustic_mode_indices()

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

        return {"phonon_calculation": self._ph_calc.to_dict()}

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
            PolarGammaPhonons.from_dict(d["phonon_calculation"])
        )
