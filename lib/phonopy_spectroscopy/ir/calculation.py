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


import numpy as np

from .dielectric_function import InfraredDielectricFunction

from .spectrum_core import (
    OpticalEigenmodeSpectrum,
    EigenmodeAverageOpticalSpectrum,
    InputPolarisedOpticalSpectrum,
)

from ..constants import ZERO_TOLERANCE

from ..phonon import (
    mode_effective_charges,
    mode_oscillator_strengths,
    PolarGammaPhonons,
)

from ..utility.geometry import (
    parse_direction,
    rotation_matrix_from_vectors,
    rotate_tensors,
)


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

    @property
    def structure(self):
        """Structure : Underlying `Structure` object."""
        return self._ph_calc.structure

    @property
    def phonon_calculation(self):
        """PolarGammaPhonons : Underlying `PolarGammaPhonons` object."""
        return self._ph_calc

    def dielectric_function(
        self,
        lw=None,
        hkl=None,
        rot=None,
        q_nac=None,
        active_only=True,
        **kwargs
    ):
        r"""Simulate the tensor infrared dielectric function.

        Parameters
        ----------
        lw : float or None, optional
            Uniform linewidth or scale factor for calculated linewidths
            (defaults: 0.5 THz uniform linewidth or scale factor of
            1.0, depending on whether calculation has linewidths).
        hkl : array_like of int or None, optional
            Reorient crystal so the normal of the surface with the
            specified Miller index is oriented antiparallel to the
            incident direction (default: `None`).
        rot : array_like or None, optional
            Optional rotation to reorient the crystal - if `hkl` is also
            set, the rotation is applied after the `hkl` reorientation.
        q_nac : array_like or None, optional
            Optional "approach direction" for applying a non-analytical
            correction (NAC) to the phonon frequencies and eigenvectors.
        active_only : bool, optional
            If `True`, and if the underlying Gamma-point phonon
            calculation has irreps, calculate the POP frequency using
            only the infrared-active modes (default: `True`).
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

        Notes
        -----
        `q_nac` is specified in real-space coordinates. If `hkl` and/or
        `rot` are set, `q_nac` is specified in the rotated frame.
        """

        # Determine rotation matrix.

        if hkl is not None or rot is not None:
            r = rotation_matrix_from_vectors(
                self._ph_calc.structure.real_space_normal(hkl, conv=True), "z"
            )

            rot = r if rot is None else np.matmul(rot, r)

        # If q_nac is set, generate a new phonon calculation with LO/TO
        # splitting.

        ph_calc = self._ph_calc

        if q_nac is not None:
            # If we are applying a rotation matrix, we need to apply
            # the reverse rotation to q to convert it to the crystal
            # reference frame.

            q_nac = parse_direction(q_nac)

            if rot is not None:
                q_nac = np.matmul(rot.T, q_nac)

            ph_calc = ph_calc.gamma_phonons_with_nac(q_nac)

        # Band indices to include in calculation.

        band_inds = np.arange(0, ph_calc.num_modes, dtype=int)

        if active_only:
            # Exclude acoustic modes.

            mask = np.isin(
                band_inds, ph_calc.acoustic_mode_indices, invert=True
            )

            if ph_calc.has_irreps:
                active_inds = ph_calc.irreps.get_subset(
                    "ir"
                ).band_indices_flat()

                mask = np.logical_and(mask, np.isin(band_inds, active_inds))

            band_inds = band_inds[mask]

        # Oscillator strengths.

        osc_strs = mode_oscillator_strengths(
            mode_effective_charges(
                ph_calc.eigendisplacements[band_inds],
                self._ph_calc.born_effective_charges,
            )
        )

        # Linewidths.

        if lw is not None and lw < ZERO_TOLERANCE:
            raise ValueError("lw cannot be zero or negative.")

        lws = None

        if ph_calc.has_linewidths:
            lws = ph_calc.linewidths[band_inds]

            if lw is not None:
                lws = lw * lws
        else:
            if lw is None:
                lw = 0.5

            lws = lw * np.ones((len(band_inds),), dtype=np.float64)

        # \eps_inf.

        eps_inf = self._ph_calc.epsilon_inf

        # Irreps.

        irreps = None

        if ph_calc.has_irreps:
            irreps = ph_calc.irreps.get_subset(band_inds, reset_inds=True)

        # Rotate oscillator strengths and \eps_inf if required.

        if rot is not None:
            osc_strs = rotate_tensors(osc_strs, rot)
            eps_inf = rotate_tensors(eps_inf, rot)

        return InfraredDielectricFunction(
            ph_calc.frequencies[band_inds],
            osc_strs,
            lws,
            ph_calc.structure.volume(),
            eps_inf,
            irreps=irreps,
            **kwargs,
        )

    def powder_optical_spectrum_ema(
        self, t=1.0, p_vol_frac=1.0, p_binder_eps=1.0, p_den=1.0, **kwargs
    ):
        """Simulate the optical spectrum of a powder using the
        effective-medium approximation

        Parameters
        ----------
        t : float, optional
            Sample thickness in mm (default: 1 mm).
        p_vol_frac : float, optional
            Volume fraction of material in a pellet (default: 1.0).
        p_binder_eps : float or tuple of numpy.ndarray, optional
            Dielectric constant or tuple of `(x, eps_x)` specifying the
            frequency-dependent dielectric function of the pellet
            "binder" material (default: 1.0 = vaccum ~ air).
        p_den : float, optional
            Density of the pellet (default: 1.0).
        **kwargs : any
            Optional arguments to `dielectric_function`.

        Returns
        -------
        sp : EigenmodeAverageOpticalSpectrum
            Simulated spectrum.

        See Also
        --------
        ir.spectrum.EigenmodeAverageOpticalSpectrum
            Object returned by this function.
        """

        eps_ir = self.dielectric_function(**kwargs)

        eps_eff = None

        # Average of the optical eigenmode eigenvalues.

        oe_sp = OpticalEigenmodeSpectrum(
            eps_ir.x,
            eps_ir.epsilon,
            x_units=eps_ir.x_units,
            branch_tracking=False,
        )

        eps_eff = np.mean(oe_sp.mode_eigenvalues, axis=-1)

        oe_sp = OpticalEigenmodeSpectrum(
            eps_ir.x,
            eps_eff,
            x_units=eps_ir.x_units,
            t=t,
            p_vol_frac=p_vol_frac,
            p_binder_eps=p_binder_eps,
            p_den=p_den,
        )

        return EigenmodeAverageOpticalSpectrum(oe_sp)

    def powder_optical_spectrum_eigenmode_average(
        self, t=1.0, p_vol_frac=1.0, p_binder_eps=1.0, p_den=1.0, **kwargs
    ):
        """Simulate the optical spectrum of a powder by averaging the
        properties of the optical eigenmodes.

        Parameters
        ----------
        t : float, optional
            Sample thickness in mm (default: 1 mm).
        p_vol_frac : float, optional
            Volume fraction of material in a pellet (default: 1.0).
        p_binder_eps : float or tuple of numpy.ndarray, optional
            Dielectric constant or tuple of `(x, eps_x)` specifying the
            frequency-dependent dielectric function of the pellet
            "binder" material (default: 1.0 = vaccum ~ air).
        p_den : float, optional
            Density of the pellet (default: 1.0).
        **kwargs : any
            Optional arguments to `dielectric_function`.

        Returns
        -------
        sp : EigenmodeAverageOpticalSpectrum
            Simulated spectrum.

        See Also
        --------
        ir.spectrum.EigenmodeAverageOpticalSpectrum
            Object returned by this function.
        """

        oe_sp = OpticalEigenmodeSpectrum.from_infrared_dielectric_function(
            self.dielectric_function(**kwargs),
            t=t,
            branch_tracking=True,
            p_vol_frac=p_vol_frac,
            p_binder_eps=p_binder_eps,
            p_den=p_den,
        )

        return EigenmodeAverageOpticalSpectrum(oe_sp)

    def single_crystal_unpolarised_optical_spectrum(
        self, hkl, t=1.0, rot=None, nac=False, **kwargs
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
        nac : bool, optional
            If `True`, apply a non-analytical correction (LO/TO
            splitting) to the phonon frequencies and eigenvectors
            (default: `False`).
        **kwargs : any
            Optional arguments to `dielectric_function`.

        Returns
        -------
        sp : EigenmodeAverageOpticalSpectrum
            Simulated spectrum.

        See Also
        --------
        ir.spectrum.EigenmodeAverageOpticalSpectrum
            Object returned by this function.
        """

        eps_ir = self.dielectric_function(
            hkl=hkl, rot=rot, q_nac=("z" if nac else None), **kwargs
        )

        oe_sp = OpticalEigenmodeSpectrum(
            eps_ir.x, eps_ir.epsilon[:, :2, :2], x_units=eps_ir.x_units, t=t
        )

        return EigenmodeAverageOpticalSpectrum(
            oe_sp, eps_eels=eps_ir.epsilon[:, 2, 2]
        )

    def single_crystal_input_polarised_optical_spectrum(
        self, hkl, i_pol, t=1.0, rot=None, nac=False, **kwargs
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
        nac : bool, optional
            If `True`, apply a non-analytical correction (LO/TO
            splitting) to the phonon frequencies and eigenvectors
            (default: `False`).
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

        eps_ir = self.dielectric_function(
            hkl=hkl, rot=rot, q_nac=("z" if nac else None), **kwargs
        )

        oe_sp = OpticalEigenmodeSpectrum(
            eps_ir.x, eps_ir.epsilon[:, :2, :2], x_units=eps_ir.x_units, t=t
        )

        return InputPolarisedOpticalSpectrum(
            oe_sp, i_pol, eps_eels=eps_ir.epsilon[:, 2, 2]
        )

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
