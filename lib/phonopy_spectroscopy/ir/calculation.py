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

from ..optics.response_function import DielectricFunction

from ..optics.optical_eigenmodes import (
    OpticalEigenmodes,
    SingleCrystalOpticalEigenmodes,
)

from ..gamma_phonons import gamma_phonons_from_dict, FullPolarGammaPhonons

from ..utility.geometry import (
    parse_direction,
    rotation_matrix_from_vectors,
    rotate_tensors,
)

# -------------------------
# InfraredCalculation class
# -------------------------


class InfraredCalculation:
    r"""Use a `PolarGammaPhonons` or `FullPolarGammaPhonons` object to
    generate simulated infrared (IR) spectra."""

    def __init__(self, ph_calc):
        """Create a new instance of the `InfraredCalculation` class.

        Parameters
        ----------
        gamma_ph : PolarGammaPhonons or FullPolarGammaPhonons
            Phonon calculation.
        """

        default_nac = True

        if not isinstance(ph_calc, FullPolarGammaPhonons):
            warnings.warn(
                "Non-analytical corrections (NAC) require the "
                "calculator is initialised with a "
                "FullPolarGammaPhonons object. NAC will be disabled "
                "by default and nac=True will raise a RuntimeError.",
                RuntimeWarning,
            )

            default_nac = False

        self._ph_calc = ph_calc
        self._default_nac = default_nac

    @property
    def structure(self):
        """Structure : Crystal structure."""
        return self._ph_calc.structure

    @property
    def phonon_calculation(self):
        """PolarGammaPhonons or FullPolarGammaPhonons : Phonon
        calculation."""
        return self._ph_calc

    @property
    def default_nac(self):
        """bool : `True` if non-analytical corrections (NAC) will be
        applied when `nac="auto"` is passed to the calculation routines,
        otherwise `False`.
        """
        return self._default_nac

    def dielectric_function(
        self,
        lw=0.5,
        hkl=None,
        rot=None,
        nac="auto",
        **kwargs,
    ):
        r"""Calculate the tensor infrared dielectric function.

        Parameters
        ----------
        lw : float or None, optional
            Uniform linewidth (default: 0.5 THz, overridden by per-mode
            calculated linewidths if set).
        hkl : array_like of int or None, optional
            Reorient crystal so the normal of the surface with the
            specified Miller index is oriented antiparallel to the
            z axis (default: `None`).
        rot : array_like or None, optional
            Optional rotation to reorient the crystal - if `hkl` is also
            set, the rotation is applied after the `hkl` reorientation.
        nac : str or bool, optional
            Apply a non-analytical correction (NAC) to the phonon
            frequencies, oscillator strengths and calculated linewidths
            (default: "auto", sets to `auto_nac`).
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
        The calculation assumes a standard experimental geometry where
        light is incident along the z axis.
        """

        # Determine rotation matrix.

        if rot is not None:
            rot = np.asanyarray(rot)

            if rot.shape != (3, 3):
                raise ValueError(
                    "rot must be an array_like with shape (3, 3)."
                )

        if hkl is not None:
            hkl = np.asanyarray(hkl, dtype=int)

            if hkl.ndim != 1 or len(hkl) != 3:
                raise ValueError("hkl must specify three Miller indices.")

            r = rotation_matrix_from_vectors(
                self._ph_calc.structure.real_space_normal(hkl, conv=True), "-z"
            )

            rot = r if rot is None else rot @ r

        if rot is None:
            # Default to identity matrix.
            rot = np.identity(3, dtype=np.float64)

        # Default NAC.

        if isinstance(nac, (str, np.str_)):
            nac_str = str(nac).lower()

            if nac_str == "auto":
                nac = self._default_nac
            else:
                raise ValueError(f'Unknown string specifier nac="{nac}".')

        ph_calc = self._ph_calc

        # If nac is set, generate a new phonon calculation with LO/TO
        # splitting.

        if nac:
            if not isinstance(ph_calc, FullPolarGammaPhonons):
                raise RuntimeError(
                    "Non-analytical corrections (NAC) can only be "
                    "performed if the calculation was initialised with "
                    "a FullPolarGammaPhonons object."
                )

            # If we are applying a rotation matrix, we need to convert
            # the incident direction to the crystal reference frame.

            q_nac = parse_direction("z")

            if rot is not None:
                q_nac = rot.T @ q_nac

            ph_calc = ph_calc.gamma_phonons_with_nac(q_nac)

        # Band indices to include in calculation.

        band_inds = ph_calc.ir_band_indices

        # Frequencies.

        freqs = ph_calc.frequencies[band_inds]

        # Oscillator strengths and \eps_inf.

        osc_strs = rotate_tensors(
            ph_calc.mode_oscillator_strengths[band_inds], rot
        )

        eps_inf = rotate_tensors(ph_calc.epsilon_inf, rot)

        # Linewidths.

        lws = None

        if ph_calc.linewidths is not None:
            lws = ph_calc.linewidths[band_inds]
        else:
            if lw <= 0.0:
                raise ValueError("lw must be > 0.")

            lws = np.full((len(band_inds),), lw)

        # Irreps.

        irreps = None

        if ph_calc.irreps is not None:
            _, irreps = ph_calc.irreps.get_subset(band_inds)

        return InfraredDielectricFunction(
            freqs,
            osc_strs,
            lws,
            ph_calc.structure.volume,
            eps_inf,
            irreps=irreps,
            **kwargs,
        )

    def optical_eigenmodes(self, nac=False, **kwargs):
        """Compute the optical eigenmodes of the bulk dielectric
        function.

        Parameters
        ----------
        nac : bool
            Apply a non-analytical correction when generating the
            dielectric function (default: `False`).
        **kwargs : any
            Other keyword arguments to `dielectric_function`.

        Returns
        -------
        oe : OpticalEigenmodes
            Optical eigenmodes.

        Notes
        -----
        As this method is intended to generate the optical eigenmodes of
        the bulk dielectric function, non-analytical corrections (NAC)
        are disabled by default.

        See Also
        --------
        dielectric_function :
            Simulate the tensor infrared dielectric function.
        optics.optical_eigenmodes.OpticalEigenmodes :
            Object returned by this function.
        """

        return OpticalEigenmodes(
            self.dielectric_function(nac=nac, **kwargs), branch_tracking=True
        )

    def incoherent_single_crystal_on_axis_spectrum(
        self,
        hkl,
        rot=None,
        i_pol=None,
        d_pol=None,
        t=1.0,
        n_f=1.0,
        n_b=1.0,
        oe_branch_tracking=True,
        **kwargs,
    ):
        """Simulate the optical spectrum of a single crystal by
        projecting an incident and detected polarisation onto the
        optical eigenmodes of the 2x2 block of the infrared dielectric
        function accessible in a standard collinear measurement
        geometry, with the incident/detected light along +/- z.

        Parameters
        ----------
        hkl : array_like of int
            Surface to orient antiparallel to the incident light
            direction.
        rot : array_like or None, optional
            Optional rotation to reorient the crystal after the `hkl`
            reorientation.
        i_pol, d_pol : SPPolarisation or None, optional
            Polarisations of incident and detected light (default:
            unpolarised incident/detected light).
        t : float, optonal
            Thickness in mm (default: 1 mm).
        n_f, n_b : float, optional
            Refractive indices of the front (indicent) and back (exit)
            media (default: 1.0 = vacuum ~ air).
        oe_branch_tracking : bool
            Perform branch tracking in the intermediate
            `SingleCrystalOpticalEigenmodes` object used to generate the
            spectrum (default: `True`).
        **kwargs : any
            Optional arguments to `dielectric_function`.

        Returns
        -------
        sp : EigenmodeAverageOpticalSpectrum
            Simulated spectrum.

        See Also
        --------
        optics.optical_eigenmodes.SingleCrystalOpticalEigenmodes :
            Intermediate object used to generate the spectrum from the
            dielectric function.
        optics.optical_eigenmodes.EigenmodeAverageOpticalSpectrum :
            Object returned by this function.
        """

        eps_ir_bulk = self.dielectric_function(hkl=hkl, rot=rot, **kwargs)

        # Extract the upper 2x2 accessible block.

        eps_ir_sc = DielectricFunction(
            eps_ir_bulk.energy,
            eps_ir_bulk.epsilon[:, :2, :2],
            e_units=eps_ir_bulk.energy_units,
        )

        oes = SingleCrystalOpticalEigenmodes(
            eps_ir_sc, branch_tracking=oe_branch_tracking
        )

        return oes.incoherent_on_axis_optical_spectrum(
            i_pol=i_pol, d_pol=d_pol, t=t, n_f=n_f, n_b=n_b
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
            gamma_phonons_from_dict(d["phonon_calculation"])
        )
