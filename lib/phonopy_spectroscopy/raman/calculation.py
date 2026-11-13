# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""High-level `RamanCalculation` object providing an API for generating
simulated Raman spectra."""

# -------
# Imports
# -------

import warnings

from itertools import product

import numpy as np

from .intensity import (
    calculate_single_crystal_raman_intensities,
    calculate_powder_raman_intensities,
)

from .spectrum import RamanSpectrum, AngleResolvedRamanSpectrum
from .tensors import RamanTensors

from ..constants import ZERO_TOLERANCE
from ..distributions import march_dollase
from ..gamma_phonons import FullPolarGammaPhonons
from ..instrument import Geometry, Polarisation
from ..units import nm_to_ev, ev_to_nm

from ..utility.geometry import (
    cartesian_to_spherical_polar,
    rotation_matrix_from_vectors,
    rotation_matrix_from_axis_angle,
    rotate_tensors,
)

from ..utility.quadrature import (
    unit_circle_quad_rule,
    unit_sphere_lebedev_quad_rule,
)

# ----------------------
# RamanCalculation class
# ----------------------


class RamanCalculation:
    """Class for using `RamanTensors` objects to generate simulated
    Raman spectra."""

    def __init__(self, r_t, symprec=1.0e-5):
        """Create a new instance of the `RamanCalculation` class.

        Parameters
        ----------
        r_t : RamanTensors
            Raman tensors.
        symprec : float
            Tolerance for analysing the symmetry of the underlying
            structure to determine whether non-analytical corrections
            are required.
        """

        default_nac = False

        is_centrosymm, _ = r_t.structure.check_centrosymmetry(symprec=symprec)

        if not is_centrosymm:
            if isinstance(r_t.phonon_calculation, FullPolarGammaPhonons):
                default_nac = True
            else:
                warnings.warn(
                    "This Raman calculation is on a structure in a "
                    "non-centrosymmetric spacegroup but the Raman "
                    "tensors were not initialised with a "
                    "FullPolarGammaPhonons object. Non-analytical "
                    "corrections will not be automatically applied and "
                    "this may lead to less accurate results.",
                    RuntimeWarning,
                )

        self._r_t = r_t

        self._symprec = symprec
        self._default_nac = default_nac

    def _get_nac(self, nac):
        """Determine whether a non-analytical correction (NAC) should
        be applied in a Raman calculation.

        Parameters
        ----------
        nac : str or bool
            `True` or `False` to enable/disable NAC, or `"auto"` to
            set to the `default_nac` property.
        """

        if isinstance(nac, (str, np.str_)):
            nac_str = str(nac).lower()

            if nac_str == "auto":
                return self._default_nac
            else:
                raise ValueError(f'Unknown string specifier nac="{nac}".')

        return nac

    def _get_measurement_params(self, geom, i_pol, s_pol, w, t, e_rt):
        """Get measurement parameters for a Raman calculation.

        Parameters
        ----------
        geom : Geometry or None
            Measurement geometry.
        i_pol, s_pol : Polarisation, str or None
            Incident/scattered light polarisations.
        w, t : float
            Measurement wavelength and temperature.
        e_rt : float or None
            Energy at which to calculate Raman tensors.

        Returns
        -------
        params : tuple
            Tuple of `(geom, i_pol, s_pol, w, t, e_rt)` with measurement
            parameters.
        """

        geom, i_pol, s_pol = RamanCalculation.default_geometry(
            geom, i_pol, s_pol
        )

        # If the Raman tensors were calculated at a single non-zero
        # energy, this overrides w.

        if len(self._r_t.energies) == 1 and not np.isclose(
            self._r_t.energies[0], 0.0, atol=ZERO_TOLERANCE
        ):
            w_calc = ev_to_nm(self._r_t.energies[0])

            if not np.isclose(w, w_calc):
                warnings.warn(
                    f"w = {w:.1f} overridden by single non-zero Raman "
                    f"tensor energy = {w_calc:.1f}.",
                    RuntimeWarning,
                )

            w = w_calc
        else:
            if w <= 0.0:
                raise ValueError("w must be > 0.")

        # If the phonon calculation was performed at a fixed
        # temperature, this overrides t.

        t_calc = self._r_t.phonon_calculation.temperature

        if t_calc is not None:
            t = t_calc
        else:
            if t <= 0.0:
                raise ValueError("t must be > 0.")

        if e_rt is None:
            # If the Reman tensors are energy-dependent, calculate them
            # at the selected wavelength if possible.

            e_rt = nm_to_ev(w)

            if len(self._r_t.energies) > 1:
                if (
                    e_rt < self._r_t.energies.min()
                    or e_rt > self._r_t.energies.max()
                ):
                    warnings.warn(
                        f"w = {w:.1f} -> e_rt = {e_rt:.3f} is outside "
                        "the bounds of Raman tensor energy range and "
                        "will be reset to the default e_rt = 0.",
                        RuntimeWarning,
                    )

                    e_rt = 0.0
        else:
            if e_rt < 0.0:
                raise ValueError("e_rt must be >= 0.")

        return (geom, i_pol, s_pol, w, t, e_rt)

    def _get_calculation_params(
        self, geom, i_pol, s_pol, rot, nac, w, t, e_rt, lw
    ):
        """Get a complete set of parameters for a Raman calculation.

        Parameters
        ----------
        geom : Geometry or None
            Measurement geometry.
        i_pol, s_pol : Polarisation or None
            Incident and scattered light polarisation(s) (defaults:
            horizontal polarisation for x or)
        rot : array_like or None
            Rotation matrix for rotating the Raman tensors from the
            crystal reference frame into the instrument frame.
        nac : bool
            Apply a non-analytical correction (NAC) in the calculation.
        t : float or None
            Measurement temperature (defaults to 300 K if `None`;
            overridden with the calculation temperature of the
            underlying phonon calculation if set).
        w : float or None
            Measurement wavelength (defaults to 785 nm if `None`).
        e_rt : float or None
            Photon energy for evaluating Raman tensors (defaults:
            equivalent to `w` if the Raman tensors are energy-dependent,
            otherwise E = 0).
        lw : float or None
            Uniform linewidth to be used if the underlying phonon
            calculation does not have linewidths.

        Returns
        -------
        params : tuple
            Tuple of
            `(freqs, r_t, lws, irreps, geom, i_pol, s_pol, w, t)` with
            calculation parameters.
        """

        # Validate and/or set defaults for measurement parameters.

        geom, i_pol, s_pol, w, t, e_rt = self._get_measurement_params(
            geom, i_pol, s_pol, w, t, e_rt
        )

        ph_calc = self._r_t.phonon_calculation
        band_inds = self._r_t.band_indices

        r_t = None

        if self._get_nac(nac):
            # Get NAC-modified phonon calculation and Raman tensors.

            if nac and not isinstance(self._ph_calc, FullPolarGammaPhonons):
                raise RuntimeError(
                    "Non-analytical corrections (NAC) can only be "
                    "performed if the calculation was initialised with "
                    "a FullPolarGammaPhonons object."
                )

            q = geom.raman_nac_direction

            if rot is not None:
                # The q in the crystal reference frame for NAC is the
                # vector that is rotated into the instrument NAC
                # direction by rot.

                q = rot.T @ q

            ph_calc, r_t = self._r_t.get_tensors_with_nac(q, e=e_rt)

            band_inds = np.arange(0, ph_calc.num_modes, dtype=int)
        else:
            r_t = self._r_t.get_tensors(e=e_rt)

        freqs = ph_calc.frequencies[band_inds]

        if rot is not None:
            r_t = rotate_tensors(r_t, rot)

        lws = None

        if ph_calc.linewidths is not None:
            lws = ph_calc.linewidths[band_inds]
        else:
            if lw <= 0.0:
                raise ValueError("lw must be > 0.")

            lws = np.full((len(band_inds),), lw)

        irreps = None

        if ph_calc.irreps is not None:
            _, irreps = ph_calc.irreps.get_subset(band_inds)

        return (freqs, r_t, lws, irreps, geom, i_pol, s_pol, w, t)

    @property
    def structure(self):
        """Structure : Crystal structure."""
        return self._r_t.structure

    @property
    def phonon_calculation(self):
        """GammaPhonons, PolarGammaPhonons or FullPolarGammaPhonons :
        Phonon calculation."""
        return self._r_t.phonon_calculation

    @property
    def raman_tensors(self):
        """RamanTensors : Raman tensors."""
        return self._r_t

    @property
    def symprec(self):
        """float : Tolerance used for analysing the symmetry of the
        underlying structure."""
        return self._symprec

    @property
    def default_nac(self):
        """bool : `True` if non-analytical corrections (NAC) will be
        applied when `nac="auto"` is passed to the calculation routines,
        otherwise `False`.
        """
        return self._default_nac

    def single_crystal(
        self,
        hkl,
        geom=None,
        i_pol=None,
        s_pol=None,
        rot=None,
        nac="auto",
        w=785.0,
        t=300.0,
        lw=0.5,
        e_rt=None,
        **kwargs,
    ):
        """Simulate a single-crystal Raman measurement.

        Parameters
        ----------
        hkl : array_like of int
            Miller index of the crystal surface to orient antiparallel
            to the incident direction.
        geom : Geometry or None, optional
            Measurement geometry (if `None`, defaults to a conventional
            backscattering geometry).
        i_pol, s_pol : Polarisation, str or None, optional
            Polarisations of the incident and scattered light. The
            scattered polarisation may be specified by one of
            `{"parallel", "cross", "sum"}`. If `None`, `i_pol` defaults
            to the equivalent of horizontal polarisation in a
            conventional backscattering geometry, and `s_pol` defaults
            to `"sum"`.
        rot : array_like or None, optional
            Rotation matrix to realign the crystal after the `hkl`
            rotation (array shape: `(3, 3)`).
        nac : str or bool, optional
            Apply a non-analytical correction (NAC) to the phonon
            frequencies, Raman tensors and calculated linewidths if
            available (default: "auto", sets to `default_nac`).
        w, t : float, optional
            Measurement wavelength and temperature (defaults: 785 nm and
            300 K; `t` is overridden by the calculation temperature if
            set).
        e_rt : float or None, optional
            Photon energy for evaluating the Raman tensors (default:
            calculated from `w` for energy-dependent Raman tensors,
            otherwise E = 0).
        lw : float, optional
            Uniform linewidth (default: 0.5 THz, overridden by per-mode
            calculated linewidths if set).
        **kwargs : any
            Keyword arguments to `NDRamanSpectrum`.

        Returns
        -------
        spectrum : RamanSpectrum
            Simulated Raman spectrum.

        See Also
        --------
        raman.spectrum.RamanSpectrum :
            Object returned by this function.
        """

        # If the underlying structure has a conventional cell defined,
        # assume the Miller index is specified in terms of this and
        # not the primitive cell (if different).

        struct = self._r_t.phonon_calculation.structure

        r = rotation_matrix_from_vectors(
            struct.real_space_normal(hkl, conv=True),
            -1.0 * geom.incident_direction,
        )

        if rot is not None:
            if rot.shape != (3, 3):
                raise ValueError(
                    "If specified, rot must be an array_like with "
                    "shape (3, 3)."
                )

            r = rot @ r

        params = self._get_calculation_params(
            geom, i_pol, s_pol, r, nac, w, t, e_rt, lw
        )

        freqs, r_t, lws, irreps, geom, i_pol, s_pol, w, t = params

        # The Raman tensors are rotated into the instrument frame by
        # _get_calculation_params().

        ints = calculate_single_crystal_raman_intensities(
            r_t, geom, i_pol, s_pol, rot=None
        )

        return RamanSpectrum(freqs, ints, lws, w, t, irreps, **kwargs)

    def single_crystal_polarisation_rotation(
        self,
        hkl,
        geom=None,
        i_pol=None,
        s_pol=None,
        start=0.0,
        end=360.0,
        step=2.5,
        **kwargs,
    ):
        """Simulate a single-crystal polarisation rotation experiment
        where the incident or scattered polarisation is rotated about
        the incident/collection direction.

        Parameters
        ----------
        hkl : array_like of int
            Miller index of the crystal surface to orient antiparallel
            to the incident direction.
        geom : Geometry or None, optional
            Measurement geometry (if `None`, defaults to a conventional
            backscattering geometry).
        i_pol, s_pol : Polarisation, str or None, optional
            Polarisations of the incident and scattered light. The
            polarisation to be rotated is specified by `"rot"`. The
            defaults are the same as for `single_crystal`,
        start, end, step : float, optional
            Start/end angle and angle step for polarisation rotation in
            degrees (defaults: 0 -> 360 deg in 2.5 deg steps).
        **kwargs : any
            Optional arguments to `single_crystal`.

        Returns
        -------
        spectrum : AngleResolvedRamanSpectrum
            Simulated Raman spectrum.

        See Also
        --------
        single_crystal : Simulate a single-crystal measurement.

        Notes
        -----
        See `single_crystal` for optional keyword arguments.
        """

        geom, pols, angles = RamanCalculation.polarisation_rotation(
            geom, i_pol, s_pol, start, end, step
        )

        # To avoid redundant work, simluate default Stokes spectra in
        # THz for each rotation angle and convert to the user-specified
        # sp_type and sp_units when constructing the 2D spectrum.

        sp_type = kwargs.pop("sp_type", "stokes")
        sp_units = kwargs.pop("sp_units", "thz")

        # Also no point in applying these until the "final" spectrum.

        sp_range = kwargs.pop("sp_range", None)
        sp_res = kwargs.pop("sp_res", None)

        # Simulate a series of single-crystal spectra and combine into
        # an AngleResolvedRamanSpectrum object.

        spectra = [
            self.single_crystal(
                hkl, geom=geom, i_pol=i_pol, s_pol=s_pol, **kwargs
            )
            for i_pol, s_pol in pols
        ]

        sp_stack = spectra[0].stack_with(*spectra[1:])

        # Combine into an AngleResolvedRamanSpectrum object.

        return AngleResolvedRamanSpectrum(
            sp_stack.band_raman_shifts,
            sp_stack.band_intensities,
            sp_stack.band_linewidths,
            sp_stack.measurement_wavelength,
            sp_stack.measurement_temperature,
            angles,
            irreps=None,
            ir_syms=sp_stack.band_irrep_symbols,
            sp_type=sp_type,
            sp_units=sp_units,
            sp_range=sp_range,
            sp_res=sp_res,
            angle_unit_text_label=r"\chi / deg",
            angle_unit_plot_label=r"$\chi$ / $^\circ$",
        )

    def single_crystal_crystal_rotation(
        self,
        hkl,
        geom=None,
        start=0.0,
        end=360.0,
        step=2.5,
        rot_axis="incident",
        **kwargs,
    ):
        """Simulate a single-crystal crystal rotation measurement where
        the crystal is rotated about the incident or collection
        direction.

        Parameters
        ----------
        hkl : array_like of int
            Miller index of the crystal surface to orient antiparallel
            to the incident direction.
        geom : Geometry or None, optional
            Measurement geometry (if `None`, defaults to a conventional
            backscattering geometry).
        start, end, step : float, optional
            Start/end angle and angle step for crystal rotation in
            degrees (defaults: 0 -> 360 deg in 2.5 deg steps).
        rot_axis : {"incident", "collection"}, optional
            Axis to rotate around (default: "incident")
        **kwargs : any
            Optional arguments to `single_crystal`.

        Returns
        -------
        spectrum : AngleResolvedRamanSpectrum
            Simulated Raman spectrum.

        See Also
        --------
        single_crystal : Simulate a single-crystal measurement.

        Notes
        -----
        See `single_crystal` for optional keyword arguments. Note that
        specifying the `rot` keyword is invalid and will raise an error.
        """

        if "rot" in kwargs and kwargs["rot"] is not None:
            raise RuntimeError(
                "The rot keyword to single_crystal is not valid for a "
                "crystal-rotation experiment - specify the required "
                "rotation(s) with phi_start, phi_end and phi_step "
                "instead."
            )

        sp_type = kwargs.pop("sp_type", "stokes")
        sp_units = kwargs.pop("sp_units", "thz")
        sp_range = kwargs.pop("sp_range", None)
        sp_res = kwargs.pop("sp_res", None)

        geom, rots, angles = RamanCalculation.crystal_rotation(
            geom, start, end, step, rot_axis
        )

        spectra = [
            self.single_crystal(hkl, geom=geom, rot=rot, **kwargs)
            for rot in rots
        ]

        sp_stack = spectra[0].stack_with(*spectra[1:])

        return AngleResolvedRamanSpectrum(
            sp_stack.band_raman_shifts,
            sp_stack.band_intensities,
            sp_stack.band_linewidths,
            sp_stack.measurement_wavelength,
            sp_stack.measurement_temperature,
            angles,
            irreps=None,
            ir_syms=sp_stack.band_irrep_symbols,
            sp_type=sp_type,
            sp_units=sp_units,
            sp_range=sp_range,
            sp_res=sp_res,
            angle_unit_text_label=r"\psi / deg",
            angle_unit_plot_label=r"$\psi$ / $^\circ$",
        )

    def _powder_integration(
        self, geom, i_pol, s_pol, po_norm, po_r, w, t, e_rt, lw, **kwargs
    ):
        """Simulate a powder Raman spectrum using numerical integration.

        Parameters
        ----------
        geom : Geometry or None
            Measurement geometry
        i_pol, s_pol : Polarisation, str or None
            Polarisations of the incident and scattered light.
        po_norm : array_like or None
            Surface normal for preferred orientation.
        po_r : float
            r parameter for preferred orientation.
        w, t : float
            Measurement wavelength and temperature.
        e_rt : float or None
            Photon energy for evaluating the Raman tensors.
        lw : float
            Uniform linewidth.
        **kwargs : any
            Keyword arguments to `RamanSpectrum`.

        Returns
        -------
        spectrum : RamanSpectrum
            Simulated Raman spectrum.

        See Also
        --------
        raman.intensity.calculate_powder_raman_intensities :
            Lower-level API function used to calculate powder Raman
            intensities.
        raman.spectrum.RamanSpectrum :
            Objects returned by this function.
        """

        params = self._get_calculation_params(
            geom, i_pol, s_pol, None, False, w, t, e_rt, lw
        )

        freqs, r_t, lws, irreps, geom, i_pol, s_pol, w, t = params

        ints = calculate_powder_raman_intensities(
            r_t, geom, i_pol, s_pol, po_r=po_r, po_norm=po_norm
        )

        return RamanSpectrum(freqs, ints, lws, w, t, irreps=irreps, **kwargs)

    def _powder_orientation_average(
        self,
        p,
        geom,
        i_pol,
        s_pol,
        nac,
        po_norm,
        po_r,
        w,
        t,
        e_rt,
        lw,
        **kwargs,
    ):
        """Simulate a powder Raman spectrum from an explicit average
        over crystal orientations.

        Parameters
        ----------
        p : int or str
            Integer precision or string specifier for the
            Lebedev + circle quadrature rule for averaging.
        geom : Geometry or None
            Measurement geometry
        i_pol, s_pol : Polarisation, str or None
            Polarisations of the incident and scattered light.
        nac : bool
            Apply a non-analytical correction (NAC) to the phonon
            frequencies, Raman tensors and calculated linewidths.
        po_norm : array_like or None
            Surface normal for preferred orientation.
        po_r : float
            r parameter for preferred orientation.
        w, t : float
            Measurement wavelength and temperature.
        e_rt : float or None
            Photon energy for evaluating the Raman tensors.
        lw : float
            Uniform linewidth.
        **kwargs : any
            Keyword arguments to `AveragedRamanSpectrum`.

        Returns
        -------
        spectrum : RamanSpectrum
            Simulated Raman spectrum.

        See Also
        --------
        raman.spectrum.AveragedRamanSpectrum :
            Objects returned by this function.
        utility.quadrature.unit_sphere_lebedev_quad_rule_available_precs
            Get a list of available precisions `p`.
        """

        # Set the default geometry to obtain the incident direction
        # q_i and the NAC vector q_nac.

        geom, i_pol, s_pol = RamanCalculation.default_geometry(
            geom, i_pol, s_pol
        )

        q_i = geom.incident_direction
        q_nac = geom.raman_nac_direction

        # Orientation vectors and \phi/\theta from a Lebedev quadrature
        # rule.

        v_l, w_l = unit_sphere_lebedev_quad_rule(p, angles=False)
        a_l = cartesian_to_spherical_polar(v_l)[:, 1:]

        n = len(v_l)

        # \psi angles from a circle rule. If there is no preferred
        # orientation, or if q_i and q_nac are collinear, a minimal rule
        # with m = 8 points is suitable. If not, a rule with the same
        # number of \psi as the \phi in the Lebedev rule provides a
        # "balanced" sampling of orientations.

        m = 8

        po = np.abs(po_r - 1.0) > ZERO_TOLERANCE
        abs_dp = np.abs(np.dot(q_nac, q_i))

        if po and np.abs(abs_dp - 1.0) > ZERO_TOLERANCE:
            m = len(np.unique(a_l[:, 0]))

        a_c, w_c = unit_circle_quad_rule(m, angles=True)

        w_sp = np.repeat(w_l, m) * np.tile(w_c, n)

        freqs = np.zeros(
            (n * m, self._r_t.ph_calc.num_modes), dtype=np.float64
        )

        ints = np.zeros_like(freqs, dtype=np.float64)
        lws = np.zeros_like(freqs, dtype=np.float64)

        # Precompute rotation matrices for circle rule.

        psi_rots = np.array(
            [rotation_matrix_from_axis_angle(q_nac, a) for a in a_c],
            dtype=np.float64,
        )

        for i, (v, r_a) in enumerate(product(v_l, psi_rots)):
            rot = r_a @ rotation_matrix_from_vectors(v, q_nac)

            params = self._get_calculation_params(
                geom, i_pol, s_pol, rot, nac, w, t, e_rt, lw
            )

            freqs, r_t, lws, _, geom, i_pol, s_pol, w, t = params

            freqs[i, :] = freqs

            ints[i, :] = calculate_single_crystal_raman_intensities(
                r_t, geom, i_pol, s_pol, rot=None
            )

            lws[i, :] = lws

            if po:
                # The angle \alpha in the March-Dollase distribution is
                # the angle between the rotated preferred direction
                # and the incident axis.

                w_sp[i] *= march_dollase(
                    np.acos(np.dot(rot @ po_norm, q_i)), po_r
                )

        or_a = np.degrees(
            np.concatenate((np.repeat(a_l, m), np.tile(a_c, n)), axis=1)
        )

        return RamanSpectrum(
            freqs, ints, lws, w, t, or_a=or_a, or_w=w_sp, **kwargs
        )

    def powder(
        self,
        geom=None,
        i_pol=None,
        s_pol=None,
        po_hkl=None,
        po_r=1.0,
        nac="auto",
        p="medium",
        w=785.0,
        t=300.0,
        e_rt=None,
        lw=0.5,
        **kwargs,
    ):
        """Simulate a powder Raman spectrum using numerical integration.

        Parameters
        ----------
        geom : Geometry or None, optional
            Measurement geometry (if `None`, defaults to a conventional
            backscattering geometry).
        i_pol, s_pol : Polarisation, str or None, optional
            Polarisations of the incident and scattered light. The
            scattered polarisation may be specified by one of
            `{"parallel", "cross", "sum"}`. If `None`, `i_pol` defaults
            to the equivalent of horizontal polarisation in a
            conventional backscattering geometry, and `s_pol` defaults
            to `"sum"`.
        po_hkl : array_like of int or None, optional
            Miller index of the preferred orientation (default: `None`).
        po_r : float, optional
            r parameter for March-Dollase distribution used to model
            preferred orientation (default: 1.0)
        nac : str or bool, optional
            Apply a non-analytical correction (NAC) to the phonon
            frequencies, Raman tensors and calculated linewidths if
            available (default: "auto", sets to `default_nac`).
        p : int, optional
            Integer precision or one of
            {"min", "low", "medium", "high", "ultra", "max"} specifying
            the Lebedev + circle quadrature rule for powder averaging
            with NAC.
        w, t : float, optional
            Measurement wavelength and temperature (defaults: 785 nm and
            300 K; `t` is overridden by the calculation temperature if
            set).
        e_rt : float or None, optional
            Photon energy for evaluating the Raman tensors in energy-
            dependent Raman calculations (default: calculated from `w`
            if set and if energy-dependent Raman tensors are available,
            otherwise E = 0.)
        lw : float, optional
            Uniform linewidth (default: 0.5 THz, overridden by per-mode
            linewidths from the underlying phonon calculation if set).
        **kwargs : any
            Keyword arguments to `RamanSpectrum`.

        Returns
        -------
        spectrum : RamanSpectrum
            Simulated Raman spectrum.

        See Also
        --------
        raman.spectrum.RamanSpectrum :
            Objects returned by this function.
        """

        if np.abs(po_r - 1.0) > ZERO_TOLERANCE and po_hkl is None:
            raise ValueError(
                "po_hkl must be specified for calculations with "
                "preferred orientation."
            )

        po_norm = None

        if po_hkl is not None:
            po_norm = self.r_t.structure.real_space_normal(po_hkl, conv=True)

        nac = self._get_nac(nac)

        if nac:
            return self._powder_orientation_average(
                p,
                geom,
                i_pol,
                s_pol,
                nac,
                po_norm,
                po_r,
                w,
                t,
                e_rt,
                lw,
                **kwargs,
            )
        else:
            return self._powder_integration(
                geom, i_pol, s_pol, po_norm, po_r, w, t, e_rt, lw, **kwargs
            )

    def powder_polarisation_rotation(
        self,
        geom=None,
        i_pol=None,
        s_pol=None,
        start=0.0,
        end=360.0,
        step=2.5,
        **kwargs,
    ):
        """Simulate a powder polarisation rotation experiment where the
        incident or scattered polarisation is rotated about the
        incident/collection direction.

        Parameters
        ----------
        geom : Geometry or None, optional
            Measurement geometry (if `None`, defaults to a conventional
            backscattering geometry).
        i_pol, s_pol : Polarisation, str or None, optional
            Polarisations of the incident and scattered light. The
            polarisation to be rotated is specified by `"rot"`. The
            defaults are the same as for `powder`,
        start, end, step : float, optional
            Start/end angle and angle step for polarisation rotation in
            degrees (defaults: 0 -> 360 deg in 2.5 deg steps).
        **kwargs : any
            Optional arguments to `powder`.

        Returns
        -------
        spectrum : AngleResolvedRamanSpectrum
            Simulated Raman spectrum.

        See Also
        --------
        powder : Simulate a powder measurement.

        Notes
        -----
        See `powder` for optional keyword arguments.
        """

        geom, pols, angles = RamanCalculation.polarisation_rotation(
            geom, i_pol, s_pol, start, end, step
        )

        sp_type = kwargs.pop("sp_type", "stokes")
        sp_units = kwargs.pop("sp_units", "thz")
        sp_range = kwargs.pop("sp_range", None)
        sp_res = kwargs.pop("sp_res", None)

        # Simulate a series of powder spectra. NAC should be handled
        # transparently by powder() and stack_with().

        spectra = [
            self.powder(geom=geom, i_pol=i_pol, s_pol=s_pol, **kwargs)
            for i_pol, s_pol in pols
        ]

        sp_stack = spectra[0].stack_with(*spectra[1:])

        return AngleResolvedRamanSpectrum(
            sp_stack.band_raman_shifts,
            sp_stack.band_intensities,
            sp_stack.band_linewidths,
            sp_stack.measurement_wavelength,
            sp_stack.measurement_temperature,
            angles,
            irreps=None,
            ir_syms=sp_stack.band_irrep_symbols,
            sp_type=sp_type,
            sp_units=sp_units,
            sp_range=sp_range,
            sp_res=sp_res,
            angle_unit_text_label=r"\chi / deg",
            angle_unit_plot_label=r"$\chi$ / $^\circ$",
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

        return {"raman_tensors": self._r_t.to_dict(), "symprec": self._symprec}

    @staticmethod
    def default_geometry(geom=None, i_pol=None, s_pol=None):
        """Get default geometry and/or polarisations for a Raman
        calculation.

        Parameters
        ----------
        geom : Geometry or None, optional
            Measurement geometry.
        i_pol, s_pol : Polarisation, str or None, optional
            Polarisations of the incident and scattered light. The
            scattered polarisation may be specified by one of
            `{"parallel", "cross", "sum"}`.

        Returns
        -------
        params : tuple
            Tuple of `(geom, i_pol, s_pol)` with geometry and
            polarisation.

        Notes
        -----
        The defaults for unspecified parameters are as follows:
          * `geom`: conventional backscattering geometry.
          * `i_pol`: equivalent of horizontal polarisation in a
            conventional backscattering geometry.
          * `s_pol`: sum of parallel and perpendicular polarisations.

        If some of the parameters are specified these defaults are not
        guaranteed to be suitable. This is particularly likely to be the
        case if a non-standard geometry is specified but one or both
        polarisations are left to the defaults.
        """

        # Default to conventional backscattering geometry.

        if geom is None:
            geom = Geometry.conventional_backscattering()

        # Default polarisations: conventional horizontal polarisation
        # = x for a backscattering geometry, rotated for the incident
        # direction, and sum of parallel and cross polarisations.

        if i_pol is None:
            i_pol = Polarisation.conventional_horizontal(
                geom.incident_direction
            )

        if s_pol is None:
            s_pol = "sum"

        if isinstance(s_pol, (str, np.str_)):
            match str(s_pol).lower():
                case "parallel":
                    s_pol = i_pol
                case "cross":
                    s_pol = Polarisation.cross_to(
                        i_pol, geom.incident_direction
                    )
                case "sum":
                    s_pol = Polarisation.sum_parallel_cross_to(
                        i_pol, geom.incident_direction
                    )

        if not (
            i_pol.check_perpendicular(geom.incident_direction)
            and s_pol.check_perpendicular(geom.collection_direction)
        ):
            raise ValueError(
                "The specified i_pol or s_pol are not possible with "
                "the supplied or default geometry."
            )

        return (geom, i_pol, s_pol)

    @staticmethod
    def polarisation_rotation(geom, i_pol, s_pol, start, end, step):
        r"""Get a set of M pairs of incident and scattered polarisations
        for simulating a polarisation rotation experiment where one of
        the polarisations is rotated through a series of angles \chi.

        Parameters
        ----------
        geom : Geometry or None
            Measurement geometry (if `None`, defaults to a conventional
            backscattering geometry).
        i_pol, s_pol : Polarisation, str or None
            Incident/scattered light polarisations. The polarisation to
            be rotated must be set to `"rot"`. `s_pol` can be specified
            by one of `{"parallel", "cross", "sum"}`.
        start, end, step : float
            Start/end rotation angle and angle step.

        Returns
        -------
        params : tuple
            Tuple of `(geom, pols, angles)` with the measurement
            geometry, polarisation pairs (shape: `(M, 2)`) and rotation
            angles (shape: `(M,)`).
        """

        i_rot = (
            isinstance(i_pol, (str, np.str_)) and str(i_pol).lower() == "rot"
        )

        s_rot = (
            isinstance(s_pol, (str, np.str_)) and str(s_pol).lower() == "rot"
        )

        if not (i_rot or s_rot):
            raise ValueError('One of i_pol or s_pol must be set to "rot".')

        if i_rot and s_rot:
            raise ValueError('i_pol and s_pol cannot both be set to "rot".')

        angles = np.arange(start, end + step / 10.0, step)

        if len(angles) == 0:
            raise ValueError(
                "start, end and step must define at least one rotation "
                "angle."
            )

        if geom is None:
            geom = Geometry.conventional_backscattering()

        i_pols, s_pols = None, None

        if i_rot:
            i_pols = Polarisation.from_angles(geom.incident_direction, angles)

            # get_default_geometry() handles the case where s_pol is a
            # string or None.

            s_pols = []

            for i_pol in i_pols:
                _, _, s_pol = RamanCalculation.default_geometry(
                    geom, i_pol, s_pol
                )
                s_pols.append(s_pol)

        if s_rot:
            s_pols = Polarisation.from_angles(
                geom.collection_direction, angles
            )

            i_pols = []

            for s_pol in s_pols:
                _, i_pol, _ = RamanCalculation.default_geometry(
                    geom, i_pol, s_pol
                )
                i_pols.append(i_pol)

        # Transpose polarisation array to swap (2, M) -> (M, 2).

        return (geom, np.array([i_pols, s_pols], dtype=object).T, angles)

    @staticmethod
    def crystal_rotation(geom, start, end, step, rot_axis):
        r"""Get a set of M rotation matrices for simulating a crystal
        rotation experiment where the crystal is rotated about either
        the incident or collection direction by a series of angles
        \psi.

        Parameters
        ----------
        geom : Geometry or None
            Measurement geometry (if `None`, defaults to a conventional
            backscattering geometry).
        start, end, step : float
            Start/end rotation angle and angle step.
        rot_axis : {"incident", "collection"}
            Rotation axis.

        Returns
        -------
        params : tuple
            Tuple of `(geom, rots, angles)` with the measurement
            geometry, rotation matrices (shape: `(M, 3, 3)`) and
            rotation angles (shape: `(M,)`).
        """

        angles = np.radians(
            np.arange(start, end + step / 10.0, step, dtype=np.float64)
        )

        if len(angles) == 0:
            raise ValueError(
                "start, end and step must define at least one rotation "
                "angle."
            )

        if geom is None:
            geom = Geometry.conventional_backscattering()

        axis = None

        match str(rot_axis).lower():
            case "incident":
                axis = geom.incident_direction
            case "collection":
                axis = geom.collection_direction
            case _:
                raise ValueError(f'Unknown rot_axis specifier "{rot_axis}".')

        rots = [rotation_matrix_from_axis_angle(axis, a) for a in angles]

        return (geom, np.array(rots, dtype=np.float64), angles)

    @staticmethod
    def from_dict(d):
        """Create a new `RamanCalculation` instance from a dictionary
        generated by `RamanCalculation.to_dict()`.

        Parameters
        ----------
        d : dict
            Dictionary generated by `to_dict()`.

        Returns
        -------
        calc : RamanCalculation
            `RamanCalculation` object constructed from the data in `d`.
        """

        return RamanCalculation(
            RamanTensors.from_dict(d["raman_tensors"]), d["symprec"]
        )
