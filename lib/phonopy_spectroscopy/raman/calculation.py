# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""High-level `RamanCalculation` object providing a simplified API for
combining Gamma-point phonon and Raman tensor calculations to generate
simulated spectra."""


# -------
# Imports
# -------


import warnings

import numpy as np

from ..instrument import Polarisation

from .intensity import (
    calculate_single_crystal_raman_intensities,
    calculate_powder_raman_intensities,
)

from .spectrum import RamanSpectrum1D, RamanSpectrum2D
from .tensors import RamanTensors

from ..constants import ZERO_TOLERANCE
from ..phonon import GammaPhonons
from ..units import nm_to_ev

from ..utility.geometry import (
    rotation_matrix_from_vectors,
    rotation_matrix_from_axis_angle,
)

from ..utility.numpy_helper import (
    np_readonly_view,
    np_asarray_copy,
    np_check_shape,
    np_expand_dims,
)


# ----------------------
# RamanCalculation class
# ----------------------


class RamanCalculation:
    """Combine `GammaPhonons` and `RamanTensors` objects for a complete
    Raman calculation to generate simulated spectra."""

    def __init__(self, gamma_ph, r_t, band_inds=None):
        """Create a new instance of the `RamanCalculation` class.

        Parameters
        ----------
        gamma_ph : GammaPhonons
            Gamma-point phonon calculation.
        r_t : RamanTensors
            Raman tensors.
        band_inds : array_like or None, optional
            Band indices for which Raman tensors were calculated
            (default: all bands).
        """

        # Check/set band indices used in calculations and update irreps
        # for included indices if required.

        irreps = None

        if band_inds is not None:
            band_inds = np_asarray_copy(band_inds, dtype=int)

            if not np_check_shape(band_inds, (None,)):
                raise ValueError(
                    "If supplied, band_inds must be an array_like with "
                    "shape (N,)."
                )

            if len(band_inds) == 0:
                raise ValueError(
                    "If supplied, band_inds must specify at least one band."
                )

            if len(band_inds) != len(set(band_inds)):
                raise ValueError(
                    "One or more indices in band_inds are duplicates."
                )

            for idx in band_inds:
                if idx < 0 or idx >= gamma_ph.num_modes:
                    raise ValueError(
                        "One or more indices in band_inds are "
                        "incompatible with the number of modes in the "
                        "phonon calculation."
                    )

            if gamma_ph.has_irreps:
                irreps = gamma_ph.irreps.get_subset(band_inds, reset_inds=True)

        else:
            band_inds = np.array(list(range(gamma_ph.num_modes)), dtype=int)

        if len(band_inds) != len(r_t.raman_tensors):
            raise ValueError(
                "The number of Raman tensors is inconsistent with the "
                "number of band indices (or the number of modes in the "
                "phonon calculation if band_inds was not supplied)."
            )

        self._gamma_ph = gamma_ph
        self._r_t = r_t

        self._band_inds = band_inds
        self._irreps = irreps

    def _get_band_inds_from_grp_inds(self, band_grp_inds):
        """Determine the band indices for a set of band groups.

        Parameters
        ----------
        band_grp_inds : array_like or None
            Band groups.

        Returns
        -------
        band_inds : numpy.ndarray
            Band indices.
        """

        if band_grp_inds is not None:
            for idx in band_grp_inds:
                if idx < 0 or idx >= self.num_band_groups:
                    raise ValueError(
                        "One or more band group indices are not "
                        "compatible with the number of irrep groups or "
                        "the number of bands in the calculation."
                    )

            if self.irreps is not None:
                # Band group indices map to irrep groups.

                band_inds = []

                for idx in band_grp_inds:
                    band_inds.extend(self.irreps.irrep_band_indices[idx])

                return np.array(band_inds, dtype=int)
            else:
                # Band group indices map to individual bands.

                return np.array(band_grp_inds, dtype=int)

        return np.array(list(range(len(self._band_inds))), dtype=int)

    def _get_calc_params(
        self, geom, i_pols, s_pols, w, e_rt, lw, band_grp_inds
    ):
        """Determine parameters for Raman calculations.

        Parameters
        ----------
        geom : Geometry
            Measurement geometry.
        i_pols, s_pols : str, Polarisation or array_like
            Incident and scattered light polarisation(s). `s_pol` may be
            specified by one of {"parallel", "cross", "sum"}.
        lw : float or None
            Uniform linewidth or scale factor for calculated linewidths,
            depending on whether calculation has linewidths.
        w : float or None
            Measurement wavelength (default: None).
        e_rt : float or None
            Photon energy for evaluating Raman tensors (default:
            calculated from `w` if set, and if energy-dependent Raman
            tensors are available, otherwise E = 0.)
        band_grp_inds : array_like or None
            Indices of band groups to include in the calculation.

        Returns
        -------
        params : dict
            Dictionary of calculation parameters with the following: `{
                "frequencies": numpy.ndarray,
                "raman_tensors": numpy.ndarray,
                "linewidths": numpy.ndarray,
                "irreps": Irreps or None
                "geometry": Geometry,
                "incident_polarisations": numpy.ndarray,
                "scattered_polarisations": numpy.ndarray,
                "is_2d_spectrum": bool
                }`
        """

        # Determine band indices to include in the calculation.

        band_inds = self._get_band_inds_from_grp_inds(band_grp_inds)

        params = {}

        # Frequencies.

        params["frequencies"] = self.frequencies[band_inds]

        # Raman tensors.

        if e_rt is None:
            e_rt = 0.0

            if w is not None:
                if self._r_t.is_energy_dependent:
                    e_rt = nm_to_ev(w)
                else:
                    warnings.warn(
                        "A measurement wavelength was specified but the "
                        "calculation is not using energy-dependent Raman "
                        "tensors. Raman intensities will be calculated "
                        "using the far from resonance approximation and "
                        "the wavelength will only be used to apply "
                        "intensity modulation.",
                        RuntimeWarning,
                    )

        r_t = self._r_t.get_tensors_at_energy(e_rt)
        params["raman_tensors"] = r_t[band_inds]

        # Linewidths.

        if lw is None:
            # Default value.

            lw = 1.0 if self._gamma_ph.has_linewidths else 0.5
        else:
            if lw < ZERO_TOLERANCE:
                raise ValueError("lw cannot be zero or negative.")

        params["linewidths"] = (
            lw * self.linewidths[band_inds]
            if self.linewidths is not None
            else lw * np.ones((len(band_inds),), dtype=np.float64)
        )

        # Irreps.

        params["irreps"] = (
            self.irreps.get_subset(band_inds, reset_inds=True)
            if self.irreps is not None
            else None
        )

        # Geometry and incident/scattered polarisations.

        i_pols, i_pols_n_dim_add = np_expand_dims(
            np.asarray(i_pols, dtype=object), (None,)
        )

        is_2d = i_pols_n_dim_add == 0

        s_pols_str = str(s_pols).lower()

        if s_pols_str == "parallel":
            s_pols = i_pols
        elif s_pols_str == "cross":
            s_pols = Polarisation.cross_to(i_pols, geom.incident_direction)
        elif s_pols_str == "sum":
            s_pols = Polarisation.sum_parallel_cross_to(
                i_pols, geom.incident_direction
            )
        else:
            s_pols, s_pols_n_dim_add = np_expand_dims(
                np.asarray(s_pols, dtype=object), (None,)
            )

            is_2d = is_2d or s_pols_n_dim_add == 0

            i_pols, s_pols = np.broadcast_arrays(i_pols, s_pols)

        params["geometry"] = geom
        params["incident_polarisations"] = i_pols
        params["scattered_polarisations"] = s_pols

        # Flag for whether the calculation is working with a single or
        # multiple polarisations.

        params["is_2d_spectrum"] = is_2d

        return params

    @property
    def frequencies(self):
        """numpy.ndarray : Band frequencies for the subset of bands in
        the calculation."""
        return self._gamma_ph.frequencies[self._band_inds]

    @property
    def linewidths(self):
        """numpy.ndarrray : Linewidths for the subset of bands in the
        calculation."""
        lws = self._gamma_ph.linewidths
        return lws[self._band_inds] if lws is not None else None

    @property
    def irreps(self):
        """Irreps : Irreps for the subset of bands in the calculation."""

        # If the calculation is for a subset of the bands, and the
        # underlying phonon calculation has irreps, the internal
        # _irreps field will be set. If _irreps is None, then return the
        # irreps from the phonon calculation.

        return (
            self._irreps if self._irreps is not None else self._gamma_ph.irreps
        )

    @property
    def band_indices(self):
        """numpy.ndarray : Indices of the bands in the calculation."""
        return np_readonly_view(self._band_inds.view())

    @property
    def structure(self):
        """Structure : Underlying `Structure` object."""
        return self._gamma_ph.structure

    @property
    def gamma_phonons(self):
        """GammaPhonons : Underlying `GammaPhonons` object."""
        return self._gamma_ph

    @property
    def raman_tensors(self):
        """RamanTensors : Underlying `RamanTensors` object."""
        return self._r_t

    @property
    def num_band_groups(self):
        """int : Number of band groups in the calculation. Band groups
        are equivalent to irrep groups if the underlying phonon
        calculation has irreps, or to the number of bands in the
        calculation otherwise."""

        if self.irreps is not None:
            return len(self.irreps.irrep_band_indices)

        return len(self._band_inds)

    def single_crystal(
        self,
        hkl,
        geom,
        i_pol,
        s_pol,
        rot=None,
        w=None,
        t=None,
        lw=None,
        e_rt=None,
        band_grp_inds=None,
        **kwargs
    ):
        """Simulate a single-crystal Raman measurement.

        Parameters
        ----------
        hkl : array_like of int
            Miller index of the crystal surface to orient antiparallel
            to the incident direction.
        geom : Geometry
            Measurement geometry.
        i_pol, s_pol : str, Polarisation or list of Polarisation
            Polarisations of the incident and scattered light. The
            scattered polarisation may also be specified by one of
            {"parallel", "cross", "sum"}.
        rot : array_like or None, optional
            Rotation matrix or set of matrices to realign the crystal
            after the `hkl` rotation.
        w, t : float or None, optional
            Measurement wavelength and temperature.
        e_rt : float or None, optional
            Photon energy for evaluating the Raman tensors in energy-
            dependent Raman calculations (default: calculated from `w`
            if set and if energy-dependent Raman tensors are available,
            otherwise E = 0.)
        lw : float or None, optional
            Uniform linewidth or scale factor for calculated linewidths
            (defaults: 0.5 THz uniform linewidth or scale factor of
            1.0, depending on whether the calculation has linewidths).
        band_grp_inds : array_like or None, optional
            Indices of "band groups" to include in the calculation
            (default: all groups). Groups are defined by irreps if
            available, or the band indices in the calculation otherwise.
        **kwargs : any
            Keyword arguments to the `RamanSpectrum1D` and
            `RamanSpectrum2D` constructors (some of these may be
            required depending on other parameters - see notes below).

        Returns
        -------
        spectrum : RamanSpectrum1D or RamanSpectrum2D
            Simulated Raman spectrum.

        See Also
        --------
        num_band_groups
            Number of band groups in the calculation.
        raman.spectrum.RamanSpectrum1D, raman.spectrum.RamanSpectrum2D
            Objects returned by this function.

        Notes
        -----
        * If both `i_pols` and `s_pols` are single polarisations, and
          `rots` is either None or a single rotation, the function
          returns a 1D spectrum.
        * If `i_pols` and/or `s_pols`, or `rots`, specify multiple
          values, the function returns an 2D spectrum.
        * A single function call cannot combine multiple polarisations
          and multiple rotations.
        * The `e_rt` parameter is provided so that the intensity
          modulation "envelope" due to the laser wavelength and changes
          in intensity due energy-dependent polarisability can be
          modelled separately.
        * The indices in `band_grp_inds` correspond directly to the
          entries in the peak tables in the `RamanSpectrum` objects
          produced with the default `band_grp_inds=None`). The number of
          band groups can be obtained with the `num_band_groups`
          property.
        * If called with multiple polarisations to generate a 2D
          spectrum, additional keyword args to the `RamanSpectrum2D`
          constructor will need to be specified.
        """

        params = self._get_calc_params(
            geom, i_pol, s_pol, w, e_rt, lw, band_grp_inds
        )

        is_2d = params["is_2d_spectrum"]

        if rot is not None:
            rot, n_dim_add = np_expand_dims(np.asarray(rot), (None, 3, 3))

            if not rot.shape[0] > 1 and not params["single_polarisation"]:
                raise RuntimeError(
                    "Multiple incident/scattered polarisations cannot be "
                    "combined with multiple crystal rotations."
                )

            is_2d = is_2d or n_dim_add == 0

        # If the underlying structure has a conventional cell defined,
        # assume the Miller index is speficied in terms of this and
        # not the primitive cell (if different).

        r = rotation_matrix_from_vectors(
            self._gamma_ph.structure.real_space_normal(hkl, conv=True),
            -1.0 * geom.incident_direction,
        )

        rot = np.matmul(rot, r) if rot is not None else r.reshape((1, 3, 3))

        ints = None

        if len(rot) == 1:
            # Single rotation, potentially with multiple polarisations.

            ints = np.zeros(
                (
                    len(params["incident_polarisations"]),
                    params["raman_tensors"].shape[0],
                ),
                dtype=np.float64,
            )

            for i, (i_pol, s_pol) in enumerate(
                zip(
                    params["incident_polarisations"],
                    params["scattered_polarisations"],
                )
            ):
                ints[i] = calculate_single_crystal_raman_intensities(
                    params["raman_tensors"],
                    params["geometry"],
                    i_pol,
                    s_pol,
                    rot[0],
                )
        else:
            # Single polarisation, potentially with multiple rotations.

            ints = np.zeros(
                (rot.shape[0], params["raman_tensors"].shape[0]),
                dtype=np.float64,
            )

            for i, r in enumerate(rot):
                ints[i] = calculate_single_crystal_raman_intensities(
                    params["raman_tensors"],
                    params["geometry"],
                    params["incident_polarisations"][0],
                    params["scattered_polarisations"][0],
                    r,
                )

        ints = np.array(ints, dtype=np.float64).T

        if is_2d:
            if (
                "d2_axis_vals" not in kwargs
                or "d2_unit_text_label" not in kwargs
            ):
                raise RuntimeError(
                    "Multiple polarisations/rotations produce 2D "
                    "spectra and require the optional d2_axis_vals and "
                    "d2_unit_text_label keywords to be specified."
                )

            return RamanSpectrum2D(
                params["frequencies"],
                ints,
                params["linewidths"],
                kwargs.pop("d2_axis_vals"),
                kwargs.pop("d2_unit_text_label"),
                irreps=params["irreps"],
                w=w,
                t=t,
                **kwargs
            )
        else:
            return RamanSpectrum1D(
                params["frequencies"],
                ints,
                params["linewidths"],
                irreps=params["irreps"],
                w=w,
                t=t,
                **kwargs
            )

    def single_crystal_polarisation_rotation(
        self,
        hkl,
        geom,
        i_pol,
        s_pol,
        chi_start=0.0,
        chi_end=360.0,
        chi_step=2.5,
        **kwargs
    ):
        """Simulate a single-crystal polarisation (chi) rotation
        measurement where one of the incident or scattered polarisations
        are rotated about the incident or collection direction
        respectively.

        Parameters
        ----------
        hkl : array_like of int
            Miller index of the crystal surface to orient antiparallel
            to the incident direction.
        geom : Geometry
            Measurement geometry.
        i_pol, s_pol : str or Polarisation
            Polarisations of the incident and scattered light. The
            polarisation to be rotated can be specified by "rot". The
            scattered polarisation may also be specified by one of
            {"parallel", "cross", "sum"}.
        chi_start, chi_end, chi_step : float, optional
            Start/end angle and angle step for polarisation rotation in
            degrees (defaults: 0 -> 360 deg in 2.5 deg steps).
        **kwargs : any
            Optional arguments to `single_crystal`.

        Returns
        -------
        spectrum : RamanSpectrum2D
            Simulated Raman spectrum.

        See Also
        --------
        single_crystal : Simulate a single-crystal measurement.

        Notes
        -----
        See `single_crystal` for optional keyword arguments.
        """

        angles = np.arange(chi_start, chi_end + chi_step / 10.0, chi_step)

        if len(angles) < 2:
            raise ValueError(
                "Fewer than two angles between chi_start = {0:.2f} -> "
                "chi_end = {1:.2f} with chi_step = {2:.2f}."
                "".format(chi_start, chi_end, chi_step)
            )

        i_pol_str = str(i_pol).lower()
        s_pol_str = str(s_pol).lower()

        if i_pol_str == "rot":
            if s_pol_str == "rot":
                raise RuntimeError(
                    'i_pol and s_pol cannot both be set to "rot".'
                )

            i_pol = Polarisation.from_angles(geom.incident_direction, angles)

        if s_pol_str == "rot":
            s_pol = Polarisation.from_angles(geom.collection_direction, angles)

        return self.single_crystal(
            hkl,
            geom,
            i_pol,
            s_pol,
            d2_axis_vals=angles,
            d2_unit_text_label="chi / deg",
            d2_unit_plot_label=r"$\chi$ / $^\circ$",
            d2_col_hdrs=["chi_{0:.2f}".format(a) for a in angles],
            **kwargs
        )

    def single_crystal_crystal_rotation(
        self,
        hkl,
        geom,
        i_pol,
        s_pol,
        phi_start=0.0,
        phi_end=360.0,
        phi_step=2.5,
        rot_axis="incident",
        **kwargs
    ):
        """Simulate a single-crystal crystal (phi) rotation measurement
        where the crystal is rotated about the incident direction.

        Parameters
        ----------
        hkl : array_like of int
            Miller index of the crystal surface to orient antiparallel
            to the incident direction.
        geom : Geometry
            Measurement geometry.
        i_pol, s_pol : str or Polarisation
            Polarisations of the incident and scattered light. The
            scattered polarisation may also be specified by one of
            {"parallel", "cross", "sum"}.
        phi_start, phi_end, phi_step : float, optional
            Start/end angle and angle step for crystal rotation in
            degrees (defaults: 0 -> 360 deg in 2.5 deg steps).
        rot_axis : {"incident", "collection"}, optional
            Axis to rotate around (default: "incident")
        **kwargs : any
            Optional arguments to `single_crystal`.

        Returns
        -------
        spectrum : RamanSpectrum2D
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

        rot_axis = str(rot_axis).lower()

        if rot_axis == "incident":
            rot_axis = geom.incident_direction
        elif rot_axis == "collection":
            rot_axis = geom.collection_direction
        else:
            raise ValueError('Unknown rot_axis="{0}".'.format(rot_axis))

        angles = np.arange(phi_start, phi_end + phi_step / 10.0, phi_step)

        if len(angles) < 2:
            raise ValueError(
                "Fewer than two angles between phi_start = {0:.2f} -> "
                "phi_end = {1:.2f} with phi_step = {2:.2f}."
                "".format(phi_start, phi_end, phi_step)
            )

        rots = np.array(
            [rotation_matrix_from_axis_angle(rot_axis, a) for a in angles],
            dtype=np.float64,
        )

        return self.single_crystal(
            hkl,
            geom,
            i_pol,
            s_pol,
            rot=rots,
            d2_axis_vals=angles,
            d2_unit_text_label="phi / deg",
            d2_unit_plot_label=r"$\phi$ / $^\circ$",
            d2_col_hdrs=["phi_{0:.2f}".format(a) for a in angles],
            **kwargs
        )

    def powder(
        self,
        geom,
        i_pol,
        s_pol,
        po_hkl=None,
        po_eta=0.0,
        w=None,
        t=None,
        e_rt=None,
        lw=None,
        band_grp_inds=None,
        method="best",
        lc_prec=5,
        **kwargs
    ):
        """Simulate a powder Raman spectrum with optional preferred
        orientation.

        Parameters
        ----------
        geom : Geometry
            Measurement geometry.
        i_pol, s_pol : str, Polarisation or array_like of Polarisation
            Polarisations of the incident and scattered light. The
            scattered polarisation may also be specified by one of
            {"parallel", "cross", "sum"}.
        po_hkl : array_like of int or None, optional
            Miller index of the preferred orientation (default: None).
        po_eta : float, optional
            Fraction of crystallites with the preferred orientation
            (default: 0.0)
        w, t : float or None, optional
            Measurement wavelength and temperature.
        e_rt : float or None, optional
            Photon energy for evaluating the Raman tensors in energy-
            dependent Raman calculations (default: calculated from `w`
            if set and if energy-dependent Raman tensors are available,
            otherwise E = 0.)
        lw : float or None, optional
            Uniform linewidth or scale factor for calculated linewidths
            (defaults: 0.5 THz uniform linewidth or scale factor of
            1.0, depending on whether calculation has linewidths).
        band_grp_inds : array_like or None, optional
            Indices of "band groups" to include in the calculation
            (default: all groups). Groups are defined by irreps if
            available, or the band indices in the calculation otherwise.
        method : {"nquad", "leb+circ", "best"}, optional
            Method for powder averaging (default: `"best"`).
        lc_prec : int, optional
            Precision of the Lebedev + circle numerical quadrature
            scheme (default: 5).
        **kwargs : any
            Keyword arguments to the `RamanSpectrum1D` and
            `RamanSpectrum2D` constructors (some of these may be
            required depending on other parameters - see notes below).

        Returns
        -------
        spectrum : RamanSpectrum1D or RamanSpectrum2D
            Simulated Raman spectrum.

        See Also
        --------
        num_band_groups
            Number of band groups in the calculation.
        raman.intensity.calculate_powder_raman_intensities
            Lower-level API function used to calculate powder Raman
            intensities.
        raman.spectrum.RamanSpectrum1D, raman.spectrum.RamanSpectrum2D
            Objects returned by this function.

        Notes
        -----
        * The `w_rt` parameter is provided so that the intensity
          modulation "envelope" due to the laser wavelength and changes
          in intensity due energy-dependent polarisability can be
          modelled separately.
        * The indices in `band_grp_inds` correspond directly to the
          entries in the peak tables in the `RamanSpectrum` objects
          produced with the default `band_grp_inds=None`). The number of
          band groups can be obtained with the `num_band_groups`
          property.
        * If called with multiple polarisations to generate a 2D
          spectrum, additional keyword args to the `RamanSpectrum2D`
          constructor will need to be specified.
        """

        params = self._get_calc_params(
            geom, i_pol, s_pol, w, e_rt, lw, band_grp_inds
        )

        if method.lower() == "best":
            # method="best" will use an analytical formula if
            # possible. If multiple incident polarisations are
            # supplied and could result in a mix of analytical and
            # numerical methods being used, raise a warning.

            may_use_analytical = (
                np.abs(po_eta) < ZERO_TOLERANCE
                and not np.iscomplex(params["raman_tensors"]).any()
            )

            if may_use_analytical:
                count = 0

                c_dir = params["geometry"].collection_direction

                for p in params["incident_polarisations"]:
                    if p.check_perpendicular(c_dir):
                        count += 1

                if count < len(params["incident_polarisations"]):
                    warnings.warn(
                        "The given set of incident polarisatios may "
                        "result in a mix of calculations with "
                        "analytical and numerical methods with method "
                        '= "best". To avoid this warning, set method = '
                        '"lebedev+circle" instead.',
                        RuntimeWarning,
                    )

        po_surf_norm = None

        if po_hkl is not None:
            po_surf_norm = self._gamma_ph.structure.real_space_normal(
                po_hkl, conv=True
            )

        ints = np.zeros(
            (
                len(params["incident_polarisations"]),
                params["raman_tensors"].shape[0],
            ),
            dtype=np.float64,
        )

        for i, (i_pol, s_pol) in enumerate(
            zip(
                params["incident_polarisations"],
                params["scattered_polarisations"],
            )
        ):
            ints[i] = calculate_powder_raman_intensities(
                params["raman_tensors"],
                params["geometry"],
                i_pol,
                s_pol,
                po_eta=po_eta,
                po_surf_norm=po_surf_norm,
                method=method,
                lc_prec=lc_prec,
            )

        ints = np.array(ints, dtype=np.float64).T

        if params["is_2d_spectrum"]:
            if (
                "d2_axis_vals" not in kwargs
                or "d2_unit_text_label" not in kwargs
            ):
                raise RuntimeError(
                    "Multiple polarisations produce 2D spectra and "
                    "require the optional d2_axis_vals and "
                    "d2_unit_text_label keywords to be specified."
                )

            return RamanSpectrum2D(
                params["frequencies"],
                ints,
                params["linewidths"],
                kwargs.pop("d2_axis_vals"),
                kwargs.pop("d2_unit_text_label"),
                irreps=params["irreps"],
                w=w,
                t=t,
                **kwargs
            )
        else:
            return RamanSpectrum1D(
                params["frequencies"],
                ints,
                params["linewidths"],
                irreps=params["irreps"],
                w=w,
                t=t,
                **kwargs
            )

    def powder_polarisation_rotation(
        self,
        geom,
        i_pol,
        s_pol,
        chi_start=0.0,
        chi_end=360.0,
        chi_step=None,
        **kwargs
    ):
        """Simulate a powder polarisation (chi) rotation measurement
        where one of the incident or scattered polarisations are rotated
        about the incident or collection direction respectively.

        Parameters
        ----------
        geom : Geometry
            Measurement geometry.
        i_pol, s_pol : str or Polarisation
            Polarisations of the incident and scattered light. The
            polarisation to be rotated can be specified by "rot". The
            scattered polarisation may also be specified by one of
            {"parallel", "cross", "sum"}.
        chi_start, chi_end : float, optional
            Start/end angle for polarisation rotation in degrees
            (defaults: 0 -> 360 deg).
        chi_step : float or None, optional
            Angle step for polarisation rotation in degrees (defaults:
            2.5 deg, or 22.5 deg with preferred orientation).
        **kwargs : any
            Optional arguments to `powder`.

        Returns
        -------
        spectrum : RamanSpectrum2D
            Simulated Raman spectrum.

        See Also
        --------
        powder : Simulate a powder measurement.

        Notes
        -----
        See `powder` for optional keyword arguments.
        """

        if chi_step is None:
            chi_step = 2.5

            if "po_eta" in kwargs and kwargs["po_eta"] > 0.0:
                chi_step = 22.5

        angles = np.arange(chi_start, chi_end + chi_step / 10.0, chi_step)

        if len(angles) < 2:
            raise ValueError(
                "Fewer than two angles between chi_start = {0:.2f} -> "
                "chi_end = {1:.2f} with chi_step = {2:.2f}."
                "".format(chi_start, chi_end, chi_step)
            )

        i_pol_str = str(i_pol).lower()
        s_pol_str = str(s_pol).lower()

        if i_pol_str == "rot":
            if s_pol_str == "rot":
                raise RuntimeError(
                    'i_pol and s_pol cannot both be set to "rot".'
                )

            i_pol = Polarisation.from_angles(geom.incident_direction, angles)

        if s_pol_str == "rot":
            s_pol = Polarisation.from_angles(geom.collection_direction, angles)

        return self.powder(
            geom,
            i_pol,
            s_pol,
            d2_axis_vals=angles,
            d2_unit_text_label="chi / deg",
            d2_unit_plot_label=r"$\chi$ / $^\circ$",
            d2_col_hdrs=["chi_{0:.2f}".format(a) for a in angles],
            **kwargs
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

        return {
            "gamma_phonons": self._gamma_ph.to_dict(),
            "raman_tensors": self._r_t.to_dict(),
            "band_indices": self._band_inds.tolist(),
        }

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
        calc : RamanCalculator
            `RamanCalculator` object constructed from the data in `d`.
        """

        return RamanCalculation(
            GammaPhonons.from_dict(d["gamma_phonons"]),
            RamanTensors.from_dict(d["raman_tensors"]),
            d["band_indices"],
        )
