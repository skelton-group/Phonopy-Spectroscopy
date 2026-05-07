# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Class to set up and post-process Raman tensor calculations."""


# -------
# Imports
# -------


import numpy as np
import warnings

from .calculation import RamanCalculation
from .tensors import RamanTensors

from ..constants import ZERO_TOLERANCE
from ..structure import Structure
from ..gamma_phonons import GammaPhonons, PolarGammaPhonons
from ..units import nm_to_ev

from ..utility.numpy_helper import (
    np_asarray_copy,
    np_readonly_view,
    np_check_shape,
    np_expand_dims,
)

from ..utility.differentiation import central_difference_coefficients


# ---------------------------------------------
# FiniteDisplacementRamanTensorCalculator class
# ---------------------------------------------


class FiniteDisplacementRamanTensorCalculator:
    """Set up and post-process Raman tensor calculations using the
    finite-displacement method."""

    def __init__(self, ph_calc, step_size=1.0e-2, prec=2, band_inds="active"):
        r"""Create a new instance of the
        `FiniteDisplacementRamanTensorCalculator` class.

        Parameters
        ----------
        ph_calc : GammaPhonons or PolarGammaPhonons
            Gamma-point phonon calculation.
        step_size : float
            Step size for displacements in sqrt(amu) Ang.
        prec : int, optional
            Precision of the central-difference scheme (defualt: 2).
        band_inds : {"all", "active"} or array_like, optional
            Band indices to perform calculations for: `"active"` to
            exclude the acoustic modes and Raman-inactive optic modes
            (requires irreps), `"all"` to include all bands, or a list
            of band indices (shape `(N,)`, default: `"active"`).
        """

        # The calculation of the Raman tensors should be performed using
        # the frequencies/eigenvectors without a non-analytical
        # correction applied. In this case, the \Gamma-point
        # eigenvectors should be real.

        if np.iscomplex(ph_calc.eigenvectors).any():
            raise ValueError("ph_calc must have real eigenvectors.")

        # Set up finite-difference steps.

        step_size = float(step_size)

        if step_size <= 0.0:
            raise ValueError("step_size must be > 0.")

        cd_steps, step_coeffs = central_difference_coefficients(1, prec)

        disp_steps = cd_steps * step_size

        # Identify band indices to use for calculation.

        band_inds_str = str(band_inds).lower()

        if band_inds_str == "all":
            band_inds = np.array(list(range(ph_calc.num_modes)), dtype=int)

        elif band_inds_str == "active":
            if ph_calc.has_irreps:
                band_inds = ph_calc.irreps.get_subset(
                    band_inds="raman", reset_inds=False
                ).band_indices_flat()
            else:
                warnings.warn(
                    'band_inds="raman" defaults to all bands apart '
                    "from the acoustic modes when the supplied phonon "
                    "calculation does not have irreps.",
                    RuntimeWarning,
                )

                band_inds = list(range(ph_calc.num_modes))

            # Exclude acoustic modes.

            mask = np.isin(
                band_inds, ph_calc.acoustic_mode_indices, invert=True
            )

            band_inds = band_inds[mask]

        else:
            band_inds = np_asarray_copy(band_inds, dtype=int)

            if not np_check_shape(band_inds, (None,)):
                raise ValueError(
                    "An explicit set of band indices must be suplied as an "
                    "array_like of shape (N,)."
                )

            if len(band_inds) != len(set(band_inds)):
                raise ValueError(
                    "band_inds contains one or more duplicate indices."
                )

            for idx in band_inds:
                if idx < 0 or idx >= ph_calc.num_modes:
                    raise ValueError(
                        "One or more indices in band_inds are "
                        "incompatible with the number of modes in the "
                        "phonon calculation."
                    )

        if len(band_inds) == 0:
            raise Exception(
                "No Raman-active modes in phonon calculation or no "
                "band_inds specified."
            )

        self._ph_calc = ph_calc
        self._band_inds = band_inds

        self._prec = prec
        self._step_size = step_size
        self._disp_steps = disp_steps
        self._step_coeffs = step_coeffs

    @property
    def phonon_calculation(self):
        """GammaPhonons or PolarGammaPhonons : Gamma-point phonon
        calculation."""
        return self._ph_calc

    @property
    def band_indices(self):
        """numpy.ndarray : Band indices in calculation."""
        return np_readonly_view(self._band_inds)

    @property
    def step_size(self):
        """float : Displacement step size."""
        return self._step_size

    @property
    def precision(self):
        """int : Finite-difference stencil precision."""
        return self._prec

    @property
    def displacement_steps(self):
        """numpy.ndarray : Displacement steps."""
        return np_readonly_view(self._disp_steps)

    @property
    def step_coefficients(self):
        """numpy.ndarray : Step coefficients."""
        return np_readonly_view(self._step_coeffs)

    @property
    def num_bands(self):
        """int : Number of bands in calculation."""
        return len(self._band_inds)

    @property
    def num_steps(self):
        """int : Number of displacement steps in the chosen central-
        difference scheme."""
        return len(self._disp_steps)

    def get_maximum_displacements(self):
        """Calculate and return an `(N, M)` array of the maximum atomic
        displacemnts for each of the selected bands and displacement
        steps.

        Returns
        -------
        max_disps : numpy.ndarray
            Maximum atomic displacements for each band and displacement
            step (shape: `(N, M)`).
        """

        edisps = self._ph_calc.eigendisplacements
        abs_disp_steps = np.abs(self._disp_steps)

        max_disps = np.zeros(
            (self.num_bands, self.num_steps), dtype=np.float64
        )

        for i, idx in enumerate(self._band_inds):
            max_disps[i] = (
                abs_disp_steps * np.linalg.norm(edisps[idx], axis=1).max()
            )

        return max_disps

    def generate_displaced_structures(self):
        """Generate and return displaced structures for the selected
        bands and displacement steps.

        Returns
        -------
        disp_structs : numpy.ndarray
            Array of `Structure` objects containing displaced structures
            (shape: `(N, M)`).
        """

        struct = self._ph_calc.structure
        at_pos = struct.cartesian_positions()

        edisps = self._ph_calc.eigendisplacements

        disp_structs = np.zeros(
            (len(self._band_inds), len(self._disp_steps)), dtype=object
        )

        for i, idx in enumerate(self._band_inds):
            for j, step in enumerate(self._disp_steps):
                disp_structs[i, j] = Structure(
                    struct.lattice_vectors,
                    at_pos + step * edisps[idx],
                    struct.atom_types,
                    struct.atomic_masses,
                    cart_to_frac=True,
                )

        return disp_structs

    def calculate_raman_tensors(
        self, dielectrics, e=None, e_cut=None, w_cut=None
    ):
        """Calculate Raman tensors from dielectric tensors or
        energy-dependent dielectric functions evaluated for the
        displaced structures generated by
        `generate_displaced_structures()` and return a
        `RamanCalculation` object.

        Parameters
        ----------
        dielectrics : array_like
            Dielectric tensors or energy-dependent dielectric functions
            (shape: `(N, M, 3, 3)` or `(N, M, O, 3, 3)`).
        e : array_like or None, optional
            array_like with the energies at which the dielectirc
            constants/functions were evaluated (eV, shape: `(O,)`).
            `e` must be specified for energy-dependent dielectric
            functions, and defaults to E = 0 for single tensors.
        e_cut, w_cut : float or None, optional
            Energy or wavelength cutoff for energy-dependent dielectric
            functions (default: `None` = no cutoff). Requires that `e`
            is specified. If both are supplied, `w_cut` overrides
            `e_cut`.

        Returns
        -------
        calc : RamanCalculator
            A `RamanCalculator` object generated from the
            `GammaPhonons` object and band indices used to set up the
            calculation, and a `RamanTensors` object constructed from
            the dielectrics data.
        """

        dielectrics, _ = np_expand_dims(
            np.asarray(dielectrics),
            (len(self._band_inds), len(self._disp_steps), None, 3, 3),
        )

        if e is not None:
            e = np.asarray(e)

            if len(e) != dielectrics.shape[2]:
                raise ValueError(
                    "If supplied, e must be an array_like with shape (O,)."
                )

        # The paper this is based on (Porezag and Pederson, Phys. Rev. B
        # 54, 7830, 1996, 10.1103/PhysRevB.54.7830) uses "polarisability
        # volumes" in Ang^3, and the prefactor of V / (4 \pi) is used in
        # the vasp_raman.py code (Fonari and Stauffer,
        # https://github.com/raman-sc/VASP/).

        # With the dielectrics in units of relative permittivity,
        # scaling by V should convert to polarisability volume, and
        # the factor of 4 \pi comes from the relation:

        # \alpha = 4 \pi \eps_0 \alpha' ->
        #   \alpha' = (1 / 4 \pi) * (\alpha / \eps_0).

        r_t = (
            (
                dielectrics
                * self._step_coeffs[
                    np.newaxis, :, np.newaxis, np.newaxis, np.newaxis
                ]
            ).sum(axis=1)
            * (1.0 / self._step_size)
            * (self._ph_calc._struct.volume() / (4.0 * np.pi))
        )

        # If required, impose an energy cutoff on energy-dependent
        # dielectric functions.

        if e is not None and (e_cut is not None or w_cut is not None):
            if w_cut is not None:
                e_cut = nm_to_ev(w_cut)

            mask = e <= e_cut

            if mask.sum() == 0:
                raise RuntimeError(
                    "The supplied e_cut/w_cut removes all the energies "
                    "in e."
                )

            e = e[mask]
            r_t = r_t[:, mask, :, :]

        return RamanCalculation(
            self._ph_calc, RamanTensors(r_t, e), self._band_inds
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
            "phonon_calculation": self._ph_calc.to_dict(),
            "band_indices": self._band_inds.tolist(),
            "step_size": self._step_size,
            "precision": self._prec,
            "displacement_steps": self._disp_steps.tolist(),
            "step_coefficients": self._step_coeffs.tolist(),
        }

    @staticmethod
    def from_dict(d):
        """Create a new `FiniteDisplacementRamanTensorCalculator`
        instance from a dictionary generated by
        `FiniteDisplacementRamanTensorCalculator.to_dict()`.

        Parameters
        ----------
        d : dict
            Dictionary generated by `to_dict()`.

        Returns
        -------
        calc : FiniteDisplacementRamanTensorCalculator
            `FiniteDisplacementRamanTensorCalculator` object constructed
            from the data in `d`.
        """

        data = d["phonon_calculation"]

        ph_calc = None

        if "eps_inf" in data and "born_charges" in data:
            ph_calc = PolarGammaPhonons.from_dict(data)
        else:
            ph_calc = GammaPhonons.from_dict(data)

        fd_calc = FiniteDisplacementRamanTensorCalculator(
            ph_calc,
            prec=d["precision"],
            step_size=d["step_size"],
            band_inds=d["band_indices"],
        )

        # The displacement_steps and step_coefficients properties are
        # generated based on the prec and stpe_size parameters, so check
        # for consistency with the data in the dictionary.

        diff = (
            np.array(d["displacement_steps"], dtype=np.float64)
            - fd_calc.displacement_steps
        )

        if (np.abs(diff) > ZERO_TOLERANCE).any():
            raise Exception(
                "displacement_steps property inconsistent with the"
                '"displacement_steps" key in the supplied dictionary '
                "(this is likely a bug)."
            )

        diff = (
            np.array(d["step_coefficients"], dtype=np.float64)
            - fd_calc.step_coefficients
        )

        if (np.abs(diff) > ZERO_TOLERANCE).any():
            raise Exception(
                "step_coefficients property inconsistent with the"
                '"step_coefficients" key in the supplied dictionary '
                "(this is likely a bug)."
            )

        return fd_calc
