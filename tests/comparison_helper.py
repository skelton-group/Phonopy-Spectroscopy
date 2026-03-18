# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Routines for comparing objects as part of unit tests."""


# -------
# Imports
# -------


import numpy as np

from phonopy_spectroscopy.phonon import PolarGammaPhonons


# ---------
# Functions
# ---------


def compare_structures(struct_cmp, struct_ref):
    """Compare two `Structure` objects and determine whether they hold
    identical data.

    Parameters
    ----------
    struct_cmp, struct_ref : Structure
        `Structure` objects to compare.

    Returns
    -------
    equiv : bool
        `True` if `struct_cmp` is equivalent to `struct_ref`, otherwise
        `False`.
    """

    return (
        np.allclose(struct_cmp.lattice_vectors, struct_ref.lattice_vectors)
        and np.allclose(struct_cmp.atom_positions, struct_ref.atom_positions)
        and np.equal(struct_cmp.atom_types, struct_ref.atom_types).all()
        and np.allclose(struct_cmp.atomic_masses, struct_ref.atomic_masses)
        and np.allclose(
            struct_cmp.conventional_transformation_matrix,
            struct_ref.conventional_transformation_matrix,
        )
    )


def compare_irreps(irreps_cmp, irreps_ref):
    """Compare two `Irreps` objects and determine
    whether they hold identical data.

    Parameters
    ----------
    irreps_cmp, irreps_ref : Irreps
        `Irreps` objects to compare.

    Returns
    -------
    equiv : bool
        `True` if `irreps_cmp` is equivalent to `irreps_ref`,
        otherwise `False`.
    """

    if irreps_cmp.point_group != irreps_ref.point_group:
        return False

    if len(irreps_cmp.irrep_symbols) != len(irreps_ref.irrep_symbols):
        return False

    # IrreducibleRepresentations constructor should enforce that
    # the irrep symbols and irrep_band_indices properties are the same
    # length.

    for sym_cmp, sym_ref in zip(
        irreps_cmp.irrep_symbols, irreps_ref.irrep_symbols
    ):
        if sym_ref != sym_cmp:
            return False

    for band_inds_cmp, band_inds_ref in zip(
        irreps_cmp.irrep_band_indices, irreps_ref.irrep_band_indices
    ):
        if (band_inds_cmp != band_inds_ref).any():
            return False

    return True


def compare_gamma_phonons(gamma_ph_cmp, gamma_ph_ref):
    """Compare two `GammaPhonons` objects and determine whether they
    hold identical data.

    Parameters
    ----------
    gamma_ph_cmp, gamma_ph_ref : GammaPhonons
        `GammaPhonons` objects to compare.

    Returns
    -------
    equiv : bool
        `True` if `gamma_ph_cmp` is equivalent to `gamma_ph_ref`,
        otherwise `False`.
    """

    if not compare_structures(gamma_ph_cmp.structure, gamma_ph_ref.structure):
        return False

    if not np.allclose(gamma_ph_cmp.frequencies, gamma_ph_ref.frequencies):
        return False

    if not np.allclose(gamma_ph_cmp.eigenvectors, gamma_ph_ref.eigenvectors):
        return False

    if gamma_ph_cmp.has_linewidths:
        if not gamma_ph_ref.has_linewidths or not np.allclose(
            gamma_ph_cmp.linewidths, gamma_ph_ref.linewidths
        ):
            return False

    if gamma_ph_cmp.has_irreps:
        if not gamma_ph_ref.has_irreps or not compare_irreps(
            gamma_ph_cmp.irreps, gamma_ph_ref.irreps
        ):
            return False
    return True


def compare_polar_gamma_phonons(gamma_ph_cmp, gamma_ph_ref):
    """Compare two `PolarGammaPhonons` objects and determine whether
    they hold identical data.

    Parameters
    ----------
    gamma_ph_cmp, gamma_ph_ref : PolarGammaPhonons
        `PolarGammaPhonons` objects to compare.

    Returns
    -------
    equiv : bool
        `True` if `gamma_ph_cmp` is equivalent to `gamma_ph_ref`,
        otherwise `False`.
    """

    if not compare_gamma_phonons(gamma_ph_cmp, gamma_ph_ref):
        return False

    return np.allclose(
        gamma_ph_cmp.epsilon_inf, gamma_ph_ref.epsilon_inf
    ) and np.allclose(
        gamma_ph_cmp.born_effective_charges,
        gamma_ph_ref.born_effective_charges,
    )


def compare_infrared_calculations(calc_cmp, calc_ref):
    """Compare two `InfraredCalculation` objects and determine whether
    they hold identical data.

    Parameters
    ----------
    calc_cmp, calc_ref : InfraredCalculation
        `InfraredCalculation` objects to compare.

    Returns
    -------
    equiv : bool
        `True` if `calc_cmp` is equivalant to `calc_ref`, otherwise
        `False`.
    """

    return compare_gamma_phonons(
        calc_ref.phonon_calculation, calc_cmp.phonon_calculation
    )


def compare_finite_displacement_raman_tensor_calculators(
    fd_calc_cmp, fd_calc_ref
):
    """Compare two `FiniteDisplacementRamanTensorCalculator` objects and
    determine whether they hold identical data.

    Parameters
    ----------
    fd_calc_cmp, fd_calc_ref : FiniteDisplacementRamanTensorCalculator
        `FiniteDisplacementRamanTensorCalculator` objects to compare.

    Returns
    -------
    equiv : bool
        `True` if `fd_calc_cmp` is equivalant to `fd_calc_ref`,
        otherwise `False`.
    """

    if not compare_gamma_phonons(
        fd_calc_ref.phonon_calculation, fd_calc_cmp.phonon_calculation
    ):
        return False

    return (
        np.equal(fd_calc_ref.band_indices, fd_calc_cmp.band_indices).all()
        and np.allclose(
            fd_calc_ref.displacement_steps, fd_calc_cmp.displacement_steps
        )
        and np.allclose(
            fd_calc_ref.step_coefficients, fd_calc_cmp.step_coefficients
        )
    )


def compare_raman_tensors(r_t_cmp, r_t_ref):
    """Compare two `RamanTensors` objects and determine whether they
    hold identical data.

    Parameters
    ----------
    r_t_cmp, r_t_ref : RamanTensors
        `RamanTensors` objects to compare.

    Returns
    -------
    equiv : bool
        `True` if `r_t_cmp` is equivalant to `r_t_ref`, otherwise
        `False`.
    """

    return (
        np.allclose(r_t_cmp.energies, r_t_ref.energies),
        np.allclose(r_t_cmp.raman_tensors, r_t_ref.raman_tensors),
    )


def compare_raman_calculations(calc_cmp, calc_ref):
    """Compare two `RamanCalculation` objects and determine whether they
    hold identical data.

    Parameters
    ----------
    calc_cmp, calc_ref : RamanCalculation
        `RamanCalculation` objects to compare.

    Returns
    -------
    equiv : bool
        `True` if `calc_cmp` is equivalant to `calc_ref`, otherwise
        `False`.
    """

    # The RamanCalculation object can be initialised with either of a
    # GammaPhonons or PolarGammaPhonons object.

    if isinstance(calc_cmp.phonon_calculation, PolarGammaPhonons):
        if isinstance(calc_ref.phonon_calculation, PolarGammaPhonons):
            if not compare_polar_gamma_phonons(
                calc_cmp.phonon_calculation, calc_ref.phonon_calculation
            ):
                return False
        else:
            return False
    else:
        if isinstance(calc_ref.phonon_calculation, PolarGammaPhonons):
            return False

        if not compare_gamma_phonons(
            calc_cmp.phonon_calculation, calc_ref.phonon_calculation
        ):
            return False

    return (
        compare_raman_tensors(calc_cmp.raman_tensors, calc_ref.raman_tensors)
        and np.equal(calc_cmp.band_indices, calc_ref.band_indices).all()
        and compare_irreps(calc_cmp.irreps, calc_ref.irreps)
    )
