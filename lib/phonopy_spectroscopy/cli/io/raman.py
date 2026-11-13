# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Helper routines for input/output for Raman calculations."""

# -------
# Imports
# -------

import numpy as np

from scipy.interpolate import interp1d

from ...constants import ZERO_TOLERANCE
from ...interfaces.vasp_interface import dielectric_from_vasprun_xml

# ----------------
# Dielectric input
# ----------------


def fd_read_dielectrics_vasp(file_list, num_bands, num_steps, e_tol=0.1):
    """Read dielectric constants or functions for a finite-differences
    Raman calculation, using linear interpolation to obtain a common
    energy axis if required.

    Parameters
    ----------
    file_list : array_like
        Files to read input data from.
    num_bands, num_steps : int
        Number of bands and displacement steps in the calculation.
    e_tol : float, optional
        Allowed tolerance on differences in average energy axis spacing.

    Returns
    -------
    epss : numpy.ndarray
        Array of `DielectricFunction` objects (shape: `(N, M)`) for N
        bands and M displacement steps.

    Notes
    -----
    The files in `file_list` should be ordered as:
        <Band00001-Step01>
        <Band00001-Step02>
        ...
        <Band00001-StepMM>
        ...
        <BandNNNNN-Step01>
        <BandMMMMM-StepMM>
    """

    if len(file_list) != num_bands * num_steps:
        raise ValueError(
            "Number of files in file_list is not consistent with the "
            "number of bands and displacement steps."
        )

    epss = np.array(
        [dielectric_from_vasprun_xml(f) for f in file_list], dtype=object
    )

    e_ref = epss[0].energy

    for f, eps in zip(file_list, epss):
        e = eps.energy

        if len(e) != len(e_ref):
            # This almost certainly indicates that the individual
            # calculations were performed with different parameters and
            # is likely a (user) mistake.

            raise RuntimeError(
                f"Input file {f} : Number of energies {len(e):,} "
                f"differs from reference {len(e_ref):,}."
            )

    if len(e_ref) > 1:
        d_e_ref = np.mean(e_ref[1:] - e_ref[:-1])

        for f, eps in zip(file_list, epss):
            d_e = np.mean(eps.energy[1:] - eps.energy[:-1])

            if np.abs(d_e - d_e_ref) > e_tol:
                raise RuntimeError(
                    f"Input file {f} : Difference in energy spacing "
                    f"from reference is {d_e:.3f} eV > {d_e_ref:.3f} "
                    "eV."
                )

    return epss.reshape(num_bands, num_steps)
