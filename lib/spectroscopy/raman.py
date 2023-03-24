# spectroscopy/raman.py


# ---------
# Docstring
# ---------

""" Core routines for setting up and post-processing Raman-activity
calculations. """


# -------
# Imports
# -------

import math

import numpy as np

from spectroscopy import constants, utilities


# -----------
# Preparation
# -----------

def generate_displaced_structures(
        structure, eigendisplacement, step_size, num_steps,
        step_type='delta_q', omit_zero=False
        ):
    """ Generate a set of structures displaced along a mode
    eigendisplacement for a Raman activity calculation.

    Arguments:
        structure -- a tuple of (lattice_vectors, atomic_symbols,
            atom_positions) specifying the structure to displace.
        eigendisplacement -- mode eigendisplacement (N x three-component
            vectors).
        step_size -- size of the displacement step.
        num_steps -- number of displacement steps to make (minimum 2).

    Keyword arguments:
        step_type -- step specified by step_size argument ('delta_q' -
            normal-mode coordinate, default; 'max_r' - max. Cartesian
            displacement; 'norm_x' - units of normalised
            eigendisplacement).
        omit_zero -- if True, and if num_steps is odd, omit the zero
            step.

    Notes:
        If num_steps is even, displacements will be generated for
        +/- step_size, +/- 2 * step_size, etc.
        If num_steps is odd, the set of displacements will include a
        zero step (overridden by the omit_zero argument).
        If step_type = 'max_r' is used, the max. displacement applies to
        +/- step_size (i.e. if num_steps > 3, some steps will be larger
        than this).
    """

    lattice_vectors, atomic_symbols, atom_positions = structure

    # Using np.shape() should allow the parameters to be supplied as
    # (correctly-dimensioned) lists.

    dim_1, dim_2 = np.shape(lattice_vectors)

    if dim_1 != 3 or dim_2 != 3:
        raise Exception(
            "Error: The lattice vectors provided in structure "
            "must be a 3x3 matrix."
            )

    num_atoms = len(atomic_symbols)

    if len(atom_positions) != num_atoms:
        raise Exception(
            "Error: The lengths of the atomic symbol and atom position "
            "lists provided in structure are inconsistent.")

    dim_1, dim_2 = np.shape(eigendisplacement)

    if dim_1 != num_atoms or dim_2 != 3:
        raise Exception(
            "Error: eigendisplacement must be an Nx3 matrix.")

    if step_size <= 0.0:
        raise Exception("Error: step_size must be positive.")

    if num_steps < 2:
        raise Exception("Error: num_steps must be at least 2.")

    if step_type not in ['delta_q', 'max_r', 'norm_x']:
        raise Exception(
            "Error: Unsupported stepType '{0}'.".format(step_type))

    # Displacements are created at 0, +/- stepsize, +/- (2 * stepsize),
    # ... depending on numsteps and whether omitzero is set.

    disp_steps = []

    i_max = num_steps // 2

    for i in range(-i_max, i_max + 1):
        if i == 0 and (omit_zero or num_steps % 2 == 0):
            continue

        disp_steps.append(i * step_size)

    # Convert atom positions to Cartesian coordinates.

    atom_positions = utilities.fractional_to_cartesian_coordinates(
        atom_positions, lattice_vectors)

    # Displace structures along eigenvectors.

    disp_structures, max_disps = [], []

    # Calculate the maximum (Cartesian) offset in the eigendisplacement.

    max_disp = max(np.linalg.norm(v) for v in eigendisplacement)

    # Calculate the norm of the eigendisplacement.

    eigendisplacement_norm = np.linalg.norm(eigendisplacement)

    # Calculate the norm of the eigendisplacement.

    for i, disp_step in enumerate(disp_steps):
        if step_type == 'max_r':
            # If the step is defined as a maximum Cartesian
            # displacement, scale it by max_disp.

            disp_step = disp_step / max_disp
            disp_steps[i] = disp_step
        elif step_type == 'norm_x':
            # If the step is defined in terms of the normalised
            # eigendisplacement, scale it by eigendisplacementnorm.

            disp_step = disp_step / eigendisplacement_norm
            disp_steps[i] = disp_step

        # Generate displaced positions.

        new_atom_positions = [
            np.add(pos, np.multiply(disp_step, v))
                for pos, v in zip(atom_positions, eigendisplacement)
            ]

        # Convert new positions back to fractional coordinates when
        # storing.

        disp_structures.append(
            (lattice_vectors, atomic_symbols,
             utilities.cartesian_to_fractional_coordinates(
                new_atom_positions, lattice_vectors))
            )

        max_disps.append(disp_step * max_disp)

    # Return displacement steps, displaced structures and max.
    # Cartesian displacements.

    return disp_steps, disp_structures, max_disps


# --------
# Analysis
# --------

def calculate_raman_activity_tensor(
        eps_tensors, q_disps, cell_volume, fit_type
        ):
    """ Calculate the Raman activity from a set of dielectric tensors
    obtained at different displacements along the normal-mode coordinate
    Q.

    Arguments:
        eps_tensors -- a set of calculated dielectric tensors at
            displacements along the normal-mode coordinate.
        q_disps -- normal-mode displacements at which the dielectric
            tensors were calculated.
        cell_volume -- volume of the unit cell (required to compute a
            normalisation factor).
        fit_type -- fitting type to use to compute the derivative
            (currently only 'finite_difference' is supported).

    Notes:
        At present, this function only supports two-point finite
        difference stencils for computing the derivative.
    """

    if fit_type not in ['finite_difference']:
        raise Exception(
            "Error: Unsupported fit_type '{0}'.".format(fit_type)
            )

    if len(eps_tensors) != len(q_disps):
        raise Exception(
            "Error: The numbers of supplied dielectric tensors "
            "and displacement steps are inconsistent."
            )

    for eps_tensor in eps_tensors:
        dim_1, dim_2 = np.shape(eps_tensor)

        if dim_1 != 3 or dim_2 != 3:
            raise Exception(
                "Error: Dielectric tensors must be 3x3 matrices.")

    # Make sure the dielectric tensors and displacements are ordered by
    # step.

    temp = sorted(
        [item for item in zip(q_disps, eps_tensors)],
        key=lambda item: item[0]
        )

    q_disps = [q_disp for q_disp, _ in temp]
    eps_tensors = [eps_tensor for _, eps_tensor in temp]

    if fit_type == 'finite_difference':
        if len(eps_tensors) != 2:
            raise Exception(
                "Error: This function currently only supports "
                "two-point finite-difference stencils.")

        # Calculate and verify the step size.

        disps_temp = sorted(
            [0.0] + [disp for disp in q_disps]
            )

        step_size = np.average(
            np.subtract(disps_temp[1:], disps_temp[:-1])
            )

        for q_disp in q_disps:
            if math.fabs(q_disp / step_size) < constants.ZERO_TOLERANCE:
                raise Exception(
                    "Error: Finite-dfference stencils require a "
                    "uniform step size.")

        # Calculate the Raman tensor.

        eps_1, eps_2 = eps_tensors

        # The two-point central-difference formula is
        # f'(x) ~= [f(x + h) - f(x - h)] / 2h.

        return (
            (cell_volume / (4.0 * math.pi))
            * np.subtract(eps_2, eps_1)
            / (2.0 * step_size)
            )

def scalar_average_raman_activity_tensor(raman_tensor):
    """ Average a Raman-activity tensor to obtain a scalar
    intensity. """

    # This formula came from D. Porezag and M. R. Pederson, Phys. Rev.
    # B: Condens. Matter Mater. Phys., 1996, 54, 7830.

    alpha = (
        (raman_tensor[0][0] + raman_tensor[1][1] + raman_tensor[2][2])
        / 3.0)

    beta_squared = 0.5 * (
        (raman_tensor[0][0] - raman_tensor[1][1]) ** 2
        + (raman_tensor[0][0] - raman_tensor[2][2]) ** 2
        + (raman_tensor[1][1] - raman_tensor[2][2]) ** 2
        + 6.0 * (raman_tensor[0][1] ** 2 + raman_tensor[0][2] ** 2 +
            raman_tensor[1][2] ** 2)
        )

    return 45.0 * alpha ** 2 + 7.0 * beta_squared
