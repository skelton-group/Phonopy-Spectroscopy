# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for modelling the dielectric functions of multiphase
mixtures."""

# -------
# Imports
# -------

import warnings

import numpy as np

from .response_function import DielectricFunction

from ..constants import ZERO_TOLERANCE

# ---------
# Functions
# ---------


def effective_dielectric_function(
    epss, fracs, shapes=None, alpha=0.5, tol=1.0e-5, max_iter=50
):
    r"""Determine the effective dielectric function or a multiphase
    mixture by solving the self-consistent Bruggeman equation.

    Parameters
    ----------
    epss : array_like of DielectricFunctions
        Dielectric functions of components (shape: `(N,)`; must be 3D).
    fracs : array_like
        Volume fractions (must sum to 1).
    shapes : array_like, optional
        Shape tensors (shape: `(3, 3)` or `(N, 3, 3)`).
    alpha : float, optional
        Relaxation weight for updating \eps_eff (default: 0.5).
    tol : float, optional
        Convergence threshold (default: 1.0e-5).
    max_iter : int, optional
        Maximum number of iterations (default: 50).

    Returns
    -------
    eps_eff : DielectricFunction
        Effective dielectric function.

    Notes
    -----
    The effective dielectric tensor :math:`\varepsilon_\mathrm{eff}`
    obeys the Bruggeman equation:

    .. math::

        \sum_i f_i
          (\varepsilon_i - \varepsilon_\mathrm{eff})
          [ I + L_i \cdot \varepsilon_\mathrm{eff}^{-1} \cdot (\varepsilon_i - \varepsilon_\mathrm{eff}) ]^{-1}
        = \sum_i f_i (\varepsilon_i - \varepsilon_\mathrm{eff}) A_i
        = 0

    This routine finds :math:`\varepsilon_\mathrm{eff}` from the
    self-consistent solution:

    .. math::

        \varepsilon_\mathrm{eff}^{n + 1} =
          \left[ f_i \varepsilon_i A_i(\varepsilon_\mathrm{eff}^n) \right]
          \left[ f_i A_i(\varepsilon_\mathrm{eff}^n) \right]^{-1}

    The optional `shapes` parameter is used to specify the shape tensors
    :math:`L_i`, which define the aspect ratio of the inclusions of each
    phase :math:`i`. Shape tensors must be "transpose" symmetric
    (:math:`L_i = L_i^\mathrm{T}`) with :math:`\mathrm{Tr}(L_i) = 1`),
    and all three eigenvalues must be between 0 and 1. If `shapes=None`
    the routine assumes spherical inclusions for all phases. If `shapes`
    specifies a single tensor it is used for all phases.

    The parameter `alpha` controls the update of :math:`\epsilon_eff`
    according to:

    .. math::

        \epsilon_\mathrm{eff}
          = (1 - \alpha) \epsilon_\mathrm{eff}^{n - 1} + \alpha \epsilon_\mathrm{eff}^{n}

    The self-consistent loop is assumed to be converged once the maximum
    element-wise change in :math:`\varepsilon_\mathrm{eff}` at all
    energies is below the threshold specified in `tol` or after at most
    `max_iter` updates, whichever occurs sooner.
    """

    epss = np.atleast_1d(np.asanyarray(epss, dtype=object))

    if epss.ndim > 1 or len(epss) == 0:
        raise ValueError("epss must be an array_like with shape (N,).")

    for eps in epss:
        if eps.num_dims != 3:
            raise ValueError("All dielectric functions must be 3D.")

    fracs = np.atleast_1d(np.asanyarray(fracs, dtype=np.float64))

    if len(fracs) != len(epss):
        raise ValueError("fracs must be an array_like with shape (N,).")

    if not np.isclose(fracs.sum(), 1.0, atol=ZERO_TOLERANCE):
        raise ValueError("fracs must sum to 1.")

    if shapes is not None:
        shapes = np.atleast_2d(np.asanyarray(shapes, dtype=np.float64))

        match shapes.shape:
            case (3, 3):
                # Broadcast to all dielectric functions.
                shapes = np.broadcast_to(shapes, (len(epss), 3, 3))
            case (len(epss), 3, 3):
                pass
            case _:
                raise ValueError(
                    "shapes must be an array_like with shape (3, 3) or "
                    "(N, 3, 3)."
                )

        if not np.allclose(
            shapes, shapes.swapaxes(-2, -1), atol=ZERO_TOLERANCE
        ):
            raise ValueError("shape tensors must be symmetric.")

        if not np.allclose(
            np.trace(shapes, axis1=-2, axis2=-1), 1.0, atol=ZERO_TOLERANCE
        ):
            raise ValueError("shape tensors must have a trace of 1.")

        evals = np.linalg.eigvalsh(shapes)

        if (evals < -1.0 * ZERO_TOLERANCE).any() or (
            (evals > 1.0 + ZERO_TOLERANCE).any()
        ):
            raise ValueError(
                "shape tensors must have eigenvalues between 0 and 1."
            )
    else:
        # Default to spherical shape tensors.
        shapes = np.broadcast_to(np.eye(3, 3) / 3.0, (len(epss), 3, 3))

    if tol <= 0.0:
        raise ValueError("tol must be > 0.")

    if max_iter <= 0:
        raise ValueError("max_iter must be > 0.")

    if len(epss) == 1:
        # Single dielectric function -> nothing to do.
        return epss[0]

    # Align energies of dielectric functions.

    epss = DielectricFunction.align_energies(epss, e=None)

    # All dielectric functions should now have the same energies.

    e = epss[0].energy

    fracs = fracs[:, np.newaxis, np.newaxis]
    ident = np.identity(3, dtype=np.float64)[np.newaxis, :, :]

    eps_eff = np.zeros((len(e), 3, 3), dtype=np.complex128)

    global_max_diff = 0.0
    global_num_non_conv = 0

    # Temporary array to store \eps tenosrs of components at current
    # energy.

    eps_cmp = np.zeros((len(epss), 3, 3), dtype=np.complex128)

    for idx in range(len(e)):
        # Update eps_temp with component tensors at current energy.

        for i, eps in enumerate(epss):
            eps_cmp[i] = eps.epsilon[idx]

        # For first energy, set initial guess for eps_eff to weighted
        # average of component tensors. For subsequent energies, use the
        # previous converged eps_eff as the initial guess.

        if idx == 0:
            eps_eff[idx] = (fracs * eps_cmp).sum(axis=0)
        else:
            eps_eff[idx] = eps_eff[idx - 1]

        # Track convergence at current energy.

        max_diff = None

        for n in range(max_iter):
            eps_eff_old = eps_eff[idx]

            # Compute A_i.

            a = np.linalg.inv(
                ident
                + (
                    shapes
                    @ np.linalg.inv(eps_eff_old)[np.newaxis, :, :]
                    @ (eps_cmp - eps_eff_old[np.newaxis, :, :])
                )
            )

            # Compute updated eps_eff.

            eps_eff_new = (
                (fracs * (eps_cmp @ a)).sum(axis=0)
            ) @ np.linalg.inv((fracs * a).sum(axis=0))

            # Mix old and new eps_eff.

            eps_eff[idx] = (1.0 - alpha) * eps_eff_old + (alpha * eps_eff_new)

            # Update max_diff and check convergence against tol.

            max_diff = np.abs(eps_eff_new - eps_eff_old).max()

            if max_diff < tol:
                break

        global_max_diff = max(max_diff, global_max_diff)

        if max_diff >= tol:
            global_num_non_conv += 1

    if global_max_diff > tol:
        warnings.warn(
            f"Convergence at {global_num_non_conv} energies not "
            f"achieved in max_iter = {max_iter}: max_diff = "
            f"{max_diff:.3e} > tol = {tol:.3e}.",
            RuntimeWarning,
        )

    return DielectricFunction(eps_eff, e, e_units=epss[0].energy_units)
