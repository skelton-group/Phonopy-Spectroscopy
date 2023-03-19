# SpectroscoPy/Raman.py


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

from SpectroscoPy import Utilities


# -----------
# Preparation
# -----------

def generatedisplacedstructures(structure, eigendisplacement, stepsize,
                                numsteps, steptype='delta_q',
                                omitzero=False
                                ):
    """
    Generate a set of structures displaced along a mode eigendisplacement for
    a Raman activity calculation.

    Arguments:
        structure -- a tuple of (lattice_vectors, atomic_symbols,
        atom_positions) specifying the structure to displace.
        eigendisplacement -- mode eigendisplacement (N x three-component
        vectors).
        stepSize -- size of the displacement step.
        numSteps -- number of displacement steps to make (minimum 2).

    Keyword arguments:
        stepType -- step specified by stepSize argument ('delta_q' -
        normal-mode coordinate, default; 'max_r' - max. Cartesian displacement;
        'norm_x' - units of normalised eigendisplacement).
        omitZero -- if True, and if numSteps is odd, omit the zero step.

    Notes:
        If numSteps is even, displacements will be generated for +/- stepSize,
        +/- 2 * stepSize, etc.
        If numSteps is odd, the set of displacements will include a zero step
        (overridden by the omitZero argument).
        If stepType = 'max_r' is used, the max. displacement applies to +/-
        stepSize (i.e. if numSteps > 3, some steps will be larger than this).
    """

    latticevectors, atomicsymbols, atompositions = structure

    # Using np.shape() should allow the parameters to be supplied as
    # (correctly-dimensioned) lists.

    dim1, dim2 = np.shape(latticevectors)

    if dim1 != 3 or dim2 != 3:
        raise Exception("Error: The lattice vectors provided in structure "
                        "must be a 3x3 matrix."
                        )

    numatoms = len(atomicsymbols)

    if len(atompositions) != numatoms:
        raise Exception("Error: The lengths of the atomic symbol and atom "
                        "position lists provided in structure are "
                        "inconsistent."
                        )

    dim1, dim2 = np.shape(eigendisplacement)

    if dim1 != numatoms or dim2 != 3:
        raise Exception("Error: eigendisplacement must be an Nx3 matrix.")

    if stepsize <= 0.0:
        raise Exception("Error: stepSize must be positive.")

    if numsteps < 2:
        raise Exception("Error: numSteps must be at least 2.")

    if steptype not in ['delta_q', 'max_r', 'norm_x']:
        raise Exception("Error: Unsupported stepType '{0}'.".format(steptype))

    # Displacements are created at 0, +/- stepsize, +/- (2 * stepsize), ...
    # depending on numsteps and whether omitzero is set.

    dispsteps = []

    imax = numsteps // 2

    for i in range(-imax, imax + 1):
        if i == 0 and (omitzero or numsteps % 2 == 0):
            continue

        dispsteps.append(i * stepsize)

    # Convert atom positions to Cartesian coordinates.

    atompositions = Utilities.fractionaltocartesiancoordinates(
        atompositions, latticevectors)

    # Displace structures along eigenvectors.

    dispstructures, maxdisps = [], []

    # Calculate the maximum (Cartesian) offset in the eigendisplacement.

    maxdisp = max(
        np.linalg.norm(v) for v in eigendisplacement
        )

    # Calculate the norm of the eigendisplacement.

    eigendisplacementnorm = np.linalg.norm(eigendisplacement)

    # Calculate the norm of the eigendisplacement.

    for i, dispstep in enumerate(dispsteps):
        if steptype == 'max_r':
            # If the step is defined as a maximum Cartesian displacement,
            # scale it by maxdisp.

            dispstep = dispstep / maxdisp
            dispsteps[i] = dispstep
        elif steptype == 'norm_x':
            # If the step is defined in terms of the normalised
            # eigendisplacement, scale it by eigendisplacementnorm.

            dispstep = dispstep / eigendisplacementnorm
            dispsteps[i] = dispstep

        # Generate displaced positions.

        newatompositions = [
            np.add(position, np.multiply(dispstep, v))
            for position, v in zip(atompositions, eigendisplacement)
            ]

        # Convert new positions back to fractional coordinates when storing.

        dispstructures.append(
            (latticevectors, atomicsymbols,
             Utilities.cartesiantofractionalcoordinates(newatompositions,
                                                        latticevectors))
            )

        maxdisps.append(dispstep * maxdisp)

    # Return displacement steps, displaced structures and max.
    # Cartesian displacements.

    return dispsteps, dispstructures, maxdisps


# --------
# Analysis
# --------

def calculateramanactivitytensor(epstensors, qdisps, cellvolume, fittype):
    """
    Calculate the Raman activity from a set of dielectric tensors obtained at
    different displacements along the normal-mode coordinate Q.

    Arguments:
        epsTensors -- a set of calculated dielectric tensors at displacements
        along the normal-mode coordinate.
        qDisps -- normal-mode displacements at which the dielectric tensors
        were calculated.
        cellVolume -- volume of the unit cell (required to compute a
        normalisation factor).
        fitType -- fitting type to use to compute the derivative (currently
        only 'finite_difference' is supported).

    Notes:
        At present, this function only supports two-point finite difference
        stencils for computing the derivative.
    """

    if fittype not in ['finite_difference']:
        raise Exception("Error: Unsupported fitType '{0}'.".format(fittype))

    if len(epstensors) != len(qdisps):
        raise Exception("Error: The numbers of supplied dielectric tensors "
                        "and displacement steps are inconsistent."
                        )

    for epstensor in epstensors:
        dim1, dim2 = np.shape(epstensor)

        if dim1 != 3 or dim2 != 3:
            raise Exception("Error: Dielectric tensors must be 3x3 matrices.")

    # Make sure the dielectric tensors and displacements are ordered by step.

    temp = [item for item in zip(qdisps, epstensors)]

    temp.sort(key=lambda item: item[0])

    qdisps = [qdisp for qdisp, _ in temp]
    epstensors = [epstensor for _, epstensor in temp]

    if fittype == 'finite_difference':
        if len(epstensors) != 2:
            raise Exception("Error: This function currently only supports "
                            "two-point finite-difference stencils."
                            )

        # Calculate and verify the step size.

        dispstemp = sorted(
            [0.0] + [disp for disp in qdisps]
            )

        stepsize = np.average(
            np.subtract(dispstemp[1:], dispstemp[:-1])
            )

        for qdisp in qdisps:
            if math.fabs(qdisp / stepsize) < 1.0e-8:
                raise Exception("Error: Finite-dfference stencils require a "
                                "uniform step size."
                                )

        # Calculate the Raman tensor.

        eps1, eps2 = epstensors

        # The two-point central-difference formula is f'(x) ~= [f(x + h)
        # - f(x - h)] / 2h.

        return (cellvolume / (4.0 * math.pi)) * np.subtract(eps2, eps1) / (
                2.0 * stepsize)


def scalaraverageramanactivitytensor(ramantensor):
    """ Average a Raman-activity tensor to obtain a scalar intensity. """

    # This formula came from D. Porezag and M. R. Pederson, Phys. Rev.
    # B: Condens. Matter Mater. Phys., 1996, 54, 7830.

    alpha = (ramantensor[0][0] + ramantensor[1][1] + ramantensor[2][2]) / 3.0

    betasquared = 0.5 * (
        (ramantensor[0][0] - ramantensor[1][1]) ** 2 +
        (ramantensor[0][0] - ramantensor[2][2]) ** 2 +
        (ramantensor[1][1] - ramantensor[2][2]) ** 2 +
        6.0 * (ramantensor[0][1] ** 2 + ramantensor[0][2] ** 2 +
               ramantensor[1][2] ** 2)
        )

    return 45.0 * alpha ** 2 + 7.0 * betasquared
