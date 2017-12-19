# SpectroscoPy/Raman.py


# ---------
# Docstring
# ---------

""" Core routines for setting up and post-processing Raman-activity calculations. """


# -------
# Imports
# -------

import math;

import numpy as np;

from SpectroscoPy import Utilities;


# -----------
# Preparation
# -----------

def GenerateDisplacedStructures(structure, eigendisplacement, stepSize, numSteps, stepType = 'delta_q', omitZero = False):
    """
    Generate a set of structures displaced along a mode eigendisplacement for a Raman activity calculation.

    Arguments:
        structure -- a tuple of (lattice_vectors, atomic_symbols, atom_positions) specifying the structure to displace.
        eigendisplacement -- mode eigendisplacement (N x three-component vectors).
        stepSize -- size of the displacement step.
        numSteps -- number of displacement steps to make (minimum 2).

    Keyword arguments:
        stepType -- step specified by stepSize argument ('delta_q' - normal-mode coordinate, default; 'max_r' - max. Cartesian displacement; 'norm_x' - units of normalised eigendisplacement).
        omitZero -- if True, and if numSteps is odd, omit the zero step.

    Notes:
        If numSteps is even, displacements will be generated for +/- stepSize, +/- 2 * stepSize, etc.
        If numSteps is odd, the set of displacements will include a zero step (overridden by the omitZero argument).
        If stepType = 'max_r' is used, the max. displacement applies to +/- stepSize (i.e. if numSteps > 3, some steps will be larger than this).
    """

    latticeVectors, atomicSymbols, atomPositions = structure;

    # Using np.shape() should allow the parameters to be supplied as (correctly-dimensioned) lists.

    dim1, dim2 = np.shape(latticeVectors);

    if dim1 != 3 or dim2 != 3:
        raise Exception("Error: The lattice vectors provided in structure must be a 3x3 matrix.");

    numAtoms = len(atomicSymbols);

    if len(atomPositions) != numAtoms:
        raise Exception("Error: The lengths of the atomic symbol and atom position lists provided in structure are inconsistent.");

    dim1, dim2 = np.shape(eigendisplacement);

    if dim1 != numAtoms or dim2 != 3:
        raise Exception("Error: eigendisplacement must be an Nx3 matrix.");

    if stepSize <= 0.0:
        raise Exception("Error: stepSize must be positive.");

    if numSteps < 2:
        raise Exception("Error: numSteps must be at least 2.");

    stepType = stepType.lower();

    if stepType not in ['delta_q', 'max_r', 'norm_x']:
        raise Exception("Error: Unsupported stepType '{0}'.".format(stepType));

    # Displacements are created at 0, +/- stepSize, +/- (2 * stepSize), ... depending on numSteps and whether omitZero is set.

    dispSteps = [];

    iMax = numSteps // 2;

    for i in range(-iMax, iMax + 1):
        if i == 0 and (omitZero or numSteps % 2 == 0):
            continue;

        dispSteps.append(i * stepSize);

    # Convert atom positions to Cartesian coordinates.

    atomPositions = Utilities.FractionalToCartesianCoordinates(atomPositions, latticeVectors);

    # Displace structures along eigenvectors.

    dispStructures, maxDisps = [], [];

    # Calculate the maximum (Cartesian) offset in the eigendisplacement.

    maxDisp = max(
        np.linalg.norm(v) for v in eigendisplacement
        );

    # Calculate the norm of the eigendisplacement.

    eigendisplacementNorm = np.linalg.norm(eigendisplacement);

    # Calculate the norm of the eigendisplacement.

    for i, dispStep in enumerate(dispSteps):
        if stepType == 'max_r':
            # If the step is defined as a maximum Cartesian displacement, scale it by maxDisp.

            dispStep = dispStep / maxDisp;
            dispSteps[i] = dispStep;
        elif stepType == 'norm_x':
            # If the step is defined in terms of the normalised eigendisplacement, scale it by eigendisplacementNorm.

            dispStep = dispStep / eigendisplacementNorm;
            dispSteps[i] = dispStep;

        # Generate displaced positions.

        newAtomPositions = [
            np.add(position, np.multiply(dispStep, v))
                for position, v in zip(atomPositions, eigendisplacement)
            ];

        # Convert new positions back to fractional coordinates when storing.

        dispStructures.append(
            (latticeVectors, atomicSymbols, Utilities.CartesianToFractionalCoordinates(newAtomPositions, latticeVectors))
            );

        maxDisps.append(dispStep * maxDisp);

    # Return displacement steps, displaced structures and max. Cartesian displacements.

    return (dispSteps, dispStructures, maxDisps);


# --------
# Analysis
# --------

def CalculateRamanActivityTensor(epsTensors, qDisps, cellVolume, fitType):
    """
    Calculate the Raman activity from a set of dielectric tensors obtained at different displacements along the normal-mode coordinate Q.

    Arguments:
        epsTensors -- a set of calculated dielectric tensors at displacements along the normal-mode coordinate.
        qDisps -- normal-mode displacements at which the dielectric tensors were calculated.
        cellVolume -- volume of the unit cell (required to compute a normalisation factor).
        fitType -- fitting type to use to compute the derivative (currently only 'finite_difference' is supported).

    Notes:
        At present, this function only supports two-point finite difference stencils for computing the derivative.
    """

    if fitType not in ['finite_difference']:
        raise Exception("Error: Unsupported fitType '{0}'.".format(fitType));

    if len(epsTensors) != len(qDisps):
        raise Exception("Error: The numbers of supplied dielectric tensors and displacement steps are inconsistent.");

    for epsTensor in epsTensors:
        dim1, dim2 = np.shape(epsTensor);

        if dim1 != 3 or dim2 != 3:
            raise Exception("Error: Dielectric tensors must be 3x3 matrices.");

    # Make sure the dielectric tensors and displacements are ordered by step.

    temp = [item for item in zip(qDisps, epsTensors)];

    temp.sort(key = lambda item : item[0]);

    qDisps = [qDisp for qDisp, _ in temp];
    epsTensors = [epsTensor for _, epsTensor in temp];

    if fitType == 'finite_difference':
        if len(epsTensors) != 2:
            raise Exception("Error: This function currently only supports two-point finite-difference stencils.");

        # Calculate and verify the step size.

        dispsTemp = sorted(
            [0.0] + [disp for disp in qDisps]
            );

        stepSize = np.average(
            np.subtract(dispsTemp[1:], dispsTemp[:-1])
            );

        for qDisp in qDisps:
            if math.fabs(qDisp / stepSize) < 1.0e-8:
                raise Exception("Error: Finite-dfference stencils require a uniform step size.");

        # Calculate the Raman tensor.

        eps1, eps2 = epsTensors;

        # The two-point central-difference formula is f'(x) ~= [f(x + h) - f(x - h)] / 2h.

        return (cellVolume / (4.0 * math.pi)) * np.subtract(eps2, eps1) / (2.0 * stepSize);

def ScalarAverageRamanActivityTensor(ramanTensor):
    """ Average a Raman-activity tensor to obtain a scalar intensity. """

    # This formula came from D. Porezag and M. R. Pederson, Phys. Rev. B: Condens. Matter Mater. Phys., 1996, 54, 7830.

    alpha = (ramanTensor[0][0] + ramanTensor[1][1] + ramanTensor[2][2]) / 3.0;

    betaSquared = 0.5 * (
        (ramanTensor[0][0] - ramanTensor[1][1]) ** 2 +
        (ramanTensor[0][0] - ramanTensor[2][2]) ** 2 +
        (ramanTensor[1][1] - ramanTensor[2][2]) ** 2 +
        6.0 * (ramanTensor[0][1] ** 2 + ramanTensor[0][2] ** 2 + ramanTensor[1][2] ** 2)
        );

    return 45.0 * alpha ** 2 + 7.0 * betaSquared;
