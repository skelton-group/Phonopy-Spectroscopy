# SpectroscoPy/Utilities.py


# ---------
# Docstring
# ---------

""" Utility functions and constants. """


# -------
# Imports
# -------

import numpy as np;

from SpectroscoPy import Constants;


# -------------
# Unit Handling
# -------------

def GetFrequencyUnitLabel(frequencyUnits):
    """ Get a label for units frequencyUnits to be used in plots and output files. """

    # If frequencyUnits is one of the units supported by ConvertFrequencyUnits(), FrequencyUnitLabels will have a label for it.

    key = frequencyUnits.lower();

    if key in Constants.FrequencyUnitLabels:
        return Constants.FrequencyUnitLabels[key];
    else:
        return frequencyUnits;

def ConvertFrequencyUnits(frequencies, unitsFrom, unitsTo):
    """
    Convert frequencies in unitsFrom to unitsTo.
    Supported values of unitsFrom/unitsTo are: 'thz', 'inv_cm' and 'mev'.
    """

    unitsFrom, unitsTo = unitsFrom.lower(), unitsTo.lower();

    # Check unit conversion is supported.

    for units in unitsFrom, unitsTo:
        if units != 'thz' and units not in Constants.FrequencyConversionFactors:
            raise Exception("Error: Unsupported frequency unit '{0}'.".format(units));

    # No point in doing any work if we don't have to...

    if unitsFrom == unitsTo:
        return frequencies;

    conversionFactor = 1.0;

    if unitsFrom != 'thz':
        # Convert frequencies to THz.

        conversionFactor /= Constants.FrequencyConversionFactors[unitsFrom];

    if unitsTo != 'thz':
        conversionFactor *= Constants.FrequencyConversionFactors[unitsTo];

    return [frequency * conversionFactor for frequency in frequencies];


# ----------------------
# Peak-Table Preparation
# ----------------------

def GroupForPeakTable(frequencies, intensities, irRepData, linewidths = None):
    """
    Average mode frequencies, intensities and linewidths (if supplied) into groups defined by the (symbol, band_indices) tuples in irRepData.
    Returns a tuple of (new_frequencies, new_intensities, irrep_symbols, new_linewidths) lists, each with the same length as irRepData.
    """

    numModes = len(frequencies);

    if len(intensities) != numModes:
        raise Exception("Error: The lengths of frequencies and intensities are inconsistent.");

    if linewidths is not None and len(linewidths) != numModes:
        raise Exception("Error: If supplied, linewidths must have the same number of elements as frequencies and intensities.");

    # Check band indices.

    includedIndices = set();

    for _, bandIndices in irRepData:
        for index in bandIndices:
            if index in includedIndices:
                raise Exception("Error: Band index {0} assigned to multiple ir. rep. groups.".format(index));

            if index < 1 or index > numModes:
                raise Exception("Error: Band index {0} is out of bounds for # modes = {1}.".format(index, numModes));

            includedIndices.add(index);

    if len(includedIndices) != len(frequencies):
        raise Exception("Error: The number of bands references in the ir. rep. groups is more than # modes = {0}.".format(numModes));

    frequenciesNew, intensitiesNew = [], [];
    linewidthsNew = [] if linewidths != None else None;

    for _, bandIndices in irRepData:
        if len(bandIndices) == 1:
            index = bandIndices[0] - 1;

            frequenciesNew.append(frequencies[index]);
            intensitiesNew.append(intensities[index]);

            if linewidths != None:
                linewidthsNew.append(linewidths[index]);

        else:
            # Average the frequencies.

            frequenciesNew.append(
                np.average([frequencies[index - 1] for index in bandIndices])
                );

            # Sum the intensities.

            intensitiesNew.append(
                np.sum([intensities[index - 1] for index in bandIndices])
                );

            if linewidths != None:
                # Average the linewidths.

                linewidthsNew.append(
                    np.average([linewidths[index - 1] for index in bandIndices])
                    );

    return (frequenciesNew, intensitiesNew, [symbol for symbol, _ in irRepData], linewidthsNew);


# ---------------
# Helper Routines
# ---------------

def EigenvectorsToEigendisplacements(eigenvectors, atomicMasses):
    """
    Return the eigenvectors after division of each component by sqrt(mass).

    Arguments:
        eigenvectors -- eigenvectors as 3N x N three-component vectors.
        atomicMasses -- set of N atomic masses.
    """

    numModes, numAtoms, eigDim3 = len(eigenvectors), len(eigenvectors[0]), len(eigenvectors[0][0]);

    if numModes != 3 * numAtoms or eigDim3 != 3:
        raise Exception("Error: eigenvectors should be a 3N x N x 3 matrix.");

    massesDim1 = len(atomicMasses);

    if massesDim1 != numAtoms:
        raise Exception("Error: The number of supplied atomic masses is inconsistent with the number of displacements in the eigenvectors.");

    sqrtMasses = np.sqrt(atomicMasses);

    eigendisplacements = np.zeros(
        (numModes, numAtoms, 3), dtype = np.float64
        );

    for i in range(0, numModes):
        # This should avoid requiring external code to supply eigenvectors as NumPy arrays.

        eigenvector = eigenvectors[i];

        for j in range(0, numAtoms):
            eigendisplacements[i, j, :] = np.divide(eigenvector[j], sqrtMasses[j]);

    return eigendisplacements;
