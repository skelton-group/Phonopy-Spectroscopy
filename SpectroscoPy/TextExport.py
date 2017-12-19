# SpectroscoPy/TextExport.py


# ---------
# Docstring
# ---------

""" Routines for exporting data in plain-text formats. """


# -------
# Imports
# -------

import csv;
import sys;

import numpy as np;

from SpectroscoPy import Constants;


# ---------
# Constants
# ---------

""" Format code for converting floating-point values to text when writing plain-text DAT files. """

_DATFileFloatFormat = "{0:.6f}";


# ------------------------
# Generic Writer Functions
# ------------------------

def WriteData(dataRows, filePath, fileFormat):
    """ Write row-wise data to an output file in the specified file format. """

    fileFormat = fileFormat.lower();

    # Dispatch to the correct writer function for the chosen file format.

    if fileFormat == 'csv':
        WriteDataCSV(dataRows, filePath);
    elif fileFormat == 'dat':
        WriteDataDAT(dataRows, filePath);
    else:
        raise Exception("Error: Unknown file format '{0}'.".format(fileFormat));

def WriteDataCSV(dataRows, filePath):
    """ Write row-wise data to an Excel-compatible CSV file. """

    outputWriter = None;

    # Workaround for a bug in the csv module, where using the csv.writer on Windows inserts extra blank lines between rows.

    if sys.platform.startswith("win"):
        # Try to open the file with newline = '' set (Python >= 3).
        # If this is not possible, issue a RuntimeWarning.

        if sys.version_info.major >= 3:
            outputWriter = open(filePath, 'w', newline = '');
        else:
            warnings.warn("CSV files output from Python < 3 on Windows platforms may have blank lines between rows.", RuntimeWarning);

    if outputWriter == None:
        outputWriter = open(filePath, 'w');

    outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

    for row in dataRows:
        outputWriterCSV.writerow(row);

    outputWriter.close();

def WriteDataDAT(dataRows, filePath):
    """ Write row-wise data to a plain-text DAT file. """

    # To generate a nicely-formatted DAT file, we need to know the column widths in advance.
    # The brute-force way to do this is to simply call str() on each data item and find the maximum length of the strings to be written in each column.

    dataRowsString = [];

    for dataRow in dataRows:
        dataRowString = [];

        for item in dataRow:
            # Figure out if item should be formatted as a float.

            floatFormat = False;

            try:
                if isinstance(item, float):
                    floatFormat = True;
                else:
                    item = float(item);

                    if not item.is_integer():
                        floatFormat = True;
            except (TypeError, ValueError):
                pass;

            if floatFormat:
                dataRowString.append(
                    _DATFileFloatFormat.format(item)
                    );
            else:
                dataRowString.append(
                    str(item)
                    );

        dataRowsString.append(dataRowString);

    # Calculate column widths and generate format codes.

    columnWidths = [len(item) for item in dataRowsString[0]];

    for dataRowString in dataRowsString[1:]:
        for i, item in enumerate(dataRowString):
            columnWidths[i] = max(columnWidths[i], len(item));

    columnFormatCodes = [
        "{{0: >{0}}}".format(columnWidth)
            for columnWidth in columnWidths
        ];

    # Write data.

    with open(filePath, 'w') as outputWriter:
        for dataRowString in dataRowsString:
            # Use four blank spaces as a column separator.

            outputWriter.write(
                '    '.join([formatCode.format(item) for formatCode, item in zip(columnFormatCodes, dataRowString)]) + '\n'
                );


# -----------
# Peak Tables
# -----------

def GroupForPeakTable(frequencies, intensities, irRepData, linewidths = None):
    """
    Average mode frequencies, intensities and linewidths (if supplied) into groups defined by the (symbol, band_indices) tuples in irRepData.
    Returns a tuple of (new_frequencies, new_intensities, irrep_symbols, new_linewidths) lists, each with the same length as irRepData.

    Notes:
        Band indices should be zero-based.
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

            if index < 0 or index >= numModes:
                raise Exception("Error: Band index {0} is out of bounds for # modes = {1}.".format(index, numModes));

            includedIndices.add(index);

    if len(includedIndices) != len(frequencies):
        raise Exception("Error: The number of bands references in the ir. rep. groups is more than # modes = {0}.".format(numModes));

    frequenciesNew, intensitiesNew = [], [];
    linewidthsNew = [] if linewidths != None else None;

    for _, bandIndices in irRepData:
        if len(bandIndices) == 1:
            index = bandIndices[0];

            frequenciesNew.append(frequencies[index]);
            intensitiesNew.append(intensities[index]);

            if linewidths != None:
                linewidthsNew.append(linewidths[index]);

        else:
            # Average the frequencies.

            frequenciesNew.append(
                np.average([frequencies[index] for index in bandIndices])
                );

            # Sum the intensities.

            intensitiesNew.append(
                np.sum([intensities[index] for index in bandIndices])
                );

            if linewidths != None:
                # Average the linewidths.

                linewidthsNew.append(
                    np.average([linewidths[index] for index in bandIndices])
                    );

    return (frequenciesNew, intensitiesNew, [symbol for symbol, _ in irRepData], linewidthsNew);

def SavePeakTable(frequencies, intensities, filePath, irRepSymbols = None, linewidths = None, frequencyUnits = None, intensityUnits = None, fileFormat = 'dat'):
    """
    Save a set of peak parameters (peak table) to a plain-text file.

    Arguments:
        frequencies -- list of mode frequencies.
        intensities -- list of band intensities.
        filePath -- path to output file.

    Keyword arguments:
        irRepSymbols -- (optional) list of irreducible representation symbols (ir. reps.).
        linewidths -- (optional) list of mode linewidths.
        frequencyUnits -- units of frequencies and linewidths, if supplied (default: Constants.DefaultFrequencyUnits).
        intensityUnits -- units of spectral intensity (default: Constants.DefaultIntensityUnits).
        fileFormat -- file format to export to ('csv' or 'dat'; default 'dat');
    """

    numModes = len(frequencies);

    if intensities == None or len(intensities) != numModes:
        raise Exception("Error: intensities cannot be None and must be the same length as frequencies.");

    if irRepSymbols != None and len(irRepSymbols) != numModes:
        raise Exception("Error: If supplied, irReps must be the same length as frequencies.");

    if linewidths != None and len(linewidths) != numModes:
        raise Exception("Error: If supplied, linewidths must be the same length as frequencies.");

    if frequencyUnits == None:
        frequencyUnits = Constants.DefaultFrequencyUnits;

    if intensityUnits == None:
        intensityUnits = Constants.DefaultIntensityUnits;

    # Compile data rows.

    headerRow = ["v [{0}]".format(Constants.GetFrequencyUnitLabel(frequencyUnits))];

    if irRepSymbols != None:
        headerRow = headerRow + ["Ir. Rep."];

    headerRow = headerRow + ["I [{0}]".format(intensityUnits)];

    if linewidths != None:
        headerRow = headerRow + [r"\Gamma [{0}]".format(Constants.GetFrequencyUnitLabel(frequencyUnits))];

    dataItems = [frequencies];

    if irRepSymbols != None:
        dataItems.append(irRepSymbols);

    dataItems.append(intensities);

    if linewidths != None:
        dataItems.append(linewidths);

    dataRows = [headerRow] + [dataRow for dataRow in zip(*dataItems)];

    # Write output file.

    WriteData(dataRows, filePath, fileFormat);


# -----------------
# Simulated Spectra
# -----------------

def SaveSpectrum(spectrum, filePath, frequencyUnits = None, intensityUnits = None, fileFormat = 'dat'):
    """
    Save a simulated spectrum to a plain-text file.

    Arguments:
        spectrum -- a tuple of (frequencies, intensities, intensities_norm) data sets.
        filePath -- path to output file.

    Keyword arguments:
        frequencyUnits -- units of the frequency axis (defaults to Constants.DefaultFrequencyUnits).
        intensityUnits -- units of the intensity axis (defaults to Constants.DefaultIntensityUnits).
        fileFormat -- file format to export to ('csv' or 'dat'; default: 'dat').
    """

    spectrumX, spectrumY, spectrumYNorm = spectrum;

    numDataPoints = len(spectrumX);

    if len(spectrumY) != numDataPoints or len(spectrumYNorm) != numDataPoints:
        raise Exception("Error: The number of frequency and intensity points in spectra are inconsistent.");

    if frequencyUnits == None:
        frequencyUnits = Constants.DefaultFrequencyUnits;

    if intensityUnits == None:
        intensityUnits = Constants.DefaultIntensityUnits;

    # Compile data rows.

    dataRows = [];

    dataRows.append(
        ["v [{0}]".format(Constants.GetFrequencyUnitLabel(frequencyUnits)), "I(v) [{0}]".format(intensityUnits), "I_Norm(v) [AU]"]
        );

    for dataRow in zip(*spectrum):
        dataRows.append(dataRow);

    # Write output file.

    WriteData(dataRows, filePath, fileFormat);
