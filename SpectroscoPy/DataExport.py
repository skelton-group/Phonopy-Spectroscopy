# SpectroscoPy/DataExport.py


# ---------
# Docstring
# ---------

""" Core routines for exporting data in plain-text formats. """


# -------
# Imports
# -------

import csv;
import sys;


# ---------
# Constants
# ---------

""" Default frequency unit. """

DefaultFrequencyUnit = "?";

""" Default intensity unit. """

DefaultIntensityUnit = "AU";

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
                item = float(item);

                if not item.is_integer():
                    floatFormat = True;
            except ValueError:
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


# ---------------------
# Text Export Functions
# ---------------------

def SavePeakTable(frequencies, intensities, filePath, irReps = None, linewidths = None, frequencyUnit = None, intensityUnit = None, fileFormat = 'dat'):
    """
    Save a set of peak parameters (peak table) to a plain-text file.

    Arguments:
        frequencies -- list of mode frequencies.
        intensities -- list of band intensities.
        filePath -- path to output file.

    Keyword arguments:
        irReps -- (optional) list of irreducible representation symbols (ir. reps.).
        linewidths -- (optional) list of mode linewidths.
        frequencyUnit -- units of frequencies and linewidths, if supplied (default: DefaultFrequencyUnit constant).
        intensityUnit -- units of spectral intensity (default: DefaultIntensityUnit constant).
        fileFormat -- file format to export to ('csv' or 'dat'; default 'dat');
    """

    numModes = len(frequencies);

    if intensities == None or len(intensities) != numModes:
        raise Exception("Error: intensities cannot be None and must be the same length as frequencies.");

    if irReps != None and len(irReps) != numModes:
        raise Exception("Error: If supplied, irReps must be the same length as frequencies.");

    if linewidths != None and len(irReps) != numModes:
        raise Exception("Error: If supplied, linewidths must be the same length as frequencies.");

    if frequencyUnit == None:
        frequencyUnit = DefaultFrequencyUnit;

    if intensityUnit == None:
        intensityUnit = DefaultIntensityUnit;

    # Compile data rows.

    headerRow = ["v [{0}]".format(frequencyUnit)];

    if irReps != None:
        headerRow = headerRow + ["Ir. Rep."];

    headerRow = headerRow + ["I [{0}]".format(intensityUnit)];

    if linewidths != None:
        headerRow = headerRow + [r"\Gamma [{0}]".format(frequencyUnit)];

    dataItems = [frequencies];

    if irReps != None:
        dataItems.append(irReps);

    dataItems.append(intensities);

    if linewidths != None:
        dataItems.append(linewidths);

    dataRows = [headerRow] + [dataRow for dataRow in zip(*dataItems)];

    # Write output file.

    WriteData(dataRows, filePath, fileFormat);

def SaveSpectrum(spectrum, filePath, frequencyUnit = None, intensityUnit = None, fileFormat = 'dat'):
    """
    Save a simulated spectrum to a plain-text file.

    Arguments:
        spectrum -- a tuple of (frequencies, intensities) data sets.
        filePath -- path to output file.

    Keyword arguments:
        frequencyUnit -- units of the frequency axis (defaults to the DefaultFrequencyUnit constant).
        intensityUnit -- units of the intensity axis (defaults to the DefaultIntensityUnit constant).
        fileFormat -- file format to export to ('csv' or 'dat'; default: 'dat').
    """

    spectrumX, spectrumY = spectrum;

    if len(spectrumX) != len(spectrumY):
        raise Exception("Error: the number of frequency and intensity points in spectrum are inconsistent.");

    if frequencyUnit == None:
        frequencyUnit = DefaultFrequencyUnit;

    if intensityUnit == None:
        intensityUnit = DefaultIntensityUnit;

    # Compile data rows.

    dataRows = [];

    dataRows.append(
        ["v [{0}]".format(frequencyUnit), "I(v) [{0}]".format(intensityUnit)]
        );

    for dataRow in zip(*spectrum):
        dataRows.append(dataRow);

    # Write output file.

    WriteData(dataRows, filePath, fileFormat);
