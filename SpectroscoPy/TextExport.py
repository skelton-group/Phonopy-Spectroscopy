# SpectroscoPy/TextExport.py


# ---------
# Docstring
# ---------

""" Routines for exporting data in plain-text formats. """


# -------
# Imports
# -------

import csv
import sys

import numpy as np

from SpectroscoPy import Constants


# ---------
# Constants
# ---------

""" Format code for converting floating-point values to text when writing 
plain-text DAT files. """

_DATFileFloatFormat = "{0:.6f}"


# ------------------------
# Generic Writer Functions
# ------------------------

def writedata(datarows, filepath, fileformat):
    """ Write row-wise data to an output file in the specified file format. """

    # Dispatch to the correct writer function for the chosen file format.

    if fileformat == 'csv':
        writedatacsv(datarows, filepath)
    elif fileformat == 'dat':
        writedatadat(datarows, filepath)
    else:
        raise Exception("Error: Unknown file format '{0}'.".format(fileformat))


def writedatacsv(datarows, filepath):
    """ Write row-wise data to an Excel-compatible CSV file. """

    outputwriter = None

    # Workaround for a bug in the csv module, where using the csv.writer on
    # Windows inserts extra blank lines between rows.

    if sys.platform.startswith("win"):
        # Try to open the file with newline = '' set (Python >= 3).
        # If this is not possible, issue a RuntimeWarning.

        if sys.version_info.major >= 3:
            outputwriter = open(filePath, 'w', newline='')
        else:
            warnings.warn("CSV files output from Python < 3 on Windows "
                          "platforms may have blank lines between rows.",
                          RuntimeWarning
                          )

    if outputwriter is None:
        outputwriter = open(filepath, 'w')

    outputwritercsv = csv.writer(outputwriter, delimiter=',',
                                 quotechar='\"', quoting=csv.QUOTE_ALL)

    for row in datarows:
        outputwritercsv.writerow(row)

    outputwriter.close()


def writedatadat(datarows, filepath):
    """ Write row-wise data to a plain-text DAT file. """

    # To generate a nicely-formatted DAT file, we need to know the column
    # widths in advance.
    # The brute-force way to do this is to simply call str() on each data item
    # and find the maximum length of the strings to be written in each column.

    datarowsstring = []

    for datarow in datarows:
        datarowstring = []

        for item in datarow:
            # Figure out if item should be formatted as a float.

            floatformat = False

            try:
                if isinstance(item, float):
                    floatformat = True
                else:
                    item = float(item)

                    if not item.is_integer():
                        floatformat = True
            except (TypeError, ValueError):
                pass

            if floatformat:
                datarowstring.append(
                    _DATFileFloatFormat.format(item)
                    )
            else:
                datarowstring.append(
                    str(item)
                    )

        datarowsstring.append(datarowstring)

    # Calculate column widths and generate format codes.

    columnwidths = [len(item) for item in datarowsstring[0]]

    for datarowstring in datarowsstring[1:]:
        for i, item in enumerate(datarowstring):
            columnwidths[i] = max(columnwidths[i], len(item))

    columnformatcodes = [
        "{{0: >{0}}}".format(columnwidth)
        for columnwidth in columnwidths
        ]

    # Write data.

    with open(filepath, 'w') as outputwriter:
        for datarowstring in datarowsstring:
            # Use four blank spaces as a column separator.

            outputwriter.write(
                '    '.join([formatCode.format(item) for formatCode, item in
                             zip(columnformatcodes, datarowstring)]) + '\n'
                )


# -----------
# Peak Tables
# -----------

def groupforpeaktable(frequencies, intensities, irrepdata, linewidths=None):
    """
    Average mode frequencies, intensities and linewidths (if supplied) into
    groups defined by the (symbol, band_indices) tuples in irRepData.
    Returns a tuple of (new_frequencies, new_intensities, irrep_symbols,
    new_linewidths) lists, each with the same length as irRepData.

    Notes:
        Band indices should be zero-based.
    """

    nummodes = len(frequencies)

    if len(intensities) != nummodes:
        raise Exception("Error: The lengths of frequencies and intensities "
                        "are inconsistent."
                        )

    if linewidths is not None and len(linewidths) != nummodes:
        raise Exception("Error: If supplied, linewidths must have the same "
                        "number of elements as frequencies and intensities."
                        )

    # Check band indices.

    includedindices = set()

    for _, bandindices in irrepdata:
        for index in bandindices:
            if index in includedindices:
                raise Exception("Error: Band index {0} assigned to multiple "
                                "ir. rep. groups.".format(index)
                                )

            if index < 0 or index >= nummodes:
                raise Exception("Error: Band index {0} is out of bounds for "
                                "# modes = {1}.".format(index, nummodes)
                                )

            includedindices.add(index)

    if len(includedindices) != len(frequencies):
        raise Exception("Error: The number of bands references in the ir. "
                        "rep. groups is more than # modes = {0}."
                        .format(nummodes)
                        )

    frequenciesnew, intensitiesnew = [], []
    linewidthsnew = [] if linewidths is not None else None

    for _, bandindices in irrepdata:
        if len(bandindices) == 1:
            index = bandindices[0]

            frequenciesnew.append(frequencies[index])
            intensitiesnew.append(intensities[index])

            if linewidths is not None:
                linewidthsnew.append(linewidths[index])

        else:
            # Average the frequencies.

            frequenciesnew.append(
                np.average([frequencies[index] for index in bandindices])
                )

            # Sum the intensities.

            intensitiesnew.append(
                np.sum([intensities[index] for index in bandindices])
                )

            if linewidths is not None:
                # Average the linewidths.

                linewidthsnew.append(
                    np.average([linewidths[index] for index in bandindices])
                    )

    return (frequenciesnew, intensitiesnew, [symbol for symbol, _
                                             in irrepdata], linewidthsnew)


def savepeaktable(frequencies, intensities, filepath, irrepsymbols=None,
                  linewidths=None, frequencyunits=None, intensityunits=None,
                  fileformat='dat'):
    """
    Save a set of peak parameters (peak table) to a plain-text file.

    Arguments:
        frequencies -- list of mode frequencies.
        intensities -- list of band intensities.
        filePath -- path to output file.

    Keyword arguments:
        irRepSymbols -- (optional) list of irreducible representation symbols
        (ir. reps.).
        linewidths -- (optional) list of mode linewidths.
        frequencyUnits -- units of frequencies and linewidths, if supplied
        (default: Constants.DefaultFrequencyUnits).
        intensityUnits -- units of spectral intensity (default:
        Constants.DefaultIntensityUnits).
        fileFormat -- file format to export to ('csv' or 'dat'; default 'dat');
    """

    nummodes = len(frequencies)

    if intensities is None or len(intensities) != nummodes:
        raise Exception("Error: intensities cannot be None and must be the "
                        "same length as frequencies."
                        )

    if irrepsymbols is not None and len(irrepsymbols) != nummodes:
        raise Exception("Error: If supplied, irReps must be the same length "
                        "as frequencies."
                        )

    if linewidths is not None and len(linewidths) != nummodes:
        raise Exception("Error: If supplied, linewidths must be the same "
                        "length as frequencies."
                        )

    if frequencyunits is None:
        frequencyunits = Constants.DefaultFrequencyUnits

    if intensityunits is None:
        intensityunits = Constants.DefaultIntensityUnits

    # Compile data rows.

    headerrow = ["v [{0}]".format(Constants.GetFrequencyUnitLabel(
        frequencyunits))]

    if irrepsymbols is not None:
        headerrow = headerrow + ["Ir. Rep."]

    headerrow = headerrow + ["I [{0}]".format(intensityunits)]

    if linewidths is not None:
        headerrow = headerrow + [r"\Gamma [{0}]".format(
            Constants.GetFrequencyUnitLabel(frequencyunits))]

    dataitems = [frequencies]

    if irrepsymbols is not None:
        dataitems.append(irrepsymbols)

    dataitems.append(intensities)

    if linewidths is not None:
        dataitems.append(linewidths)

    datarows = [headerrow] + [datarow for datarow in zip(*dataitems)]

    # Write output file.

    writedata(datarows, filepath, fileformat)


# -----------------
# Simulated Spectra
# -----------------

def savespectrum(spectrum, filepath, frequencyunits=None, intensityunits=None,
                 fileformat='dat'):
    """
    Save a simulated spectrum to a plain-text file.

    Arguments:
        spectrum -- a tuple of (frequencies, intensities, intensities_norm)
        data sets.
        filePath -- path to output file.

    Keyword arguments:
        frequencyUnits -- units of the frequency axis (defaults to
        Constants.DefaultFrequencyUnits).
        intensityUnits -- units of the intensity axis (defaults to
        Constants.DefaultIntensityUnits).
        fileFormat -- file format to export to ('csv' or 'dat'; default:
        'dat').
    """

    spectrumx, spectrumy, spectrumynorm = spectrum

    numdatapoints = len(spectrumx)

    if len(spectrumy) != numdatapoints or len(spectrumynorm) != numdatapoints:
        raise Exception("Error: The number of frequency and intensity points "
                        "in spectra are inconsistent."
                        )

    if frequencyunits is None:
        frequencyunits = Constants.DefaultFrequencyUnits

    if intensityunits is None:
        intensityunits = Constants.DefaultIntensityUnits

    # Compile data rows.

    datarows = []

    datarows.append(
        ["v [{0}]".format(Constants.GetFrequencyUnitLabel(frequencyunits)),
         "I(v) [{0}]".format(intensityunits), "I_Norm(v) [AU]"]
        )

    for datarow in zip(*spectrum):
        datarows.append(datarow)

    # Write output file.

    writedata(datarows, filepath, fileformat)
