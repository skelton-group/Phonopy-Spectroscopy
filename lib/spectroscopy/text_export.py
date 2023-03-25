# spectroscopy/text_export.py


# ---------
# Docstring
# ---------

""" Routines for exporting data in plain-text formats. """


# -------
# Imports
# -------

import csv
import sys
import warnings

import numpy as np

from spectroscopy import constants


# ---------
# Constants
# ---------

""" Format code for converting floating-point values to text when
writing plain-text DAT files. """

_DAT_FILE_FLOAT_FORMAT = "{0:.6f}"


# ------------------------
# Generic Writer Functions
# ------------------------

def write_data(data_rows, file_path, file_format):
    """ Write row-wise data to an output file in the specified file
    format. """

    # Dispatch to the correct writer function for the chosen file
    # format.

    if file_format == 'csv':
        write_data_csv(data_rows, file_path)
    elif file_format == 'dat':
        write_data_dat(data_rows, file_path)
    else:
        raise Exception(
            "Error: Unknown file format '{0}'.".format(file_format)
            )

def write_data_csv(data_rows, file_path):
    """ Write row-wise data to an Excel-compatible CSV file. """

    output_writer = None

    # Workaround for a bug (?) in the csv module, where using the
    # csv.writer on Windows inserts extra blank lines between rows.

    if sys.platform.startswith("win"):
        # Try to open the file with newline = '' set (Python >= 3). If
        # this is not possible, issue a RuntimeWarning.

        if sys.version_info.major >= 3:
            output_writer = open(file_path, 'w', newline='')
        else:
            warnings.warn(
                "CSV files output from Python < 3 on Windows platforms "
                "may have blank lines between rows.", RuntimeWarning)

    if output_writer is None:
        output_writer = open(file_path, 'w')

    output_writer_csv = csv.writer(
        output_writer, delimiter=',', quotechar='\"',
        quoting=csv.QUOTE_ALL)

    for row in data_rows:
        output_writer_csv.writerow(row)

    output_writer.close()

def write_data_dat(data_rows, file_path):
    """ Write row-wise data to a plain-text DAT file. """

    # To generate a nicely-formatted DAT file, we need to know the
    # column widths in advance. The brute-force way to do this is to
    # simply call str() on each data item and find the maximum length of
    # the strings to be written in each column.

    data_rows_string = []

    for data_row in data_rows:
        data_row_string = []

        for item in data_row:
            # Figure out if item should be formatted as a float.

            float_format = False

            try:
                if isinstance(item, float):
                    float_format = True
                else:
                    item = float(item)

                    if not item.is_integer():
                        float_format = True
            except (TypeError, ValueError):
                pass

            if float_format:
                data_row_string.append(_DAT_FILE_FLOAT_FORMAT.format(item))
            else:
                data_row_string.append(str(item))

        data_rows_string.append(data_row_string)

    # Calculate column widths and generate format codes.

    column_widths = [
        len(item) for item in data_rows_string[0]
        ]

    for data_row_string in data_rows_string[1:]:
        for i, item in enumerate(data_row_string):
            column_widths[i] = max(column_widths[i], len(item))

    column_format_codes = [
        "{{0: >{0}}}".format(column_width)
            for column_width in column_widths
        ]

    # Write data.

    with open(file_path, 'w') as output_writer:
        for data_row_string in data_rows_string:
            items = [
                format_code.format(item) for format_code, item in
                    zip(column_format_codes, data_row_string)
                ]
            
            # Use four blank spaces as a column separator.

            output_writer.write('    '.join(items) + '\n')


# -----------
# Peak Tables
# -----------

def group_for_peak_table(
        frequencies, intensities, irrep_data, linewidths=None
        ):
    """ Average mode frequencies, intensities and linewidths (if
    supplied) into groups defined by the (symbol, band_indices) tuples
    in irrep_data.
    
    Returns a tuple of (new_frequencies, new_intensities, irrep_symbols,
    new_linewidths) lists, each with the same length as irrep_data.

    Notes:
        Band indices should be zero-based.
    """

    num_modes = len(frequencies)

    if len(intensities) != num_modes:
        raise Exception(
            "Error: The lengths of frequencies and intensities are "
            "inconsistent.")

    if linewidths is not None and len(linewidths) != num_modes:
        raise Exception(
            "Error: If supplied, linewidths must have the same number "
            "of elements as frequencies and intensities.")

    # Check band indices.

    included_indices = set()

    for _, band_indices in irrep_data:
        for index in band_indices:
            if index in included_indices:
                raise Exception(
                    "Error: Band index {0} assigned to multiple irrep "
                    "groups.".format(index))

            if index < 0 or index >= num_modes:
                raise Exception(
                    "Error: Band index {0} is out of bounds for # "
                    "modes = {1}.".format(index, num_modes))

            included_indices.add(index)

    if len(included_indices) != len(frequencies):
        raise Exception(
            "Error: The number of bands references in the ir. rep. "
            "groups is more than # modes = {0}.".format(num_modes))

    frequencies_new, intensities_new = [], []
    linewidths_new = [] if linewidths is not None else None

    for _, band_indices in irrep_data:
        if len(band_indices) == 1:
            index = band_indices[0]

            frequencies_new.append(frequencies[index])
            intensities_new.append(intensities[index])

            if linewidths is not None:
                linewidths_new.append(linewidths[index])

        else:
            # Average the frequencies.

            frequencies_new.append(
                np.average([frequencies[index] for index in band_indices])
                )

            # Sum the intensities.

            intensities_new.append(
                np.sum([intensities[index] for index in band_indices])
                )

            if linewidths is not None:
                # Average the linewidths.

                linewidths_new.append(
                    np.average([linewidths[index] for index in band_indices])
                    )

    return (
        frequencies_new, intensities_new, [symbol for symbol, _ in irrep_data],
        linewidths_new
        )

def save_peak_table(
        frequencies, intensities, file_path,
        irrep_symbols=None, linewidths=None,
        frequency_units=constants.DEFAULT_FREQUENCY_UNITS,
        intensity_units=constants.DEFAULT_INTENSITY_UNITS, file_format='dat'
        ):
    """ Save a set of peak parameters (peak table) to a plain-text file.

    Arguments:
        frequencies -- list of mode frequencies.
        intensities -- list of band intensities.
        file_path -- path to output file.

    Keyword arguments:
        irrep_symbols -- (optional) list of irreducible representation
            symbols (irreps).
        linewidths -- (optional) list of mode linewidths.
        frequency_units -- units of frequencies and linewidths, if
            supplied (default: constants.DEFAULT_FREQUENCY_UNITS).
        intensity_units -- units of spectral intensity (default:
            constants.DEFAULT_INTENSITY_UNITS).
        file_format -- file format to export to ('csv' or 'dat';
            default 'dat').
    """

    num_modes = len(frequencies)

    if intensities is None or len(intensities) != num_modes:
        raise Exception(
            "Error: intensities cannot be None and must be the same "
            "length as frequencies.")

    if irrep_symbols is not None and len(irrep_symbols) != num_modes:
        raise Exception(
            "Error: If supplied, irrep_symbols must be the same length "
            "as frequencies.")
    
    if linewidths is not None and len(linewidths) != num_modes:
        raise Exception(
            "Error: If supplied, linewidths must be the same length as"
            "frequencies.")

    if frequency_units is None:
        frequency_units = constants.DEFAULT_FREQUENCY_UNITS

    if intensity_units is None:
        intensity_units = constants.DEFAULT_INTENSITY_UNITS

    # Compile data rows.

    header_row = [
        "v [{0}]".format(constants.get_frequency_unit_label(frequency_units))
        ]

    if irrep_symbols is not None:
        header_row = header_row + ["Ir. Rep."]

    header_row = header_row + ["I [{0}]".format(intensity_units)]

    if linewidths is not None:
        header_row = header_row + [
            r"\Gamma [{0}]".format(
                constants.get_frequency_unit_label(frequency_units)
                )
            ]

    data_items = [frequencies]

    if irrep_symbols is not None:
        data_items.append(irrep_symbols)

    data_items.append(intensities)

    if linewidths is not None:
        data_items.append(linewidths)

    data_rows = [header_row] + [datarow for datarow in zip(*data_items)]

    # Write output file.

    write_data(data_rows, file_path, file_format)


# -----------------
# Simulated Spectra
# -----------------

def save_spectrum(
        spectrum, file_path,
        frequency_units=None, intensity_units=None, file_format='dat'
        ):
    """ Save a simulated spectrum to a plain-text file.

    Arguments:
        spectrum -- a tuple of
            (frequencies, intensities, intensities_norm) data sets.
        file_path -- path to output file.

    Keyword arguments:
        frequency_units -- units of the frequency axis (defaults to
            constants.DEFAULT_FREQUENCY_UNITS).
        intensity_units -- units of the intensity axis (defaults to
            constants.DEFAULT_INTENSITY_UNITS).
        file_format -- file format to export to ('csv' or 'dat';
            default: 'dat').
    """

    spectrum_x, spectrum_y, spectrum_y_norm = spectrum

    num_data_points = len(spectrum_x)

    if (len(spectrum_y) != num_data_points
            or len(spectrum_y_norm) != num_data_points):
        raise Exception(
            "Error: The number of frequency and intensity points in "
            "spectra are inconsistent.")

    # Compile data rows.

    header_row = [
        "v [{0}]".format(constants.get_frequency_unit_label(frequency_units)),
        "I(v) [{0}]".format(intensity_units), "I_Norm(v) [AU]"
        ]

    data_rows = [header_row]

    for data_row in zip(*spectrum):
        data_rows.append(data_row)

    # Write output file.

    write_data(data_rows, file_path, file_format)
