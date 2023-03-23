# spectroscopy/cli/runtime.py


# ---------
# Docstring
# ---------

""" Implements standard runtime routines. """


# -------
# Imports
# -------

import math
import os

from spectroscopy.ir import calculate_ir_intensity

from spectroscopy.raman import (
    generate_displaced_structures, calculate_raman_activity_tensor,
    scalar_average_raman_activity_tensor)

from spectroscopy.spectrum import simulate_spectrum

from spectroscopy.text_export import (
    group_for_peak_table, save_peak_table, save_spectrum)

from spectroscopy.yaml_export import save_raman_tensors

from spectroscopy.plotting import plot_spectrum

from spectroscopy.constants import get_frequency_unit_label

from spectroscopy.utilities import (
    calculate_cell_volume, convert_frequency_units)

from spectroscopy.cli.intermediate_io import (
    read_raman_dataset, write_raman_dataset)


# ---------
# Constants
# ---------

""" Frequency units used for the peak table if the user requests a 
wavelength unit for the simulated spectrum. """

_DEFAULT_PEAK_TABLE_FREQUENCY_UNITS = 'inv_cm'


""" Units of Raman displacement step (normal-mode coordinate) for
plain-text intermediate YAML data sets. """

_RAMAN_STEP_UNITS = "sqrt(amu) * Ang"

""" Units of maximum displacement distance in Raman calculations for 
plain-text intermediate YAML data sets. """

_RAMAN_DISTANCE_UNITS = "Ang"

""" Units of the cell volume for plain-text intermediate YAML data
sets. """

_RAMAN_VOLUME_UNITS = "Ang ^ 3"


""" Units of (tensor) Raman activity for plain-text YAML output
file. """

_RAMAN_ACTIVITY_UNITS = "sqrt(amu) * Ang ^ 2"


""" Units of IR intensity for plots/text output files. """

_IR_INTENSITY_UNITS = "$e^{2}$ amu$^{-1}$"

""" Units of (scalar) Raman intensity for plots/text output files. """

_RAMAN_INTENSITY_UNITS = r"$\AA^{4}$ amu$^{-1}$"


# ---------------
# IR Spectroscopy
# ---------------

def run_mode_ir(
        frequencies, frequency_units, eigendisplacements, bec_tensors,
        args, linewidths=None, irrep_data=None
        ):
    """ Calculates IR intensities for spectrum simulation.

    Main arguments:
        eigendisplacements -- Gamma-point eigendisplacements.
        bec_tensors -- Born effective-charge tensors.

    All other arguments are passed through to the run_mode_output()
    function.
    """

    # Calculate IR intensities.

    ir_intensities = [
        calculate_ir_intensity(eigendisplacement, bec_tensors)
            for eigendisplacement in eigendisplacements
        ]

    # Hand over to the RunMode_Output() routine.

    run_mode_output(
        'ir', frequencies, frequency_units, ir_intensities,
        _IR_INTENSITY_UNITS, args, linewidths=linewidths,
        irrep_data=irrep_data)


# ------------------
# Raman Spectroscopy
# ------------------

def run_mode_raman_disp(
        structure, frequencies, frequency_units, eigendisplacements,
        args):
    """ Prepare and return a sequence of displaced structures for a
    Raman activity calculation. Metadata about the displacements
    performed is written to "[<output_prefix>_]Raman.yaml".

    Arguments:
        structure -- structure to displace as a tuple of
            (lattice_vectors, atomic_symbols, atom_positions_frac).
        frequencies -- Gamma-point mode frequencies.
        frequencyUnits -- units of mode frequencies.
        eigendisplacements -- Gamma-point eigendisplacements.
        args -- parsed command-line arguments with the RamanStep,
            RamanStepType, RamanBandIndices and, optionally,
            OutputPrefix attributes.

    Return value:
        A list of tuples of (band_index, band_frequency, structures,
        disp_steps, max_disps) data.
        Band indices are zero-based.
    """

    num_modes = len(frequencies)

    band_indices = args.RamanBandIndices

    if band_indices is None:
        # Exclude the acoustic modes by identifying the three modes
        # whose absolute frequency is closest to zero. This is fairly
        # crude compared to e.g. checking the eigenvector coefficients,
        # but the user can always override it using command-line
        # arguments.

        temp = [
            (i, math.fabs(frequency))
                for i, frequency in enumerate(frequencies)
            ]

        temp.sort(key=lambda item: item[1])

        exclude_indices = [i for i, _ in temp[:3]]

        print(
            "INFO: Automatically excluding bands {0}, {1}, {2} from "
            "displacement set (to override use --raman-bands)."
            .format(*[index + 1 for index in exclude_indices]))

        band_indices = [
            i for i in range(0, num_modes)
                if i not in exclude_indices
            ]

    for index in band_indices:
        if index < 0 or index >= num_modes:
            raise Exception(
                "Error: Band index {0} out of range for # modes = "
                "{1}.".format(index + 1, num_modes))

    # Build a dictionary of displacement sets: key -> band index,
    # value -> (frequency, [(step, structure) tuples]). At present, we
    # use a default of two finite-difference steps per mode.

    disp_sets = {}

    for index in band_indices:
        disp_structures = generate_displaced_structures(
            structure, eigendisplacements[index], args.RamanStep, 2,
            step_type=args.RamanStepType, omit_zero=True
            )

        disp_sets[index] = (frequencies[index], disp_structures)

    # Write data set to a YAML file.

    lattice_vectors, _, _ = structure
    cell_volume = calculate_cell_volume(lattice_vectors)

    frequencies, mode_disp_steps, mode_max_disps = [], [], []

    for index in band_indices:
        frequency, (disp_steps, _, max_disps) = disp_sets[index]

        frequencies.append(frequency)
        mode_disp_steps.append(disp_steps)
        mode_max_disps.append(max_disps)

    file_name = "Raman.yaml"

    if args.OutputPrefix is not None:
        file_name = "{0}_{1}".format(args.OutputPrefix, file_name)

    write_raman_dataset(
        file_name, cell_volume, band_indices, frequencies,
        mode_disp_steps, mode_max_disps, eps_tensor_sets=None,
        volume_units=_RAMAN_VOLUME_UNITS, frequency_units=frequency_units,
        step_units=_RAMAN_STEP_UNITS, distance_units=_RAMAN_DISTANCE_UNITS)

    # Return the band indices, frequencies and displacement sets.

    disp_sets_list = []

    for index in band_indices:
        frequency, (disp_steps, structures, max_disps) = disp_sets[index]

        disp_sets_list.append(
            (index, frequency, structures, disp_steps, max_disps))

    return disp_sets_list

def run_mode_raman_read(eps_tensors, args):
    """ Collect dielectric tensors calculated for displaced structures
    prepared for a Raman calculation.

    Metadata on the displaced structures is read from
    "[<output_prefix>_]Raman.yaml". An updated file is then written
    which contains the displacement metadata and calculated dielectric
    tensors.

    Arguments:
        eps_tensors -- list of dielectric tensors calculated for
            displaced structures.
        args -- parsed command-line arguments with the optional
            OutputPrefix attribute.

    Notes:
        The dielectric tensors are assumed to be in the same order as
        the displacements in the YAML data set written out by
        run_mode_raman_disp().
    """
    # Get the name of the intermediate YAML file.

    dataset_file = "Raman.yaml"

    if args.OutputPrefix is not None:
        dataset_file = "{0}_{1}".format(args.OutputPrefix, dataset_file)

    if not os.path.isfile(dataset_file):
        raise Exception(
            "Error: Displacement data set file \"{0}\" not found."
            .format(dataset_file))

    # Read displacement data set: tuple of (cell_volume, volume_units,
    # band_indices, frequencies, frequency_units, disp_step_sets,
    # step_units, max_disps_sets, distance_units, eps_tensor_sets).

    cell_volume, volume_units, band_indices, frequencies, frequency_units, \
        disp_step_sets, step_units, max_disps_sets, distance_units, _ \
        = read_raman_dataset(dataset_file)
    
    # Verify units.

    if volume_units != _RAMAN_VOLUME_UNITS:
        raise Exception(
            "Error: Unexpected volume units '{0}' read from "
            "displacement dataset file \"{1}\"."
            .format(volume_units, dataset_file))

    if step_units != _RAMAN_STEP_UNITS:
        raise Exception(
            "Error: Unexpected step units '{0}' read from "
            "displacement dataset file \"{1}\"."
            .format(step_units, dataset_file))

    if distance_units != _RAMAN_DISTANCE_UNITS:
        raise Exception(
            "Error: Unexpected distance units '{0}' read from "
            "displacement dataset file \"{1}\"."
            .format(distance_units, dataset_file))

    # Check the expected number of dielectric tensors has been supplied.

    num_eps_tensors_expected = sum(
        len(disp_steps) for disp_steps in disp_step_sets
        )

    if len(eps_tensors) != num_eps_tensors_expected:
        raise Exception(
            "Error: Incorrect number of dielectric tensors supplied "
            "- {0} expected, {1} supplied."
            .format(num_eps_tensors_expected, len(eps_tensors)))

    # Sort the dielectric tensors into sets.

    eps_tensor_sets = []

    input_pointer = 0

    for disp_step_set in disp_step_sets:
        eps_tensor_set = []

        for i in range(0, len(disp_step_set)):
            eps_tensor_set.append(eps_tensors[input_pointer])

            input_pointer += 1

        eps_tensor_sets.append(eps_tensor_set)

    # Write an updated data set containing the extracted dielectric
    # tensors.

    write_raman_dataset(
        dataset_file, cell_volume, band_indices, frequencies,
        disp_step_sets, max_disps_sets, eps_tensor_sets, volume_units,
        frequency_units, step_units, distance_units)

def run_mode_raman_postproc(args, linewidths=None, irrep_data=None):
    """ Post process a data set prepared for a Raman calculation.

    A data set containing frequencies, displacement steps and dielectric
    tensors is read from "[<output_prefix>_]Raman.yaml". The Raman
    activity tensors are calculated from the derivatives of the static
    dielectric tensor with respect to the normal-mode coordinates and
    written to "[<output_prefix>_]RamanActivity.yaml". The activity
    tensors are averaged to obtai scalar intensities that are then used
    to plot a spectrum.

    Arguments:
        args -- parsed command-line arguments with the optional
            OutputPrefix attribute, plus settings to be passed to the
            RunMode_Output() routine.

    Keyword arguments:
        linewidths -- linewidths used to simulate the spectrum.
        irrep_data -- irreps used to group modes in the peak table.
    """
    # Work out the name of the intermediate YAML file and check the file
    # exists.

    dataset_file = "Raman.yaml"

    if args.OutputPrefix is not None:
        dataset_file = "{0}_{1}".format(args.OutputPrefix, dataset_file)

    if not os.path.isfile(dataset_file):
        raise Exception(
            "Error: Displacement data set file \"{0}\" not found."
            .format(dataset_file))

    # Read in the Raman data set and verify hard-coded units.

    cell_volume, volume_units, band_indices, frequencies, frequency_units, \
        disp_step_sets, step_units, max_disps_sets, distance_units, \
        eps_tensor_sets = read_raman_dataset(dataset_file)

    if volume_units != _RAMAN_VOLUME_UNITS:
        raise Exception(
            "Error: Unexpected volume units '{0}' read from Raman "
            "dataset file \"{1}\".".format(volume_units, dataset_file))

    if step_units != _RAMAN_STEP_UNITS:
        raise Exception(
            "Error: Unexpected step units '{0}' read from Raman "
            "dataset file \"{1}\".".format(step_units, dataset_file))

    if distance_units != _RAMAN_DISTANCE_UNITS:
        raise Exception(
            "Error: Unexpected distance units '{0}' read from Raman "
            "dataset file \"{1}\".".format(distance_units, dataset_file))

    # Calculate Raman activity tensors and averaged scalar intensities.

    raman_tensors, raman_intensities = [], []

    for disp_steps, eps_tensors in zip(disp_step_sets, eps_tensor_sets):
        raman_tensor = calculate_raman_activity_tensor(
            eps_tensors, disp_steps, cell_volume, 'finite_difference')

        raman_tensors.append(raman_tensor)

        raman_intensities.append(
            scalar_average_raman_activity_tensor(raman_tensor))

    # Write out Raman tensors.

    file_name = "Raman-Tensors.yaml"

    if args.OutputPrefix is not None:
        file_name = "{0}_{1}".format(args.OutputPrefix, file_name)

    # Convert frequency units, if required.

    frequencies_output = None

    if frequency_units != args.Units:
        frequencies_output = convert_frequency_units(
            frequencies, frequency_units, args.Units)
    else:
        frequencies_output = frequencies

    save_raman_tensors(
        band_indices, frequencies_output, raman_tensors, file_name,
        frequency_units=args.Units,
        activity_units=_RAMAN_ACTIVITY_UNITS)

    # If linewidths are supplied, discard data for modes for which Raman
    # activities were not calculated.

    if linewidths is not None:
        linewidths = [
            linewidths[index] for index in band_indices
            ]

    # If ir. rep. data is supplied, we need to: (1) re-map the sets of
    # band indices to account for modes that have been discarded; and
    # (2) remove sets of data for which all modes have been discarded.

    index_mapping = {
        index: i for i, index in enumerate(band_indices)
        }

    if irrep_data is not None:
        irrep_data_new = []

        for symbol, indices in irrep_data:
            indices_new = []

            for index in indices:
                if index in index_mapping:
                    indices_new.append(index_mapping[index])

            if len(indices_new) > 0:
                irrep_data_new.append((symbol, indices_new))

        irrep_data = irrep_data_new

    # Hand over to run_mode_output() routine with the scalar Raman
    # intensities.

    run_mode_output(
        'raman', frequencies, frequency_units, raman_intensities,
        _RAMAN_INTENSITY_UNITS, args, linewidths=linewidths,
        irrep_data=irrep_data)


# -----------
# Data Output
# -----------

def run_mode_output(
        spectrum_type, frequencies, frequency_units, intensities,
        intensity_units, args, linewidths=None, irrep_data=None
        ):
    """ Processes input data and arguments and saves a peak table and simulated
    spectrum.

    Arguments:
        spectrum_type -- type of spectrum being output
            ('ir' or 'raman').
        frequencies -- mode frequencies.
        frequency_units -- units of frequencies (must be one of the
            units supported by the convert_frequency_units() function).
        intensities -- band intensities.
        intensity_units -- units of intensities.
        args -- parsed command-line arguments with the options added by
            the update_parser() method in the Parser module.

    Keyword arguments:
        linewidths -- per-mode linewidths.
        irrep_data -- set of (symbol, band_indices) tuples assigning
            irrep symbols to modes.
    """
    # Most of the functions called from the other SpectroscoPy modules
    # have their own error-handling, and we defer to those where
    # possible to reduce "boilerplate" code.

    num_modes = len(frequencies)

    if len(intensities) != num_modes:
        raise Exception(
            "Error: The numbers of frequencies and intensities are "
            "inconsistent.")

    if linewidths is None:
        if args.Linewidth < 0.0:
            raise Exception(
                "Error: Negative mode linewidths are unphysical (!?).")

    if args.Irreps and irrep_data is None:
        raise Exception(
            "Error: If the Irreps argument is set, irrep. data must be "
            "available.")

    # Convert frequency units if required.

    frequency_units_convert = args.Units

    if frequency_units_convert == 'um':
        # Outputting spectra in wavelength units (e.g. um) is supported,
        # but it does not make sense to output peak parameters like
        # this (central frequencies and linewidths are energies). If the
        # user requests a wavelength unit, we work in a default energy
        # unit, output the peak table, simulate the spectrum, then
        # convert just before output.

        frequency_units_convert = _DEFAULT_PEAK_TABLE_FREQUENCY_UNITS

        print(
            "INFO: Peak table frequency units reset to '{0}'."
            .format(frequency_units_convert))

    frequencies = convert_frequency_units(
        frequencies, frequency_units, frequency_units_convert)

    if linewidths is not None:
        linewidths = convert_frequency_units(
            linewidths, frequency_units, frequency_units_convert)
    else:
        # Convert user-supplied nominal linewidth (in inv. cm) to
        # frequency units.

        args.Linewidth, = convert_frequency_units(
            [args.Linewidth], 'inv_cm', frequency_units_convert)

    frequency_units = frequency_units_convert

    file_format = args.DataFormat.lower()

    # Set prefix for output files.

    outputprefix = None

    if spectrum_type == 'ir':
        output_prefix = "IR"
    elif spectrum_type == 'raman':
        output_prefix = "Raman"
    else:
        # We _could_ work around this, but if it happens it's almost
        # certainly a programming error.

        raise Exception(
            "Error: Unknown spectrum_type '{0}' - if you are running "
            "one of the command-line scripts this is probably a bug."
            .format(spectrum_type))

    if args.OutputPrefix is not None:
        output_prefix = "{0}_{1}".format(args.OutputPrefix, output_prefix)

    # Save peak table.

    file_name = "{0}-PeakTable.{1}".format(output_prefix, file_format)

    if args.Irreps:
        pt_frequencies, pt_intensities, pt_irrep_symbols, pt_linewidths = \
            group_for_peak_table(
                frequencies, intensities, irrep_data, linewidths=linewidths)

        save_peak_table(
            pt_frequencies, pt_intensities, file_name,
            irrep_symbols=pt_irrep_symbols, linewidths=pt_linewidths,
            frequency_units=frequency_units, intensity_units=intensity_units,
            file_format=file_format)
    else:
        save_peak_table(
            frequencies, intensities, file_name,
            linewidths=linewidths, frequency_units=frequency_units,
            intensity_units=intensity_units, file_format=file_format)

    # Simulate spectrum and convert units if required.

    spectrum_range = args.SpectrumRange

    if spectrum_range is not None and args.Units != frequency_units:
        spectrum_min, spectrum_max = args.SpectrumRange

        if args.Units == 'um':
            # If the user specifies a spectrum min and/or max <= 0 when
            # using wavelength units, this will cause problems when
            # converting units.

            if spectrum_min <= 0.0 or spectrum_max <= 0.0:
                raise Exception(
                    "Error: If a spectrum range is specified for "
                    "wavelength units, both min and max must be > 0.")

        spectrum_range = convert_frequency_units(
            spectrum_range, args.Units, frequency_units)

    spectrum = simulate_spectrum(
        frequencies, intensities,
        linewidths if linewidths is not None else args.Linewidth,
        spectrum_range=spectrum_range,
        spectrum_resolution=args.SpectrumResolution,
        instrument_broadening_width=args.InstrumentBroadeningWidth,
        instrument_broadening_shape=args.InstrumentBroadeningShape)

    if args.Units != frequency_units:
        spectrum_x, spectrum_y, spectrum_y_norm = spectrum

        if args.Units == 'um':
            # Reverse arrays.

            spectrum_x = spectrum_x[::-1]
            spectrum_y = spectrum_y[::-1]
            spectrum_y_norm = spectrum_y_norm[::-1]

            # Frequencies <= 0 cause problems when converting to
            # wavelength units. Since wavelength units are not
            # anticipated to be a routine use of this program, for now 
            # we simply mask them out and print a warning message for
            # the user.

            mask = spectrum_x > 0.0

            clip = not mask.all()

            if clip:
                spectrum_x = spectrum_x[mask]
                spectrum_y = spectrum_y[mask]
                spectrum_y_norm = spectrum_y_norm[mask]

            spectrum_x = convert_frequency_units(
                spectrum_x, frequency_units, args.Units)

            if clip:
                print(
                    "INFO: Spectrum clipped to ({0:.2f}, {1:.2f})."
                    .format(spectrum_x[0], spectrum_x[-1]))

        spectrum = (spectrum_x, spectrum_y, spectrum_y_norm)

    # Save spectrum.

    file_name = "{0}-Spectrum.{1}".format(outputprefix, file_format)

    save_spectrum(
        spectrum, file_name, frequency_units=args.Units,
        intensity_units=intensity_units, file_format=file_format)

    # Plot spectrum.

    file_name = "{0}-Spectrum.png".format(output_prefix)

    plot_spectrum(
        spectrum, file_name, frequency_units=args.Units,
        plot_type=spectrum_type)
