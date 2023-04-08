# spectroscopy/cli/runtime.py


# TODO: Sort out unit handling.
# TODO:     3: Move code to "rebase" band indices to ireps to another
# TODO:         module.
# TODO:     4. If using "energy" units, convert at earliest opporrunity
# TODO:         so that Raman tensors/peak tables are in requested
# TODO:         units.
# TODO:     5. For "wavelength" units, simulate spectrum in "energy"
# TODO:         units then "clip" before output.


# ---------
# Docstring
# ---------

"""
Implements standard runtime routines.

Notes:
    * The following input units are assumed...:
        frequencies/linewidths -- THz.
        distances and volumes -- Ang/Ang^3.
        eigendisplacements -- sqrt(amu)^-1.
        charges -- |e|.
    * ... resulting in the following units for derived quantities:
        normal-mode amplitudes -- sqrt(amu) * Ang
        infrared intensities -- e^2 amu^-1
        Raman activity tensors -- sqrt(amu) * Ang^2
        Raman intensities -- Ang^4 amu^-1
"""


# -------
# Imports
# -------

import math
import os

import numpy as np

from spectroscopy.ir import calculate_ir_intensity

from spectroscopy.raman import (
    generate_displaced_structures, calculate_raman_activity_tensor,
    scalar_average_raman_activity_tensor)

from spectroscopy.spectrum import simulate_spectrum

from spectroscopy.text_export import (
    group_for_peak_table, save_peak_table, save_scalar_spectrum_or_spectra,
    save_raman_intensity_theta)

from spectroscopy.yaml_export import save_raman_tensors

from spectroscopy.plotting import (
    plot_scalar_spectrum_or_spectra, plot_2d_polarised_raman_spectrum)

from spectroscopy.constants import get_irrep_activities

from spectroscopy.units import (
    get_distance_unit_text_label, get_frequency_unit_text_label,
    get_frequency_unit_plot_label, convert_frequency_units)

from spectroscopy.utilities import (
    hkl_to_real_space_normal, rotation_matrix_from_vectors, rotation_matrix_xy,
    rotate_tensors)

from spectroscopy.cli.intermediate_io import (
    read_raman_dataset, write_raman_dataset)


# ---------
# Constants
# ---------

""" Unit label for normal-mode amplitudes for text output. """

_MODE_AMPLITUDE_UNITS_TEXT_LABEL = "sqrt(amu) * Ang"

""" Unit label for IR intensities for text output. """

_IR_INTENSITY_UNITS_TEXT_LABEL = "e^2 amu^-1"

""" Unit label for IR intensities for plots. """

_IR_INTENSITY_UNITS_PLOT_LABEL = "$e^{2}$ amu$^{-1}$"

""" Unit label for (tensor) Raman activities for text output. """

_RAMAN_ACTIVITY_UNITS_TEXT_LABEL = "sqrt(amu) * Ang^2"

""" Unit label for (scalar) Raman intensities for text output. """

_RAMAN_INTENSITY_UNITS_TEXT_LABEL = r"Ang^4 amu^-1"

""" Unit label for (scalar) Raman intensities for plots. """

_RAMAN_INTENSITY_UNITS_PLOT_LABEL = r"$\mathrm{\AA}^{4}$ amu$^{-1}$"

""" Fallback measurement units used when the user requests a wavelength
unit and it does not make sense to use one. """

_WAVELENGTH_FALLBACK_MEASUREMENT_UNITS = 'inv_cm'


# ---------------
# IR Spectroscopy
# ---------------

def run_mode_ir(
        frequencies, eigendisplacements, bec_tensors, args,
        linewidths=None, irrep_data=None
        ):
    """ Calculates IR intensities and passes data to _output_spectrum()
    to generate and output peak table and simulated spectrum.

    Arguments:
        frequencies -- Gamma-point mode frequencies.
        eigendisplacements -- Gamma-point eigendisplacements.
        bec_tensors -- Born effective-charge tensors.
            assigning mode irreps.
        args -- parsed command-line arguments with parameters
            for _output_spectrum().
    
    Keyword arguments:
        linewidths -- (optional) mode linewidths.
        irrep_data -- (optional) list of (symbol, band_indices) tuples.
    """

    # Calculate IR intensities.

    ir_intensities = [
        calculate_ir_intensity(eigendisplacement, bec_tensors)
            for eigendisplacement in eigendisplacements
        ]

    # Hand over to the _output_scalar_spectrum() routine.

    _output_scalar_spectrum_or_spectra(
        'ir', frequencies, ir_intensities, args,
        linewidths=linewidths, irrep_data=irrep_data)


# ------------------
# Raman Spectroscopy
# ------------------

def run_mode_raman_disp(
        structure, frequencies, eigendisplacements, args,
        point_group = None, irrep_data = None):
    """ Prepare and return a sequence of displaced structures for a
    Raman activity calculation, and write metadata to an intermediate
    Raman.yaml file.

    Arguments:
        structure -- tuple of (lattice_vectors, atomic_symbols,
            atom_positions_frac) specifying the structure.
        frequencies -- Gamma-point mode frequencies.
        eigendisplacements -- Gamma-point eigendisplacements.
        args -- parsed command-line arguments with the RamanStep,
            RamanStepType, RamanBandIndices and OutputPrefix
            attributes.
        point_group -- (optional) point group of the structure.
        irrep_data -- (optional) irreducible representations of the
            modes.

    Return value:
        A set of (band_index, band_frequency, structures, disp_steps,
        max_disps) tuples containing the displaced structures for each
        mode for the calling code to output.
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
            "displacement set as these appear to be acoustic modes "
            "(to override use --raman-bands)."
            .format(*[index + 1 for index in exclude_indices]))
        
        # If irrep data is supplied, exclude modes with irreps that
        # are not Raman active.

        if point_group is not None and irrep_data is not None:
            irrep_activities = get_irrep_activities(point_group, 'raman')

            if irrep_activities is not None:
                active_irreps, inactive_irreps = irrep_activities

                for symbol, band_indices in irrep_data:
                    if symbol is not None:
                        if symbol in inactive_irreps:
                            exclude_indices.extend(band_indices)
                        elif symbol not in active_irreps:
                            # This is (hopefully!) mostly for debugging
                            # purposes...

                            print(
                                "WARNING: Unknown irrep '{0}' for "
                                "point group '{1}'."
                                .format(symbol, point_group))

                if len(exclude_indices) > 3:
                    band_list_str = ", ".join(
                        "{0}".format(index + 1)
                            for index in exclude_indices[3:])

                    print(
                        "INFO: Automatically excluding band(s) {0} "
                        "from displacement set as these appear to be "
                        "Raman inactive (to override use "
                        "--raman-bands).".format(band_list_str))
            else:
                # ... as is this.

                print(
                    "WARNING: Unknown point group '{0}'."
                    .format(point_group))

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
        lattice_vectors, None, band_indices, frequencies,
        mode_disp_steps, mode_max_disps, file_name, eps_tensor_sets=None,
        distance_unit=get_distance_unit_text_label('ang'),
        frequency_unit=get_frequency_unit_text_label('thz'),
        step_unit=_MODE_AMPLITUDE_UNITS_TEXT_LABEL
        )

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

    Metadata on the displaced structures is read from the intermediate
    Raman.yaml file and updated with the calculated dielectric tensors.

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

    # Read displacement data set from intermediate Raman.yaml file and
    # verify hard-coded units.

    lattice_vectors, cell_volume, band_indices, frequencies, disp_step_sets, \
        max_disps_sets, _, distance_unit, frequency_unit, step_unit \
        = read_raman_dataset(dataset_file)
    
    _raman_check_dataset_units(
        distance_unit, frequency_unit, step_unit, dataset_file)

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

    # For backwards compatibility.

    if lattice_vectors is not None:
        cell_volume = None

    write_raman_dataset(
        lattice_vectors, cell_volume, band_indices, frequencies,
        disp_step_sets, max_disps_sets, dataset_file, eps_tensor_sets,
        distance_unit, frequency_unit, step_unit)

def run_mode_raman_postproc(args, linewidths=None, irrep_data=None):
    """ Post process a data set prepared for a Raman calculation.

    A data set containing frequencies, displacement steps and dielectric
    tensors is read from the intermediate Raman.yaml file, and the Raman
    activity tensors are calculated from the derivatives of the static
    dielectric tensor with respect to the normal-mode coordinates.
    
    Depending on args, the activity tensors are orientationally averaged
    to obtain a powder/gas spectrum, or are combined with user-input
    parameters for a polarised Raman calculation.

    Arguments:
        args -- parsed command-line arguments with the optional
            OutputPrefix attribute and parameters to be passed to
            the _raman_postproc_*() functions.

    Keyword arguments:
        linewidths -- (optional) mode linewidths.
        irrep_data -- (optional) list of (symbol, band_indices) tuples.
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

    lattice_vectors, cell_volume, band_indices, frequencies, disp_step_sets, \
        _, eps_tensor_sets, distance_unit, frequency_unit, step_unit \
        = read_raman_dataset(dataset_file)

    _raman_check_dataset_units(
        distance_unit, frequency_unit, step_unit, dataset_file)

    # Calculate Raman activity tensors and write to file.

    raman_tensors = []

    for disp_steps, eps_tensors in zip(disp_step_sets, eps_tensor_sets):
        raman_tensor = calculate_raman_activity_tensor(
            eps_tensors, disp_steps, cell_volume, 'finite_difference')

        raman_tensors.append(raman_tensor)

    # If linewidths are supplied, discard data for modes for which Raman
    # activities were not calculated.

    if linewidths is not None:
        linewidths = [linewidths[index] for index in band_indices]

    # If irrep data is supplied, we need to: (1) re-map the sets of band
    # indices to account for modes that have been discarded; and (2)
    # remove sets of data for which all modes have been discarded.

    index_mapping = {index: i for i, index in enumerate(band_indices)}

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

    # Depending on args, hand over to another _raman_* routines to
    # process the Raman tensors.

    if args.SurfaceHKL is not None:
        _raman_postproc_pol(
            band_indices, frequencies, raman_tensors, linewidths, irrep_data,
            lattice_vectors, args)
    else:
        _raman_postproc_scalar(
            band_indices, frequencies, raman_tensors, linewidths, irrep_data,
            args)

def _raman_check_dataset_units(
        distance_unit, frequency_unit, step_unit, file_path):
    """ Validates units in the intermediate Raman.yaml file produced
     and read by the run_mode_raman_* routines. """

    if distance_unit != get_distance_unit_text_label('ang'):
        raise Exception(
            "Error: Unexpected distance unit '{0}' read from "
            "displacement dataset file \"{1}\"."
            .format(distance_unit, file_path))

    if frequency_unit != get_frequency_unit_text_label('thz'):
        raise Exception(
            "Error: Unexpected frequency unit '{0}' read from "
            "displacement dataset file \"{1}\"."
            .format(frequency_unit, file_path))
    
    if step_unit != _MODE_AMPLITUDE_UNITS_TEXT_LABEL:
        raise Exception(
            "Error: Unexpected step unit '{0}' read from "
            "displacement dataset file \"{1}\"."
            .format(step_unit, file_path))

def _raman_postproc_scalar(
        band_indices, frequencies, raman_tensors, linewidths,
        irrep_data, args):
    """ Post process the Raman activity tensors generated by
    run_mode_raman_postproc() to generate an orientationally-averaged
    (powder) spectrum.

    Main arguments:
        raman_tensors -- Raman activity tensors.
    
    All other arguments are passed through to the output_spectrum()
    function.
    """
    # Write out Raman tensors.

    file_name = "Raman-Tensors.yaml"

    if args.OutputPrefix is not None:
        file_name = "{0}_{1}".format(args.OutputPrefix, file_name)

    save_raman_tensors(
        band_indices, convert_frequency_units(frequencies, 'thz', args.Units),
        raman_tensors, args.Units, _RAMAN_ACTIVITY_UNITS_TEXT_LABEL, file_name)

    # Calculate orientationally-averaged Raman intensities.

    raman_intensities = [
        scalar_average_raman_activity_tensor(t)
            for t in raman_tensors
        ]

    # Hand over to output_spectrum() routine with the scalar Raman
    # intensities.

    _output_scalar_spectrum_or_spectra(
        'raman', frequencies, raman_intensities, args,
        linewidths=linewidths, irrep_data=irrep_data)

def _raman_postproc_pol(
        band_indices, frequencies, raman_tensors, linewidths, irrep_data,
        lattice_vectors, args
        ):
    """ Post process the Raman activity tensors generated by
    _run_mode_raman_postproc() to perform a polarised Raman simulation.

    Main arguments:
        raman_tensors -- Raman activity tensors.
        args -- parsed command-line arguments including parameters
            for polarised Raman simulation.
    
    All other arguments are passed through to the output_spectrum()
    function.
    """

    if lattice_vectors is None:
        raise Exception(
            "Error: Polarised Raman calculations cannot be performed "
            "using old-format Raman.yaml files - please regenerate the "
            "Raman.yaml file by repeating the -d and -r steps (you "
            "should _not_ need to re-run the first-principles "
            "calculations again.")

    # Calculate (real-space) surface normal for supplied (h, k, l) and
    # rotation matrix to orient surface along the z direction.

    surface_normal = hkl_to_real_space_normal(args.SurfaceHKL, lattice_vectors)
    
    surface_rot_matrix = rotation_matrix_from_vectors(
        surface_normal, [0.0, 0.0, 1.0])

    # Base file name for Raman tensors.

    file_name_base = "Raman-Pol-Tensors-{0}{1}{2}".format(*args.SurfaceHKL)

    if args.OutputPrefix is not None:
        file_name_base = "{0}_{1}".format(args.OutputPrefix, file_name_base)

    # Output the Raman tensors with the surface rotation applied.

    surface_raman_tensors = rotate_tensors(surface_rot_matrix, raman_tensors)

    file_name = "{0}.yaml".format(file_name_base)

    save_raman_tensors(
        band_indices,
        convert_frequency_units(frequencies, 'thz', args.Units),
        surface_raman_tensors, args.Units,
        _RAMAN_ACTIVITY_UNITS_TEXT_LABEL, file_name,
        surface_hkl = args.SurfaceHKL, surface_rotation = surface_rot_matrix)

    # If the user has specified one or more theta, generate intensities
    # and spectra for those; otherwise, generate for a full 360 deg
    # rotation.

    theta_vals = args.Theta

    if theta_vals is None:
        theta_vals = np.arange(
            0.0, 360.0 + args.ThetaStep / 10.0, args.ThetaStep)

    intensity_sets = []

    for theta in theta_vals:
        theta_rot_matrix = rotation_matrix_xy(theta)
        combined_rot_matrix = np.dot(theta_rot_matrix, surface_rot_matrix)

        theta_raman_tensors = rotate_tensors(
            combined_rot_matrix, raman_tensors)
        
        # If this is a user-specified theta, output the Raman tensors.

        if args.Theta is not None:
            file_name = "{0}-{1:.1f}deg.yaml".format(file_name_base, theta)

            save_raman_tensors(
                band_indices,
                convert_frequency_units(frequencies, 'thz', args.Units),
                theta_raman_tensors, args.Units,
                _RAMAN_ACTIVITY_UNITS_TEXT_LABEL, file_name,
                surface_hkl = args.SurfaceHKL, theta = theta,
                surface_rotation = surface_rot_matrix,
                theta_rotation = theta_rot_matrix,
                combined_rotation = combined_rot_matrix)

        theta_intensities = [
            np.dot(args.IncidentPol, np.dot(t, args.ScatteredPol)) ** 2
                for t in theta_raman_tensors
            ]

        intensity_sets.append(theta_intensities)
    
    # Depending on whether or not the user specified theta values, hand
    # over to one of the output routines.

    if args.Theta is not None:
        theta_text_labels = [
            "theta = {0:.1f} deg".format(theta) for theta in theta_vals]
        
        theta_plot_labels = [
            "{0:.1f} $^\circ$".format(theta) for theta in theta_vals]
        
        _output_scalar_spectrum_or_spectra(
            'raman-pol', frequencies, intensity_sets, args, linewidths=linewidths,
            irrep_data=irrep_data,intensity_set_text_labels=theta_text_labels,
            intensity_set_plot_labels=theta_plot_labels)
    else:
        _output_2d_polarised_raman_spectrum(
            frequencies, theta_vals, intensity_sets, args,
            linewidths=linewidths, irrep_data=irrep_data)

# -----------
# Data Output
# -----------

def _output_scalar_spectrum_or_spectra(
        spectrum_type, frequencies, intensities_or_intensity_sets, args,
        linewidths=None, irrep_data=None, intensity_set_text_labels=None,
        intensity_set_plot_labels=None
        ):
    """ Processes input data and arguments and saves a peak table and
    simulated single (scalar) spectrum.

    Arguments:
        spectrum_type -- type of spectrum being output ('ir', 'raman'
            or 'raman-pol').
        frequencies -- mode frequencies.
        intensities_or_intensity_sets -- band intensities or list of
            sets of band intensities.
        args -- parsed command-line arguments with the options added by
            the update_parser() method in the parser module.

    Keyword arguments:
        linewidths -- per-mode linewidths.
        irrep_data -- set of (symbol, band_indices) tuples assigning
            irrep symbols to modes.
        intensity_set_text_labels, intensity_set_plot_labels -- labels
            for each intensity set to use in text output/plots.
    """
    # Most of the functions called from the other SpectroscoPy modules
    # have their own error-handling, and we defer to those where
    # possible to reduce "boilerplate" code.

    frequencies = convert_frequency_units(frequencies, 'thz', args.Units)

    if linewidths is not None:
        linewidths = convert_frequency_units(linewidths, 'thz', args.Units)

    dim = np.ndim(intensities_or_intensity_sets)

    if dim == 1:
        intensities_or_intensity_sets = [intensities_or_intensity_sets]

    file_format = args.DataFormat.lower()

    # Set variables based on plot type for output files.

    output_prefix = None
    intensity_unit_text_label = None
    intensity_unit_plot_label = None
    plot_line_colour_or_colours = None

    if spectrum_type == 'ir':
        output_prefix = "IR"
        intensity_unit_text_label = _IR_INTENSITY_UNITS_TEXT_LABEL
        intensity_unit_plot_label = _IR_INTENSITY_UNITS_PLOT_LABEL
        plot_line_colour_or_colours = 'r'
    elif spectrum_type == 'raman':
        output_prefix = "Raman"
        intensity_unit_text_label = _RAMAN_INTENSITY_UNITS_TEXT_LABEL
        intensity_unit_plot_label = _RAMAN_INTENSITY_UNITS_PLOT_LABEL
        plot_line_colour_or_colours = 'b'
    elif spectrum_type == 'raman-pol':
        output_prefix = "Raman-Pol"
        intensity_unit_text_label = _RAMAN_INTENSITY_UNITS_TEXT_LABEL
        intensity_unit_plot_label = _RAMAN_INTENSITY_UNITS_PLOT_LABEL
    else:
        # We _could_ work around this, but if it happens it's almost
        # certainly a programming error.

        raise Exception(
            "Error: Unknown spectrum_type '{0}' - if you are running "
            "one of the command-line scripts this is probably a bug."
            .format(spectrum_type))

    if len(intensities_or_intensity_sets) > 1:
        pass # TODO

    if args.OutputPrefix is not None:
        output_prefix = "{0}_{1}".format(args.OutputPrefix, output_prefix)

    # If irrep data is provided, group modes (this will sum the
    # intensities of degenerate modes, so will not affect the spectrum).

    file_name = "{0}-PeakTable.{1}".format(output_prefix, file_format)

    irrep_symbols = None

    if irrep_data is not None:
        frequencies_new, linewidths_new = None, None
        intensity_sets_new = []

        for intensities in intensities_or_intensity_sets:
            frequencies_new, intensities, irrep_symbols, linewidths_new = \
                group_for_peak_table(
                    frequencies, intensities, irrep_data, linewidths=linewidths)
            
            intensity_sets_new.append(intensities)
        
        frequencies, linewidths = frequencies_new, linewidths_new
        intensities_or_intensity_sets = intensity_sets_new

    # Save peak table.

    save_peak_table(
        frequencies, intensities_or_intensity_sets,
        get_frequency_unit_text_label(args.Units), intensity_unit_text_label,
        file_name, intensity_set_labels=intensity_set_text_labels,
        irrep_symbols=irrep_symbols, linewidths=linewidths,
        file_format=file_format)

    # Simulate spectrum.

    v, i_v_sets = None, []

    for intensities in intensities_or_intensity_sets:
        v, i_v = simulate_spectrum(
            frequencies, intensities,
            linewidths if linewidths is not None else args.Linewidth,
            spectrum_range=args.SpectrumRange,
            spectrum_resolution=args.SpectrumResolution,
            instrument_broadening_width=args.InstrumentBroadeningWidth,
            instrument_broadening_shape=args.InstrumentBroadeningShape)
        
        i_v_sets.append(i_v)
    
    # Save spectrum.

    file_name = "{0}-Spectrum.{1}".format(output_prefix, file_format)

    save_scalar_spectrum_or_spectra(
        v, i_v_sets, get_frequency_unit_text_label(args.Units),
        intensity_unit_text_label, file_name,
        intensity_set_labels=intensity_set_text_labels,
        file_format=file_format)

    # Plot spectrum.

    file_name = "{0}-Spectrum.png".format(output_prefix)

    plot_scalar_spectrum_or_spectra(
        v, i_v_sets, get_frequency_unit_plot_label(args.Units),
        intensity_unit_plot_label, file_name, normalise=True,
        legend_labels=intensity_set_plot_labels,
        line_colour_or_colours=plot_line_colour_or_colours)

def _output_2d_polarised_raman_spectrum(
        frequencies, theta_vals, intensity_sets, args, linewidths=None,
        irrep_data=None):
    """ Processes a set of frequencies, intensities as a function of the
    rotation angle theta and arguments and saves a polar plot for each
    mode plus a 2D spectrum showing the intensity profile as a function
    of frequency and theta.

    Arguments:
        frequencies -- mode frequencies.
        theta_vals -- values of theta.
        intensity_sets -- a set of band intensities for each value of
            theta
        args -- parsed command-line arguments with the options added by
            the update_parser() method in the parser module.

    Keyword arguments:
        linewidths -- per-mode linewidths.
        irrep_data -- set of (symbol, band_indices) tuples assigning
            irrep symbols to modes.
    """

    # If irrep data is supplied, combine intensities.

    irrep_symbols = None

    if irrep_data is not None:
        frequencies_new, linewidths_new = None, None
        
        for i, intensities in enumerate(intensity_sets):
            frequencies_new, intensities, irrep_symbols, linewidths_new = \
                group_for_peak_table(
                    frequencies, intensities, irrep_data,
                    linewidths=linewidths)
            
            intensity_sets[i] = intensities
        
        frequencies, linewidths = frequencies_new, linewidths_new

    # Convert frequency units if required.

    frequencies = convert_frequency_units(frequencies, 'thz', args.Units)

    if linewidths is not None:
        linewidths = convert_frequency_units(linewidths, 'thz', args.Units)

    # Data format and prefix for output files.

    file_format = args.DataFormat.lower()

    output_prefix = "Raman-Pol"

    if args.OutputPrefix is not None:
        output_prefix = "{0}_{1}".format(args.OutputPrefix, output_prefix)

    # Outut a polar plot for each mode.

    file_name = "{0}-IntensityTheta.{1}".format(output_prefix, file_format)

    save_raman_intensity_theta(
        theta_vals, intensity_sets, "deg", _RAMAN_INTENSITY_UNITS_TEXT_LABEL,
        file_name, file_format=file_format)

    # TODO: Polar plots.

    # Simulate (1D) spectrum for each set of intensities (i.e. each
    # value of theta).

    spectrum_x_ref, spectra_set = None, []
    
    for intensities in intensity_sets:
        spectrum_x, spectrum_y = simulate_spectrum(
            frequencies, intensities,
            linewidths if linewidths is not None else args.Linewidth,
            spectrum_range=args.SpectrumRange,
            spectrum_resolution=args.SpectrumResolution,
            instrument_broadening_width=args.InstrumentBroadeningWidth,
            instrument_broadening_shape=args.InstrumentBroadeningShape)
    
        if spectrum_x_ref is not None:
            if not (spectrum_x == spectrum_x_ref).all():
                raise Exception(
                    "Error: Inconsistent frequencies returned by "
                    "simulate_spectrum() - this is a bug.")
        else:
            spectrum_x_ref = spectrum_x

        spectra_set.append(spectrum_y)
    
    # Plot 2D spectrum.

    file_name = "{0}-2DSpectrum.{1}".format(output_prefix, file_format)

    spectra_labels = [
        "theta = {0:.1f} deg".format(theta) for theta in theta_vals]

    save_scalar_spectrum_or_spectra(
        spectrum_x_ref, spectra_set, get_frequency_unit_text_label(args.Units),
        _RAMAN_INTENSITY_UNITS_TEXT_LABEL, file_name,
        intensity_set_labels=spectra_labels, file_format=file_format)

    file_name = "{0}-2DSpectrum.png".format(output_prefix)

    plot_2d_polarised_raman_spectrum(
        spectrum_x_ref, theta_vals, spectra_set,
        get_frequency_unit_plot_label(args.Units),  r"$^\circ$", file_name,
        normalise=True)