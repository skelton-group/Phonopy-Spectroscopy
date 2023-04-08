# SpectroscoPy/CLI/IntermediateIO.py


# ---------
# Docstring
# ---------

""" Routines for storing and retrieving intermediate data in multi-step
calculations. """


# -------
# Imports
# -------

import warnings

import yaml

# Try to use the C YAML Loader if possible; if not, fall back to the
# standard routines.

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

import numpy as np

from spectroscopy.utilities import calculate_cell_volume


# ---------
# Constants
# ---------

""" Default step unit written to intermediate YAML files. """

_DEFAULT_STEP_UNIT = "Unspecified"

""" Default distance unit written to intermediate YAML files. """

_DEFAULT_DISTANCE_UNIT = "Unspecified"


# ------------------
# Raman Calculations
# ------------------

""" Default string written to intermediate YAML files when a unit is 
unspecified. """

_RAMAN_DATASET_DEFAULT_UNIT = "Unspecified"


def read_raman_dataset(file_path):
    """ Read data for a Raman calculation from a YAML-format file
    produced by write_raman_dataset().

    Return value:
        A tuple of (lattice_vectors, cell_volume, band_indices,
        frequencies, disp_step_sets, max_disps_sets, eps_tensor_sets,
        distance_unit, frequency_unit, step_unit) data read from
        the file.
    
    Notes:
        * The order of the elements in the return value matches the
            arguments passed to write_raman_dataset() except for
            cell_volume.
        * For compatibility reasons, the lattice_vectors may be None,
            in which case cell_volume should be set.
    """
    # Read and parse data set.

    dataset = None

    with open(file_path, 'r') as input_reader:
        dataset = yaml.load(input_reader, Loader=Loader)

    # Validate keys.
    
    for required_key in [
            'displacement_sets', 'distance_unit', 'frequency_unit', 'step_unit'
            ]:
        if required_key not in dataset:
            # For backwards compatibility.

            error = True

            if required_key.endswith('unit'):
                old_key = required_key + 's'

                if old_key in dataset:
                    dataset[required_key] = dataset[old_key]
                    error = False

            if error:
                raise Exception(
                    "Error: Missing key '{0}' in displacement dataset."
                    .format(required_key))

    if len(dataset['displacement_sets']) == 0:
        raise Exception(
            "Error: Displacement dataset contains no displacement sets.")

    has_eps_tensors = False

    for item_1 in dataset['displacement_sets']:
        for required_key in ['band_index', 'frequency', 'displacements']:
            if required_key not in item_1:
                raise Exception(
                    "Error: Missing key '{0}' in one or more "
                    "displacement sets.".format(required_key))

        for item_2 in item_1['displacements']:
            for required_key in [
                    'displacement_step', 'max_cartesian_displacement'
                ]   :
                if required_key not in item_2:
                    raise Exception(
                        "Error: Missing key '{0}' in one or more "
                        "displacement entries.".format(required_key))

            if 'epsilon_static' in item_2:
                has_eps_tensors = True

    if has_eps_tensors:
        for item_1 in dataset['displacement_sets']:
            for item_2 in item_1['displacements']:
                if 'epsilon_static' not in item_2:
                    raise Exception(
                        "Error: If dielectric tensors are supplied, "
                        "they must be present for all displacements.")
                else:
                    dim_ok = True

                    if len(item_2['epsilon_static']) != 3:
                        dim_ok = False
                    else:
                        for row in item_2['epsilon_static']:
                            if len(row) != 3:
                                dim_ok = False

                    if not dim_ok:
                        raise Exception(
                            "Error: One or more dielectric tensors in "
                            "the data set does not have the correct "
                            "dimension.")

    # For compatibility with older versions of the code, allow
    # the returned lattice_vectors to be None, and calculate the cell
    # volume if the 'volume' key is not present in the YAML file.

    lattice_vectors, cell_volume = None, None

    if 'lattice_vectors' in dataset:
        lattice_vectors = dataset['lattice_vectors']
        cell_volume = calculate_cell_volume(lattice_vectors)
    else:
        if 'cell_volume' in dataset:
            cell_volume = dataset['cell_volume']
        else:
            raise Exception(
                "Error: Dataset must contain one of the "
                "'lattice_vectors' or 'cell_volume' keys.")

    # Extract, reformat and return data.

    band_indices = []
    frequencies = []
    disp_step_sets, max_disps_sets = [], []

    eps_tensor_sets = [] if has_eps_tensors else None

    for item_1 in dataset['displacement_sets']:
        band_indices.append(item_1['band_index'] - 1)

        frequencies.append(item_1['frequency'])

        disp_step_sets.append(
            [item_2['displacement_step'] for item_2 in item_1['displacements']]
            )

        max_disps_sets.append(
            [item_2['max_cartesian_displacement']
                for item_2 in item_1['displacements']])

        if has_eps_tensors:
            eps_tensor_set = [
                np.array(item_2['epsilon_static'], dtype=np.float64)
                    for item_2 in item_1['displacements']
                ]

            eps_tensor_sets.append(eps_tensor_set)

    return (
        lattice_vectors, cell_volume, band_indices, frequencies,
        disp_step_sets, max_disps_sets, eps_tensor_sets,
        dataset['distance_unit'], dataset['frequency_unit'],
        dataset['step_unit']
        )

def write_raman_dataset(
        lattice_vectors, cell_volume, band_indices, frequencies,
        disp_step_sets, max_disps_sets, file_path, eps_tensor_sets=None,
        distance_unit = _RAMAN_DATASET_DEFAULT_UNIT,
        frequency_unit = _RAMAN_DATASET_DEFAULT_UNIT,
        step_unit = _RAMAN_DATASET_DEFAULT_UNIT
        ):
    """ Write band indices, frequencies, displacement data and,
    optionally, calculated dielectric tensors to a YAML-format file.

    Arguments:
        lattice_vectors -- lattice vectors (stored for polarised Raman
            calculations during post processing).
        cell_volume -- unit-cell volume (included for compatibilty and
            will only be written when lattice_vectors is not set).
        band_indices -- indices of the bands (modes) for which data is
            to be output.
        frequencies -- mode frequencies.
        disp_step_sets -- sets of displacements performed along each
            mode in normal-mode coordinates.
        max_disps_dets -- sets of maximum Cartesian atomic displacements
            across all atoms for each displacement.
        file_path -- path to the output file.

    Keyword arguments:
        eps_tensor_sets -- sets of 3x3 dielectric tensors calculated at
            each displacement.
        distance_unit -- unit of distance (default:
            _RAMAN_DATASET_DEFAULT_UNIT).
        frequency_unit -- unit of frequency (default:
            _RAMAN_DATASET_DEFAULT_UNIT).
        step_unit -- unit of the normal-mode amplitude (defaut:
            _RAMAN_DATASET_DEFAULT_UNIT).
    """
    
    with open(file_path, 'w') as output_writer:
        # Write units for physical quantities.

        output_writer.write("distance_units: {0}\n".format(distance_unit))
        output_writer.write("frequency_units: {0}\n".format(frequency_unit))
        output_writer.write("step_units: {0}\n".format(step_unit))

        output_writer.write("\n")

        # Write lattice vectors or cell volume.

        if lattice_vectors is not None:
            output_writer.write("lattice_vectors:\n")

            for v in lattice_vectors:
                output_writer.write(
                    "     - [ {0: 15.8f}, {1: 15.8f}, {2: 15.8f} ]\n".format(*v)
                    )
            
            if cell_volume is not None:
                warnings.warn(
                    "When lattice vectors are specified the unit-cell "
                    "volume is not written to the intermediate YAML "
                    "files.", RuntimeWarning
                    )
        else:
            if cell_volume is None:
                raise Exception(
                    "cell_volume must be specified when "
                    "lattice_vectors is set to None.")
            
            output_writer.write("cell_volume: {0:.4f}\n".format(cell_volume))
        
        output_writer.write("\n")

        # Write displacement data for each mode.

        output_writer.write("displacement_sets:\n")

        for i, band_index in enumerate(band_indices):
            frequency = frequencies[i]
            disp_steps, max_disps = disp_step_sets[i], max_disps_sets[i]

            eps_tensors = (
                eps_tensor_sets[i] if eps_tensor_sets is not None else None)

            output_writer.write("- # {0}\n".format(i + 1))

            # Convert band indices to "one-based".

            output_writer.write("  band_index: {0}\n".format(band_index + 1))

            output_writer.write("  frequency: {0}\n".format(frequency))

            output_writer.write("  displacements:\n")

            # Write displacement steps, max. Cartesian displacements,
            # and calculated dielectric tensors, if available.

            for j in range(0, len(disp_steps)):
                output_writer.write("   - # Step {0}\n".format(j + 1))
                
                output_writer.write(
                    "     displacement_step: {0:.8f}\n"
                    .format(disp_steps[j]))

                output_writer.write(
                    "     max_cartesian_displacement: {0:.8f}\n"
                    .format(max_disps[j]))

                if eps_tensors is not None:
                    output_writer.write("     epsilon_static:\n")

                    for row in eps_tensors[j]:
                        output_writer.write(
                            "     - [ {0: 15.8f}, {1: 15.8f}, {2: 15.8f} ]\n"
                            .format(*row)
                            )
