# SpectroscoPy/CLI/IntermediateIO.py


# ---------
# Docstring
# ---------

""" Routines for storing and retrieving intermediate data in multi-step
calculations. """


# -------
# Imports
# -------

import yaml

# Try to use the C YAML Loader if possible; if not, fall back to the
# standard routines.

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

import numpy as np


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

_RAMAN_DATASET_DEFAULT_UNITS = "Unspecified"


def read_raman_dataset(file_path):
    """ Read data for a Raman calculation from a YAML-format file
    produced by write_raman_dataset().

    Return value:
        A tuple of (cell_volume, volume_units, band_indices,
        frequencies, frequency_units, disp_step_sets, step_units,
        max_disps_sets, distance_units, eps_tensor_sets) data read from
        the file. (The order of the elements matches the arguments
        passed to write_raman_dataset()).
    """
    # Read and parse data set.

    dataset = None

    with open(file_path, 'r') as input_reader:
        dataset = yaml.load(input_reader, Loader=Loader)

    # Validate.

    for required_key in [
            'frequency_units', 'step_units', 'distance_units',
            'volume_units', 'cell_volume', 'displacement_sets'
            ]:
        if required_key not in dataset:
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
        dataset['cell_volume'], dataset['volume_units'],
        band_indices, frequencies, dataset['frequency_units'],
        disp_step_sets, dataset['step_units'], max_disps_sets,
        dataset['distance_units'], eps_tensor_sets
        )

def write_raman_dataset(
        file_path, cell_volume, band_indices, frequencies,
        disp_step_sets, max_disps_sets,  eps_tensor_sets=None,
        volume_units = _RAMAN_DATASET_DEFAULT_UNITS,
        frequency_units = _RAMAN_DATASET_DEFAULT_UNITS,
        step_units = _RAMAN_DATASET_DEFAULT_UNITS,
        distance_units = _RAMAN_DATASET_DEFAULT_UNITS
        ):
    """ Write band indices, frequencies, displacement data and,
    optionally, calculated dielectric tensors to a YAML-format file.

    Arguments:
        file_path -- path to the output file.
        cell_volume -- volume of the unit cell (stored for normalising
            the Raman activity during post processing).
        band_indices -- indices of the bands (modes) for which data is
            to be output.
        frequencies -- mode frequencies.
        disp_step_sets -- sets of displacements performed along each
            mode in normal-mode coordinates.
        max_disps_dets -- sets of maximum Cartesian atomic displacements
            across all atoms for each displacement.

    Keyword arguments:
        eps_tensor_sets -- sets of 3x3 dielectric tensors calculated at
            each displacement.
        volume_units -- units of the cell volume (default:
            _RAMAN_DATASET_DEFAULT_UNITS)
        frequency_units -- units of the frequencies (default:
            _RAMAN_DATASET_DEFAULT_UNITS).
        step_units -- units of the normal-mode coordinate (defaut:
            _RAMAN_DATASET_DEFAULT_UNITS).
        distance_units -- units of the maximum distance (default:
            _RAMAN_DATASET_DEFAULT_UNITS).
    """
    
    with open(file_path, 'w') as output_writer:
        # Write units for physical quantities.

        output_writer.write("frequency_units: {0}\n".format(frequency_units))
        output_writer.write("step_units: {0}\n".format(step_units))
        output_writer.write("distance_units: {0}\n".format(distance_units))
        output_writer.write("volume_units: {0}\n".format(volume_units))

        output_writer.write("\n")

        # Write cell volume.

        output_writer.write("cell_volume: {0: 12.4f}\n".format(cell_volume))
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
                    "     displacement_step: {0: 12.8f}\n"
                    .format(disp_steps[j]))

                output_writer.write(
                    "     max_cartesian_displacement: {0: 12.8f}\n"
                    .format(max_disps[j]))

                if eps_tensors is not None:
                    output_writer.write("     epsilon_static:\n")

                    for k in range(0, 3):
                        output_writer.write(
                            "     - [  {0: 15.8f},  {1: 15.8f},  {2: 15.8f}  "
                            "]\n".format(*[item for item in eps_tensors[j][k]])
                            )
