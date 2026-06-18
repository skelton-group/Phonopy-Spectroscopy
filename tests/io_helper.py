# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Toutines for input/output as part of unit tests."""

# -------
# Imports
# -------

import numpy as np

# ---------
# Functions
# ---------


def generate_fd_raman_dielectric_input_file_list(
    band_inds, num_steps, file_path_template
):
    """Generate a list of dielectric constanst/function input files for
    a finite-displacement Raman calculation.

    Parameters
    ----------
    band_inds : array_like of int
        Band indices in calculation.
    num_steps : int
        Number of displacement steps in calculation.
    file_path_template : str
        Template for file paths, used as
        `file_path_template.format(band_idx, step_idx)`.

    Returns
    -------
    file_list : list of str
        List of input files for the calculation.
    """

    file_list = []

    for band_idx in band_inds:
        for j in range(num_steps):
            file_list.append(file_path_template.format(band_idx + 1, j + 1))

    return file_list
