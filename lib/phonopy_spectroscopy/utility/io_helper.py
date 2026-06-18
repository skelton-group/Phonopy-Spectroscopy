# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines for input/output."""

# -------
# Imports
# -------

import json

import yaml

# This code tries to use the C version of the YAML Loader and falls back
# to the standard Loader if this is not available. The code was taken
# from https://pyyaml.org/wiki/PyYAMLDocumentation.

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

# ---------
# Functions
# ---------


def load_json(file_path):
    """Load a JSON-format file and return a Python data structure.

    Parameters
    ----------
    file_path : str
        File path.

    Returns
    -------
    data : any
        Python object containing parsed data.
    """

    with open(file_path, "r") as f:
        return json.load(f)


def save_json(obj, file_path):
    """Save a Python data structure to a JSON-format file.

    Parameters
    ----------
    obj : any
        Python object to save.
    file_path : str
        File path.
    """

    with open(file_path, "w") as f:
        json.dump(obj, f)


def load_yaml(file_path):
    """Load a YAML-format file and return a Python data structure.

    Parameters
    ----------
    file_path : str
        File path.

    Returns
    -------
    data : any
        Python object containing parsed data.
    """

    with open(file_path, "rb") as f:
        return yaml.load(f, Loader=Loader)
