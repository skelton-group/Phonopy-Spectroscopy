# spectroscopy/interfaces/misc.py


# ---------
# Docstring
# ---------

""" Routines for reading and writing various types of file not
associated with particular codes. """


# -------
# Imports
# -------

import numpy as np


# ---------
# XYZ Files
# ---------

_XYZ_DEFAULT_TITLE_LINE = "Unknown Molecule"


def read_xyz(file_path):
    """ Read a molecular structure from an XYZ-format file.

    Return value:
        A (title_line, atomic_symbols, atom_positions) tuple containing
        the data.

    Notes:
        For XYZ files containing multiple structures, only the first one
        is read.
        If the title line in the file is blank, it will be set to the
        _XYZ_DEFAULT_TITLE_LINE constant.
    """

    title_line = None
    atomic_symbols, atom_positions = None, None

    with open(file_path, 'r') as input_reader:
        # Read atom count.

        atom_count = int(next(input_reader).strip())

        # Read title line.

        title_line = next(input_reader).strip()

        if title_line == "":
            # If the title line is empty, replace it with the default.

            title_line = _XYZ_DEFAULT_TITLE_LINE

        # Read atom data.

        atomic_symbols, atom_positions = [], []

        for i in range(0, atom_count):
            atomic_symbol, x, y, z = next(input_reader).strip().split()

            atomic_symbols.append(atomic_symbol)

            atom_positions.append(
                np.array([float(x), float(y), float(z)], dtype=np.float64))

    return title_line, atomic_symbols, atom_positions
