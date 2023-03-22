# SpectroscoPy/Interfaces/Misc.py


# ---------
# Docstring
# ---------

""" Routines for reading and writing various types of file not associated
with particular codes. """


# -------
# Imports
# -------

import numpy as np


# ---------
# XYZ Files
# ---------

_XYZ_DefaultTitleLine = "Unknown Molecule"


def readxyz(filepath):
    """
    Read a molecular structure from an XYZ-format file.

    Return value:
        A (title_line, atomic_symbols, atom_positions) tuple containing the
        data.

    Notes:
        For XYZ files containing multiple structures, only the first one is
        read.
        If the title line in the file is blank, it will be set to the
        _XYZ_DefaultTitleLine constant.
    """

    titleline = None
    atomicsymbols, atompositions = None, None

    with open(filepath, 'r') as inputReader:
        # Read atom count.

        atomcount = int(
            next(inputReader).strip()
            )

        # Read title line.

        titleline = next(inputReader).strip()

        if titleline == "":
            # If the title line is empty, replace it with the default.

            titleline = _XYZ_DefaultTitleLine

        # Read atom data.

        atomicsymbols, atompositions = [], []

        for i in range(0, atomcount):
            atomicsymbol, x, y, z = next(inputReader).strip().split()

            atomicsymbols.append(atomicsymbol)

            atompositions.append(
                np.array([float(x), float(y), float(z)], dtype=np.float64)
                )

    return titleline, atomicsymbols, atompositions
