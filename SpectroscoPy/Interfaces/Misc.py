# SpectroscoPy/Interfaces/Misc.py


# ---------
# Docstring
# ---------

""" Routines for reading and writing various types of file not associated with particular codes. """


# -------
# Imports
# -------

import numpy as np;


# ---------
# XYZ Files
# ---------

_XYZ_DefaultTitleLine = "Unknown Molecule";

def ReadXYZ(filePath):
    """
    Read a molecular structure from an XYZ-format file.

    Return value:
        A (title_line, atomic_symbols, atom_positions) tuple containing the data.

    Notes:
        For XYZ files containing multiple structures, only the first one is read.
        If the title line in the file is blank, it will be set to the _XYZ_DefaultTitleLine constant.
    """

    titleLine = None;
    atomicSymbols, atomPositions = None, None;

    with open(filePath, 'r') as inputReader:
        # Read atom count.

        atomCount = int(
            next(inputReader).strip()
            );

        # Read title line.

        titleLine = next(inputReader).strip();

        if titleLine == "":
            # If the title line is empty, replace it with the default.

            titleLine = _XYZ_DefaultTitleLine;

        # Read atom data.

        atomicSymbols, atomPositions = [], [];

        for i in range(0, atomCount):
            atomicSymbol, x, y, z = next(inputReader).strip().split();

            atomicSymbols.append(atomicSymbol);

            atomPositions.append(
                np.array([float(x), float(y), float(z)], dtype = np.float64)
                );

    return (titleLine, atomicSymbols, atomPositions);
