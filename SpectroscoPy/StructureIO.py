# SpectroscoPy/StructureIO.py


# ---------
# Docstring
# ---------

""" Routines for importing and exporting structural data. """


# -------
# Imports
# -------

import numpy as np;


# ---------
# Constants
# ---------

""" Default title line for POSCAR files written by the WritePOSCAR() function. """

_DefaultTitleLine = "Unknown Structure";


# ------------
# POSCAR Files
# ------------

def WritePOSCAR(structure, filePath, titleLine = None):
    """
    Write a structure to a VASP POSCAR file (VASP 5 format).

    Arguments:
        structure -- a tuple of (lattice_vectors, atomic_symbols, atom_positions_frac) data.
        filePath -- file path to save structure to.

    Keyword arguments:
        titleLine -- title line to include on line 1 of the POSCAR file (default: _DefaultTitleLine constant).

    Notes:
        The atom positions should be in fractional coordinates.
    """

    latticeVectors, atomicSymbols, atomPositions = structure;

    # Group atom positions by symbol.
    # To avoid reordering atom types, we keep track of the order in which new symbols were added.

    atomGroupsSymbols, atomGroups = [], { };

    for symbol, position in zip(atomicSymbols, atomPositions):
        if symbol not in atomGroupsSymbols:
            atomGroups[symbol] = [];
            atomGroupsSymbols.append(symbol);

        atomGroups[symbol].append(position);

    with open(filePath, 'w') as outputWriter:
        # Comment line.

        outputWriter.write(
            "{0}\n".format(titleLine if titleLine != None else _DefaultTitleLine)
            );

        # Write scale factor.

        outputWriter.write("  {0: >19.16f}\n".format(1.0));

        # Write lattice vectors.

        for i in range(0, 3):
            vx, vy, vz = latticeVectors[i];
            outputWriter.write("  {0: >21.16f}  {1: >21.16f}  {2: >21.16f}\n".format(vx, vy, vz));

        # Write atomic symbols and atom counts.

        outputWriter.write(
            "  {0}\n".format(
                "  ".join(["{0: >3}".format(symbol) for symbol in atomGroupsSymbols])
                )
            );

        outputWriter.write(
            "  {0}\n".format(
                "  ".join(["{0: >3}".format(len(atomGroups[symbol])) for symbol in atomGroupsSymbols])
                )
            );

        # Write positions.

        outputWriter.write("Direct\n");

        for symbol in atomGroupsSymbols:
            for f1, f2, f3 in atomGroups[symbol]:
                outputWriter.write("  {0: >21.16f}  {1: >21.16f}  {2: >21.16f}\n".format(f1, f2, f3));


# ---------
# XYZ Files
# ---------

def ReadXYZ(filePath):
    """
    Read a molecular structure from an XYZ-format file.

    Return value:
        A (title_line, atomic_symbols, atom_positions) tuple containing the data.

    Notes:
        For XYZ files containing multiple structures, only the first one is read.
        If the title line in the file is blank, it will be set to the _DefaultTitleLine constant.
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

            titleLine = _DefaultTitleLine;

        # Read atom data.

        atomicSymbols, atomPositions = [], [];

        for i in range(0, atomCount):
            atomicSymbol, x, y, z = next(inputReader).strip().split();

            atomicSymbols.append(atomicSymbol);

            atomPositions.append(
                np.array([float(x), float(y), float(z)], dtype = np.float64)
                );

    return (titleLine, atomicSymbols, atomPositions);
