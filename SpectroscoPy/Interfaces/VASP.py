# SpectroscoPy/Interfaces/VASP.py


# ---------
# Docstring
# ---------

""" Routines for interfacing with the Vienna Ab-initio Simulation Package
(VASP) code. """


# -------
# Imports
# -------

import math
import re

import numpy as np


# ------------
# POSCAR Files
# ------------

""" Default title line for POSCAR files written by the WritePOSCAR() 
function. """

_POSCAR_DefaultTitleLine = "Unknown Structure"


def readposcar(filepath):
    """
    Read a crystal structure from a VASP POSCAR file.

    Return value:
        A tuple of (title_line, lattice_vectors, atomic_symbols,
        atom_positions_frac) data.

    Notes:
        If the title line is blank, it will be set to the
        _POSCAR_DefaultTitleLine constant.
        This function only supports VASP 5.x-format POSCAR files (i.e. with
        atomic symbols as well as counts given in the header).
        Negative scale factors on line 2, interpreted by VASP as the volume
        of the cell, are not supported.
    """

    titleline = None
    latticevectors = None
    atomicsymbols, atompositions = None, None

    with open(filepath, 'r') as inputReader:
        # Read title from line 1.

        titleline = next(inputReader).strip()

        if titleline == "":
            titleline = _POSCAR_DefaultTitleLine

        # Read scale factor.

        scalefactor = float(
            next(inputReader).strip()
            )

        if scalefactor < 0.0:
            raise Exception("Error: Unsupported negative scale factor in "
                            "input file \"{0}\".".format(filepath)
                            )

        # Read lattice vectors.

        latticevectors = []

        for i in range(0, 3):
            latticevectors.append(
                [float(item) for item in next(inputReader).strip().split()[:3]]
                )

        latticevectors = [
            scalefactor * np.array(v, dtype=np.float64)
            for v in latticevectors
            ]

        # Read atom types and positions.

        atomtypes = next(inputReader).strip().split()

        atomcounts = None

        try:
            atomcounts = [
                int(item) for item in next(inputReader).strip().split()
                ]
        except ValueError:
            raise Exception("Error: Failed to read atom counts from input "
                            "file \"{0}\".".format(filepath)
                            )

        if len(atomtypes) != len(atomcounts):
            raise Exception("Error: Inconsistent number of atom types and "
                            "atom counts in input file \"{0}\".".format(
                                                            filepath)
                            )

        atomicsymbols = []

        for atomtype, atomcount in zip(atomtypes, atomcounts):
            atomicsymbols = atomicsymbols + [atomtype] * atomcount

        # Check for direct or Cartesian coordinates.

        key = next(inputReader).strip()[0].lower()

        if key == 's':
            # Selective dynamics keyword -> read first letter of the
            # following line.

            key = next(inputReader).strip()[0].lower()

        coordinatetype = None

        if key == 'c' or key == 'k':
            coordinatetype = 'cartesian'
        elif key == 'd':
            coordinatetype = 'fractional'
        else:
            raise Exception("Error: Unknown coordinate type in input file "
                            "\"{0}\".".format(filepath)
                            )

        # Read atomic positions.

        atompositions = []

        for i in range(0, sum(atomcounts)):
            atompositions.append(
                [float(item) for item in next(inputReader).strip().split()[:3]]
                )

        atompositions = [
            np.array(v, dtype=np.float64) for v in atompositions
            ]

        if coordinatetype == 'cartesian':
            # Scale and convert to fractional coordinates.

            atompositions = Utilities.cartesiantofractionalcoordinates(
                [scalefactor * v for v in atompositions], latticevectors
                )

    return titleline, latticevectors, atomicsymbols, atompositions


def writeposcar(structure, filepath, titleline=None):
    """
    Write a structure to a VASP POSCAR file (VASP 5 format).

    Arguments:
        structure -- a tuple of (lattice_vectors, atomic_symbols,
        atom_positions_frac) data.
        filePath -- file path to save structure to.

    Keyword arguments:
        titleLine -- title line to include on line 1 of the POSCAR file
        (default: _POSCSAR_DefaultTitleLine constant).

    Notes:
        The atom positions should be in fractional coordinates.
    """

    latticevectors, atomicsymbols, atompositions = structure

    # Group atom positions by symbol.
    # To avoid reordering atom types, we keep track of the order in which
    # new symbols were added.

    atomgroupssymbols, atomgroups = [], {}

    for symbol, position in zip(atomicsymbols, atompositions):
        if symbol not in atomgroupssymbols:
            atomgroups[symbol] = []
            atomgroupssymbols.append(symbol)

        atomgroups[symbol].append(position)

    with open(filepath, 'w') as outputWriter:

        outputWriter.write(
            "{0}\n".format(titleLine if titleline is not None else
                           _POSCAR_DefaultTitleLine)
            )

        # Write scale factor.

        outputWriter.write("  {0: >19.16f}\n".format(1.0))

        # Write lattice vectors.

        for i in range(0, 3):
            vx, vy, vz = latticevectors[i]
            outputWriter.write("  {0: >21.16f}  {1: >21.16f}  {2: "
                               ">21.16f}\n".format(vx, vy, vz)
                               )

        # Write atomic symbols and atom counts.

        outputWriter.write(
            "  {0}\n".format(
                "  ".join(["{0: >3}".format(symbol) for symbol in
                           atomgroupssymbols])
                )
            )

        outputWriter.write(
            "  {0}\n".format(
                "  ".join(["{0: >3}".format(len(atomgroups[symbol])) for
                           symbol in atomgroupssymbols])
                )
            )

        # Write positions.

        outputWriter.write("Direct\n")

        for symbol in atomgroupssymbols:
            for f1, f2, f3 in atomgroups[symbol]:
                outputWriter.write("  {0: >21.16f}  {1: >21.16f}  {2: "
                                   ">21.16f}\n".format(f1, f2, f3))


# ------------
# OUTCAR files
# ------------

""" Generic regex for capturing integer values. """

_OUTCAR_GenericIntegerRegex = re.compile(r"(?P<value>\d+)")

""" Regex for extrating atomic masses from the POMASS tags in re-printed 
POTCAR headers at the top of the OUTCAR file. """

_OUTCAR_POMASSRegex = re.compile(r"POMASS =\s*(?P<pomass>\d+\.\d+)")

""" Regex for exctracting the number of atoms in the cell. """

_OUTCAR_NIONSRegex = re.compile(r"NIONS\s*=\s*(?P<nions>\d+)")

""" Regex for extracting phonon frequencies in THz, rad. THz, inv. Cm and 
meV. """

_OUTCAR_ModeFrequencyRegex = re.compile(
    r"^\s*(?P<mode_index>\d+) (?P<freq_sign>f|f/i)\s*=\s*("
    r"?P<freq_thz>\d+\.\d+) THz\s*(?P<freq_radthz>\d+\.\d+) 2PiTHz\s*("
    r"?P<freq_invcm>\d+\.\d+) cm-1\s*(?P<freq_mev>\d+\.\d+) meV\s*$"
    )


def parseoutcar(filepath, extractlist):
    """
    Parse a VASP OUTCAR file and extract the data specified by extractList.

    extractList may contain any of the following items:
        'atomic_masses' -- atomic masses (taken from the POMASS tags in the
        re-printed POTCAR headers at the top of the OUTCAR file).
        'born_charges' -- Born effective-charge tensors.
        'epsilon_static' -- macroscopic dielectric tensor (_not_ the ionic
        contribution).
        'phonon_modes' -- phonon frequencies (in THz, rad. THz, inv. cm and
        meV) and eigenvectors.

    If one or more requested data items cannot be extracted, an error is
    raised.

    Return value:
        A dictionary with the requested data stored under the keys specified
        in extractList.

    Notes:
        If the atomic masses are overridden in the INCAR file, this will not
        be picked up by this function; this is because the fixed column width
        used to write out the POMASS INCAR tag makes it awkward to parse.
    """

    extractlist = [tag.lower() for tag in extractlist]

    for tag in extractlist:
        if tag not in ['atomic_masses', 'born_charges', 'epsilon_static',
                       'phonon_modes']:
            raise Exception("Error: Unrecognised extract tag '{0}'"
                            ".".format(tag)
                            )

    outputdata = {}

    with open(filepath, 'r') as inputReader:
        # Regardless of what we're returning, we parse the first part of the
        # OUTCAR file to obtain:
        # (1) the POMASS values written in the pseudopotential files;
        # (2) the number of atoms; and
        # (3) the number of atoms of each type.

        numatoms = None
        typemasses, typecounts = None, None

        for line in inputReader:
            if "POMASS" in line:
                match = _OUTCAR_POMASSRegex.search(line)

                if match:
                    if typemasses is None:
                        typemasses = []

                    typemasses.append(
                        float(match.group('pomass'))
                        )

            elif "NIONS" in line:
                match = _OUTCAR_NIONSRegex.search(line)
                numatoms = int(match.group('nions'))

            elif "ions per type =" in line:
                matches = _OUTCAR_GenericIntegerRegex.findall(line)

                typecounts = [
                    int(match) for match in matches
                    ]

            if numatoms is not None and typemasses is not None and \
                    typecounts is not None:
                if sum(typecounts) is not numatoms or len(typemasses) is not\
                        len(typecounts):
                    raise Exception("Error: Failed to read atom information "
                                    "from OUTCAR file \"{0}\".".format(
                                                            filepath)
                                    )

                break

        if 'atomic_masses' in extractlist:
            atomicmasses = []

            for mass, count in zip(typemasses, typecounts):
                atomicmasses = atomicmasses + [mass] * count

            outputdata['atomic_masses'] = atomicmasses

        # Loop over the remaining lines and check for other things to extract.

        for line in inputReader:
            line = line.strip()

            if 'born_charges' in extractlist:
                if line == "BORN EFFECTIVE CHARGES (in e, cummulative output)":
                    next(inputReader)

                    bectensors = []

                    for _ in range(0, numatoms):
                        next(inputReader)

                        bectensor = []

                        for i in range(0, 3):
                            rowindex, i1, i2, i3 = next(inputReader).strip(
                            ).split()

                            bectensor.append(
                                [float(i1), float(i2), float(i3)]
                                )

                        bectensors.append(
                            np.array(bectensor, dtype=np.float64)
                            )

                    outputdata['born_charges'] = bectensors

            if 'epsilon_static' in extractlist:
                # Make sure we don't capture the ionic part computed with
                # certain combinations of tags -- this is not what we want
                # for Raman calculations.

                if line[:36] == "MACROSCOPIC STATIC DIELECTRIC TENSOR" and \
                        "IONIC CONTRIBUTION" not in line:
                    next(inputReader)

                    dielectrictensor = []

                    for i in range(0, 3):
                        i1, i2, i3 = next(inputReader).strip().split()

                        dielectrictensor.append(
                            [float(i1), float(i2), float(i3)]
                            )

                    outputdata['epsilon_static'] = np.array(
                        dielectrictensor, dtype=np.float64)

            if 'phonon_modes' in extractlist:
                if line == "Eigenvectors and eigenvalues of the dynamical " \
                           "matrix":
                    if 'phonon_modes' in outputdata:
                        # If NWRITE = 3 was set in the INCAR file, a second
                        # set of frequencies and eigenvectors after division
                        # by sqrt(mass) are written to the OUTCAR file.
                        # Therefore, if we encounter multiple sets of phonon
                        # data, we skip it.

                        continue

                    frequencies = [[], [], [], []]
                    eigenvectors = []

                    for i in range(0, 3):
                        next(inputReader)

                    for i in range(0, 3 * numatoms):
                        # Read frequencies in various units.

                        match = _OUTCAR_ModeFrequencyRegex.match(next(
                            inputReader)
                        )

                        if int(match.group('mode_index')) != i + 1:
                            raise Exception("Error: Unexpected mode ordering "
                                            "in a frequencies/eigenvectors "
                                            "block in OUTCAR file \""
                                            "{0}\".".format(filepath)
                                            )

                        freqsign = -1.0 if match.group('freq_sign') == 'f/i'\
                            else 1.0

                        for e, groupName in enumerate(['freq_thz',
                                                       'freq_radthz',
                                                       'freq_invcm',
                                                       'freq_mev']):
                            frequencies[e].append(
                                math.copysign(float(match.group(groupName)),
                                              freqsign)
                                )

                        next(inputReader)

                        # Read eigenvector.

                        eigenvector = []

                        for j in range(0, numatoms):
                            eigenvector.append(
                                [float(element) for element in next(
                                    inputReader).strip().split()[-3:]]
                                )

                        eigenvectors.append(
                            np.array(eigenvector, dtype=np.float64)
                            )

                        next(inputReader)

                    outputdata['phonon_modes'] = (frequencies, eigenvectors)

    # Check we were able to extract everything that was asked for.

    for tag in extractlist:
        if tag not in outputdata:
            raise Exception("Error: Failed to extract data for tag '{0}' "
                            "from OUTCAR file \"{1}\".".format(tag, filepath)
                            )

    return outputdata
