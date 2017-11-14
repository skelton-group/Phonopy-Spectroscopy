# SpectroscoPy/Interfaces/VASP.py


# ---------
# Docstring
# ---------

""" Contains routines for interfacing with the Vienna Ab-initio Simulation Package (VASP) code. """


# -------
# Imports
# -------

import math;
import re;

import numpy as np;


# ---------
# Constants
# ---------

""" Generic regex for identifying integer values. """

_GenericIntegerRegex = re.compile(r"(?P<value>\d+)");

""" Regex for extrating atomic masses from the POMASS tags in re-printed POTCAR headers at the top of the OUTCAR file. """

_OUTCAR_POMASSRegex = re.compile(r"POMASS =\s*(?P<pomass>\d+\.\d+);");

""" Regex for exctracting the number of atoms in the cell. """

_OUTCAR_NIONSRegex = re.compile("NIONS\s*=\s*(?P<nions>\d+)");

""" Regex for extracting phonon frequencies in THz, rad. THz, inv. Cm and meV. """

_OUTCAR_ModeFrequencyRegex = re.compile(
    r"^\s*(?P<mode_index>\d+) (?P<freq_sign>f|f/i)\s*=\s*(?P<freq_thz>\d+\.\d+) THz\s*(?P<freq_radthz>\d+\.\d+) 2PiTHz\s*(?P<freq_invcm>\d+\.\d+) cm-1\s*(?P<freq_mev>\d+\.\d+) meV\s*$"
    );


# ---------
# Functions
# ---------

def ParseOUTCAR(extractList, filePath = r"OUTCAR"):
    """
    Parse a VASP OUTCAR file and extract the data specified by extractList.

    extractList may contain any of the following items:
        'atomic_masses' -- atomic masses (taken from the POMASS tags in the re-printed POTCAR headers at the top of the OUTCAR file).
        'phonon_modes' -- phonon frequencies (in THz, rad. THz, inv. cm and meV) and eigenvectors.
        'born_charges' -- Born effective-charge tensors.

    If one or more requested data items cannot be extracted, an error is raised.

    Return value:
        A dictionary with the requested data stored under the keys specified in extractList.

    Notes:
        If the atomic masses are overridden in the INCAR file, this will not be picked up by this function; this is because the fixed column width used to write out the POMASS INCAR tag makes it awkward to parse.
    """

    extractList = [tag.lower() for tag in extractList];

    for tag in extractList:
        if tag not in ['atomic_masses', 'phonon_modes', 'born_charges']:
            raise Exception("Error: Unrecognised extract tag '{0}'.".format(tag));

    outputData = { };

    with open(filePath, 'r') as inputReader:
        # Regardless of what we're returning, we parse the first part of the OUTCAR file to obtain:
        # (1) the POMASS values written in the pseudopotential files;
        # (2) the number of atoms; and
        # (3) the number of atoms of each type.

        numAtoms = None;
        typeMasses, typeCounts = None, None;

        for line in inputReader:
            if "POMASS" in line:
                match = _OUTCAR_POMASSRegex.search(line);

                if match:
                    if typeMasses == None:
                        typeMasses = [];

                    typeMasses.append(
                        float(match.group('pomass'))
                        );

            elif "NIONS" in line:
                match = _OUTCAR_NIONSRegex.search(line);
                numAtoms = int(match.group('nions'));

            elif "ions per type =" in line:
                matches = _GenericIntegerRegex.findall(line);

                typeCounts = [
                    int(match) for match in matches
                    ];

            if numAtoms != None and typeMasses != None and typeCounts != None:
                if sum(typeCounts) != numAtoms or len(typeMasses) != len(typeCounts):
                    raise Exception("Error: Failed to read atom information from OUTCAR file \"{0}\".".format(filePath));

                break;

        if 'atomic_masses' in extractList:
            atomicMasses = [];

            for mass, count in zip(typeMasses, typeCounts):
                atomicMasses = atomicMasses + [mass] * count;

            outputData['atomic_masses'] = atomicMasses;

        if 'born_charges' in extractList:
            # Search for the header of the section where the Born charges are output.

            for line in inputReader:
                if line.strip() == "BORN EFFECTIVE CHARGES (in e, cummulative output)":
                    break;

            next(inputReader);

            becTensors = [];

            for _ in range(0, numAtoms):
                next(inputReader);

                becTensor = [];

                for i in range(0, 3):
                    rowIndex, i1, i2, i3 = next(inputReader).strip().split();

                    becTensor.append(
                        [float(i1), float(i2), float(i3)]
                        );

                becTensors.append(
                    np.array(becTensor, dtype = np.float64)
                    );

            outputData['born_charges'] = becTensors;

        if 'phonon_modes' in extractList:
            # Search for the header of the section where the (raw) eigenvectors are output.

            for line in inputReader:
                if line.strip() == "Eigenvectors and eigenvalues of the dynamical matrix":
                    break;

            frequencies = [[], [], [], []];
            eigenvectors = [];

            for i in range(0, 3):
                next(inputReader);

            for i in range(0, 3 * numAtoms):
                # Read frequencies in various units.

                match = _OUTCAR_ModeFrequencyRegex.match(next(inputReader));

                if int(match.group('mode_index')) != i + 1:
                    raise Exception("Error: Unexpected mode ordering in a frequencies/eigenvectors block in OUTCAR file \"{0}\".".format(filePath));

                freqSign = -1.0 if match.group('freq_sign') == 'f/i' else 1.0;

                for i, groupName in enumerate(['freq_thz', 'freq_radthz', 'freq_invcm', 'freq_mev']):
                    frequencies[i].append(
                        math.copysign(float(match.group(groupName)), freqSign)
                        );

                next(inputReader);

                # Read eigenvector.

                eigenvector = [];

                for j in range(0, numAtoms):
                    eigenvector.append(
                        [float(element) for element in next(inputReader).strip().split()[-3:]]
                        );

                eigenvectors.append(
                    np.array(eigenvector, dtype = np.float64)
                    );

                next(inputReader);

            outputData['phonon_modes'] = (frequencies, eigenvectors);
    
    # Check we were able to extract everything that was asked for.

    for tag in extractList:
        if tag not in outputData:
            raise Exception("Error: Failed to extract data for tag '{0}' from OUTCAR file \"{1}\".".format(tag, filePath));

    return outputData;
