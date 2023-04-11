# spectroscopy/cli/parser.py


# ---------
# Docstring
# ---------

""" Implements standardised argparse aruguments for the Runtime module. """


# -------
# Imports
# -------

import math

import numpy as np

from spectroscopy.units import convert_frequency_units
from spectroscopy.utilities import rotation_matrix_xy


# ---------
# Functions
# ---------

def update_parser(parser, spectrum_type):
    """ Updates parser with groups of arguments for the supplied
    spectrum_type, read by the run_* methods in the runtime module."""

    # If spectrum_type is 'raman', add Raman-specific arguments.

    if spectrum_type == 'raman':
        parser.set_defaults(RunMode=None)

        group = parser.add_argument_group("Raman calculations")

        group.add_argument(
            "-d", "--create-disp",
            action='store_const', const="raman_disp", dest="RunMode",
            help="Generate displaced structures for a Raman "
                 "calculation")

        group.add_argument(
            "-r", "--read",
            metavar="<file>", nargs='+', dest="RamanInputFiles",
            help="Collect dielectric constants from a sequence of "
                 "input files")

        group.add_argument(
            "-p", "--post-process",
            action='store_const', const='raman_postproc', dest="RunMode",
            help="Post-process a Raman calculation")

        group.add_argument(
            "--amplitude",
            metavar="<step>", type=float, dest="RamanStep", default=1.0e-2,
            help="Step size for generating displaced structures "
                 "(default: 0.01)")

        group.add_argument(
            "--step-type",
            metavar="< delta_q | max_r | norm_x >",
            type=str, dest="RamanStepType", default='norm_x',
            help="Step size is given in normal-mode amplitude in "
                 "sqrt(amu) Ang ('delta_q'), max. Cartesian "
                 "displacement in Ang ('max_r') or eigendisplacement "
                 "norm ('norm_x'); default: 'norm_x'")

        group.add_argument(
            "--bands",
            metavar="<band_indices>", type=str, dest="RamanBandIndices",
            help="Band indices to generate displacements for (1-3n_a, "
                 "default: all except acoustic modes)")

    # Add spectrum-simulation parameters.

    group = parser.add_argument_group("Spectrum simulation")

    group.add_argument(
        "--spectrum-range",
        metavar="<min> <max>", type=str, dest="SpectrumRange", default=None,
        help="Frequency range of simulated spectrum in --output-units "
            "(default: automatically determined)")

    group.add_argument(
        "--spectrum-resolution",
        metavar="<resolution>",
        type=float, dest="SpectrumResolution", default=None,
        help="Frequency resolution of simulated spectrum in "
            "--output-units (default: >1,000 points or minimum spacing "
            "of 1)")

    group.add_argument(
        "--linewidth",
        metavar="<linewidth>", type=float, dest="Linewidth", default=0.5,
        help="Uniform mode linewidth in --output-units (default: "
             "0.5 THz ~= 16.7 cm^-1 ~= 2.1 meV)")

    group.add_argument(
        "--instrument-broadening-width",
        metavar="<broadening>",
        type=float, dest="InstrumentBroadeningWidth", default=None,
        help="Instrument broadening width in --output-units (default: "
             "no broadening).")

    group.add_argument(
        "--instrument-broadening-shape",
        metavar="<lineshape>",
        type=str, dest="InstrumentBroadeningShape", default='gaussian',
        help="Instrument broadening shape ('gaussian' or 'lorentzian'; "
             "default: 'gaussian')")

    # If spectrum_type is 'raman', add additional parameters for
    # polarised Raman.

    if spectrum_type == 'raman':
        group = parser.add_argument_group("Polarised Raman")

        group.add_argument(
            "--surface",
            metavar="<h> <k> <l>", type=str, dest="SurfaceHKL", default=None,
            help="Miller indices of the crystal surface to rotate to "
                 "the z axis of the instrument frame")
        
        group.add_argument(
            "--theta", metavar="<theta> [<theta_2>, ...]",
            type=str, dest="Theta", default=None,
            help="Optionally, specify one or more rotation angles in "
                 "the xy plane of the instrument frame (i.e. about the "
                 "z axis) **in degrees**")
        
        group.add_argument(
            "--theta-step",
            metavar="<step>", type=float, dest="ThetaStep", default=2.5,
            help="If not specifying a theta with the --theta argument "
                "the spectrum will be simulated between 0 -> 360 deg "
                "and --theta-step specifies the step (default: 2.5 "
                "deg)")
        
        group.add_argument(
            "--incident-pol",
            metavar="<pol>", type = str, dest="IncidentPol", default='x',
            help="Polarisation of incident light (possible values: "
                "'x', 'y', or an (x, y) vector; default: 'x'")

        group.add_argument(
            "--scattered-pol", metavar="<pol>",
            type = str, dest="ScatteredPol", default='parallel',
            help="Polarisation of scattered light (possible values: "
                "'parallel', 'cross', 'x', 'y', or an (x, y) vector; "
                "default: 'parallel')")

    # Add data-output parameters.

    group = parser.add_argument_group("Data output")

    group.add_argument(
        "-o", "--output-prefix",
        metavar="<prefix>", type=str, dest="OutputPrefix", default=None,
        help="Prefix for output files (default: None)")

    group.add_argument(
        "--data-format",
        metavar="dat | csv", type=str, dest="DataFormat", default='dat',
        help="Format for plain-text output files ('dat' or 'csv'; "
             "default: 'dat')")

    group.add_argument(
        "--output-units",
        metavar="< thz | inv_cm | mev | um >",
        type=str, dest="Units",  default='inv_cm',
        help="Output measurement units ('thz', 'inv_cm' or 'mev'; "
             "default: 'inv_cm')")

def _parse_polarisation_string(polarisation_str):
    """ Parse a polarisation specifier into a normalised vector.

    polarisation_str can be one of 'x' or 'y' or two or three values
    specifying a two- or three component vector.
    """

    vals = polarisation_str.strip().split()

    if len(vals) == 1:
        val = vals[0].lower()

        if val == 'x':
            return np.array([1.0, 0.0, 0.0], dtype = np.float64)
        
        if val == 'y':
            return np.array([0.0, 1.0, 0.0], dtype = np.float64)
        
        raise Exception(
            "Error: Invalid polarisation specifier '{0}'.".format(val))
    else:
        v_1, v_2 = None, None

        if len(vals) == 2:
            v_1, v_2 = vals
        else:
            raise Exception(
                "Error: Polarisations must be specified as a "
                "two-component vector.")
        
        v = np.array(
            [float(v_1), float(v_2), 0.0], dtype = np.float64)
        
        return v / np.linalg.norm(v)

def post_process_args(args, spectrum_type):
    """ Post-process arguments added by UpdateParser() after
    collecting with ArgumentParser.parse_args(). """

    if args.SpectrumRange is not None:
        # Convert to a (min, max) tuple.

        spectrum_min, spectrum_max = args.SpectrumRange.strip().split()
        args.SpectrumRange = (float(spectrum_min), float(spectrum_max))
    
    # Convert linewidth to selected output units.

    args.Linewidth = convert_frequency_units(args.Linewidth, 'thz', args.Units)

    if spectrum_type == 'raman':
        if args.RamanInputFiles is not None:
            # The -r/--read flag was set -> set RunMode to 'raman_read'.

            args.RunMode = 'raman_read'

        if args.RamanBandIndices is not None:
            # Convert to a list of zero-based indices.

            args.RamanBandIndices = [
                int(item) - 1 for item in args.RamanBandIndices.strip().split()
                ]

        if args.SurfaceHKL is not None:
            h, k, l = args.SurfaceHKL.strip().split()
            args.SurfaceHKL = (int(h), int(k), int(l))
        
        if args.Theta is not None:
            vals = args.Theta.strip().split()

            if len(vals) > 0:
                args.Theta = [float(val) for val in vals]
            else:
                args.Theta = None
        
        # Parse incident/scattering polarisation vectors.

        i_pol = _parse_polarisation_string(args.IncidentPol)

        args.IncidentPol = i_pol

        s_pol = args.ScatteredPol.strip().lower()

        if s_pol == 'parallel':
            args.ScatteredPol = i_pol
        
        elif s_pol == 'cross':
            # Detection in cross polarisation corresponds to a 90 degree
            # rotation of the incident light polarisation.
            
            args.ScatteredPol = np.dot(rotation_matrix_xy(90.0), i_pol)
        
        else:
            args.ScatteredPol = _parse_polarisation_string(args.IncidentPol)
