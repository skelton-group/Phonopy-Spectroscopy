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
            help="Generate displaced structures for a Raman calculation")

        group.add_argument(
            "-r", "--read",
            metavar="<file>", nargs='+', dest="RamanInputFiles",
            help="Collect dielectric constants from a sequence of input files")

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
            help="Step size is given in normal-mode amplitude (dQ; "
                 "'delta_q'), max. Cartesian displacement ('max_r') or "
                 "eigendisplacement norm ('norm_x'); default: 'norm_x'")

        group.add_argument(
            "--bands",
            metavar="<band_indices>", type=str, dest="RamanBandIndices",
            help="Band indices to generate displacements for (1-N, "
                 "default: all except acoustic modes)")

    # Add units group.

    group = parser.add_argument_group("Units")

    group.add_argument(
        "--units",
        metavar="< thz | inv_cm | mev | um >",
        type=str, dest="Units",  default='inv_cm',
        help="Units for simulated spectrum ('thz', 'inv_cm', 'mev' or "
             "'um'; default: 'inv_cm')")

    # Add spectrum-simulation parameters.

    group = parser.add_argument_group("Spectrum simulation")

    group.add_argument(
        "--linewidth",
        metavar="<linewidth>", type=float, dest="Linewidth", default=16.5,
        help="Uniform mode linewidth ***in inverse cm*** (default: "
             "16.5 cm^-1 ~= 0.5 THz)")

    group.add_argument(
        "--spectrum-range",
        metavar="<min> <max>", type=str, dest="SpectrumRange", default=None,
        help="Frequency range of simulated spectrum (default: "
             "automatically determined)")

    group.add_argument(
        "--spectrum-resolution",
        metavar="<resolution>",
        type=float, dest="SpectrumResolution", default=None,
        help="Frequency resolution of simulated spectrum (default: "
             ">1,000 points or minimum spacing of 1)")

    group.add_argument(
        "--instrument-broadening-width",
        metavar="<broadening>",
        type=float, dest="InstrumentBroadeningWidth", default=None,
        help="Instrument broadening width (default: None)")

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
            help="Miller indices of the crystal surface to rotate to the z "
                "axis of the instrument frame.")
        
        group.add_argument(
            "--theta", metavar="<theta> [<theta_2>]",
            type=str, dest="Theta", default=None,
            help="Optionally, specify a rotation angle or (min, max) pair "
                "of angles in the xy plane of the instrument frame (i.e. "
                "about the z axis (default: 0 -> 360 deg)")
        
        group.add_argument(
            "--theta-step",
            metavar="<step>", type=float, dest="ThetaStep", default=None,
            help="If specifying a range of theta using the --theta "
                "argument, specify a step within this range (default: "
                "1 deg or 101 points, whichever is smaller)")
        
        group.add_argument(
            "--incident-pol",
            metavar="<pol>", type = str, dest="IncidentPol", default='x',
            help="Polarisation of incident light (possible values: 'x', "
                "'y', an (x, y) vector, or a (x, y, z) vector; default: "
                "'x'")

        group.add_argument(
            "--scattered-pol", metavar="<pol>",
            type = str, dest="ScatteredPol", default='none',
            help="Polarisation of scattered light (possible values: "
                "'none', 'parallel', 'cross', 'x', 'y', an (x, y) vector, "
                "or a (x, y, z) vector; default: 'none')")

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

def _parse_polarisation_string(polarisation_str):
    """ Parse a polarisation specifier into a normalised vector.

    polarisation_str can be one of 'x' or 'y' or two or three values
    specifying a two- or three-componetn vector.
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
        v_1, v_2, v_3 = None, None, None

        if len(vals) == 2:
            v_1, v_2 = vals
            v_3 = 0.0
        elif len(vals) == 3:
            v_1, v_2, v_3 = vals
        else:
            raise Exception(
                "Error: Polarisations must be specified as two- or "
                "three-component vectors.")
        
        v = np.array(
            [float(v_1), float(v_2), float(v_3)], dtype = np.float64)
        
        return v / np.linalg.norm(v)

def post_process_args(args, spectrum_type):
    """ Post-process arguments added by UpdateParser() after
    collecting with ArgumentParser.parse_args(). """

    if spectrum_type == 'raman':
        if args.RamanInputFiles is not None:
            # The -r/--read flag was set -> set RunMode to 'raman_collect'.

            args.RunMode = 'raman_read'

        if args.RamanBandIndices is not None:
            # Convert to a list of zero-based indices.

            args.RamanBandIndices = [
                int(item) - 1 for item in args.RamanBandIndices.strip().split()
                ]

    if args.SpectrumRange is not None:
        # Convert to a (min, max) tuple.

        spectrum_min, spectrum_max = args.SpectrumRange.strip().split()
        args.SpectrumRange = (float(spectrum_min), float(spectrum_max))

    if args.SurfaceHKL is not None:
        h, k, l = args.SurfaceHKL.strip().split()
        args.SurfaceHKL = (int(h), int(k), int(l))
    
    if args.Theta is not None:
        vals = args.Theta.strip().split()

        if len(vals) == 0 or len(vals) > 2:
            raise Exception(
                "Error: If supplied, --theta must specify one or two "
                "angles")
        
        args.Theta = tuple(float(val) for val in vals)
    else:
        args.Theta = (0.0, 360.0)
    
    if len(args.Theta) == 2:
        theta_1, theta_2 = args.Theta
        
        if args.ThetaStep is not None:
            if int(math.floor((theta_2 - theta_1) / args.ThetaStep)) < 1:
                raise Exception(
                    "If supplied, --theta-step must produce more than one "
                    "value between the range specified by --theta.")
        else:
            args.ThetaStep = min(1.0, (theta_2 - theta_1) / 100.0)

    # Parse incident/scattering polarisation vectors.

    i_pol = _parse_polarisation_string(args.IncidentPol)

    args.IncidentPol = i_pol

    s_pol = args.ScatteredPol.strip().lower()

    if s_pol == 'none':
        # Normalised vector v_x = v_y to treat all in-plane elements
        # of the Raman tensor with equal weighting.

        args.ScatteredPol = (
            np.array([1.0, 1.0, 0.0], dtype = np.float64) / math.sqrt(2.0))
    
    elif s_pol == 'parallel':
        args.ScatteredPol = i_pol
    
    elif s_pol == 'cross':
        # Detection in cross polarisation corresponds to a 90 degree
        # rotation of the incident light polarisation.
        
        args.ScatteredPol = np.dot(rotation_matrix_xy(90.0), i_pol)
    
    else:
        args.ScatteredPol = _parse_polarisation_string(args.IncidentPol)
