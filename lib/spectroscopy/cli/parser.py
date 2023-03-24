# spectroscopy/cli/parser.py


# ---------
# Docstring
# ---------

""" Implements standardised argparse aruguments for the Runtime module. """


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
        metavar="(<min>, <max>)", type=str, dest="SpectrumRange", default=None,
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
