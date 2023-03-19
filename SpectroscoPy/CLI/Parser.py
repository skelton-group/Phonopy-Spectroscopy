# SpectroscoPy/CLI/Parser.py


# ---------
# Docstring
# ---------

""" Implements standardised argparse aruguments for the Runtime module. """


# ---------
# Functions
# ---------

def update_parser(parser, spectrumtype, supportedfeatures=None):
    """
    Updates parser with groups of arguments for the supplied spectrumType,
    read by the Run* methods in the Runtime module.

    The supportedFeatures parameter specifies a list of optional features
    supported by the CLI, which may include the following:
        'ir_reps' : the CLI can provide mode irreducible representations
        for building the peak table.
    """

    if supportedfeatures is not None:
        supportedfeatures = [
            feature.lower() for feature in supportedfeatures
            ]
    else:
        supportedfeatures = []

    # If runMode is 'raman', add Raman-specific arguments.

    if spectrumtype == 'raman':
        parser.set_defaults(
            RunMode=None
            )

        group = parser.add_argument_group("Raman calculations")

        group.add_argument(
            "-d", "--create_disp",
            action='store_const', const="raman_disp", dest="RunMode",
            help="Generate displaced structures for a Raman calculation"
            )

        group.add_argument(
            "-r", "--read",
            metavar="<file>",
            nargs='+', dest="RamanInputFiles",
            help="Collect dielectric constants from a sequence of input files"
            )

        group.add_argument(
            "-p", "--post_process",
            action='store_const', const='raman_postproc', dest="RunMode",
            help="Post-process a Raman calculation"
            )

        group.add_argument(
            "--amplitude",
            metavar="<step>",
            type=float, dest="RamanStep",
            default=1.0e-2,
            help="Step size for generating displaced structures (default: "
                 "0.01)"
            )

        group.add_argument(
            "--step_type",
            metavar="< delta_q | max_r | norm_x >",
            type=str, dest="RamanStepType",
            default='norm_x',
            help="Step size is given in normal-mode amplitude (dQ; "
                 "'delta_q') max. Cartesian displacement ('max_r') or "
                 "eigendisplacement norm ('norm_x'); default: 'norm_x'"
            )

        group.add_argument(
            "--bands",
            metavar="<band_indices>",
            type=str, dest="RamanBandIndices",
            help="Band indices to generate displacements for (1-N, default: "
                 "all except acoustic modes)"
            )

    # Add frequency units group.

    group = parser.add_argument_group("Frequency units")

    group.add_argument(
        "--freq_units",
        metavar="< thz | inv_cm | mev | um >",
        type=str, dest="FrequencyUnits",
        default='inv_cm',
        help="Frequency units for simulated spectrum ('thz', 'inv_cm', "
             "'mev' or 'um'; default: 'inv_cm')"
        )

    # If the CLI supports irreducible representations, add a "Peak table"
    # argument group.

    parser.set_defaults(
        IrReps=False
        )

    if 'ir_reps' in supportedfeatures:
        group = parser.add_argument_group("Peak table")

        group.add_argument(
            "--ir_reps",
            action='store_true', dest='IrReps',
            help="Include irreducible representation symbols in peak table"
            )

    # Add spectrum-simulation parameters.

    group = parser.add_argument_group("Spectrum simulation")

    group.add_argument(
        "--linewidth",
        metavar="<linewidth>",
        type=float, dest="Linewidth",
        default=16.5,
        help="Uniform mode linewidth ***in inverse cm*** (default: 16.5 "
             "cm^-1 ~= 0.5 THz)"
        )

    group.add_argument(
        "--spectrum_range",
        metavar="(<min>, <max>)",
        type=str, dest="SpectrumRange",
        default=None,
        help="Frequency range of simulated spectrum (default: "
             "automatically determined)"
        )

    group.add_argument(
        "--spectrum_resolution",
        metavar="<resolution>",
        type=float, dest="SpectrumResolution",
        default=None,
        help="Frequency resolution of simulated spectrum (default: >1,"
             "000 points or minimum spacing of 1)"
        )

    group.add_argument(
        "--instrument_broadening",
        metavar="<broadening>",
        type=float, dest="InstrumentBroadening",
        default=None,
        help="Instrument broadening width (default: None)"
        )

    group.add_argument(
        "--instrument_broadening_shape",
        metavar="<lineshape>",
        type=str, dest="InstrumentBroadeningShape",
        default='gaussian',
        help="Instrument broadening shape ('gaussian' or 'lorentzian'; "
             "default: 'gaussian')"
        )

    # Add data-output parameters.

    group = parser.add_argument_group("Data output")

    group.add_argument(
        "-o", "--output",
        metavar="<prefix>",
        type=str, dest="OutputPrefix",
        default=None,
        help="Prefix for output files (default: None)"
        )

    group.add_argument(
        "--data_format",
        metavar="dat | csv",
        type=str, dest="DataFormat",
        default='dat',
        help="Format for plain-text output files ('dat' or 'csv'; default: "
             "'dat')"
        )


def post_process_args(args, spectrumtype):
    """ Post-process arguments added by UpdateParser() after
    collecting with ArgumentParser.parse_args(). """

    if spectrumtype == 'raman':
        if args.RamanInputFiles is not None:
            # The -r/--read flag was set -> set RunMode to 'raman_collect'.

            args.RunMode = 'raman_read'

        if args.ramanBandIndices is not None:
            # Convert to a list of zero-based indices.

            args.RamanBandIndices = [
                int(item) - 1 for item in args.RamanBandIndices.strip().split()
                ]

    if args.SpectrumRange is not None:
        # Convert to a (min, max) tuple.

        spectrummin, spectrummax = args.SpectrumRange.strip().split()
        args.SpectrumRange = (float(spectrummin), float(spectrummax))
