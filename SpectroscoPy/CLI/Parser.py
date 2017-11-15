# SpectroscoPy/CLI/Parser.py


# ---------
# Docstring
# ---------

""" Implements standardised argparse aruguments for the Runtime module. """


# ---------
# Functions
# ---------

def UpdateParser(parser, supportedFeatures = None):
    """
    Updates parser with groups of arguments read by the Run* methods in the Runtime module.

    The supportedFeatures parameter specifies a list of optional features supported by the CLI, which may include the following:
        'ir_reps' : the CLI can provide mode irreducible representations for building the peak table.
    """

    if supportedFeatures != None:
        supportedFeatures = [
            feature.lower() for feature in supportedFeatures
            ];
    else:
        supportedFeatures = [];

    # Add run-mode argument group.

    group = parser.add_argument_group("Run mode");

    group.add_argument(
        "--ir",
        action = 'store_const', const = 'ir',
        dest = 'RunMode',
        help = "Simulate IR spectrum"
        );

    # Add frequency units group.

    group = parser.add_argument_group("Frequency units");

    group.add_argument(
        "--freq_units",
        metavar = "<units>",
        type = str, dest = "FrequencyUnits",
        help = "Frequency units ('thz', 'inv_cm' or 'mev'; default: 'inv_cm')"
        );

    # If the CLI supports irreducible representations, add a "Peak table" argument group.

    parser.set_defaults(
        IrReps = False
        );

    if 'ir_reps' in supportedFeatures:
        group = parser.add_argument_group("Peak table");

        group.add_argument(
            "--ir_reps",
            action = 'store_true', dest = 'IrReps',
            help = "Include irreducible representation symbols in peak table"
            );

    # Add spectrum-simulation parameters.

    group = parser.add_argument_group("Spectrum simulation");

    group.add_argument(
        "--linewidth",
        metavar = "<linewidth>",
        type = float, dest = "Linewidth",
        default = None,
        help = "Uniform mode linewidth (required if per-mode linewidths are not supplied)"
        );

    group.add_argument(
        "--spectrum_range",
        metavar = "(<min>, <max>)",
        type = str, dest = "SpectrumRange",
        default = None,
        help = "Frequency range of simulated spectrum (default: automatically determined)"
        );

    group.add_argument(
        "--spectrum_resolution",
        metavar = "<resolution>",
        type = float, dest = "SpectrumResolution",
        default = None,
        help = "Frequency resolution of simulated spectrum (default: >1,000 points or minimum spacing of 1)"
        );

    group.add_argument(
        "--instrument_broadening_width",
        metavar = "<broadening>",
        type = float, dest = "InstrumentBroadeningWidth",
        default = None,
        help = "Instrument broadening width (default: None)"
        );

    group.add_argument(
        "--instrument_broadening_shape",
        metavar = "<lineshape>",
        type = str, dest = "InstrumentBroadeningShape",
        default = 'gaussian',
        help = "Instrument broadening shape ('gaussian' or 'lorentzian'; default: 'gaussian')"
        );

    # Add data-output parameters.

    group = parser.add_argument_group("Data output");

    group.add_argument(
        "-o", "--output",
        metavar = "<prefix>",
        type = str, dest = "OutputPrefix",
        default = None,
        help = "Prefix for output files (default: None)"
        );

    group.add_argument(
        "--data_format",
        metavar = "<format>",
        type = str, dest = "DataFormat",
        default = 'dat',
        help = "Format for plain-text output files ('dat' or 'csv'; default: 'dat')"
        );


def PostProcessArgs(args):
    """ Post-process arguments added by UpdateParser() after collecting with ArgumentParser.parse_args(). """

    if args.SpectrumRange != None:
        # Convert to a (min, max) tuple.

        spectrumMin, spectrumMax = args.SpectrumRange.strip().split();
        args.SpectrumRange = (float(spectrumMin), float(spectrumMax));
