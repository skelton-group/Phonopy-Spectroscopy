# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines implementing command-line argument handling."""

# -------
# Imports
# -------

from argparse import ArgumentParser

from .code_interfaces import get_supported_code_interfaces
from .utility.parsing import parse_3x3_matrix

# ------------
# Parser setup
# ------------


def init_parser_for_setup():
    """Initialise and return an `ArgumentParser` object to parse
    command-line arguments, initialised with common setup options.

    Returns
    -------
    parser : tuple
        Tuple of `(parser, setup_group)` with the initialised
        `ArgumentParser` and argument group with setup options.
    """

    parser = ArgumentParser()

    parser.add_argument(
        "--setup",
        dest="run_mode_setup",
        action="store_true",
        default=False,
        help="perform calculation setup",
    )

    setup_group = parser.add_argument_group("setup")

    setup_group.add_argument(
        "--cell-file",
        dest="cell_file",
        type=str,
        default=None,
        help="crystal structure (cell) file",
    )

    setup_group.add_argument(
        "--freqs-evecs-file",
        dest="freqs_evecs_file",
        type=str,
        default=None,
        help="frequencies and eigenvectors file",
    )

    setup_group.add_argument(
        "--lws-file",
        dest="lws_file",
        type=str,
        default=None,
        help="linewidths file (optional)",
    )

    setup_group.add_argument(
        "--lws-temp",
        dest="lws_temp",
        type=float,
        default=300.0,
        help="linewidths temperature (optional, default: 300 K)",
    )

    setup_group.add_argument(
        "--irreps-file",
        dest="irreps_file",
        type=str,
        default=None,
        help="irreps file (optional)",
    )

    setup_group.add_argument(
        "--eps-born-file",
        dest="eps_born_file",
        type=str,
        default=None,
        help="high-frequency dielectric constant/Born charges file",
    )

    setup_group.add_argument(
        "--pa",
        "--primitive-axis",
        dest="primitive_matrix",
        type=str,
        default=None,
        help="transformation matrix from conventional to primitive cell",
    )

    setup_group.add_argument(
        "--interface",
        dest="code_interface",
        type=str,
        choices=get_supported_code_interfaces(),
        default="vasp",
        help="code interface to use",
    )

    return (parser, setup_group)


def update_parser_for_ir(parser, setup_group):
    """Update an initialised `ArgumentParser` object with additional
    options for IR calculations.

    Parameters
    ----------
    parser : ArgumentParser
        Argument parser initialised with setup argument group.
    setup_group : _ArgumentGroup
        Argument group containing setup arguments.
    """

    parser.add_argument(
        "--props",
        dest="run_mode_props",
        action="store_true",
        default=False,
        help="Calculate dielectric properties",
    )


def update_parser_for_raman(parser, setup_group):
    """Update an initialised `ArgumentParser` object with additional
    options for Raman calculations.

    Parameters
    ----------
    parser : ArgumentParser
        Argument parser initialised with setup argument group.
    setup_group : _ArgumentGroup
        Argument group containing setup arguments.
    """

    parser.add_argument(
        "-r",
        "--read",
        dest="run_mode_read",
        action="store_true",
        default=False,
        help="Run in read mode",
    )

    parser.add_argument(
        "input-files",
        dest="input_files",
        nargs="*",
        help="Input files for creating Raman calculations.",
    )

    setup_group.add_argument(
        "--calc-type",
        dest="raman_calc_type",
        type=str,
        choices=["fd"],
        default="fd",
        help="Type of Raman calculation",
    )

    setup_group.add_argument(
        "--fd-step",
        dest="raman_fd_step",
        type=float,
        default="1.0e-2",
        help="Step size for finite-difference calculations",
    )

    setup_group.add_argument(
        "--fd-prec",
        dest="raman_fd_prec",
        type=int,
        default=2,
        help="Precision of finite-difference stencil in FD calculations",
    )

    parser.add_argument_group(
        "--sort-files",
        dsest="raman_sort_files",
        action="store_true",
        default=False,
        help="Sort input files by name.",
    )


# ---------------
# Post processing
# ---------------


def post_process_args(args):
    """Post-process parsed command-line arguments.

    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments.
    """

    if args.primitive_matrix is not None:
        args.primitive_matrix = parse_3x3_matrix(args.primitive_matrix)


def post_process_args_for_ir(args):
    """Post-process parsed command-line arguments for IR workflow.

    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments.
    """

    pass


def post_process_args_for_raman(args):
    """Post-process parsed command-line arguments for Raman workflow.

    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments.
    """

    pass
