# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Implementation of the phonopy-ir script."""

# -------
# Imports
# -------

import os

import numpy as np

from .common import load_gamma_phonons

from ..args import (
    init_parser_for_setup,
    update_parser_for_ir,
    post_process_args,
    post_process_args_for_ir,
)

from ..code_interfaces import get_code_interface
from ..utility.printing import print_3x3_matrix

from ...gamma_phonons import PolarGammaPhonons, FullPolarGammaPhonons
from ...ir.calculation import InfraredCalculation
from ...utility.io_helper import load_json, save_json

# ---------
# Functions
# ---------


def run_mode_setup(proc_args):
    """Implementation of --setup mode for phonopy-ir script.

    Parameters
    ----------
    proc_args : Namespace
        Post-processed command-line arguments.
    """

    print("-----------------")
    print("Calculation setup")
    print("-----------------")
    print("")

    loader = get_code_interface(proc_args.code_interface)

    ph_calc = load_gamma_phonons(loader, proc_args)

    if not isinstance(ph_calc, (PolarGammaPhonons, FullPolarGammaPhonons)):
        raise RuntimeError(
            "Failed to load a high-frequency dielectric constant and "
            "Born charges - please specify --eps-born-file."
        )

    ir_calc = InfraredCalculation(ph_calc)

    save_json(ir_calc.to_dict(), r"ir_calc.json")

    print("Infrared calculation written to: ir_calc.json")
    print("")


def run_mode_props(proc_args):
    """Implementation of --props mode for phonopy-ir script.

    Parameters
    ----------
    proc_args : Namespace
        Post-processed command-line arguments.
    """

    print("---------------------")
    print("Dielectric properties")
    print("---------------------")
    print("")

    if not os.path.isfile(r"ir_calc.json"):
        raise RuntimeError(
            "ir_calc.json not found - please run in --setup mode."
        )

    ir_calc = InfraredCalculation.from_dict(load_json(r"ir_calc.json"))

    print("Infrared calculation read from: ir_calc.json")
    print("")

    ph_calc = ir_calc.phonon_calculation

    for eps, name, sym in [
        (ph_calc.epsilon_inf, "High-frequency", r"\eps_inf"),
        (ph_calc.epsilon_ionic, "Ionic", r"\eps_ionic"),
        (ph_calc.epsilon_static, "Static", r"\eps_static"),
    ]:
        print("{0} dielectric constant {1}:".format(name, sym))
        print("")

        print_3x3_matrix(eps)
        print("")

    pop_freq = ir_calc.phonon_calculation.pop_frequency()

    print("POP frequency w_po = {0:.3f} THz".format(pop_freq))

    print("")


def main(args=None):
    """Entry-point for phonopy-ir script.

    Parameters
    ----------
    args : list of str or None, optional
        Command-line arguments for testing (default: None).
    """

    # Collect, parse and post process and parse command-line arguments.

    parser, setup_group = init_parser_for_setup()
    update_parser_for_ir(parser, setup_group)

    args = parser.parse_args(args=args)

    post_process_args(args)
    post_process_args_for_ir(args)

    if args.run_mode_setup:
        run_mode_setup(args)

    if args.run_mode_props:
        run_mode_props(args)
