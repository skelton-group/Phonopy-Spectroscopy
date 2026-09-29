# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Routines implementing command-line argument handling."""

# -------
# Imports
# -------

from argparse import ArgumentParser

# ------------
# Parser setup
# ------------


def parser_init():
    """_summary_"""

    parser = ArgumentParser()

    parser.add_argument(
        "--cell",
        dest="cell_file",
        type=str,
        default="POSCAR",
        help="Crystal structure (POSCAR or phonopy.yaml)",
    )

    parser.add_argument(
        "--freqs-evecs",
        dest="freqs_evecs_file",
        type=str,
        default=None,
        help=(
            "Frequencies and eigenvectors (mesh.yaml, mesh.hdf5, "
            "band.yaml, or band.hdf5)"
        ),
    )

    parser.add_argument(
        "--lw",
        "--linewidth",
        dest="linewidth",
        type=float,
        default=None,
        help=(
            "Uniform linewidth or scale factor for calculated "
            "linewidths (default: 0.5 THz or 1.0)"
        ),
    )

    parser.add_argument(
        "--lws-file",
        dest="linewidths_file",
        type=str,
        default=None,
        help="Linewidths (kappa-m*.hdf5 or kappa-m*-g*.hdf5)",
    )

    parser.add_argument(
        "--lws-temp",
        dest="linewidths_temp",
        type=float,
        default=300.0,
    )

    parser.add_argument(
        "--irreps",
        dest="irreps_file",
        type=str,
        default="irreps.yaml",
        help="Irreps (irreps.yaml)",
    )


def parser_update_ir(parser):
    """_summary_

    Parameters
    ----------
    parser : _type_
        _description_
    """

    parser.add_argument(
        "--born",
        dest="born_file",
        type=str,
        default="BORN",
        help="Born charges and high-frequency dielectric constant (BORN)",
    )

    parser.add_argument(
        "--eps_hf",
        dest="epsilon_inf",
        type=str,
        default=None,
        help="High-frequency dielectric constant",
    )


def parser_update_raman(parser):
    """_summary_

    Parameters
    ----------
    parser : _type_
        _description_
    """

    pass


# ------------------
# Default parameters
# ------------------

# ---------------
# Post processing
# ---------------


def args_post_proc(args):
    """_summary_

    Parameters
    ----------
    args : _type_
        _description_
    """

    pass


def args_post_proc_ir(args):
    """_summary_

    Parameters
    ----------
    args : _type_
        _description_
    """

    pass


def args_post_proc_raman(args):
    """_summary_

    Parameters
    ----------
    args : _type_
        _description_
    """

    pass
