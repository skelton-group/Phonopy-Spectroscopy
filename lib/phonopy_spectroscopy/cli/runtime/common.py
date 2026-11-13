# -*- coding: utf-8 -*-


# -------
# Imports
# -------

import numpy as np

from phonopy_spectroscopy.gamma_phonons import get_gamma_phonons

# ---------
# Functions
# ---------


def load_gamma_phonons(loader, proc_args):
    """Load a `GammaPhonons` or `PolarGammaPhonons` for the phonopy-ir
    and phonopy-raman scripts.

    Parameters
    ----------
    interface : LoaderBase
        Loader object for code interface.
    proc_args : Namespace
        Post-processed command-line arguments.

    Returns
    -------
    gamma_ph : GammaPhonons or PolarGammaPhonons
        `GammaPhonons` or `PolarGammaPhonons` object, depending on
        whether a high-frequency dielectric constant and Born charges
        are available.
    """

    conv_trans = None

    if proc_args.primitive_matrix is not None:
        conv_trans = np.linalg.inv(proc_args)

    res = loader.try_load_structure(proc_args.cell_file, conv_trans=conv_trans)

    if res is None:
        raise RuntimeError(
            "Failed to load crystal structure - please specify --cell-file."
        )

    struct, f = res

    print("Structure read from: {0}".format(f))

    res = loader.try_load_gamma_freqs_evecs(proc_args.freqs_evecs_file)

    if res is None:
        raise RuntimeError(
            "Failed to load frequencies/eigenvectors - please specify "
            "--freqs-evecs-file."
        )

    (freqs, evecs), f = res

    print("Frequencies and eigenvectors read from: {0}".format(f))

    lws = loader.try_load_linewidths(proc_args.lws_file, t=proc_args.lws_temp)

    if lws is not None:
        lws, f = lws

        print(
            "Linewidths read from: {0} (T = {1:.0f} K)"
            "".format(f, proc_args.lws_temp)
        )

    irreps = loader.try_load_irreps(proc_args.irreps_file)

    if irreps is not None:
        irreps, f = irreps
        print("Irreps read from: {0}".format(f))

    eps_inf, born_charges = None, None

    res = loader.try_load_eps_inf_and_born(struct, proc_args.eps_born_file)

    if res is not None:
        (eps_inf, born_charges), f = res
        print(r"\eps_inf and Born charges read from: {0}".format(f))

    print("")

    return get_gamma_phonons(
        struct,
        freqs,
        evecs,
        eps_inf=eps_inf,
        born_charges=born_charges,
        lws=lws,
        irreps=irreps,
        t=(proc_args.lws_temp if lws is not None else None),
    )
