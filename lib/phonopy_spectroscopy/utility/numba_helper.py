# -*- coding: utf-8 -*-


# ---------
# Docstring
# ---------


"""Helper routines for working with Numba."""


# ------------------
# "Dummy" attributes
# ------------------


def dummy_njit(*args, **kwargs):
    """Dummy attribute to substitute for @njit when Numba is not
    available."""

    if len(args) == 1 and callable(args[0]):
        return args[0]
    else:
        return lambda f: f
