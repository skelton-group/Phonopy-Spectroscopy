# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Helper routines for handling plot styles."""

# -------
# Imports
# -------

import matplotlib.pyplot as plt

# ---------
# Functions
# ---------


def set_plot_style(style="phonopy_spectroscopy.cli.plotting.default"):
    plt.style.use(style)
