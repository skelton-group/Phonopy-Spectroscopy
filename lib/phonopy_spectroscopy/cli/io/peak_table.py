# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Helper routines for outputting peak tables and spectra to plain-text
files."""

# -------
# Imports
# -------

import numpy as np

# -------
# General
# -------


def df_to_txt(df, file_path, col_fmts, preamble_lines=None):
    """Write a Pandas `DataFrame` to a text file.

    Parameters
    ----------
    df : pandas.DataFrame
        `DataFrame` to write.
    preamble_lines : array_like
        Line(s) to write to the start of the output file (will be
        prepended with "# ").
    col_fmts : array_like
        Format specifiers for columns.
    file_path : str
        File path.
    """

    col_fmts = np.atleast_1d(np.asanyarray(col_fmts, dtype=str))

    if col_fmts.ndim != 1:
        raise ValueError("col_fmts must be 1D.")

    if preamble_lines is not None:
        preamble_lines = np.atleast_1d(np.asanyarray(preamble_lines))

        if preamble_lines.ndim != 1:
            raise ValueError("preamble_lines must be 1D.")

    if len(col_fmts) != len(df.columns):
        raise ValueError(
            "col_fmts must contain a format specifier for each column in df."
        )

    # Determine column widths as the maximum of the lengths of the
    # DataFrame headers and the lengths of the values in the first row
    # after formatting according to col_fmts.

    row_test = [fmt.format(v) for fmt, v in zip(col_fmts, df.iloc[0].values)]
    col_widths = [max(len(s), len(h)) for s, h in zip(row_test, df.columns)]

    with open(file_path, "w") as f:
        if preamble_lines is not None:
            for l in preamble_lines:
                f.write(f"# {l}\n")

        headers = [
            "{{0: >{0}}}".format(w).format(h)
            for w, h in zip(col_widths, df.columns)
        ]

        f.write("# " + "  ".join(headers) + "\n")

        for row in df.itertuples(index=False, name=None):
            row = [
                "{{0: >{0}}}".format(w).format(fmt.format(v))
                for w, fmt, v in zip(col_widths, col_fmts, row)
            ]

            f.write("  " + "  ".join(row) + "\n")


# ----------------
# Infrared spectra
# ----------------


def infrared_dielectric_function_peak_table_to_txt(eps_ir, file_path):
    """Write a peak table from an `InfraredDielectricFunction` object to
    a plain-text file.

    Parameters
    ----------
    eps_ir : InfraredDielectricFunction
        Dielectric function.
    file_path : str
        File path.
    """

    df = eps_ir.peak_table()

    preamble_lines = [
        f"energy : {eps_ir.energy_unit_text_label}",
        f"osc_str : {eps_ir.mode_oscillator_strength_unit_text_label}",
    ]

    col_fmts = ["{0: >12.5f}"]
    col_fmts += ["{0: >12.5e}"] * 9
    col_fmts += ["{0: >10.5f}", "{0: >5}"]

    df_to_txt(df, file_path, col_fmts, preamble_lines=preamble_lines)


def dielectric_function_to_txt(eps, file_path):
    """Write the dielectric function from a `DielectricFunction` object
    to a plain-text file.

    Parameters
    ----------
    eps : DielectricFunction
        Dielectric function.
    file_path : str
        File path.
    """

    df = eps.dielectric_function()

    preamble_lines = [
        f"energy : {eps.energy_unit_text_label}",
        f"epsilon : {eps.epsilon_unit_text_label}",
    ]

    col_fmts = ["{0: >12.5f}"] + ["{0: >12.5e}"] * (len(df.columns) - 1)

    df_to_txt(df, file_path, col_fmts, preamble_lines=preamble_lines)


def optical_eigenmodes_to_txt(oes, file_path):
    """Write the quantities from an `OpticalEigenmodes` object to a
    plain-text file.

    Parameters
    ----------
    oes : OpticalEigenmodes
        Optical eigenmodes.
    file_path : str
        File path.
    """

    df = oes.eigenmode_optical_properties()

    preamble_lines = [
        "energy : {0}".format(oes.energy_unit_text_label),
        "epsilon : {0}".format(oes.eigenvalue_unit_text_label),
        "absorption coefficient : {0}".format(
            oes.absorption_coefficient_unit_text_label
        ),
        "optical conductivity : {0}".format(
            oes.optical_conductivity_unit_text_label
        ),
    ]

    col_fmts = ["{0: >12.5f}"] + ["{0: >12.5e}"] * (len(df.columns) - 1)

    df_to_txt(df, file_path, col_fmts, preamble_lines=preamble_lines)


def optical_spectrum_to_txt(sp, file_path):
    """Write the optical spectra from an object derived from the
    `OpticalSpectrumBase` class to a plain-text file.

    Parameters
    ----------
    sp : OpticalSpectrumBase
        Optical spectrum.
    file_path : str
        File path.
    """

    df = sp.spectrum()

    preamble_lines = ["energy : {0}".format(sp.energy_unit_text_label)]
    col_fmts = ["{0: >12.5f}"] + ["{0: >12.5e}"] * (len(df.columns) - 1)

    df_to_txt(df, file_path, col_fmts, preamble_lines=preamble_lines)


# -------------
# Raman spectra
# -------------


def raman_peak_table_to_txt(sp, file_path):
    """Write the peak table from an `NDRamanSpectrum object or derived
    class to a plain-text file.

    Parameters
    ----------
    sp : NDRamanSpectrum
        Raman spectrum.
    file_path : str
        File path.
    """

    df = sp.peak_table(or_stack=True)

    num_sp = np.prod(
        sp.band_intensities.shape[: sp.num_spectral_dims], dtype=int
    )

    col_fmts = [
        "{0: >12.5f}",
        "{0: >12.5e}",
        "{0: >12.5e}",
        "{0: >10.5f}",
        "{0: >5}",
    ] * num_sp

    preamble_lines = [
        f"raman_shift : {sp.raman_shift_unit_text_label}",
        f"intensity : {sp.intensity_unit_text_label}",
        f"cross_section : {sp.cross_section_unit_text_label}",
    ]

    df_to_txt(df, file_path, col_fmts, preamble_lines=preamble_lines)


def raman_spectrum_to_txt(sp, file_path):
    """Write the spectra from an `NDRamanSpectrum` object or derived
    class to a plain-text file.

    Parameters
    ----------
    sp : NDRamanSpectrum
        Raman spectrum.
    file_path : str
        File path.
    """

    df = sp.spectrum()

    num_sp = np.prod(sp.cross_section.shape[: sp.num_spectral_dims], dtype=int)
    col_fmts = ["{0: >12.5f}"] + ["{0: >10.5e}"] * num_sp

    preamble_lines = [
        f"raman_shift : {sp.raman_shift_unit_text_label}",
        f"cross_section : {sp.cross_section_unit_text_label}",
    ]

    df_to_txt(df, file_path, col_fmts, preamble_lines=preamble_lines)
