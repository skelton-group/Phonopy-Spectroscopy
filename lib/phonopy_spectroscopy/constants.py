# -*- coding: utf-8 -*-

# ---------
# Docstring
# ---------

"""Defines constants."""

# -------
# Imports
# -------

import math

# ------------------
# Physical constants
# ------------------

BOLTZMANN_CONSTANT_EV = 8.617333262e-5

"""float : Boltzmann constant in eV/K."""

ELEMENTARY_CHARGE = 1.602176634e-19

"""float: Elementary charge in C."""

PLANCK_CONSTANT_J = 6.62607015e-34

"""float : Planck constant in J/s."""

PLANCK_CONSTANT_EV = 4.135667696e-15

"""float : Planck constant in eV/s."""

SPEED_OF_LIGHT = 299792458.0

"""float : Speed of light in m/s."""

VACUUM_PERMITIVITY = 8.8541878188e-12

"""float : Vacuum permittivity in F/m."""

# ---------------------
# Conversions: distance
# ---------------------

BOHR_TO_ANG = 0.529177249

"""float : 1 Bohr radius a_0 in Angstroms."""

ANG_TO_M = 1e-10

"""float : 1 Angstrom in m."""

# -----------------
# Conversions: mass
# -----------------

AMU_TO_KG = 1.660539066e-27

"""float : 1 atomic mass unit (amu) in kg."""

# ---------------------------------
# Conversions: Frequency and energy
# ---------------------------------

EV_TO_J = 1.602176634e-19

"""float : 1 eV in J."""

THZ_TO_INV_CM = 33.35641

"""float : 1 THz in inverse cm."""

THZ_TO_EV = 4.13567 / 1000.0

"""float : 1 THz in eV."""

# ---------------------------
# Conversions: Internal units
# ---------------------------

VASP_TO_THZ = (
    (math.sqrt(EV_TO_J / AMU_TO_KG))
    * (1.0 / ANG_TO_M)
    * (1.0 / (2.0 * math.pi * 1.0e12))
)

"""float : Conversion factor for sqrt(eigenvalues) of a dynamical matrix
in VASP internal units to ordinal frequency in THz."""

# Born charges in e, high-frequency dielectric constant \eps_inf in
# \eps_0, cell volume in Ang^3 and masses in AMU.

VASP_NAC_PREFACTOR = (ELEMENTARY_CHARGE**2 / (VACUUM_PERMITIVITY)) / (
    EV_TO_J * ANG_TO_M
)

r"""float : Conversion factor for the non-analytical correction to the
dynamical matrix (NAC) to VASP internal units."""

# Born chrages in units of E, inverse Hessian in Ang^2 / eV and volume
# in Ang^3. This is equivalent to the EDEPS constant in the VASP source
# code.

VASP_DIELECTRIC_TO_RELATIVE_PERMITTIVITY = (
    ELEMENTARY_CHARGE**2 * (1.0 / EV_TO_J) * (1.0 / ANG_TO_M)
) / VACUUM_PERMITIVITY

"""float : Conversion factor for a dielectric constant in VASP internal
units to relative permittivity."""

# Mode oscillator strengths in e^2 / amu, volume in Ang^3 and
# ordinal frequencies in THz.

VASP_INFRARED_DIELECTIC_TO_RELATIVE_PERMITTIVITY = (
    ELEMENTARY_CHARGE**2
    * (1.0 / AMU_TO_KG)
    * (1.0 / ANG_TO_M**3)
    * (1.0 / ((1e12 * 2.0 * math.pi) ** 2))
) / VACUUM_PERMITIVITY

"""float : Conversion factor for an infrared dielectric function
calculated with oscillator strengths in VASP internal units and
ordinal frequencies in THz to relative permittivity."""

# ----
# Misc
# ----

ZERO_TOLERANCE = 1.0e-8

"""float : Zero tolerance."""
