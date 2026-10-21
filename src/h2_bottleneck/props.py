"""Thermophysical properties of pure gases via CoolProp.

All inputs and outputs are SI: pressure in Pa, temperature in K.
CoolProp uses reference-quality equations of state:
  Hydrogen - Leachman et al. (2009); Methane - Setzmann and Wagner (1991),
the same ones used by Scafidi et al. (2023).
"""
from CoolProp.CoolProp import PropsSI

# Standard reference conditions for "standard cubic metres" (Sm3).
# 15 degC and 1 atm, the convention used in UK gas metering.
P_SC = 101_325.0   # Pa
T_SC = 288.15      # K

GASES = ("Hydrogen", "Methane")


def _check(gas):
    if gas not in GASES:
        raise ValueError(f"gas must be one of {GASES}, got {gas!r}")


def rho(p, T, gas):
    """Density (kg/m3) at pressure p (Pa) and temperature T (K)."""
    _check(gas)
    return PropsSI("D", "P", p, "T", T, gas)


def mu(p, T, gas):
    """Dynamic viscosity (Pa.s) at p (Pa), T (K)."""
    _check(gas)
    return PropsSI("V", "P", p, "T", T, gas)


def Z(p, T, gas):
    """Compressibility factor (-): Z = p / (rho * R_specific * T)."""
    _check(gas)
    return PropsSI("Z", "P", p, "T", T, gas)


def rho_sc(gas):
    """Density at standard conditions (kg/Sm3)."""
    return rho(P_SC, T_SC, gas)


def q_sc_to_mdot(q_sc, gas):
    """Convert standard volumetric flow (Sm3/s) to mass flow (kg/s)."""
    return q_sc * rho_sc(gas)
