"""Reservoir (node 5): closed-tank material balance.

The reservoir is one closed tank of fixed gas pore volume V at constant
temperature T (decisions A2 and A8). The state is the mass of gas in the tank.
Pressure follows from the reference equation of state in CoolProp:

    inventory -> density (m / V) -> pressure p(rho, T)

This is the p/Z material balance (Lee and Wattenbarger, 1996) solved exactly:
p/Z = rho * R_specific * T, so p/Z is linear in inventory for a closed tank,
and no iteration on Z is needed.

The pore volume is set once, from Rough's methane GIIP at the initial pressure
(Scafidi et al., 2023, Table 1), and is then reused for hydrogen. Rock and
water compressibility are neglected, and the tank is not compartmentalised;
Heinemann et al. (2025) show pressure support from the Carboniferous
underburden, so a single tank is a stated limitation.

All inputs and outputs are SI: Pa, K, kg, m3, Sm3 (15 degC, 1 atm), s, J.
Flow sign convention: q_sc > 0 is injection, q_sc < 0 is withdrawal.
"""
from dataclasses import dataclass, replace

from CoolProp.CoolProp import PropsSI

from . import props


def pressure_from_density(rho, T, gas):
    """Pressure (Pa) of a pure gas at density rho (kg/m3) and temperature T (K)."""
    if gas not in props.GASES:
        raise ValueError(f"gas must be one of {props.GASES}, got {gas!r}")
    return PropsSI("P", "D", rho, "T", T, gas)


def pore_volume_from_giip(giip_sc, p_i, T, gas="Methane"):
    """Gas pore volume (m3) that holds giip_sc (Sm3) of gas at p_i (Pa), T (K)."""
    mass = giip_sc * props.rho_sc(gas)
    return mass / props.rho(p_i, T, gas)


@dataclass(frozen=True)
class Tank:
    """Closed reservoir tank. Immutable: step() returns a new Tank.

    V    gas pore volume, m3
    T    reservoir temperature, K
    gas  "Hydrogen" or "Methane" (pure gas, decision A1)
    m    gas in place, kg
    """
    V: float
    T: float
    gas: str
    m: float

    # --- construction -----------------------------------------------------
    @classmethod
    def at_pressure(cls, V, p, T, gas):
        """Tank of volume V (m3) filled with `gas` to pressure p (Pa)."""
        return cls(V=V, T=T, gas=gas, m=props.rho(p, T, gas) * V)

    @classmethod
    def from_inputs(cls, x, gas, p=None):
        """Rough tank from inputs.load() output.

        The volume comes from the methane GIIP at the initial pressure. The tank
        is then filled with `gas` to pressure p (Pa); p defaults to p_i.
        """
        T = x["reservoir_temperature"]
        p_i = x["reservoir_pressure_initial"]
        V = pore_volume_from_giip(x["GIIP"], p_i, T, "Methane")
        return cls.at_pressure(V, p_i if p is None else p, T, gas)

    # --- state ------------------------------------------------------------
    @property
    def rho(self):
        """Average gas density in the tank, kg/m3."""
        return self.m / self.V

    @property
    def p(self):
        """Average reservoir pressure, Pa."""
        return pressure_from_density(self.rho, self.T, self.gas)

    @property
    def Z(self):
        """Compressibility factor at the current state, -."""
        return props.Z(self.p, self.T, self.gas)

    @property
    def inventory_sc(self):
        """Gas in place, Sm3."""
        return self.m / props.rho_sc(self.gas)

    def energy(self, hhv):
        """Energy of the gas in place (J) for a heating value hhv (J/kg)."""
        return self.m * hhv

    # --- time stepping ----------------------------------------------------
    def step(self, q_sc, dt):
        """Tank after a constant rate q_sc (Sm3/s; + inject, - withdraw) for dt (s).

        A closed tank's mass balance is linear, so this update is exact for a
        constant rate; only the rate itself depends on the step size.
        """
        m_new = self.m + q_sc * props.rho_sc(self.gas) * dt
        if m_new <= 0.0:
            raise ValueError(
                f"withdrawal of {-q_sc * dt:.3e} Sm3 exceeds the "
                f"{self.inventory_sc:.3e} Sm3 in place")
        return replace(self, m=m_new)


def working_gas(V, T, gas, p_low, p_high):
    """Working and cushion gas (kg) for a pressure window p_low..p_high (Pa).

    Returns (working_mass, cushion_mass). This is a static estimate: it assumes
    the wells and surface can move the reservoir across the whole window.
    """
    m_low = props.rho(p_low, T, gas) * V
    m_high = props.rho(p_high, T, gas) * V
    if m_high <= m_low:
        raise ValueError("p_high must be greater than p_low")
    return m_high - m_low, m_low