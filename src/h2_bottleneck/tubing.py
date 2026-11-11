"""Tubing model: pressure between bottomhole and wellhead for a pure gas. SI units.

Method as published by Scafidi et al. (2023, eqs. 8-13): the average temperature and Z
method (Lee and Wattenbarger, 1996, ch. 4), written in SI from the physics, with their
Jain (1976) friction factor (eq. 11). This follows the printed equations; their
supplementary code differs from them in four places (see validation/scafidi2023.py).

Energy balance along the tubing (z = depth along it, positive down):
    dp/dz = rho*g*cos(theta) + sign(q) * f * mdot^2 / (2*d*rho*A^2)          (1)
q > 0 withdrawal (flow up), q < 0 injection (flow down), as in well.py.

Holding Z and T at their tubing averages makes rho = (rho_bar/p_bar)*p, and (1)
integrates in closed form (Scafidi et al. eq. 12):
    p_wf^2 = p_wh^2 * e^s + (c/a) * (e^s - 1),   s = 2*a*L                    (2)
    a = rho_bar*g*cos(theta)/p_bar,   c = sign(q)*f*mdot^2*p_bar/(2*d*A^2*rho_bar)
iterated because p_bar depends on the unknown end pressure.

RK4 integration of (1) with local properties (eq. 3) is kept as the reference
solution for verification only (V2.3).

Assumptions: A1 pure gas, steady flow; A14 temperature averaged between wellhead
and reservoir; A15 one tubing size, cos(theta) = TVD/MD; kinetic energy neglected.
"""
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from h2_bottleneck import props
from h2_bottleneck import well as wl

G = 9.80665     # m/s2
P_MAX = 1000e5  # above this the tubing is treated as unable to carry the rate (and CoolProp stays in range)


def friction_factor(Re, rel_rough):
    """Darcy (Moody) friction factor, Jain (1976) as Scafidi et al. eq. 11, with log10:
    1/sqrt(f) = 1.14 - 2*log10(eps/d + 21.25/Re^0.9); 64/Re if laminar (Re <= 2000)."""
    if Re <= 0:
        return 0.0
    if Re <= 2000:
        return 64.0 / Re
    return (1.14 - 2 * math.log10(rel_rough + 21.25 / Re**0.9)) ** -2


@dataclass(frozen=True)
class Tubing:
    md: float     # measured depth to perforations, m
    tvd: float    # true vertical depth, m
    d: float      # inner diameter, m
    eps: float    # roughness, m
    T_wh: float   # wellhead temperature, K
    T_bh: float   # bottomhole (reservoir) temperature, K

    def __post_init__(self):
        if not 0 < self.tvd <= self.md:
            raise ValueError(f"need 0 < TVD <= MD, got TVD={self.tvd}, MD={self.md}")

    @property
    def area(self):
        return math.pi * self.d**2 / 4

    @property
    def cos_theta(self):
        return self.tvd / self.md

    @property
    def T_avg(self):
        return (self.T_wh + self.T_bh) / 2

    @classmethod
    def from_inputs(cls, x, f):
        """Same tubing in every well, vertical (MD = TVD = depth) unless inputs.csv has tubing_md."""
        return cls(md=x.get("tubing_md", x["depth"]), tvd=x["depth"], d=x["tubing_diameter"],
                   eps=x["tubing_roughness"], T_wh=x["wellhead_temperature"], T_bh=f.T)


# ---- Average temperature and Z method (main method) ------------------------
def _coefficients(p_bar, q, tub, gas):
    """Gas-column exponent s and friction/gravity ratio c/a (Pa^2) at the average conditions."""
    T = tub.T_avg
    rho, mu = props.rho(p_bar, T, gas), props.mu(p_bar, T, gas)
    mdot = props.rho_sc(gas) * abs(q)
    a = rho * G * tub.cos_theta / p_bar
    f = friction_factor(4 * mdot / (math.pi * tub.d * mu), tub.eps / tub.d)
    c = np.sign(q) * f * mdot**2 * p_bar / (2 * tub.d * tub.area**2 * rho)
    return 2 * a * tub.md, c / a


def _solve(p_known, q, tub, gas, downward, tol=1e-6):
    """Iterate eq. (2) from the known end; NaN if the tubing cannot carry this rate
    (pressure would fall to zero, or rise above P_MAX)."""
    p_other = p_known
    for _ in range(100):
        if p_other > P_MAX:
            return math.nan
        s, c_a = _coefficients((p_known + p_other) / 2, q, tub, gas)
        if downward:
            p2 = p_known**2 * math.exp(s) + c_a * (math.exp(s) - 1)
        else:
            p2 = (p_known**2 - c_a * (math.exp(s) - 1)) * math.exp(-s)
        if p2 <= 0:
            return math.nan
        new = math.sqrt(p2)
        if abs(new - p_other) < tol * new:
            return new if new <= P_MAX else math.nan
        p_other = (p_other + new) / 2 if _ > 20 else new       # damp if slow to converge
    return new


def pwf_from_pwh(p_wh, q, tub, gas):
    """Bottomhole pressure (Pa) from wellhead pressure (Pa) and rate q (Sm3/s)."""
    return _solve(p_wh, q, tub, gas, downward=True)


def pwh_from_pwf(p_wf, q, tub, gas):
    """Wellhead pressure (Pa) from bottomhole pressure (Pa) and rate q (Sm3/s)."""
    return _solve(p_wf, q, tub, gas, downward=False)


# ---- RK4 reference solution (verification only) ------------------------------
def pwf_rk4(p_wh, q, tub, gas, n=60):
    """Integrate eq. (1) down the tubing with RK4, properties at each local p and T (eq. 3)."""
    mdot = props.rho_sc(gas) * abs(q)

    def dpdz(z, p):
        T = tub.T_wh + (tub.T_bh - tub.T_wh) * z / tub.md
        rho, mu = props.rho(p, T, gas), props.mu(p, T, gas)
        f = friction_factor(4 * mdot / (math.pi * tub.d * mu), tub.eps / tub.d)
        return rho * G * tub.cos_theta + np.sign(q) * f * mdot**2 / (2 * tub.d * rho * tub.area**2)

    h, p = tub.md / n, p_wh
    for i in range(n):
        z = i * h
        k1 = dpdz(z, p)
        k2 = dpdz(z + h / 2, p + h / 2 * k1)
        k3 = dpdz(z + h / 2, p + h / 2 * k2)
        k4 = dpdz(z + h, p + h * k3)
        p += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return p


# ---- Tubing + inflow for one well -------------------------------------------
def whp_at_bhp_limit(p_res, w, f, tub, gas, m, mode):
    """WHP needed to hold a well at BHP_max ("injection") or BHP_min ("withdrawal")."""
    p_wf = w.p_bhp_max if mode == "injection" else w.p_bhp_min
    if (mode == "injection") != (p_wf > p_res):          # well cannot flow that way
        return dict(q=0.0, p_wf=p_wf, p_wh=math.nan)
    q = wl.rate_from_pwf(p_wf, p_res, w, f, gas, m)
    return dict(q=q, p_wf=p_wf, p_wh=pwh_from_pwf(p_wf, q, tub, gas))


def operating_point(p_wh, p_res, w, f, tub, gas, m):
    """Rate at which inflow and tubing agree for a given WHP, within the BHP limits.

    Find p_wf where the tubing, at the inflow rate for that p_wf, gives back p_wf.
    The answer lies between p_res and the static BHP (no flow). Injection is cut at
    BHP_max and withdrawal at BHP_min; beyond that limit the well is held at it.
    Returns dict(q, p_wf, limit), limit in {"none", "wellhead", "bhp_max", "bhp_min"}.
    """
    p_static = pwf_from_pwh(p_wh, 0.0, tub, gas)
    if abs(p_static - p_res) < 1.0:
        return dict(q=0.0, p_wf=p_res, limit="none")
    injecting = p_static > p_res
    lo, hi = sorted((p_static, p_res))
    if injecting:
        hi, limit = min(hi, w.p_bhp_max), "bhp_max"
    else:
        lo, limit = max(lo, w.p_bhp_min), "bhp_min"
    if lo >= hi:                                           # p_res already at or past the limit
        return dict(q=0.0, p_wf=p_res, limit=limit)

    def resid(p_wf):
        q = wl.rate_from_pwf(p_wf, p_res, w, f, gas, m)
        r = pwf_from_pwh(p_wh, q, tub, gas) - p_wf
        return r if np.isfinite(r) else math.copysign(1e12, q)   # tubing cannot carry q

    if resid(lo) * resid(hi) > 0:                          # match lies beyond the limit: hold it there
        p_wf = hi if injecting else lo
    else:
        p_wf, limit = brentq(resid, lo, hi, xtol=100.0), "wellhead"
    return dict(q=wl.rate_from_pwf(p_wf, p_res, w, f, gas, m), p_wf=p_wf, limit=limit)


# ---- All wells at one wellhead (manifold) pressure, node 3 -------------------
def manifold_rate(p_wh, p_res, wells, f, tubs, gas, m):
    """Total field rate (Sm3/s) with every well at wellhead pressure p_wh; tubs maps name -> Tubing."""
    ops = [operating_point(p_wh, p_res, w, f, tubs[w.name], gas, m) for w in wells]
    lim = [o["limit"] for o in ops]
    return dict(q=sum(o["q"] for o in ops), n_wellhead_limited=lim.count("wellhead"),
                n_bhp_limited=lim.count("bhp_max") + lim.count("bhp_min"))


def manifold_whp_for_rate(q_target, p_res, wells, f, tubs, gas, m, p_lo=10e5, p_hi=450e5):
    """Common WHP at which the wells move q_target (Sm3/s). Field rate falls as WHP rises."""
    def g(p):
        return manifold_rate(p, p_res, wells, f, tubs, gas, m)["q"] - q_target
    if not g(p_lo) >= 0 >= g(p_hi):                        # out of reach: wells at their limits
        return dict(p_wh=math.nan, reachable=False)
    return dict(p_wh=brentq(g, p_lo, p_hi, xtol=1e3), reachable=True)