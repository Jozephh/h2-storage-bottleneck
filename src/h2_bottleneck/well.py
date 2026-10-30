"""Well inflow performance relationship (IPR) for a real gas. SI units throughout.

Pseudosteady-state radial inflow with non-Darcy (turbulent) flow near the well
(Houpeurt, 1959; Lee and Wattenbarger, 1996, ch. 4), the same equation used by
Scafidi et al. (2023, eq. 1) but written in SI units:

    m(p_res) - m(p_wf) = C * q * (L + s + D * |q|)

    m(p)   real-gas pseudopressure, 2 * integral of p / (mu * Z) dp      [Pa/s]
    C      p_sc * T / (pi * k * h * T_sc)                                [Pa/s per Sm3/s]
    L      ln(r_e / r_w) - 3/4   (circular drainage area, well at centre) [-]
    s      skin factor                                                   [-]
    D      non-Darcy coefficient                                         [s/Sm3]
    q      standard volumetric rate per well, Sm3/s: + withdrawal, - injection

The same equation covers injection: a negative q makes the right-hand side
negative, so m(p_wf) > m(p_res), i.e. p_wf > p_res. D * |q| is always added,
because turbulence resists flow in both directions.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid

from h2_bottleneck import props

# Molar masses, kg/mol (NIST / CoolProp values)
MOLAR_MASS = {"Hydrogen": 2.01588e-3, "Methane": 16.0428e-3}


class Pseudopressure:
    """Lookup table of real-gas pseudopressure m(p) for one gas at one temperature.

    m(p) = 2 * integral from p_min to p of p' / (mu(p') * Z(p')) dp'
    Built once (CoolProp is slow), then evaluated by linear interpolation.
    m increases monotonically with p, so the inverse p(m) is also an interpolation.
    """

    def __init__(self, T, gas, p_min=1e5, p_max=4e7, n=2000):
        self.T, self.gas = T, gas
        p = np.linspace(p_min, p_max, n)
        integrand = np.array([2 * x / (props.mu(x, T, gas) * props.Z(x, T, gas)) for x in p])
        self.p = p
        self.m = np.concatenate([[0.0], cumulative_trapezoid(integrand, p)])

    def __call__(self, p):
        p = np.asarray(p, dtype=float)
        if np.any(p < self.p[0]) or np.any(p > self.p[-1]):
            raise ValueError(f"pressure outside pseudopressure table ({self.p[0]:.3g}-{self.p[-1]:.3g} Pa)")
        return np.interp(p, self.p, self.m)

    def inverse(self, m):
        m = np.asarray(m, dtype=float)
        if np.any(m < self.m[0]) or np.any(m > self.m[-1]):
            raise ValueError("pseudopressure outside table: rate too high for this reservoir pressure")
        return np.interp(m, self.m, self.p)


@dataclass
class Formation:
    """Reservoir properties seen by every well (field averages, A11)."""
    k: float    # permeability to gas, m2
    h: float    # net thickness, m
    T: float    # reservoir temperature, K
    re: float   # drainage radius, m

    @classmethod
    def from_inputs(cls, x):
        return cls(k=x["permeability"], h=x["net_thickness"],
                   T=x["reservoir_temperature"], re=x["drainage_radius"])


@dataclass
class Well:
    """One well: geometry, skin, reference D-factor and bottomhole pressure limits."""
    name: str
    rw: float          # wellbore radius, m
    skin: float        # -
    D_ref: float       # non-Darcy coefficient measured for the reference gas, s/Sm3
    p_bhp_max: float   # Pa
    p_bhp_min: float   # Pa

    @classmethod
    def from_row(cls, r):
        """Build from one row of inputs.load_wells()."""
        return cls(name=r["well"], rw=r["d_wellbore"] / 2, skin=r["skin"], D_ref=r["D"],
                   p_bhp_max=r["p_bhp_max"], p_bhp_min=r["p_bhp_min"])


def darcy_constant(f):
    """C = p_sc * T / (pi * k * h * T_sc), Pa/s per Sm3/s."""
    return props.P_SC * f.T / (np.pi * f.k * f.h * props.T_SC)


def geometry_term(f, w):
    """L = ln(r_e / r_w) - 3/4 for a circular drainage area with the well at its centre."""
    return np.log(f.re / w.rw) - 0.75


def non_darcy_D(w, p, T, gas, ref_gas="Methane"):
    """Non-Darcy coefficient for `gas`, scaled from the measured reference-gas value.

    D is proportional to M / mu(p_wf) for the same rock and well (Lee and
    Wattenbarger, 1996; Scafidi et al., 2023, eq. 5), so
        D_gas = D_ref * (M_gas / M_ref) * (mu_ref / mu_gas), both viscosities at p, T.
    Heinemann et al. (2025) kept the methane value for hydrogen (decision A6).
    """
    if gas == ref_gas:
        return w.D_ref
    return w.D_ref * (MOLAR_MASS[gas] / MOLAR_MASS[ref_gas]) * (props.mu(p, T, ref_gas) / props.mu(p, T, gas))


def pwf_from_rate(q, p_res, w, f, gas, mfun, tol=1.0, max_iter=50):
    """Bottomhole flowing pressure (Pa) for a well rate q (Sm3/s; + withdrawal, - injection).

    D depends on viscosity at p_wf, which is unknown, so iterate: start at p_res,
    compute D, solve for p_wf through the inverse pseudopressure, repeat to tol (Pa).
    """
    C, L, m_res = darcy_constant(f), geometry_term(f, w), mfun(p_res)
    pwf = p_res
    for _ in range(max_iter):
        D = non_darcy_D(w, pwf, f.T, gas)
        new = float(mfun.inverse(m_res - C * q * (L + w.skin + D * abs(q))))
        if abs(new - pwf) < tol:
            return new
        pwf = new
    raise RuntimeError(f"p_wf did not converge for well {w.name}")


def rate_from_pwf(pwf, p_res, w, f, gas, mfun):
    """Well rate (Sm3/s; + withdrawal, - injection) for a given bottomhole pressure.

    Rearranges dm = a*q + b*q*|q|, with a = C(L + s) and b = C*D, as a quadratic in |q|:
        |q| = (-a + sqrt(a^2 + 4*b*|dm|)) / (2b),  sign of q = sign of dm.
    """
    C, L = darcy_constant(f), geometry_term(f, w)
    a = C * (L + w.skin)
    if a <= 0:
        raise ValueError(f"L + skin must be positive (well {w.name})")
    b = C * non_darcy_D(w, pwf, f.T, gas)
    dm = float(mfun(p_res) - mfun(pwf))
    if b == 0:
        return dm / a
    return np.sign(dm) * (-a + np.sqrt(a * a + 4 * b * abs(dm))) / (2 * b)
