"""Validation against Scafidi et al. (2023), Geol. Soc. Spec. Publ. 528, 267-283 (V2.1, V2.2).

Two parts:

1. ``run``: a faithful rebuild of their supplementary code (SP528-2022-71_si_001), in their
   field units, with a switch for each place where the code departs from their printed
   equations. All switches on reproduces their published results (si_002).

     T_KELVIN  : eq. 12 friction term receives the average temperature in kelvin, not degR
     FF_LN     : eq. 11 friction factor (Jain, 1976) uses ln instead of log10
     A_LN      : eq. 3 inflow coefficient, 1.151*ln(...) instead of 1.151*log10(...)
     SG_100BAR : specific gravity = density ratio to air at 0 degC and 100 bar, not M/M_air

Their published Rough results (30 wells) are in PUBLISHED_ROUGH_30.

2. ``this_model``: this project's model (well.py inflow + tubing.py) run under Scafidi's
   conditions with their generic well, for a like-for-like comparison with (1) corrected.

Everything else in ``run`` follows their code: 1 psi pseudopressure grid, D at the viscosity
of each p_wf (eq. 5), their tubing iteration, WHP = 50 % of p_res, surface 15 degC, GIIP
scaled for H2 by the expansion-factor ratio (eq. 15), CGR from eqs. 16-18.
"""
import math

import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq

from h2_bottleneck import props
from h2_bottleneck import tubing as tb
from h2_bottleneck import well as wl

ISSUES = {"T_KELVIN": "temperature in kelvin in eq. 12", "FF_LN": "natural log in eq. 11 (friction factor)",
          "A_LN": "natural log in eq. 3 (inflow coefficient)", "SG_100BAR": "specific gravity at 100 bar"}
CORRECTED = {k: False for k in ISSUES}

# Their Table 1: p_res psia, k mD, porosity, h ft, L ft, GIIP MMscf, T_res degR
FIELDS = {
    "Rough": (4533.0, 75.0, 0.125, 95.0, 9000.0, 366000.0, 656.67),
    "Grijpskerk": (5700.0, 25.0, 0.15, 590.0, 11220.5, 381000.0, 702.27),
    "Hatfield Moors": (650.0, 143.2, 0.191, 57.5, 1400.0, 6100.0, 527.67),
    "Cousland": (660.0, 70.0, 0.15, 48.0, 1582.0, 870.0, 527.67),
}
GAS = {"Hydrogen": "H2", "Methane": "CH4"}            # this project's names -> their CoolProp names
_RE, _S, _EPS, _DD, _RF, _CA, _TSC, _PSC = 4920.0, 1.0, 0.0006, 0.5, 0.9, 31.62, 491.67, 14.5038
_T_SURF_R, _PSI, _R2K = 15 * 9 / 5 + 491.67, 6894.76, 5 / 9


# ---- 1. Their code ------------------------------------------------------------------
def run(field, gas, d_wb, n_wells, T_KELVIN=True, FF_LN=True, A_LN=True, SG_100BAR=True, tp=90.0):
    """Their model for one case. gas 'CH4' or 'H2'; d_wb in inches.
    Returns dict(q_well [MMscf/d], pwf [psia], cgr [%])."""
    PR, K, PORO, H, L, GIIP, TRES_R = FIELDS[field]
    rw, tres_k, pwh = d_wb / 12 / 2, TRES_R * _R2K, PR * _DD
    M = PropsSI("M", "T", tres_k, "P", _PSC * _PSI, gas) * 1000
    sg = (PropsSI("D", "T", 273.15, "P", 100e5, gas) / PropsSI("D", "T", 273.15, "P", 100e5, "AIR")
          if SG_100BAR else M / 28.9647)
    log = np.log if A_LN else np.log10
    A = (1.422e6 * TRES_R / (K * H)) * (1.151 * log(10.06 * math.pi * _RE**2 / (_CA * rw * rw)) - 0.75 + _S)
    beta = 1.88e10 * K**-1.47 * PORO**-0.53

    p = np.arange(1.0, PR + 1.0, 1.0)                                # pseudopressure, eq. 2
    z = np.array([PropsSI("Z", "T", tres_k, "P", x * _PSI, gas) for x in p])
    mu = np.array([PropsSI("V", "T", tres_k, "P", x * _PSI, gas) * 1000 for x in p])
    pp = np.cumsum(2 * p / (mu * z))

    def q_ipr(pwf):                                                  # eqs. 1, 4, 5, 7
        i = min(max(int(round(pwf)) - 1, 0), len(p) - 1)
        b = 1.422e6 * TRES_R * (2.715e-12 * beta * K * M * _PSC / (H * mu[i] * rw * _TSC)) / (K * H)
        return (-A + (A * A + 4 * b * (pp[-1] - pp[i])) ** 0.5) / (2 * b)

    def tpr(q_mscf):                                                 # eqs. 8-13, their iteration
        lg = np.log if FF_LN else np.log10
        tav = (TRES_R + _T_SURF_R) / 2
        pav = (pwh + 0.25 * (pwh / 100) * (L / 100) + pwh) / 2
        zav = PropsSI("Z", "T", tav * _R2K, "P", pav * _PSI, gas)
        vav = PropsSI("V", "T", tav * _R2K, "P", pav * _PSI, gas) * 1000
        nre = max((20 * sg * q_mscf) / (vav * d_wb), 100)
        f = 4 * (2.28 - 4 * lg(_EPS / d_wb + 21.25 / nre**0.9)) ** -2 if nre > 2000 else 64 / nre
        zz = vv = 1.0
        while round(zz, 3) > 0.005 and round(vv, 3) > 0.005:
            s = 0.0375 * sg * L / (zav * tav)
            t_f = tav * _R2K if T_KELVIN else tav
            pwf = (pwh**2 * math.exp(s) + 6.67e-4 * q_mscf**2 * f * t_f**2 * zav**2 / d_wb**5 * (math.exp(s) - 1)) ** 0.5
            pav2 = (pwf + pwh) / 2
            z2 = PropsSI("Z", "T", tav * _R2K, "P", pav2 * _PSI, gas)
            v2 = PropsSI("V", "T", tav * _R2K, "P", pav2 * _PSI, gas) * 1000
            zz, vv, zav, vav = (z2 - zav) / z2, (v2 - vav) / v2, z2, v2
        return pwf

    pwf = brentq(lambda pw: tpr(q_ipr(pw) * 1000) - pw, 2.0, PR - 1.0)   # IPR meets TPR
    q = q_ipr(pwf)
    if gas == "H2":                                                  # eq. 15
        e = lambda g: PR / (PropsSI("Z", "T", tres_k, "P", PR * _PSI, g) * TRES_R)
        GIIP *= e("H2") / e("CH4")
    q_plat = 1 / (1 / (q * n_wells) + tp / GIIP)                     # eq. 16
    wgv = q_plat * tp                                                # eq. 17
    return dict(q_well=q, pwf=pwf, cgr=100 * (GIIP * _RF - wgv) / (GIIP * _RF))   # eq. 18


# Their published CGR (%) for Rough, 30 wells (medium scenario), from the supplementary results
# file SP528-2022-71_si_002.csv (Scafidi et al. 2023, https://doi.org/10.6084/m9.figshare.c.6192546).
# Table 4 of the paper quotes the CH4 values for 2 in (69.6 %) and 8 in (35.0 %).
PUBLISHED_ROUGH_30 = {
    ("CH4", 2): 69.59, ("CH4", 4): 42.70, ("CH4", 6): 37.02, ("CH4", 8): 34.95,
    ("H2", 2): 35.52, ("H2", 4): 17.02, ("H2", 6): 14.59, ("H2", 8): 13.59,
}


# ---- 2. This project's model under their conditions (Rough, from inputs.csv) -------------
_P_SC, _T_SC = 1e5, 273.15                       # their standard conditions: 0 degC, 1 bar


def _to_their_m3(gas):
    return props.rho_sc(gas) / props.rho(_P_SC, _T_SC, gas)


def generic_well(d, x, f, gas):
    """Their generic well: r_w = d/2 (m), skin from inputs.csv (1), D from eq. 5 at p_res, no BHP limits."""
    k_md, h_ft, rw_ft = f.k / 9.869233e-16, f.h / 0.3048, d / 2 / 0.3048
    mu_cp = props.mu(x["reservoir_pressure_initial"], f.T, gas) * 1e3
    beta = 1.88e10 * k_md**-1.47 * x["porosity"] ** -0.53
    D_field = 2.715e-12 * beta * k_md * wl.MOLAR_MASS[gas] * 1e3 * 14.5038 / (h_ft * mu_cp * rw_ft * 491.67)
    return wl.Well(name="generic", rw=d / 2, skin=x["skin"], D_ref=D_field * 86400 * _to_their_m3(gas) / 28316.85,
                   p_bhp_max=1e9, p_bhp_min=1e5)


def this_model(d_in, gas, x, f, tub, m, wells=None, n_wells=30):
    """This model at their conditions: WHP = 50 % of p_res, tubing d_in inches.
    wells=None uses their generic well; otherwise the mean rate of the given wells.
    Returns dict(q_well [MMscf/d, their standard conditions], cgr [%])."""
    p_i = x["reservoir_pressure_initial"]
    t_d = tb.Tubing(md=tub.md, tvd=tub.tvd, d=d_in * 0.0254, eps=tub.eps, T_wh=tub.T_wh, T_bh=tub.T_bh)
    ws = wells or [generic_well(d_in * 0.0254, x, f, gas)]
    q = np.mean([tb.operating_point(0.5 * p_i, p_i, w, f, t_d, gas, m)["q"] for w in ws])   # Sm3/s
    q_day = q * 86400 * _to_their_m3(gas)                                                    # their m3/d
    giip = x["GIIP"] * (props.Z(p_i, f.T, "Methane") / props.Z(p_i, f.T, gas))               # eq. 15, m3
    wgv = 90 / (1 / (n_wells * q_day) + 90 / giip)
    return dict(q_well=q_day / 28316.85, cgr=100 * (1 - wgv / (0.9 * giip)))
