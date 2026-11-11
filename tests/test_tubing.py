"""Tests for h2_bottleneck.tubing (average temperature and Z method, Scafidi et al., 2023, eqs. 8-13)."""
import math

import pytest

from h2_bottleneck import inputs
from h2_bottleneck import tubing as tb
from h2_bottleneck import well as wl

GASES = ("Hydrogen", "Methane")


@pytest.fixture
def tub():
    return tb.Tubing(md=2743.2, tvd=2743.2, d=7 * 0.0254, eps=1.524e-5, T_wh=288.15, T_bh=364.65)


@pytest.fixture(scope="module")
def field():
    x = inputs.load()
    f = wl.Formation.from_inputs(x)
    return x, f, inputs.load_wells().set_index("well", drop=False)


# ---- Friction factor (eq. 11) -------------------------------------------------
def test_friction_factor():
    assert tb.friction_factor(1000.0, 0.0) == pytest.approx(0.064)                          # laminar 64/Re
    assert tb.friction_factor(1e5, 0.0) == pytest.approx(0.316 / 1e5**0.25, rel=0.02)       # Blasius, smooth
    assert tb.friction_factor(1e9, 1e-3) == pytest.approx((-2 * math.log10(1e-3 / 3.7)) ** -2, rel=0.01)  # fully rough
    Re, rr = 3e6, 3e-4                                                                       # eq. 11 with log10
    assert tb.friction_factor(Re, rr) == pytest.approx(4 * (2.28 - 4 * math.log10(rr + 21.25 / Re**0.9)) ** -2)


# ---- Gas column -----------------------------------------------------------------
def test_static_column_ideal_gas_limit():
    # Low-pressure isothermal H2: p_bh = p_wh * exp(M g L / (R T)), Z ~ 1 (hand calculation)
    t = tb.Tubing(md=2800.0, tvd=2800.0, d=0.15, eps=1.5e-5, T_wh=350.0, T_bh=350.0)
    expected = 3e5 * math.exp(wl.MOLAR_MASS["Hydrogen"] * tb.G * t.tvd / (8.314462618 * 350.0)) - 3e5
    assert tb.pwf_from_pwh(3e5, 0.0, t, "Hydrogen") - 3e5 == pytest.approx(expected, rel=0.01)


def test_hydrogen_column_much_lighter(tub):
    dp = {g: tb.pwf_from_pwh(150e5, 0.0, tub, g) - 150e5 for g in GASES}
    assert 6.0 < dp["Methane"] / dp["Hydrogen"] < 14.0


# ---- Flow -----------------------------------------------------------------------
@pytest.mark.parametrize("gas", GASES)
def test_down_then_up_returns_start(tub, gas):
    for q in (-30.0, 0.0, 30.0):                                  # Sm3/s, ~2.6 million Sm3/d
        p_wf = tb.pwf_from_pwh(120e5, q, tub, gas)
        assert tb.pwh_from_pwf(p_wf, q, tub, gas) == pytest.approx(120e5, rel=1e-5)


@pytest.mark.parametrize("gas", GASES)
def test_friction_sign(tub, gas):
    static = tb.pwf_from_pwh(120e5, 0.0, tub, gas)
    assert tb.pwf_from_pwh(120e5, -30.0, tub, gas) < static     # injection: friction lowers BHP
    assert tb.pwf_from_pwh(120e5, +30.0, tub, gas) > static     # withdrawal: friction raises BHP


def test_tubing_limited_gives_nan():
    t = tb.Tubing(md=2743.2, tvd=2743.2, d=2 * 0.0254, eps=1.524e-5, T_wh=288.15, T_bh=364.65)
    assert math.isnan(tb.pwh_from_pwf(300e5, 20.0, t, "Methane"))    # withdrawal: pressure runs out
    assert math.isnan(tb.pwh_from_pwf(300e5, -40.0, t, "Methane"))   # injection: needs > 1000 bar


# ---- V2.3: average method against RK4 ----------------------------------------
@pytest.mark.parametrize("gas", GASES)
@pytest.mark.parametrize("q", [-30.0, 0.0, 30.0])
def test_average_method_against_rk4(tub, gas, q):
    column = tb.pwf_from_pwh(150e5, 0.0, tub, gas) - 150e5
    assert abs(tb.pwf_from_pwh(150e5, q, tub, gas) - tb.pwf_rk4(150e5, q, tub, gas)) < 0.04 * column


def test_rk4_converged(tub):
    assert tb.pwf_rk4(150e5, -30.0, tub, "Methane", n=10) == pytest.approx(
        tb.pwf_rk4(150e5, -30.0, tub, "Methane", n=60), abs=10.0)                # within 0.1 mbar


# ---- Inputs and coupling with the inflow model ------------------------------
def test_from_inputs(field):
    x, f, _ = field
    t = tb.Tubing.from_inputs(x, f)
    assert t.tvd == pytest.approx(9000 * 0.3048) and t.T_bh == pytest.approx(91.5 + 273.15)
    with pytest.raises(ValueError):
        tb.Tubing(md=2000.0, tvd=2500.0, d=0.15, eps=0.0, T_wh=300.0, T_bh=360.0)


def test_operating_point_satisfies_both_models(field):
    x, f, wells = field
    w = wl.Well.from_row(wells.loc["IA04"])
    t, m = tb.Tubing.from_inputs(x, f), wl.Pseudopressure(f.T, "Hydrogen")
    op = tb.operating_point(230e5, 220e5, w, f, t, "Hydrogen", m)
    assert op["q"] < 0 and op["limit"] == "wellhead"
    assert tb.pwf_from_pwh(230e5, op["q"], t, "Hydrogen") == pytest.approx(op["p_wf"], abs=200.0)


def test_injection_allowed_below_bhp_min(field):
    # BHP_min limits withdrawal only: a reservoir below it must still accept gas
    x, f, wells = field
    w = wl.Well.from_row(wells.loc["IA01"])                    # BHP_min 168.5 bar
    t, m = tb.Tubing.from_inputs(x, f), wl.Pseudopressure(f.T, "Hydrogen")
    op = tb.operating_point(160e5, 150e5, w, f, t, "Hydrogen", m)
    assert op["q"] < 0 and op["limit"] == "wellhead"


def test_manifold_whp_for_rate_round_trip(field):
    x, f, wells = field
    ws = [wl.Well.from_row(r) for _, r in wells.iterrows()][:5]
    t, m = tb.Tubing.from_inputs(x, f), wl.Pseudopressure(f.T, "Hydrogen")
    tubs = {w.name: t for w in ws}
    r = tb.manifold_whp_for_rate(-40.0, 220e5, ws, f, tubs, "Hydrogen", m)
    assert r["reachable"]
    assert tb.manifold_rate(r["p_wh"], 220e5, ws, f, tubs, "Hydrogen", m)["q"] == pytest.approx(-40.0, rel=1e-3)