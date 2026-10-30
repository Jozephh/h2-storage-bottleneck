import numpy as np
import pytest

from h2_bottleneck import inputs, props
from h2_bottleneck import well as wl


@pytest.fixture(scope="module")
def setup():
    x = inputs.load()
    f = wl.Formation.from_inputs(x)
    wells = inputs.load_wells().set_index("well", drop=False)
    m = {g: wl.Pseudopressure(f.T, g) for g in ("Hydrogen", "Methane")}
    return f, wells, m


def test_pseudopressure_matches_ideal_gas_at_low_pressure(setup):
    f, _, m = setup
    p1, p2 = 2e5, 5e5                      # near-ideal: Z ~ 1, mu ~ constant
    mu = props.mu(3.5e5, f.T, "Hydrogen")
    expected = (p2**2 - p1**2) / mu        # m = p^2 / (mu Z) for an ideal gas
    assert m["Hydrogen"](p2) - m["Hydrogen"](p1) == pytest.approx(expected, rel=0.01)


def test_zero_rate_gives_reservoir_pressure(setup):
    f, wells, m = setup
    w = wl.Well.from_row(wells.loc["IA02"])
    assert wl.pwf_from_rate(0.0, 200e5, w, f, "Hydrogen", m["Hydrogen"]) == pytest.approx(200e5)


@pytest.mark.parametrize("gas", ["Hydrogen", "Methane"])
def test_withdrawal_lowers_and_injection_raises_pwf(setup, gas):
    f, wells, m = setup
    w = wl.Well.from_row(wells.loc["IB04"])
    q = 10.0                               # Sm3/s ~ 0.86 million Sm3/d
    assert wl.pwf_from_rate(+q, 200e5, w, f, gas, m[gas]) < 200e5
    assert wl.pwf_from_rate(-q, 200e5, w, f, gas, m[gas]) > 200e5


@pytest.mark.parametrize("gas", ["Hydrogen", "Methane"])
@pytest.mark.parametrize("q", [-20.0, -5.0, 5.0, 20.0])
def test_rate_pwf_round_trip(setup, gas, q):
    f, wells, m = setup
    w = wl.Well.from_row(wells.loc["IA02"])
    pwf = wl.pwf_from_rate(q, 220e5, w, f, gas, m[gas])
    assert wl.rate_from_pwf(pwf, 220e5, w, f, gas, m[gas]) == pytest.approx(q, rel=1e-3)


def test_hydrogen_non_darcy_coefficient_smaller(setup):
    f, wells, _ = setup
    w = wl.Well.from_row(wells.loc["IA02"])
    ratio = wl.non_darcy_D(w, 200e5, f.T, "Hydrogen") / w.D_ref
    assert 0.1 < ratio < 0.4               # (M_H2/M_CH4) * (mu_CH4/mu_H2) ~ 0.13 * 1.7


def test_hydrogen_flows_more_than_methane_for_same_pressures(setup):
    f, wells, m = setup
    w = wl.Well.from_row(wells.loc["IA02"])
    q_h2 = wl.rate_from_pwf(180e5, 220e5, w, f, "Hydrogen", m["Hydrogen"])
    q_ch4 = wl.rate_from_pwf(180e5, 220e5, w, f, "Methane", m["Methane"])
    assert q_h2 > q_ch4 > 0                # lower viscosity -> larger pseudopressure drop


def test_rate_too_high_raises(setup):
    f, wells, m = setup
    w = wl.Well.from_row(wells.loc["IA02"])
    with pytest.raises(ValueError):
        wl.pwf_from_rate(1e5, 50e5, w, f, "Methane", m["Methane"])
