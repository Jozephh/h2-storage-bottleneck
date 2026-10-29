import pytest

from h2_bottleneck import inputs, props
from h2_bottleneck.reservoir import Tank, pore_volume_from_giip, working_gas

DAY = 86_400.0


@pytest.fixture(scope="module")
def x():
    return inputs.load()


def test_methane_tank_returns_giip_and_initial_pressure(x):
    # Round trip: the volume is built from GIIP at p_i, so it must give them back.
    tank = Tank.from_inputs(x, "Methane")
    assert tank.inventory_sc == pytest.approx(x["GIIP"], rel=1e-6)
    assert tank.p == pytest.approx(x["reservoir_pressure_initial"], rel=1e-6)


def test_pore_volume_order_of_magnitude(x):
    # Hand estimate: ~4.3e7 m3 (10.36e9 Sm3 CH4 at 312.5 bar, 91.5 degC, Z ~ 1.0)
    V = pore_volume_from_giip(x["GIIP"], x["reservoir_pressure_initial"],
                              x["reservoir_temperature"])
    assert 3.5e7 < V < 5.0e7


def test_p_over_z_linear_in_inventory(x):
    # Closed tank: p/Z is proportional to gas in place (classic p/Z plot).
    full = Tank.from_inputs(x, "Methane")
    ratios = []
    for frac in (1.0, 0.7, 0.4, 0.1):
        t = Tank(full.V, full.T, "Methane", full.m * frac)
        ratios.append((t.p / t.Z) / t.m)
    assert max(ratios) == pytest.approx(min(ratios), rel=1e-6)


def test_inject_then_withdraw_returns_to_start(x):
    t0 = Tank.from_inputs(x, "Hydrogen", p=200e5)
    q = 40e6 / DAY                      # 40 million Sm3/d
    t1 = t0.step(q, 30 * DAY)
    assert t1.p > t0.p
    t2 = t1.step(-q, 30 * DAY)
    assert t2.p == pytest.approx(t0.p, rel=1e-9)


def test_over_withdrawal_rejected(x):
    t = Tank.from_inputs(x, "Hydrogen", p=50e5)
    with pytest.raises(ValueError):
        t.step(-t.inventory_sc, 2.0)    # twice the gas in place


def test_hydrogen_stores_about_28_percent_of_methane_energy(x):
    # Notebook 01: 27.4-31.4 % over 20-300 bar, minimum near 250 bar.
    V, T = 1.0, x["reservoir_temperature"]
    h2 = Tank.at_pressure(V, 250e5, T, "Hydrogen").energy(x["hhv_hydrogen"])
    ch4 = Tank.at_pressure(V, 250e5, T, "Methane").energy(x["hhv_methane"])
    assert 0.26 < h2 / ch4 < 0.30


def test_static_working_gas_same_order_as_heinemann(x):
    # Heinemann et al. (2025): 7.5 TWh working, 16.0 TWh cushion (ratio 2.1),
    # average pressure cycling ~180-275 bar. A static tank ignores well limits,
    # so it should land close but somewhat higher. Hand estimate: ~8.5 TWh, 2.2.
    tank = Tank.from_inputs(x, "Hydrogen")
    w, c = working_gas(tank.V, tank.T, "Hydrogen", 180e5, 275e5)
    twh = 3.6e15
    assert 7.0 < w * x["hhv_hydrogen"] / twh < 10.0
    assert 1.9 < c / w < 2.5