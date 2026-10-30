from h2_bottleneck import inputs
import pytest

def test_inputs_load_and_convert():
    x = inputs.load()
    assert abs(x["reservoir_pressure_initial"] / 1e5 - 312.5) < 1   # 4533 psia ~ 312.5 bar
    assert x["number_of_wells"] == 24


def test_giip_converted_to_standard_m3():
    x = inputs.load()
    assert abs(x["GIIP"] / 1e9 - 10.36) < 0.05     # 366,000 MMSCF ~ 10.4 billion Sm3


def test_wells_loaded_and_converted():
    w = inputs.load_wells()
    assert len(w) == 29                                    # Heinemann Table 1
    assert (w["p_bhp_max"] > w["p_bhp_min"]).all()
    assert w["p_bhp_max"].max() / 1e5 == pytest.approx(306.4)   # IB04
    assert w["p_bhp_min"].min() / 1e5 == pytest.approx(98.2)    # IB08
    ia02 = w.set_index("well").loc["IA02"]
    assert ia02["D"] == pytest.approx(0.000364 * 86400 / 28.3168)  # ≈ 1.111 s/Sm3
    assert ia02["skin"] == 7