import pytest
from h2_bottleneck import props


def test_hydrogen_density_reference_value():
    # NIST: H2 at 300 K, 100 bar ~ 7.63 kg/m3
    assert props.rho(100e5, 300.0, "Hydrogen") == pytest.approx(7.8, abs=0.3)


def test_standard_density_hydrogen():
    # ~0.085 kg/Sm3 at 15 degC, 1 atm
    assert props.rho_sc("Hydrogen") == pytest.approx(0.085, abs=0.002)


def test_hydrogen_Z_above_one_methane_below_one():
    T = 364.65  # 91.5 degC, Rough reservoir (Heinemann et al., 2025)
    assert props.Z(250e5, T, "Hydrogen") > 1.0
    assert props.Z(250e5, T, "Methane") < 1.0


def test_unknown_gas_rejected():
    with pytest.raises(ValueError):
        props.rho(1e5, 300.0, "Nitrogen")
