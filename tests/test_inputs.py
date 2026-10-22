from h2_bottleneck import inputs


def test_inputs_load_and_convert():
    x = inputs.load()
    assert abs(x["reservoir_pressure_initial"] / 1e5 - 312.5) < 1   # 4533 psia ~ 312.5 bar
    assert x["number_of_wells"] == 24


def test_giip_converted_to_standard_m3():
    x = inputs.load()
    assert abs(x["GIIP"] / 1e9 - 10.36) < 0.05     # 366,000 MMSCF ~ 10.4 billion Sm3
