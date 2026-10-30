"""Load model inputs from data/inputs.csv and convert to SI."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TO_SI = {"psia": 6894.757, "bar": 1e5, "ft": 0.3048, "in": 0.0254, "km": 1e3,
         "mD": 9.869233e-16, "MMSCF": 28316.85, "kWh/kg": 3.6e6,
         "MW": 1e6, "m": 1.0, "-": 1.0}


def load(case="base"):
    """Return {parameter: value in SI} for case 'base', 'low' or 'high'."""
    df = pd.read_csv(ROOT / "data" / "inputs.csv")
    out = {}
    for _, r in df.iterrows():
        v = r[case] if pd.notna(r[case]) else r["base"]
        if pd.isna(v):
            continue                                   # data gap: no value yet
        units = str(r["units"]).strip()                # tolerate stray spaces in the CSV
        if units == "degC":
            out[r["parameter"]] = float(v) + 273.15    # temperature is an offset, not a factor
            continue
        if units not in TO_SI:
            raise KeyError(f"Unknown units '{units}' for {r['parameter']} in inputs.csv")
        out[r["parameter"]] = float(v) * TO_SI[units]
    return out


KSCF_TO_SM3 = 28.3168          # 1000 standard cubic feet (60 °F, 1 atm) in Sm3


def load_wells():
    """Per-well data from Heinemann et al. (2025) Table 1, converted to SI.

    Columns: well, d_wellbore (m), skin (-), D (s/Sm3), p_bhp_max (Pa), p_bhp_min (Pa).
    D is the non-Darcy coefficient: D * q is dimensionless when q is in Sm3/s.
    Source units are day/kscf, i.e. D * q with q in kscf/day.
    """
    df = pd.read_csv(ROOT / "data" / "wells.csv")
    return pd.DataFrame({
        "well": df["well"],
        "d_wellbore": df["wellbore_diameter_m"],
        "skin": df["skin"],
        "D": df["D_day_per_kscf"] * 86400 / KSCF_TO_SM3,
        "p_bhp_max": df["bhp_max_bar"] * 1e5,
        "p_bhp_min": df["bhp_min_bar"] * 1e5,
    })