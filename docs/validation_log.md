# Validation log

Every check of the model against an independent value: what was checked, what it should give and where that value comes from, what the model gave, and the commit of the code that gave it.

- **This file** records **checks**, meaning evidence that the code is right.
- **`decisions.md`** records **judgement calls**: assumptions and choices, each with its reason and source.

Fill in the **Commit** column with `git rev-parse --short HEAD` for the commit holding the code that produced the result. Never delete a failed check: add a row for the fix and the re-run.

**Result codes:**

- **Pass:** within the tolerance.
- **Consistent:** a comparison with a different method, where the difference is explained in the notes.
- **Fail:** outside the tolerance.
- **Pending:** not run yet.

IDs follow the validation ladder in `nodal_model_design.md` (V1–V6).

## V1: Properties and material balance

| ID | Date | Commit | Check | Reference value (source) | Obtained | Tolerance | Result | Notes |
|---|---|---|---|---|---|---|---|---|
| V1.1 | 2026-10-22 | 62b3661 | H2 density, 300 K, 100 bar | 7.63 kg/m³ (NIST Chemistry WebBook; Leachman et al., 2009) | | ±0.3 kg/m³ | Pass | `test_props.py` |
| V1.2 | 2026-10-22 | 62b3661 | H2 density at standard conditions | 0.085 kg/Sm³ at 15 °C, 1 atm | | ±0.002 | Pass | `test_props.py` |
| V1.3 | 2026-10-22 | 62b3661 | Z at 250 bar, 91.5 °C | Z(H2) > 1, Z(CH4) < 1 | | Sign | Pass | `test_props.py` |
| V1.4 | 2026-10-22 | 62b3661 | GIIP unit conversion | 366,000 MMSCF ≈ 10.36 × 10⁹ Sm³ | | ±0.05 × 10⁹ | Pass | `test_inputs.py` |
| V1.5 | 2026-10-25 | 62b3661 | Energy volume multiplier, H2 vs CH4 | About 3× (DESNZ, 2023, p.13) | 3.19–3.65× over 20–300 bar | — | Consistent | Notebook 01. Higher than 3× at reservoir pressure because Z(H2) > 1 and Z(CH4) < 1 |
| V1.6 | 2026-10-29 | e35c617 | Material balance round trip: CH4 tank at p_i | 366 bcf at 312.5 bar (Scafidi et al., 2023, Table 1) | 366.0 bcf at 312.5 bar | 10⁻⁶ relative | Pass | Notebook 02; `test_reservoir.py` |
| V1.7 | 2026-10-29 | e35c617 | p/Z linear in inventory (closed tank) | Constant ratio (Lee and Wattenbarger, 1996) | Varies by 1.2 × 10⁻⁷ % (CH4), 1.3 × 10⁻¹³ % (H2) | 10⁻⁶ relative | Pass | Notebook 02, Figure 02a |
| V1.8 | 2026-10-29 | e35c617 | Mass conservation over a cycle | Returns to 180.00 bar | 180.00 bar | 10⁻⁹ relative | Pass | Notebook 02, Figure 02 cycle |
| V1.9 | 2026-10-29 | e35c617 | Pore volume against hand estimate | About 4.3 × 10⁷ m³ (ideal-gas hand calculation) | 4.28 × 10⁷ m³ | 3.5–5.0 × 10⁷ | Pass | `test_reservoir.py` |
| V1.10 | 2026-10-29 | e35c617 | Historic CH4 window against BHP limits | Inside per-well BHP limits of 98–185 and 282–306 bar (Heinemann et al., 2025, Table 1); inventories from Centrica (2006, p.6) | 147–248 bar | Inside limits | Consistent | Average reservoir pressure must lie between the injection and withdrawal BHPs |
| V1.11 | 2026-10-29 | e35c617 | Static H2 working gas, 180–275 bar | 7.5 TWh working, 16.0 TWh cushion, ratio 2.1 (Heinemann et al., 2025) | 8.49 TWh, 18.44 TWh, ratio 2.17 | Same order | Consistent | Both about 13–15 % high with the ratio within 3 %, so the difference is a larger connected volume (decision A9), not deliverability. Part of V4 |
| V1.12 | 2026-10-29 | e35c617 | Total stored H2 energy against Amid et al. (2016) | 25.2 TWh, lower heating value (Amid et al., 2016, s.4) | About 25.4 TWh LHV (30.07 TWh HHV at p_i) | — | Consistent | Their 42 % uses a deliverability-set window, so only the totals are compared. Check against the published version |

## V2–V6: To come

| ID | Check | Reference value (source) | Result |
|---|---|---|---|
| V2.1 | Rough CH4, 30 wells, 2 in tubing: cushion gas | About 70 % (Scafidi et al., 2023, Table 4); check standard conditions first (decision C3) | Pending |
| V2.2 | Rough CH4, 30 wells, 8 in tubing: cushion gas | About 35 % (Scafidi et al., 2023, Table 4) | Pending |
| V2.3 | Gas column, CH4 and H2, 2.7 km | Hand calculation, ∫ρg dz | Pending |
| — | Pipeline against a worked example | Menon (2005), general flow equation example | Pending |
| — | Compressor single stage against ideal-gas polytropic formula | Hand calculation | Pending |
| V3 | Field CH4 rates | 44.8 million Sm³/d withdrawal; about 15 average injection (Centrica, 2006) | Pending |
| V4 | H2 with surface constraints off, against Heinemann et al. (2025) | 7.5 TWh (29 wells); 6.8 TWh (20 wells) | Pending (static tank: V1.11) |
| V6 | Compression energy against fuel and CO2 data | Centrica Storage (2015) | Pending |
