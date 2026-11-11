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

## V2: Well model (inflow and tubing)

| ID | Date | Commit | Check | Reference value (source) | Obtained | Tolerance | Result | Notes |
|---|---|---|---|---|---|---|---|---|
| V2.0 | 2026-11-11 | cb53c93 | Rebuild of Scafidi et al. (2023) supplementary code reproduces their published Rough CGRs (30 wells, CH4 and H2, 2–8 in) | Published values (Scafidi et al., 2023, si_002; held in `PUBLISHED_ROUGH_30`) | Maximum difference 0.02 percentage points over 8 cases | ±0.5 points | Pass | `validation/scafidi2023.py`; notebook 04b. Confirms the rebuild is faithful, so the differences below are attributable |
| V2.1 | 2026-11-11 | cb53c93 | Rough, 30 wells, 2 in: cushion gas requirement | Scafidi et al. (2023) equations as printed (their code with the four departures corrected): CH4 86.7 %, H2 65.3 % | CH4 87.3 %, H2 66.3 % (generic well) | ±2 points | Pass | Published 69.6 % / 35.5 %. Difference due to temperature in kelvin in eq. 12, ln in eqs. 3 and 11, specific gravity at 100 bar (`04b_scafidi_attribution.csv`). With the 29 Heinemann wells: 87.3 % / 66.2 % (tubing-limited, so well data have little effect) |
| V2.2 | 2026-11-11 | cb53c93 | Rough, 30 wells, 8 in: cushion gas requirement | Corrected, as V2.1: CH4 24.3 %, H2 4.9 % | CH4 24.5 %, H2 3.7 % (generic well) | ±2 points | Pass | Published 35.0 % / 13.6 %. Largest difference over all four sizes 1.2 points (CH4, 4 in). With the 29 Heinemann wells: 33.3 % / 8.7 % (inflow-limited, so well data matter) |
| V2.3 | 2026-11-11 | cb53c93 | Average temperature and Z method (eq. 12) against RK4 integration, WHP 150 bar, both gases, 7 in and 2 in | RK4 with local properties, 60 segments (converged: 10 vs 60 segments < 0.1 mbar) | 7 in: H2 ≤ 0.46 %, CH4 ≤ 2.31 % of gas column (≤ 0.69 bar); 2 in: H2 ≤ 0.60 %, CH4 up to 4.61 % (1.37 bar, injection 0.2 million Sm³/d) | 4 % of gas column at 7 in (base case) | Pass | Also enforced by `test_average_method_against_rk4`. Average method slightly underestimates the CH4 column (−0.53 bar static). Larger CH4 error in 2 in injection is expected: friction nearly cancels the column, so small errors in each term are amplified |
| V2.4 | 2026-11-11 | cb53c93 | Scafidi et al. (2023) Table 4 field cases with their code corrected (CH4, 8 in) | Reported actual CGR: Grijpskerk 88.5 % (Juez-Larré et al., 2016); Hatfield Moors 66.7 % (Chaudry et al., 2008) | Grijpskerk (1 well): 88.6 % as published → 87.9 % corrected. Hatfield Moors (2 wells): 68.0 % → 48.5 % | — | Consistent | Grijpskerk agreement survives correction; Hatfield Moors agreement does not, suggesting compensating differences in the published code. Their model only; not this model |

## V3–V6: To come

| ID | Check | Reference value (source) | Result |
|---|---|---|---|
| — | Pipeline against a worked example | Menon (2005), general flow equation example | Pending |
| — | Compressor single stage against ideal-gas polytropic formula | Hand calculation | Pending |
| V3 | Field CH4 rates | 44.8 million Sm³/d withdrawal; about 15 average injection (Centrica, 2006) | Pending |
| V4 | H2 with surface constraints off, against Heinemann et al. (2025) | 7.5 TWh (29 wells); 6.8 TWh (20 wells) | Pending (static tank: V1.11) |
| V6 | Compression energy against fuel and CO2 data | Centrica Storage (2015) | Pending |