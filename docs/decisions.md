# Decisions and assumptions

Format: **ID – decision** – reason (alternative considered). Date added.

## Conventions
- **C1 – SI units inside all code** – avoids conversion errors; field units converted on input. 2026-10-22
- **C2 – Standard conditions 15 °C, 1 atm** – UK convention; other conventions (e.g. 60 °F, Scafidi et al., 2023) converted before comparison. 2026-10-22
- **C3 – Standard conditions differ between sources** – Heinemann uses 60 °F and 1 atm (kscf; H2 0.0851 kg/m³). D-factors are converted on that basis. 2026-10-29

## Data
- **D1 – HHVs from Scafidi et al. (2023) Table 2** – consistency with validation tool (cross-check Bossel, 2003). 2026-10-22
- **D2 – CoolProp reference equations of state** – same as Scafidi et al. (2023). 2026-10-22
- **D3 – Per-well data from Heinemann et al. (2025) Table 1** – wellbore diameter, skin, D-factor and BHP limits for the 29 wells modelled there. Diameters are wellbore, not tubing (tubing size remains a data gap). 2026-10-29

## Assumptions
- **A1 – Pure gases, no H2/CH4 mixing** – as Scafidi et al. (2023); limitation for discussion. 2026-10-22
- **A2 – Closed volumetric reservoir, no aquifer** – Rough is a closed gas field. 2026-10-22
- **A4 – BHP limits fixed at the historic natural-gas envelope** – as Heinemann et al. (2025); treated as a reservoir constraint, not infrastructure. 2026-10-29
- **A8 – Isothermal reservoir at 91.5 °C** – injected gas reaches reservoir temperature near the well (Heinemann et al., 2025, s.2.2) (alternative: thermal tank model). 2026-10-29
- **A9 – Single tank, pore volume fixed from methane GIIP at p_i, reused for hydrogen; rock and water compressibility neglected** – data available only as GIIP (Scafidi et al., 2023, Table 1). Limitation: Heinemann et al. (2025) show pressure support from the Carboniferous underburden and compartmentalisation (alternative: two-tank model as a sensitivity case). 2026-10-29
- **A10 – Material balance solved exactly from density** – p(ρ, T) from the reference equation of state, so p/Z needs no iteration (Lee and Wattenbarger, 1996). 2026-10-29

## Corrections
- **X1 – inputs.py rejects unknown units** – unconverted HHV gave zero energy values. 2026-10-22
- **X2 – MMSCF factor (28,316.85 m³) is at 60 °F, 14.696 psia, not 15 °C, 1 atm** – about 0.2 % on GIIP and the pore volume; accepted and logged, not corrected. 2026-10-29



