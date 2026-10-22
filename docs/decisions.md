# Decisions and assumptions

Format: **ID – decision** – reason (alternative considered). Date added.

## Conventions
- **C1 – SI units inside all code** – avoids conversion errors; field units converted on input. 2026-10-22
- **C2 – Standard conditions 15 °C, 1 atm** – UK convention; other conventions (e.g. 60 °F, Scafidi et al., 2023) converted before comparison. 2026-10-22

## Data
- **D1 – HHVs from Scafidi et al. (2023) Table 2** – consistency with validation tool (cross-check Bossel, 2003). 2026-10-22
- **D2 – CoolProp reference equations of state** – same as Scafidi et al. (2023). 2026-10-22

## Assumptions
- **A1 – Pure gases, no H2/CH4 mixing** – as Scafidi et al. (2023); limitation for discussion. 2026-10-22
- **A2 – Closed volumetric reservoir, no aquifer** – Rough is a closed gas field. 2026-10-22

## Corrections
- **X1 – inputs.py rejects unknown units** – unconverted HHV gave zero energy values. 2026-10-22