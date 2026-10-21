# h2-storage-bottleneck – which infrastructure component limits hydrogen storage? (Rough case study)

MEng research project (Heriot-Watt, B41FP 2026–27). Identifies whether wells, the
export pipeline or compression limits hydrogen storage service at the Rough field.

## Reproduce
```bash
conda env create -f environment.yml
conda activate b41fp
nbstripout --install          # strips notebook outputs before commit
pytest                        # all checks should pass
jupyter lab                   # run notebooks/ in numerical order
```

[![tests](https://github.com/<Jozephh>/h2-storage-bottleneck/actions/workflows/tests.yml/badge.svg)](https://github.com/<Jozephh>/h2-storage-bottleneck/actions)

Code: MIT licence (see `LICENSE`). Cite via `CITATION.cff`.

## Layout
| Folder | Contents |
|---|---|
| `data/inputs.csv` | Every model input with source, units and low/base/high range |
| `src/h2_bottleneck/` | Model code (properties, reservoir, well, pipeline, compressor, network) |
| `notebooks/` | Numbered analyses that call `src/` – no physics defined in notebooks |
| `tests/` | Unit and validation tests |
| `results/` | Generated figures/tables (regenerate from code) |
| `docs/validation_log.md` | Record of every validation check and its outcome |
| `docs/decisions.md` | Modelling decisions and assumptions, with reasons |

## Conventions
- SI units inside all code; convert only at input/output.
- No hard-coded inputs: everything is read from `data/inputs.csv`.
- Each result in the report cites the git commit (`git rev-parse --short HEAD`) that produced it.
