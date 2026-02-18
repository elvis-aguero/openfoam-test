# Task: Generate convergence plots for R1 and R3

## Status
- state: done
- owner_session: local
- created: 2026-02-18

## Goal
Generate `convergence_summary.png` for:
- `case_H0.0082_D0.0083_flat_tilt_T0.0_m0.0003_R1`
- `case_H0.0082_D0.0083_flat_tilt_T0.0_m0.0003_R3`

## Plan
- Reuse each case's local `adaptive_stop.py` plotting helper.
- Generate PNG from existing `postProcessing/adaptive_stop/adaptive_metrics.csv`.
- Verify files and explain when pipeline generates them.

## Log
- 2026-02-18: Confirmed both cases contain `adaptive_stop.py` and `postProcessing/adaptive_stop/adaptive_metrics.csv`.

## Messages
- User asked for R1/R3 convergence PNGs and when these are generated in pipeline.
- 2026-02-18: Generated both convergence PNGs via each case's `_write_convergence_summary_figure` helper.
- 2026-02-18: Verified outputs exist:
  - `case_H0.0082_D0.0083_flat_tilt_T0.0_m0.0003_R1/postProcessing/adaptive_stop/convergence_summary.png`
  - `case_H0.0082_D0.0083_flat_tilt_T0.0_m0.0003_R3/postProcessing/adaptive_stop/convergence_summary.png`
- 2026-02-18: Root-cause check: `of13 python3 -c "import matplotlib"` fails (`ModuleNotFoundError`).
- 2026-02-18: In `adaptive_stop.py`, convergence PNG generation happens at script end and is skipped silently if matplotlib import fails.

## Outcome
- Requested PNGs generated now for R1/R3.
- Missing auto-generated PNGs in pipeline are explained by missing matplotlib in `of13` Python environment.
