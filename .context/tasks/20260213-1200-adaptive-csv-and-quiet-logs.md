# Task: Add adaptive CSV metrics and reduce .out verbosity

## Status
- state: open
- owner_session: local
- created: 2026-02-13

## Goal
When running built cases, emit a CSV with timestamp, normU, maxDeltaAlpha, interface max height, and interface min height; also make Slurm `.out` logs significantly less verbose.

## Plan
- Inspect adaptive-stop and probe plumbing in `main.py` and `circularTiltingTank/adaptive_stop.py`.
- Add CSV writer in adaptive-stop loop using existing probe metrics.
- Compute interface min/max height from `alpha.water` probe columns.
- Route full solver stdout/stderr to per-case log files while keeping `.out` concise.
- Validate on a representative case and document output paths.

## Log
- 2026-02-13: Task created and code paths inspected for planning.

## Messages
- User requested plan before modifying codebase.
- 2026-02-13: Implemented probe-based adaptive CSV output in `circularTiltingTank/adaptive_stop.py` with columns timestamp, normU, maxDeltaAlpha, interfaceHeightMax, interfaceHeightMin.
- 2026-02-13: Added interface min/max estimation from alpha probe columns via linear interpolation at alpha=0.5.
- 2026-02-13: Routed foamRun stdout/stderr to `postProcessing/adaptive_stop/foamRun.stdout.log` and `postProcessing/adaptive_stop/foamRun.stderr.log` to keep Slurm `.out` concise.
- 2026-02-13: Updated `main.py` probe generation to 100 near-wall probes (20 azimuth x 5 z-levels) for better interface-height estimates in newly built cases.
- 2026-02-13: Syntax-checked modified Python files with `python3 -m py_compile`.
- 2026-02-13: User requested main.py pipeline only for runtime validation.
- 2026-02-13: Built and submitted fresh coarse case via main.py pipeline: `case_H0.0083_D0.0083_flat_tilt_T0.0_m0.001` with dt=0.05 and duration=0.3.
- 2026-02-13: Slurm job 366965 completed from t=0; `.out` showed adaptive summary lines while solver verbosity moved to `postProcessing/adaptive_stop/foamRun.stdout.log`.
- 2026-02-13: Verified CSV at `postProcessing/adaptive_stop/adaptive_metrics.csv` with requested columns and non-empty rows.
- 2026-02-13: Updated adaptive CSV interface metrics to be relative to still interface level (H/2 from case_params.json; fallback to first observed midpoint).
- 2026-02-13: Renamed CSV columns to `interfaceHeightMaxRel` and `interfaceHeightMinRel` and updated console summary to `hRel=[min,max]`.
- 2026-02-13: Syntax check passed for `circularTiltingTank/adaptive_stop.py`.
- 2026-02-13: Updated `compare_interfaces` summary output in `main.py` to include undisturbed interface area `A0=pi*(D/2)^2` and `l2_rms_over_undisturbed_area = l2_rms / A0`.
- 2026-02-13: Added `undisturbed_interface_area` and `l2_rms_over_undisturbed_area` keys to `postProcessing/interface_compare/comparison_summary.csv`.
- 2026-02-13: Syntax check passed for `main.py`.
- 2026-02-13: Updated adaptive-stop `.out` heartbeat to include solver runtime metrics parsed from solver log: `dt`, `Co_max`, `iCo_max`.
- 2026-02-13: Added incremental parser for `foamRun.stdout.log` lines (`Time`, `deltaT`, `Courant Number`, `Interface Courant Number`) to keep `.out` concise but informative.
- 2026-02-13: Syntax check passed for `circularTiltingTank/adaptive_stop.py`.
