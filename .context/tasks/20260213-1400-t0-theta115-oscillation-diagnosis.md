# Task: Diagnose non-steady oscillations in T0.0 theta0=115 run

## Status
- state: open
- owner_session: local
- created: 2026-02-13

## Goal
Diagnose likely causes for persistent velocity/interface oscillations in completed case without modifying code.

## Plan
- Identify exact case and confirm setup (tilt/contact angle/BCs).
- Extract solver-history evidence from logs (adaptive metrics, Courant, alpha bounds).
- Check for hidden forcing sources (g, motion tables, BCs).
- Provide ranked causes and practical validation checks.

## Log
- 2026-02-13: Identified target case `case_H0.0083_D0.0083_flat_tilt_T0.0_m0.0002` (`contact_angle=115`).
- 2026-02-13: Confirmed no external forcing (`g=(0 0 -9.81)`, `6DoF.dat` all zeros).
- 2026-02-13: Extracted adaptive history from slurm log: normU plateau ~2e-3 and maxDeltaAlpha plateau ~4e-3 through t=5.
- 2026-02-13: Observed persistent alpha undershoot near end (`Min(alpha.water)~-1.2e-4`), consistent with numerical interface ringing/spurious currents.
- 2026-02-13: Confirmed run ended by endTime (not adaptive steady detection).

## Messages
- User requested diagnosis only, no code changes.
- 2026-02-13: Built comparison case via main.py pipeline for falsification: `case_H0.0082_D0.0083_flat_tilt_T0.0_m0.0002` with `contact_angle=90`, `tilt_deg=0`, `mesh=0.0002`, `n_cpus=16`.
- 2026-02-13: Verified `case_params.json` and `0/alpha.water` contain `H=0.0082` and `theta0=90.0`.
