# Task: Verify Non-90 Theta0 Build Enforcement

## Status
- state: done
- owner_session: local
- created: 2026-02-12

## Goal
Build a dry-run case with `theta0 != 90`, verify `0/alpha.water` is set at build time, and minimally patch build logic if not.

## Plan
- Trigger a dry build path with `contact_angle = 110`.
- Inspect generated `0/alpha.water` and `case_params.json`.
- If mismatch exists, patch minimal build-time logic and re-run verification.

## Log
- 2026-02-12: Task created.
- 2026-02-12: Reproduced mismatch via dry build: `case_params.json` stored `contact_angle=110`, but `0/alpha.water` remained `theta0 90`.
- 2026-02-12: Identified root cause in `main.py`: wall-block parser only matched `walls {` on same line; template uses two-line block (`walls` then `{`).
- 2026-02-12: Patched `_patch_alpha_water_bc` and `_set_contact_angle` to support both styles.
- 2026-02-12: Verified fix with dry build `case_H0.0083_D0.0083_flat_tilt_T0.0_m0.00012`: `case_params.json` has `contact_angle=115.0` and `0/alpha.water` has `theta0 115.0`.

## Messages
- User requested dry build validation and minimal fix if theta is not enforced.
