# Task 20260206-1605-plot-wall-h-theta

## Summary
- Add plot of wall interface height vs theta with OpenFOAM vs analytical overlay.

## Status
- open

## Plan
- Inspect existing plot script and interface outputs.
- Add wall-height vs theta plotting support.
- Update task log with changes.

## Log
- 2026-02-06: Task created.
- 2026-02-06: Added wall-height vs theta plotting mode to plot_l2_vs_tilt.py.
- 2026-02-06: Left task open per user authorization requirement.
- 2026-02-06: Added unit tests for wall-profile helpers.
- 2026-02-06: Removed unit tests and switched to default wall-profile plotting filtered by contact angle.
- 2026-02-06: Updated wall profile extraction to use OpenFOAM interface CSV point clouds.
- 2026-02-06: Contact angle filtering now uses case_params.json when present.

## Messages
- 
