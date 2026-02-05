# Task: Investigate OpenFOAM case bug in main.py

## Status
- state: open
- owner_session: local
- created: 2026-02-02

## Goal
Diagnose and fix reported bug in `main.py` related to newly built OpenFOAM case.

## Plan
- Collect bug details from user (symptoms, logs, repro steps).
- Inspect relevant code paths in `main.py` and case files.
- Implement fix and validate with dry-run diagnostics only.

## Log
- 2026-02-03: Build now skips only if case_params.json matches; otherwise creates a suffixed case name (_1, _2, ...).
- 2026-02-03: Added capillary-length-normalized L2 RMS metric to interface comparison summary.
- 2026-02-03: Added case status classification (NEW/UNFINISHED/FINISHED) based on latest time folders vs duration; menu now reflects accurate status.
- 2026-02-03: Added case_params.json persistence and read-only param loading; run-time no longer overwrites controlDict based on name parsing. Menu status now uses persisted duration.
- 2026-02-03: Fixed resume detection to ignore bare processor directories (no time folders), ensuring `setFields` runs for fresh cases; updated local run and Slurm script logic.
- 2026-02-02: Fixed snappy parallel launch (mpirun), updated surfaceFeaturesDict usage, and corrected snappyHexMesh geometry/file paths; regeneration required for existing cases.
- 2026-02-02: Fixed snappyHexMesh pipeline: OF13 geometry file syntax, surfaceFeaturesDict, parallel snappy mesh build, and snappy case inference in meshTool.
- 2026-02-02: Performed static review of snappyHexMesh pipeline; reported issues and risks to user.
- 2026-02-02: Added mergeTolerance to generated snappyHexMeshDict to satisfy OF13 requirements.
- 2026-02-02: Added snappyHexMesh mesher option (STL generation + blockMesh/snappy dicts) and Makefile support; main.py now supports mesher selection and passes MESH_TOOL.
- 2026-02-02: Hardened resource estimator with capillary and viscous timestep limits (sigma/rho/nu from case or template).
- 2026-02-02: Fixed local run to honor decomposeParDict CPU count; added helper to read numberOfSubdomains.
- 2026-02-02: Task created; awaiting bug details from user.

## Messages
- User reported they believe they found a bug after building an OpenFOAM case, details pending.
