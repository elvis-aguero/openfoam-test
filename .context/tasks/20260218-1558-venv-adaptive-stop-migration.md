# Task: Migrate adaptive_stop Python execution to repo venv

## Status
- state: done
- owner_session: local
- created: 2026-02-18

## Goal
Ensure future runs use repo venv Python for adaptive stop CSV/PNG generation and only use `of13` for OpenFOAM binaries.

## Plan
- Update template Makefile to call `../sloshing/bin/python` for adaptive_stop.
- Add `--of-prefix` support in adaptive_stop.py so solver still runs via `of13` when needed.
- Add clear warnings for convergence plot skip reasons.
- Update docs to clarify Python/OpenFOAM env split.
- Validate via static checks/dry-run command expansion.

## Log
- 2026-02-18: Confirmed current OSCAR path executes `of13 python3 adaptive_stop.py` from template Makefile.

## Messages
- User requested avoiding `of13 python3` in favor of repo venv where possible.
- 2026-02-18: Updated `circularTiltingTank/Makefile` to run adaptive stop via `../sloshing/bin/python` and pass `--of-prefix` for OpenFOAM commands.
- 2026-02-18: Added fail-fast guard when venv python is missing.
- 2026-02-18: Updated `circularTiltingTank/adaptive_stop.py` with `--of-prefix`, command diagnostics, and explicit convergence-plot skip warnings.
- 2026-02-18: Updated docs (`USER_GUIDE.md`, `KNOWLEDGE_BASE.md`) to clarify Python/OpenFOAM environment split.
- 2026-02-18: Validation: `python3 -m py_compile circularTiltingTank/adaptive_stop.py` passed.
- 2026-02-18: Validation: `make -n` for run/resume shows adaptive stop now executes with `../sloshing/bin/python` while OF commands remain prefixed by `of13`.

## Outcome
- Future generated cases will use repo venv for adaptive-stop CSV/PNG generation.
- OpenFOAM binaries continue to run through `of13` on OSCAR.
