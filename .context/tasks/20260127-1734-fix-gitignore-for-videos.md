# 20260127-1734-fix-gitignore-for-videos Fix Gitignore for Video Tracking

## Owner + Lease
- owner_session: Google/Antigravity/2026-01-27T17:34
- lease_expires: 2026-01-27T18:34

## Goal / Acceptance Criteria
- Ensure `.mp4` videos generated in `postProcessing` directories are not ignored by Git.

## Constraints / Non-goals
- Do not track large raw mesh/field data (which is correctly ignored).

## Repo Touchpoints
- `.gitignore`: Needs un-ignore rule for `.mp4`.

## Plan
1.  Identify the specific ignore rules in `.gitignore` blocking `.mp4` files.
2.  Add an exception (negative ignore) for `**/*.mp4`.
3.  Verify with `git check-ignore`.

## Work Log
- 2026-01-27 17:34: Task started.
- 2026-01-27 17:35: Analyzed `.gitignore`. Found that `postProcessing/` and `case_*/**` were broadly ignoring the artifacts.
- 2026-01-27 17:38: Updated `.gitignore` to un-ignore `**/*.mp4`, `**/*.csv`, and specific paths while ensuring heavy simulation data (numerical fields, processor folders) remains ignored.
- 2026-01-27 17:39: Verified that `video_lateral.mp4` is no longer ignored using `git check-ignore`.

## Messages
None.

## Handoff
- `.mp4` and `.csv` files are now un-ignored globally to allow tracking of post-processing artifacts (videos and summaries).
- Large OpenFOAM data directories (e.g., `processor*`, numerical time folders like `1.5/`, and `polyMesh`) are still ignored to prevent repo bloating.
