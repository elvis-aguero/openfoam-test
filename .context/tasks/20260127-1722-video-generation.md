# 20260127-1722-video-generation Implement Headless Video Generation

## Owner + Lease
- owner_session: Google/Antigravity/2026-01-27T17:22
- lease_expires: 2026-01-27T18:22

## Goal / Acceptance Criteria
- Implement a new menu option in `main.py` (or extend the Postprocess menu) to generate simulation videos.
- capabilities:
  - Work in a headless environment (Oscar CCV).
  - Generate a lateral view (side view) of the simulation.
  - Visualize `alpha.water` (water interface).
  - Stop/cover the duration until steady state (or end of simulation).
- Output: MP4 or GIF file in the case directory.

## Constraints / Non-goals
- Must run without X11 forwarding if possible (offscreen rendering).
- Use `pyvista` if possible as it is already in the project context.
- Minimal new dependencies if possible.

## Repo Touchpoints
- `main.py`: Main entry point to add menu option.
- `postprocessing.py`: (New/Existing?) Module to handle video generation.

## Plan
1.  Analyze `main.py` to identify where to add the new post-processing option.
2.  Check for existing post-processing code (e.g. `postprocessing/` folder or scripts).
3.  Create a script/function to:
    - Load OpenFOAM case data (VTK/VTP or native).
    - Set up a PyVista plotter for offscreen rendering.
    - Camera setup: Lateral view (XZ or YZ plane depending on orientation).
    - Iterate through time steps.
    - Write frames to video.
4.  Integrate into `main.py`.
5.  Verify dependencies (pyvista, imageio, ffmpeg).

## Work Log
- 2026-01-27 17:22: Task started.
- 2026-01-27 17:25: Analyzed `main.py` and `requirements.txt`. Discovered existing but unused `generate_video` using raw VTK.
- 2026-01-27 17:35: Decided to implement `generate_lateral_video` using `pyvista` for cleaner code and better headless support (via offscreen plotting).
- 2026-01-27 17:40: Implemented `generate_lateral_video` in `main.py` which creates side-view videos.
- 2026-01-27 17:40: Updated `menu_postprocess` to include the new option.
- 2026-01-27 17:40: Updated `run_postprocess_oscar` logic (arg passing) and `argparse` to support `--action video`.
- 2026-01-27 17:40: Updated `USER_GUIDE.md` to document the new feature.
- 2026-01-27 17:42: Fixed markdown formatting in User Guide.

## Messages
None.

## Handoff
- The "Generate Lateral Video" feature is fully implemented in `main.py`.
- It relies on `pyvista` and `imageio`.
- It supports headless rendering (specifically sets `off_screen=True`).
- On Oscar, it uses the existing `xvfb-run` wrapper in the Slurm script which is perfect.
- The video is saved as `postProcessing/video_lateral.mp4` in the case folder.
