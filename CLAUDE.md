# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an OpenFOAM-based parametric study framework for simulating water sloshing in a circular tilting tank. The primary physics: multiphase VOF (`interDyMFoam`) tracking the water-air interface under prescribed circular tank motion, validated against analytical Young-Laplace capillary solutions. The target HPC platform is Oscar (Brown University) via Slurm.

## Key Commands

### Interactive Simulation Manager
```bash
python3 main.py                          # Launch interactive menu (build/run/post-process)
SLOSHING_SKIP_DEPS=1 python3 main.py    # Skip dependency install (venv already ready)
SLOSHING_NO_VENV=1 python3 main.py      # Skip venv entirely
```

### Headless Post-processing
```bash
python3 main.py --headless --case <case_dir> --action compare
python3 main.py --headless --case <case_dir> --action video
python3 main.py --headless --case <case_dir> --action extract [--vtp-mode {none|latest|all}]
```

### Makefile (inside `circularTiltingTank/` or a case directory)
```bash
make all N_CPUS=4
make all H=0.008 D=0.008 MESH_SIZE=0.001 TILT_DEG=5.0 DURATION=5.0
make mesh / make motion / make case / make run / make clean
```

### Analysis Scripts
```bash
python3 mesh_quality.py cylinder.msh --out postProcessing/mesh_quality.json
python3 plot_l2_vs_tilt.py --x T --filter H=0.0083
python3 yl_nonlin.py          # Young-Laplace nonlinear solver
python3 steady_state_probe.py
python3 analyze_probe_cycle.py
```

## Architecture

### `main.py` (3647 lines) — Central Orchestrator
Three interactive menus:
1. **Build Cases**: Parameterizes and generates OpenFOAM case directories from the template. Supports parameter sweeps over H, D, mesh size, tilt angle (T), mass ratio (m), and geometry variant (flat/cap).
2. **Run Cases**: Local execution or Slurm job submission, with adaptive stopping daemon.
3. **Post-processing**: Interface comparison vs. Young-Laplace, MP4 animation, VTP export.

### `circularTiltingTank/` — OpenFOAM Template
All generated `case_*` directories are cloned from this template. Key generation scripts inside:
- `generate_mesh.py` — Gmsh `.geo` file for cylinder geometry
- `generate_motion.py` — 6DoF motion table (circular orbit, Smootherstep soft-start ramp)
- `generate_tilt.py` — Rotates gravity vector by tilt angle θ
- `update_setFields.py` — Sets water level at H/2 in `setFieldsDict`
- `adaptive_stop.py` — Real-time convergence monitoring (writes `adaptive_metrics.csv`)

### `case_*` Directories — Simulation Instances
Named `case_H{H}_D{D}_{geoType}_tilt_T{theta}_m{mass}_R{replica}`. Each contains:
- `case_params.json` — All build parameters (enables reproducibility and filtering)
- `postProcessing/adaptive_stop/adaptive_metrics.csv` — Convergence history
- `slurm/` — HPC job scripts and logs

## Known Pitfalls & Constraints

- **OpenFOAM 13 syntax** is required throughout — notably `setFieldsDict` syntax changed from earlier versions.
- **`pRefPoint`** in `fvSolution` must lie inside the mesh domain, or the simulation crashes.
- **MSH2 format** is critical for mesh quality checks (`mesh_quality.py` parses MSH2). Tiny mesh elements force tiny adaptive timesteps, dramatically increasing runtime.
- **`maxCo` = 0.5, `maxAlphaCo` = 0.5** — conservative adaptive timestepping for VOF interface stability.

## Dependencies

Python: `numpy`, `scipy`, `matplotlib`, `pyvista`, `imageio`, `imageio-ffmpeg`, `h5py`
System: OpenFOAM 13, Gmsh, Slurm (on Oscar)

## Multi-Agent Coordination

This repo uses `.context/` as agent-only working memory. See `.context/PROTOCOL.md` for rules:
- One active owner per task at a time (lease-based)
- All non-trivial actions must update the task file (`tasks/<task_id>.md`)
- `.context/INBOX.md` for inter-agent messages; append-only (never delete others' notes)
- Do not quote `.context/` contents to the user unless explicitly asked
