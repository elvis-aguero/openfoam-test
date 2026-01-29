# Agent Knowledge Base: Circular Tilting Tank

> [!IMPORTANT]
> **AGENT PROTOCOL**: All agents MUST consult this file before running simulations. If you solve a new bug or learn a technical quirk, you are REQUIRE to update this file with a new entry in the "Technical Lessons Learned" section.

This file tracks technical discoveries, version-specific quirks, and project-wide context for the OpenFOAM tilting simulation.

## 📝 Technical Lessons Learned

### 1. OpenFOAM 13 (Foundation) - `setFieldsDict`
As of OpenFOAM 13, the `setFieldsDict` syntax has transitioned from list-based to dictionary-based. 

**Old Syntax (Deprecated/Unsupported):**
```cpp
defaultFieldValues ( volScalarFieldValue alpha.water 0 );
regions ( boxToCell { box (...); fieldValues ( volScalarFieldValue alpha.water 1 ); } );
```

**New Syntax (OpenFOAM 13):**
```cpp
defaultValues { alpha.water 0; }
zones { water { type box; box (...); values { alpha.water 1; } } }
```
*   `defaultFieldValues` -> `defaultValues`
*   `regions` -> `zones`
*   Selection (e.g., `boxToCell`) is replaced by a named selection block with a `type`.
*   `fieldValues` -> `values`

### 3. Static Tilt (Rotated Gravity)
The tilted-tank case is implemented by rotating gravity in `constant/g` about the +X axis by a fixed `tilt_deg`.
*   **Mechanism**: `g = (0, g*sin(theta), -g*cos(theta))` with `theta = tilt_deg`.
*   **Result**: The mesh remains static while the gravity vector encodes the tilt.

### 2. Python-to-OpenFOAM Dictionary Generation
When generating OpenFOAM dictionaries using Python `f-strings`, curly braces `{}` must be escaped (doubled) as `{{}}` to prevent Python from interpreting them as variable placeholders.

### 4. Adaptive Stopping via max(U) + Interface Stillness
Adaptive stopping is implemented by monitoring both:
- `max(|U|)` from `probes` (field `U`)
- Interface stillness from `probes` (field `alpha.water`) as `max(|delta alpha|)` between successive samples

When both stay below thresholds for a sliding time window, `adaptive_stop.py` requests a graceful stop via `stopAt writeNow`, then rewrites `endTime` to the final time written. This relies only on `probes` (portable across OpenFOAM builds) and requires `runTimeModifiable yes`.

### 5. Alpha Contact Angle Compatibility
Some OpenFOAM builds lack the `constantAlphaContactAngle` patch field. When available, the compatible `contactAngle` BC is used on `walls`, and `main.py` includes a patch step to replace unsupported contact-angle types in existing cases before running.

### 6. Mesh Quality Preflight (Gmsh MSH2)
Tiny elements (min edge length << target `lc`) can force extremely small adaptive `deltaT` (e.g., `1e-5`) and explode wall-clock time. The build menu runs a quick Gmsh preflight + MSH2 parser check and warns before generating cases if the mesh looks risky.

### 7. pRefPoint Must Be Inside the Mesh
`system/fvSolution` sets `pRefPoint` for pressure reference. If it lies outside the domain (common when templates are reused for smaller tanks), the pressure solve can become unstable and may crash with `sigFpe` inside GAMG/PCG. `main.py` patches `pRefPoint` per-case to `(0 0 H/2)`.

### 8. Mesh Preflight Is Non-Blocking
The build menu starts the Gmsh mesh-quality preflight in a background thread. It prints results when ready and does not delay case building; resource estimates use the analytic mesh estimate if the preflight hasn’t finished.

## 🏃 Current Workflow

**Interactive Manager** (`python3 main.py`):

1.  **Build Case Setups**: 
    - Configure parameters (H, D, mesh, geo, tilt_deg, duration, dt)
    - Support for parameter sweeps (MATLAB-style ranges: `0.1:0.05:0.2`)
    - Generates `case_*` folders with mesh, static motion, and initial conditions
    
2.  **Run Cases**:
    - Lists available cases with completion status
    - On Oscar: Submits to Slurm with dynamic resource allocation
    - Locally: Runs sequentially if OpenFOAM is installed
    
3.  **Postprocess**:
    - Generate MP4 videos (VTK + `imageio`)
    - Extract interface data (PyVista; VTP files + CSV summary)

**Default Parameters**:
- Tank: H=0.01m, D=0.0083m
- Mesh: 0.0005m characteristic length
- Tilt: 5.0 degrees about +X
- Duration: 5s
- Time step: 0.1s

## 🎯 Project Goals

## ✅ Simulation Verification Checklist
How to assess if the run was successful:

1.  **Courant Number (Co)**: Check the log for `max Courant Number`.
    -   **Good**: < 0.5 (Ideal for accuracy).
    -   **Acceptable**: < 0.9 (Standard stability limit).
    -   **Bad**: > 1.0 (Likely to crash or produce non-physical results).
2.  **Phase Boundedness**:
    -   `alpha.water` should never exceed 1.0 or drop below 0.0.
    -   Watch for "bounding alpha.water" warnings in the log.
3.  **Mass Conservation**:
    -   The total volume of water should not change significantly.
    -   Check `time step continuity errors` in the log.
4.  **Visual Inspection (ParaView)**:
    -   Is the interface sharp? (Smearing indicates excessive diffusion).
    -   Does the water surface move smoothly? (Jaggedness suggests mesh issues).
