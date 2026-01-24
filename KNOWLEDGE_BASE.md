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

### 4. Adaptive Stopping via max(U)
Adaptive stopping is implemented by monitoring `max(U)` from a `fieldMinMax` function object (`system/functions`) and requesting a graceful stop when `max(U)` stays below a threshold. The wrapper script `adaptive_stop.py` sets `stopAt writeNow` and then updates `endTime` to the final time after the run, so `runTimeModifiable yes` is required.

### 5. Alpha Contact Angle Compatibility
Some OpenFOAM builds lack the `constantAlphaContactAngle` patch field. When available, the compatible `contactAngle` BC is used on `walls`, and `main.py` includes a patch step to replace unsupported contact-angle types in existing cases before running.

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
- Duration: 10s
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
