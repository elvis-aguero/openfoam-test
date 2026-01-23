# 🧭 Circular Tilting Tank - User Guide

Welcome! This guide will help you run the tilting simulation using the **interactive Tilting Manager**.

## 🚀 Quick Start
Run the manager:
```bash
python3 main.py
```

You will be prompted:
1.  **Are you on Oscar?** (Y/N) - This determines whether jobs are submitted to Slurm.
2.  **Main Menu**:
    *   **1) Build Case Setups**: Create one or more simulation cases.
    *   **2) Run Cases**: Run or submit selected cases.
    *   **3) Postprocess**: Compare analytical vs OpenFOAM interfaces.

---

## 🛠️ Building Cases (Sweeps)

When you select "Build Case Setups", you can:
1.  **Use defaults**: Just press Enter at the prompt to generate one default case.
2.  **Override a single parameter**: Enter the parameter name (e.g., `H`) and a new value.
3.  **Sweep over a parameter**: Enter a list like `0.1, 0.15, 0.2` or a MATLAB-style range like `0.1:0.05:0.2`.

**Zip vs Cartesian**:
*   If all sweep lists have the **same length**, they are **zipped** (paired 1-to-1).
*   If lengths **differ**, a **Cartesian product** is generated (all combinations). You will be asked to confirm.

---

## 🏃 Running Cases

When you select "Run Cases":
1.  The manager scans for `case_*` folders.
2.  It shows the status: cases that are complete will have `(DONE)` next to them.
3.  Enter the indices of cases you want to run (e.g., `1, 3-5, all`).
4.  On **Oscar**: Jobs are submitted to Slurm with smart resource allocation.
5.  **Locally**: If OpenFOAM is installed, simulations run sequentially.
6.  **Adaptive stop**: Runs stop early once the interface is steady (max velocity stays below a threshold for a window). Tune in `case_*/system/adaptiveStopDict` if needed.

---

## 🎬 Postprocessing

Select **3) Postprocess** from the main menu:

**Compare Interfaces (Analytical vs OpenFOAM)**:
- Generates the analytical interface and extracts the latest OpenFOAM interface (if present)
- Saves both and computes the L2 difference between them
- Output:
  - `case_*/postProcessing/interface_compare/analytical_interface.csv`
  - `case_*/postProcessing/interface_compare/openfoam_interface_t*.vtp` (and `.csv` point cloud)
  - `case_*/postProcessing/interface_compare/comparison_summary.csv`

---

## 📊 Viewing Results (ParaView)

1.  Open **ParaView**.
2.  Open `case.foam` inside your `case_...` folder.
3.  Click **Apply**, check `alpha.water`, and press **Play** (▶️).

---

## 🧹 Cleaning Up
```bash
rm -rf case_*
```

---

## 🔧 Tilt-Specific Notes

**Tilt parameter**:
- `tilt_deg`: static tilt angle about the +X axis (degrees).

**How tilt is implemented**:
- Gravity is rotated in `constant/g` to emulate a fixed tank tilt.
- Motion is static (no translation or rotation), so the mesh remains fixed.
