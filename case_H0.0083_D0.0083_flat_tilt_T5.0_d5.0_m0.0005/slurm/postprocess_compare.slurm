#!/usr/bin/env bash
#SBATCH -J post_compare_case_H0.0083_D0.0083_flat_tilt_T5.0_d5.0_m0.0005
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH -o case_H0.0083_D0.0083_flat_tilt_T5.0_d5.0_m0.0005/slurm/slurm_postprocessing.log
#SBATCH --open-mode=append

set -euo pipefail

# --- Load Consistent Python Module ---
module load python/3.13

# --- Activate Shared Environment ---
VENV_DIR="sloshing"
if [ ! -d "$VENV_DIR" ]; then
  echo "📦 Venv not found. Creating $VENV_DIR on compute node..."
  python3 -m venv $VENV_DIR
  source $VENV_DIR/bin/activate
  pip install --upgrade pip
  pip install -r requirements.txt
else
  source $VENV_DIR/bin/activate
fi
# -----------------------------------

echo '------------------------------------------------------------'
echo 'Action: compare | Case: case_H0.0083_D0.0083_flat_tilt_T5.0_d5.0_m0.0005'
echo 'Date: $(date)'
echo 'Python: $(which python)'
export SLOSHING_OFFSCREEN=1
export VTK_DEFAULT_RENDER_WINDOW_OFFSCREEN=1
cd case_H0.0083_D0.0083_flat_tilt_T5.0_d5.0_m0.0005
if command -v xvfb-run >/dev/null 2>&1; then
  xvfb-run -s "-screen 0 1280x720x24" python ../main.py --headless --case . --action compare
else
  python ../main.py --headless --case . --action compare
fi
echo 'End: $(date)'
echo '------------------------------------------------------------'
