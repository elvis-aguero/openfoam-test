#!/usr/bin/env bash
#SBATCH -J case_H0.0083_D0.0083_flat_tilt_T15.0_m0.0005
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH -o case_H0.0083_D0.0083_flat_tilt_T15.0_m0.0005/slurm/slurm.%j.out
#SBATCH -e case_H0.0083_D0.0083_flat_tilt_T15.0_m0.0005/slurm/slurm.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=elvis_vera@brown.edu

set -euo pipefail
export OMP_NUM_THREADS=1

echo 'Case: case_H0.0083_D0.0083_flat_tilt_T15.0_m0.0005'
# Check if we are resuming (parallel processors or serial time folders existed)
cd case_H0.0083_D0.0083_flat_tilt_T15.0_m0.0005
if [ -d 'processor0' ] || (ls -d [0-9]* 2>/dev/null | grep -v '^0$' | grep -q .) ; then
    echo 'Found existing progress. Resuming simulation...'
    make resume OSCAR=1 N_CPUS=1 ADAPTIVE_STOP=1
else
    echo 'Starting fresh simulation...'
    make run OSCAR=1 N_CPUS=1 ADAPTIVE_STOP=1
fi
echo 'End: $(date)'