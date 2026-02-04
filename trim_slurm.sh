#!/bin/bash
set +H
for file in case_*/slurm/slurm.*.out; do
  if [ -f "$file" ] && [ $(stat -c%s "$file") -gt 100000000 ]; then
    echo "Processing $file..."
    
    # Create trimmed version with:
    # 1. First 1000 lines (startup info)
    # 2. All lines with "error", "Error", "ERROR", "warning", "Warning", "WARNING"
    # 3. Last 2000 lines (end of run, convergence info)
    # 4. Any line with "Time = " (timestep info)
    
    head -n 1000 "$file" > "${file}.trimmed"
    grep -i -E "error|warning|fatal|abort" "$file" >> "${file}.trimmed"
    grep "Time = " "$file" | tail -n 1000 >> "${file}.trimmed"
    tail -n 2000 "$file" >> "${file}.trimmed"
    
    # Remove duplicates
    awk '!seen[$0]++' "${file}.trimmed" > "${file}.trimmed2"
    mv "${file}.trimmed2" "${file}.trimmed"
    
    # Replace if smaller
    if [ $(stat -c%s "${file}.trimmed") -lt $(stat -c%s "$file") ]; then
      mv "${file}.trimmed" "$file"
    else
      rm -f "${file}.trimmed"
    fi
  fi
done
