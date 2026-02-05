#!/usr/bin/env python3
import csv
import glob
import os
import re
import sys

def _read_metric(path, key):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            for row in reader:
                if len(row) >= 2 and row[0].strip() == key:
                    try:
                        return float(row[1])
                    except ValueError:
                        return None
    except FileNotFoundError:
        return None
    return None

def _extract_tilt(case_dir):
    # case_H{H}_D{D}_{geo}_tilt_T{tilt}_m{mesh}[_suffix]
    m = re.search(r"_tilt_T([\d.]+)", case_dir)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None

def _load_case_params(case_dir):
    path = os.path.join(case_dir, "case_params.json")
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            import json
            return json.load(f)
    except Exception:
        return None

def _matches_sweep(params):
    if not params:
        return False
    # Sweep we just ran:
    # H=0.0083, D=0.0083, mesh=0.0001, dt=0.001, duration=4.0, mesher=snappy, n_cpus=8
    # tilt angles: 5, 10, 15, 20, 25
    return (
        float(params.get("H", -1)) == 0.0083
        and float(params.get("D", -1)) == 0.0083
        and float(params.get("mesh", -1)) == 0.0001
        and float(params.get("dt", -1)) == 0.001
        and float(params.get("duration", -1)) == 4.0
        and str(params.get("mesher", "")) == "snappy"
        and int(float(params.get("n_cpus", -1))) == 8
        and float(params.get("tilt_deg", -999)) in (5.0, 10.0, 15.0, 20.0, 25.0)
    )

def main():
    cases = sorted([d for d in glob.glob("case_*") if os.path.isdir(d)])
    rows = []
    for case in cases:
        params = _load_case_params(case)
        if not _matches_sweep(params):
            continue
        tilt = _extract_tilt(case)
        if tilt is None:
            continue
        summary = os.path.join(case, "postProcessing", "interface_compare", "comparison_summary.csv")
        l2 = _read_metric(summary, "l2_rms")
        if l2 is None:
            continue
        rows.append((tilt, l2, case))

    if not rows:
        print("No comparison_summary.csv with l2_rms found.")
        return 1

    # If multiple cases share the same tilt, keep all but sort by tilt then name.
    rows.sort(key=lambda r: (r[0], r[2]))

    out_csv = os.path.join("postProcessing", "l2_rms_vs_tilt.csv")
    os.makedirs("postProcessing", exist_ok=True)
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tilt_deg", "l2_rms", "case_dir"])
        for tilt, l2, case in rows:
            w.writerow([tilt, l2, case])

    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"Wrote {out_csv}. Matplotlib not available: {e}")
        return 0

    tilts = [r[0] for r in rows]
    l2s = [r[1] for r in rows]

    plt.figure(figsize=(6, 4))
    plt.plot(tilts, l2s, marker="o", linestyle="-")
    plt.xlabel("Tilt angle (deg)")
    plt.ylabel("L2 RMS")
    plt.title("Interface Error vs Tilt")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    out_png = os.path.join("postProcessing", "l2_rms_vs_tilt.png")
    plt.savefig(out_png, dpi=200)
    print(f"Wrote {out_csv} and {out_png}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
