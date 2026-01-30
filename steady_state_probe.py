#!/usr/bin/env python3
import argparse
import re
import sys
import numpy as np

TOKEN_RE = re.compile(r"\([^\)]*\)|[^\s]+")

def tokenize(line: str):
    return TOKEN_RE.findall(line.strip())

def parse_probe_file(path: str):
    times, rows = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            toks = tokenize(s)
            if len(toks) < 2:
                continue
            try:
                t = float(toks[0])
            except ValueError:
                continue
            times.append(t)
            rows.append(toks[1:])
    if len(times) < 20:
        raise RuntimeError(f"Too few data lines parsed ({len(times)}).")
    t = np.array(times, dtype=float)

    first = rows[0]
    is_vector = any(v.startswith("(") and v.endswith(")") for v in first)

    if not is_vector:
        data = np.array([[float(x) for x in r] for r in rows], dtype=float)  # (N, nProbes)
        return t, data, "scalar"

    vecs = []
    for r in rows:
        probe_vecs = []
        for tok in r:
            if not (tok.startswith("(") and tok.endswith(")")):
                raise RuntimeError("Vector probe file expected '(x y z)' tokens; found non-parenthesized token.")
            nums = tok[1:-1].split()
            if len(nums) != 3:
                raise RuntimeError(f"Vector token not length-3: {tok}")
            probe_vecs.append([float(nums[0]), float(nums[1]), float(nums[2])])
        vecs.append(probe_vecs)
    data = np.array(vecs, dtype=float)  # (N, nProbes, 3)
    return t, data, "vector"

def pick_series(data, kind: str, probe_idx: int, component: str):
    if kind == "scalar":
        if not (0 <= probe_idx < data.shape[1]):
            raise RuntimeError(f"--probe out of range [0, {data.shape[1]-1}]")
        return data[:, probe_idx]

    if not (0 <= probe_idx < data.shape[1]):
        raise RuntimeError(f"--probe out of range [0, {data.shape[1]-1}]")
    v = data[:, probe_idx, :]
    if component == "mag":
        return np.linalg.norm(v, axis=1)
    if component == "x":
        return v[:, 0]
    if component == "y":
        return v[:, 1]
    if component == "z":
        return v[:, 2]
    raise RuntimeError("--component must be one of: mag, x, y, z")

def robust_scale(y):
    # Robust amplitude scale: inter-quantile range
    q05, q50, q95 = np.quantile(y, [0.05, 0.5, 0.95])
    iqr = float(q95 - q05)
    fallback = max(float(abs(q95)), float(abs(q50)), 1e-12)
    return iqr if iqr > 1e-12 else fallback

def linreg_slope(t, y):
    # slope of y vs t (units: y/s)
    tt = t - t.mean()
    denom = float(np.dot(tt, tt))
    if denom <= 0:
        return 0.0, float(y.mean())
    m = float(np.dot(tt, y - y.mean()) / denom)
    b = float(y.mean() - m * t.mean())
    return m, b

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", help="Probe output file (scalar or vector).")
    ap.add_argument("--probe", type=int, default=0, help="Probe index (0-based).")
    ap.add_argument("--component", default="mag", help="For vector: mag|x|y|z (default mag).")
    ap.add_argument("--window_frac", type=float, default=0.2, help="Fraction of tail window to test (default 0.2).")
    ap.add_argument("--window_s", type=float, default=None, help="Optional tail window duration in seconds.")
    ap.add_argument("--drift_tol", type=float, default=0.01, help="Max fractional drift over window (default 1%).")
    ap.add_argument("--noise_tol", type=float, default=0.01, help="Max detrended fractional noise (default 1%).")
    args = ap.parse_args()

    t, data, kind = parse_probe_file(args.file)
    y = pick_series(data, kind, args.probe, args.component if kind == "vector" else "mag")

    # Basic time stats
    dt = np.diff(t)
    dt_med = float(np.median(dt))
    duration = float(t[-1] - t[0])

    # Define last window
    if args.window_s is not None:
        Tw = float(args.window_s)
    else:
        Tw = max(5 * dt_med, float(args.window_frac) * duration)
    t_start = float(t[-1] - Tw)
    mask = t >= t_start
    if mask.sum() < 20:
        # Ensure enough points
        mask = t >= t[-1] - min(duration, 50 * dt_med)
    tw = t[mask]
    yw = y[mask]

    # Trend in last window
    m, b = linreg_slope(tw, yw)  # y ≈ m t + b
    yw_fit = m * tw + b
    resid = yw - yw_fit

    # Robust normalization scale (global)
    s = robust_scale(y)

    # Metrics
    drift_frac = abs(m) * (tw[-1] - tw[0]) / s             # fractional change across window due to trend
    noise_frac = float(np.std(resid)) / s                  # detrended fluctuation level
    mean_last = float(np.mean(yw))
    mean_frac = abs(mean_last) / s                         # absolute level relative to scale

    # Two-window consistency (optional but helpful): compare last half vs previous half
    mid = int(len(tw) * 0.5)
    if mid >= 10:
        y1, y2 = yw[:mid], yw[mid:]
        mean_jump_frac = abs(float(np.mean(y2) - np.mean(y1))) / s
        std1, std2 = float(np.std(y1)), float(np.std(y2))
        std_ratio = (std2 + 1e-12) / (std1 + 1e-12)
    else:
        mean_jump_frac = float("nan")
        std_ratio = float("nan")

    print("=== Steady-State Probe Test (objective) ===")
    print(f"file: {args.file}")
    print(f"kind: {kind}  probe: {args.probe}" + (f"  component: {args.component}" if kind=="vector" else ""))
    print(f"samples: {len(t)}  duration: {duration:.6g} s  dt(median): {dt_med:.6g} s")
    print(f"tail window: {tw[0]:.6g} → {tw[-1]:.6g}  (Tw={tw[-1]-tw[0]:.6g} s, N={len(tw)})")
    print(f"robust scale S: {s:.6g}")
    print("\n--- Metrics (normalized) ---")
    print(f"drift_frac = |slope|*Tw/S: {drift_frac:.6g}")
    print(f"noise_frac = std(detrended)/S: {noise_frac:.6g}")
    print(f"mean_frac  = |mean_last|/S: {mean_frac:.6g}")
    if np.isfinite(mean_jump_frac):
        print(f"mean_jump_frac (2-half tail): {mean_jump_frac:.6g}")
        print(f"std_ratio (2nd/1st half tail): {std_ratio:.6g}")

    passed = (drift_frac <= args.drift_tol) and (noise_frac <= args.noise_tol)
    print("\n--- Decision ---")
    print(f"thresholds: drift_tol={args.drift_tol}  noise_tol={args.noise_tol}")
    print("RESULT: PASS (steady over tail window)" if passed else "RESULT: FAIL (not steady over tail window)")
    sys.exit(0 if passed else 1)

if __name__ == "__main__":
    main()
