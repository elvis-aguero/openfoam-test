#!/usr/bin/env python3
import argparse
import re
import sys
import numpy as np

TOKEN_RE = re.compile(r"\([^\)]*\)|[^\s]+")

def tokenize(line: str):
    return TOKEN_RE.findall(line.strip())

def parse_probe_file(path: str):
    times = []
    rows = []
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
            vals = toks[1:]
            times.append(t)
            rows.append(vals)
    if len(times) < 10:
        raise RuntimeError(f"Too few data lines parsed ({len(times)}).")

    times = np.array(times, dtype=float)

    # Determine if values are scalar-per-probe or vector-per-probe
    first = rows[0]
    is_vector = any(v.startswith("(") and v.endswith(")") for v in first)

    if not is_vector:
        # Scalars: each token is a float value for a probe
        data = np.array([[float(x) for x in r] for r in rows], dtype=float)
        return times, data, "scalar"

    # Vectors: each probe is a "(x y z)" token
    vecs = []
    for r in rows:
        probe_vecs = []
        for tok in r:
            if tok.startswith("(") and tok.endswith(")"):
                nums = tok[1:-1].split()
                if len(nums) != 3:
                    raise RuntimeError(f"Vector token not length-3: {tok}")
                probe_vecs.append([float(nums[0]), float(nums[1]), float(nums[2])])
            else:
                # Sometimes OpenFOAM prints vectors split; fail loudly
                raise RuntimeError(
                    "Detected vector-like file but found non-parenthesized token. "
                    "If this is U, ensure your probes output is in (x y z) format."
                )
        vecs.append(probe_vecs)
    data = np.array(vecs, dtype=float)  # shape: (N, nProbes, 3)
    return times, data, "vector"

def pick_series(data, probe: int, kind: str, component: str):
    if kind == "scalar":
        if probe < 0 or probe >= data.shape[1]:
            raise RuntimeError(f"--probe out of range [0, {data.shape[1]-1}]")
        return data[:, probe]

    # vector
    if probe < 0 or probe >= data.shape[1]:
        raise RuntimeError(f"--probe out of range [0, {data.shape[1]-1}]")

    v = data[:, probe, :]
    if component == "mag":
        return np.linalg.norm(v, axis=1)
    if component == "x":
        return v[:, 0]
    if component == "y":
        return v[:, 1]
    if component == "z":
        return v[:, 2]
    raise RuntimeError("--component must be one of: mag, x, y, z")

def uniform_resample(t, y):
    # Resample to uniform grid using median dt (robust to occasional adjustTimeStep)
    dt = np.diff(t)
    dt_med = float(np.median(dt))
    t0, t1 = float(t[0]), float(t[-1])
    n = int(np.floor((t1 - t0) / dt_med)) + 1
    tu = t0 + dt_med * np.arange(n)
    yu = np.interp(tu, t, y)
    return tu, yu, dt_med

def dominant_period_fft(yu, dt):
    # Detrend (remove mean and linear trend)
    x = np.arange(len(yu))
    A = np.vstack([x, np.ones_like(x)]).T
    m, b = np.linalg.lstsq(A, yu, rcond=None)[0]
    yd = yu - (m * x + b)

    # FFT
    n = len(yd)
    if n < 32:
        return None
    Y = np.fft.rfft(yd)
    freqs = np.fft.rfftfreq(n, d=dt)

    # Ignore DC
    amps = np.abs(Y)
    amps[0] = 0.0
    k = int(np.argmax(amps))
    f = float(freqs[k])
    if f <= 0:
        return None
    period = 1.0 / f

    # Peak sharpness metric (dominant amplitude vs median nonzero amplitude)
    nz = amps[amps > 0]
    sharp = float(amps[k] / (np.median(nz) + 1e-12)) if len(nz) else float("inf")
    return period, sharp

def dominant_period_autocorr(yu, dt):
    # Remove mean
    yd = yu - np.mean(yu)
    n = len(yd)
    if n < 64:
        return None
    ac = np.correlate(yd, yd, mode="full")[n-1:]
    if ac[0] <= 0:
        return None
    ac = ac / ac[0]

    # Find first strong local maximum after lag>=2
    # Limit search to lags up to n/2
    Lmax = n // 2
    best_lag = None
    best_val = -1.0
    for lag in range(2, Lmax-1):
        if ac[lag] > ac[lag-1] and ac[lag] > ac[lag+1]:
            if ac[lag] > best_val:
                best_val = float(ac[lag])
                best_lag = lag
    if best_lag is None:
        return None
    return best_lag * dt, best_val

def classify(period_s, dt_med):
    steps = period_s / dt_med
    nearest = int(np.round(steps))
    return steps, nearest

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", help="Probe output file, e.g., postProcessing/probes/0/p")
    ap.add_argument("--probe", type=int, default=0, help="Probe index (0-based)")
    ap.add_argument("--component", default="mag", help="For vector probes: mag|x|y|z (default mag)")
    ap.add_argument("--target_steps", type=float, default=3.0, help="Test against this many timesteps (default 3)")
    args = ap.parse_args()

    t, data, kind = parse_probe_file(args.file)
    y = pick_series(data, args.probe, kind, args.component if kind == "vector" else "mag")

    # Basic dt stats (before resampling)
    dt = np.diff(t)
    dt_med_raw = float(np.median(dt))
    dt_min = float(np.min(dt))
    dt_max = float(np.max(dt))

    tu, yu, dt_med = uniform_resample(t, y)

    fft_res = dominant_period_fft(yu, dt_med)
    ac_res = dominant_period_autocorr(yu, dt_med)

    print("=== Probe Oscillation Diagnostic ===")
    print(f"file: {args.file}")
    print(f"kind: {kind}   probe: {args.probe}" + (f"   component: {args.component}" if kind=="vector" else ""))
    print(f"samples(raw): {len(t)}   duration: {t[-1]-t[0]:.6g} s")
    print(f"dt(raw): median={dt_med_raw:.6g}  min={dt_min:.6g}  max={dt_max:.6g} s")
    print(f"dt(resample): {dt_med:.6g} s   samples(uniform): {len(tu)}")

    periods = []

    if fft_res is not None:
        p_fft, sharp = fft_res
        steps, nearest = classify(p_fft, dt_med)
        periods.append(("FFT", p_fft, steps, nearest, sharp))
    if ac_res is not None:
        p_ac, peak = ac_res
        steps, nearest = classify(p_ac, dt_med)
        periods.append(("AC", p_ac, steps, nearest, peak))

    if not periods:
        print("RESULT: Could not detect a stable dominant period (insufficient data or too noisy).")
        sys.exit(2)

    # Choose the method with stronger confidence metric
    periods.sort(key=lambda r: r[4], reverse=True)
    method, period_s, steps, nearest, metric = periods[0]

    print("\n--- Dominant period estimate ---")
    print(f"method: {method}   confidence_metric: {metric:.4g}")
    print(f"period: {period_s:.6g} s")
    print(f"period_in_timesteps (using median dt): {steps:.6g}")
    print(f"nearest_integer_steps: {nearest}")

    target = float(args.target_steps)
    # Objective flag: near target and very near an integer
    err_target = abs(steps - target)
    err_int = abs(steps - nearest)

    print("\n--- Objective test: 'period tied to ~3 timesteps?' ---")
    print(f"target_steps: {target}")
    print(f"|steps - target|: {err_target:.6g}")
    print(f"|steps - nearest_integer|: {err_int:.6g}")

    # Thresholds: tuned to be strict but practical for CFD noise
    # - within 0.25 steps of target
    # - within 0.15 steps of an integer
    # - plus confidence metric reasonably high
    tied = (err_target <= 0.25) and (err_int <= 0.15) and (metric >= 5.0)

    if tied:
        print("RESULT: CONSISTENT with a timestep-locked (~3-step) numerical limit cycle.")
        sys.exit(0)
    else:
        print("RESULT: NOT consistent with a clean 3-timestep-locked cycle (based on this probe/metric).")
        sys.exit(1)

if __name__ == "__main__":
    main()

