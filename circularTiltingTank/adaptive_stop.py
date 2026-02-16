#!/usr/bin/env python3
import argparse
import csv
import json
import math
import os
import re
import signal
import subprocess
import sys
import time

FUNC_NAME = "probesU"
CONTROL_DICT = os.path.join("system", "controlDict")
ADAPTIVE_DICT = os.path.join("system", "adaptiveStopDict")

DEFAULT_CONFIG = {
    "enabled": True,
    "normU": 5e-5,
    "maxDeltaAlpha": 1e-4,
    "minTime": 1.0,
    "window": 1.0,
    "minSamples": 5,
    "checkInterval": 2.0,
    "logInterval": 30.0,
    "csvEnabled": True,
    "csvPath": os.path.join("postProcessing", "adaptive_stop", "adaptive_metrics.csv"),
    "solverLogDir": os.path.join("postProcessing", "adaptive_stop"),
    "quietSolverOutput": True,
}

FLOAT_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")
TIME_LINE_RE = re.compile(r"^\s*Time\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)")
DT_LINE_RE = re.compile(r"^\s*deltaT\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)")
CO_LINE_RE = re.compile(
    r"^\s*Courant Number mean:\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)\s+max:\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)"
)
ICO_LINE_RE = re.compile(
    r"^\s*Interface Courant Number mean:\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)\s+max:\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)"
)


def parse_bool(value):
    value = value.strip().lower()
    if value in ("yes", "true", "on", "1"):
        return True
    if value in ("no", "false", "off", "0"):
        return False
    return None


def parse_value(raw):
    raw = raw.strip()
    bool_value = parse_bool(raw)
    if bool_value is not None:
        return bool_value
    try:
        if re.search(r"[.eE]", raw):
            return float(raw)
        return int(raw)
    except ValueError:
        return raw


def load_adaptive_config(path):
    config = DEFAULT_CONFIG.copy()
    if not os.path.exists(path):
        return config
    in_block = False
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                line = line.split("//", 1)[0].strip()
                if not line:
                    continue
                if not in_block:
                    if line.startswith("adaptiveStop"):
                        in_block = True
                    continue
                if "}" in line:
                    break
                match = re.match(r"([A-Za-z_][A-Za-z0-9_]*)\s+([^;]+);", line)
                if match:
                    key, value = match.groups()
                    config[key] = parse_value(value)
    except Exception:
        return config
    return config


def _parse_rows(path):
    rows = []
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                values = [float(v) for v in FLOAT_RE.findall(line)]
                if len(values) < 2:
                    continue
                rows.append((values[0], values[1:]))
    except FileNotFoundError:
        return rows
    except Exception:
        return rows
    return rows


def _series_norm_u(rows):
    out = []
    for t, data in rows:
        if not data:
            out.append((t, 0.0))
            continue
        sum_sq = sum(v * v for v in data)
        norm = sum_sq**0.5
        # Use norm(U) / numel(U) as requested
        val = norm / len(data)
        out.append((t, val))
    return out


def _series_max_delta_scalar(rows):
    out = []
    prev = None
    for t, data in rows:
        if prev is None:
            prev = data
            continue
        m = 0.0
        for a, b in zip(data, prev):
            d = abs(a - b)
            if d > m:
                m = d
        out.append((t, m))
        prev = data
    return out


def _parse_probe_locations(functions_path=os.path.join("system", "functions")):
    points = []
    if not os.path.exists(functions_path):
        return points
    in_locations = False
    try:
        with open(functions_path, "r", encoding="utf-8", errors="ignore") as handle:
            for raw in handle:
                line = raw.split("//", 1)[0].strip()
                if not line:
                    continue
                if not in_locations:
                    if line.startswith("probeLocations"):
                        in_locations = True
                    continue
                if line.startswith(");"):
                    break
                if "(" not in line or ")" not in line:
                    continue
                vals = [float(v) for v in FLOAT_RE.findall(line)]
                if len(vals) >= 3:
                    points.append((vals[0], vals[1], vals[2]))
    except Exception:
        return []
    return points


def _estimate_interface_minmax(alpha_values, probe_points):
    if not alpha_values or not probe_points:
        return None, None
    n = min(len(alpha_values), len(probe_points))
    if n == 0:
        return None, None

    columns = {}
    for i in range(n):
        x, y, z = probe_points[i]
        key = (round(x, 12), round(y, 12))
        columns.setdefault(key, []).append((z, alpha_values[i]))

    heights = []
    target = 0.5
    for samples in columns.values():
        if len(samples) < 2:
            continue
        samples.sort(key=lambda p: p[0])
        z_iface = None
        for (z0, a0), (z1, a1) in zip(samples[:-1], samples[1:]):
            lo = a0 if a0 < a1 else a1
            hi = a1 if a1 > a0 else a0
            if lo <= target <= hi:
                if abs(a1 - a0) < 1e-14:
                    z_iface = z0
                else:
                    z_iface = z0 + (target - a0) * (z1 - z0) / (a1 - a0)
                break
        if z_iface is None:
            alphas = [a for _, a in samples]
            z_min = samples[0][0]
            z_max = samples[-1][0]
            if min(alphas) >= target:
                z_iface = z_max
            elif max(alphas) <= target:
                z_iface = z_min
            else:
                z_iface = min(samples, key=lambda p: abs(p[1] - target))[0]
        heights.append(z_iface)

    if not heights:
        return None, None
    return min(heights), max(heights)


def _format_csv_value(value):
    if value is None:
        return "nan"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "nan"
    return f"{float(value):.12g}"


def _read_still_interface_height():
    path = "case_params.json"
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            params = json.load(handle)
        h = params.get("H")
        if h is None:
            return None
        return 0.5 * float(h)
    except Exception:
        return None


def _append_metrics_row(csv_path, row):
    directory = os.path.dirname(csv_path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    new_file = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        if new_file:
            writer.writerow(
                [
                    "timestamp",
                    "normU",
                    "maxDeltaAlpha",
                    "interfaceHeightMaxRel",
                    "interfaceHeightMinRel",
                ]
            )
        writer.writerow(
            [
                _format_csv_value(row.get("timestamp")),
                _format_csv_value(row.get("normU")),
                _format_csv_value(row.get("maxDeltaAlpha")),
                _format_csv_value(row.get("interfaceHeightMaxRel")),
                _format_csv_value(row.get("interfaceHeightMinRel")),
            ]
        )


def _to_float_or_nan(value):
    try:
        if value is None:
            return float("nan")
        text = str(value).strip()
        if not text or text.lower() == "nan":
            return float("nan")
        return float(text)
    except Exception:
        return float("nan")


def _write_convergence_summary_figure(csv_path, fig_path):
    """
    Create a convergence summary figure with one subplot per metric column.
    Returns True on success, False otherwise.
    """
    if not os.path.exists(csv_path):
        return False

    try:
        with open(csv_path, "r", encoding="utf-8", errors="ignore") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            rows = list(reader)
    except Exception:
        return False

    if not rows or not fieldnames:
        return False

    x_key = "timestamp" if "timestamp" in fieldnames else fieldnames[0]
    metric_keys = [k for k in fieldnames if k != x_key]
    if not metric_keys:
        return False

    x_values = []
    for i, row in enumerate(rows):
        x_val = _to_float_or_nan(row.get(x_key))
        if not math.isfinite(x_val):
            x_val = float(i)
        x_values.append(x_val)

    # Matplotlib is optional; plotting should never break adaptive stopping.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False

    n_metrics = len(metric_keys)
    n_cols = 1 if n_metrics == 1 else 2
    n_rows = int(math.ceil(n_metrics / n_cols))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(6.2 * n_cols, 2.8 * n_rows),
        squeeze=False,
        sharex=True,
        constrained_layout=True,
    )
    axes_flat = [ax for row_axes in axes for ax in row_axes]

    for i, key in enumerate(metric_keys):
        ax = axes_flat[i]
        y_values = [_to_float_or_nan(row.get(key)) for row in rows]
        xs = []
        ys = []
        for x, y in zip(x_values, y_values):
            if math.isfinite(x) and math.isfinite(y):
                xs.append(x)
                ys.append(y)

        if xs:
            ax.plot(xs, ys, "-", linewidth=1.5, color="#1f77b4")
            ax.scatter(xs, ys, s=8, color="#1f77b4", alpha=0.6)
        else:
            ax.text(0.5, 0.5, "No finite data", ha="center", va="center", transform=ax.transAxes)

        ax.set_ylabel(key)
        ax.grid(True, alpha=0.25)

    for j in range(n_metrics, len(axes_flat)):
        axes_flat[j].axis("off")

    for ax in axes[-1]:
        if ax.get_visible():
            ax.set_xlabel(x_key)

    fig.suptitle("Adaptive Stop Convergence Summary")
    out_dir = os.path.dirname(fig_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    try:
        fig.savefig(fig_path, dpi=150)
        return True
    finally:
        plt.close(fig)


def _scan_solver_log(log_path, state, metrics):
    """
    Parse newly appended solver log text and cache latest runtime metrics.
    """
    if not log_path or not os.path.exists(log_path):
        return
    try:
        with open(log_path, "r", encoding="utf-8", errors="ignore") as handle:
            handle.seek(state.get("pos", 0))
            for line in handle:
                m = TIME_LINE_RE.match(line)
                if m:
                    metrics["solver_time"] = float(m.group(1))
                    continue
                m = DT_LINE_RE.match(line)
                if m:
                    metrics["delta_t"] = float(m.group(1))
                    continue
                m = CO_LINE_RE.match(line)
                if m:
                    metrics["co_mean"] = float(m.group(1))
                    metrics["co_max"] = float(m.group(2))
                    continue
                m = ICO_LINE_RE.match(line)
                if m:
                    metrics["ico_mean"] = float(m.group(1))
                    metrics["ico_max"] = float(m.group(2))
            state["pos"] = handle.tell()
    except Exception:
        return


def should_stop_metric(samples, config, threshold_key):
    if not samples:
        return False
    last_time = samples[-1][0]
    if last_time < float(config["minTime"]):
        return False
    window = max(0.0, float(config["window"]))
    min_samples = max(1, int(config["minSamples"]))
    window_samples = [
        sample
        for sample in samples
        if sample[0] >= last_time - window and sample[0] >= float(config["minTime"])
    ]
    if len(window_samples) < min_samples:
        return False
    threshold = float(config.get(threshold_key, 0.0))
    for _, value in window_samples:
        if value > threshold:
            return False
    return True


def update_control_dict(stop_at=None, end_time=None):
    if not os.path.exists(CONTROL_DICT):
        return False
    with open(CONTROL_DICT, "r", encoding="utf-8", errors="ignore") as handle:
        content = handle.read()
    changed = False
    if stop_at is not None:
        content, count = re.subn(
            r"(^\s*stopAt\s+)[^;]+;",
            r"\1" + stop_at + ";",
            content,
            flags=re.M,
        )
        if count == 0:
            content += f"\nstopAt        {stop_at};\n"
        changed = True
    if end_time is not None:
        content, count = re.subn(
            r"(^\s*endTime\s+)[^;]+;",
            r"\1" + end_time + ";",
            content,
            flags=re.M,
        )
        if count == 0:
            content += f"\nendTime       {end_time};\n"
        changed = True
    if changed:
        with open(CONTROL_DICT, "w") as handle:
            handle.write(content)
    return changed


def _find_probe_dat(field_name):
    search_roots = ["."]
    if os.path.isdir("processor0"):
        search_roots.append("processor0")
    for name in os.listdir("."):
        if name.startswith("processor") and os.path.isdir(name):
            search_roots.append(name)
    seen = set()
    for base in search_roots:
        if base in seen:
            continue
        seen.add(base)
        func_dir = os.path.join(base, "postProcessing", FUNC_NAME)
        if not os.path.isdir(func_dir):
            continue
        candidates = []
        for root, _, files in os.walk(func_dir):
            for name in files:
                # probes writes per-field files like ".../U"
                if name == field_name:
                    candidates.append(os.path.join(root, name))
        if not candidates:
            continue
        candidates.sort(key=lambda p: os.path.getmtime(p), reverse=True)
        return candidates[0]
    return None


def find_latest_time():
    times = []
    search_roots = ["."]
    for name in os.listdir("."):
        if name.startswith("processor") and os.path.isdir(name):
            search_roots.append(name)
    for base in search_roots:
        try:
            entries = os.listdir(base)
        except OSError:
            continue
        for name in entries:
            path = os.path.join(base, name)
            if not os.path.isdir(path):
                continue
            if name in ("constant", "system"):
                continue
            try:
                times.append(float(name))
            except ValueError:
                continue
    if not times:
        return None
    return max(times)


def format_time(value):
    if value is None:
        return None
    if abs(value - int(value)) < 1e-9:
        return str(int(value))
    return f"{value:.6g}"


def build_command(args):
    if args.parallel and args.np > 1:
        return ["mpirun", "-np", str(args.np), "foamRun", "-parallel"]
    return ["foamRun"]


def main():
    parser = argparse.ArgumentParser(
        description="Run foamRun with adaptive steady-state stopping."
    )
    parser.add_argument("--parallel", action="store_true", help="Run foamRun in parallel.")
    parser.add_argument("--np", type=int, default=1, help="Number of MPI ranks.")
    args = parser.parse_args()

    config = load_adaptive_config(ADAPTIVE_DICT)
    env_toggle = os.environ.get("ADAPTIVE_STOP")
    if env_toggle is not None:
        if str(env_toggle).strip().lower() in ("0", "no", "false", "off"):
            config["enabled"] = False

    cmd = build_command(args)
    if not config.get("enabled", True):
        return subprocess.run(cmd).returncode

    print("Adaptive stop enabled: watching norm(U)/numel(U) and interface stillness (alpha.water).")
    probe_points = _parse_probe_locations()
    still_h_ref = _read_still_interface_height()
    csv_enabled = bool(config.get("csvEnabled", True))
    csv_path = str(config.get("csvPath", DEFAULT_CONFIG["csvPath"]))
    quiet_solver_output = bool(config.get("quietSolverOutput", True))
    solver_log_dir = str(config.get("solverLogDir", DEFAULT_CONFIG["solverLogDir"]))
    solver_stdout = None
    solver_stderr = None
    popen_kwargs = {}
    if quiet_solver_output:
        os.makedirs(solver_log_dir, exist_ok=True)
        stdout_path = os.path.join(solver_log_dir, "foamRun.stdout.log")
        stderr_path = os.path.join(solver_log_dir, "foamRun.stderr.log")
        solver_stdout = open(stdout_path, "a", encoding="utf-8", buffering=1)
        solver_stderr = open(stderr_path, "a", encoding="utf-8", buffering=1)
        popen_kwargs["stdout"] = solver_stdout
        popen_kwargs["stderr"] = solver_stderr
        print(f"Adaptive stop: solver stdout -> {stdout_path}")
        print(f"Adaptive stop: solver stderr -> {stderr_path}")

    proc = subprocess.Popen(cmd, **popen_kwargs)
    stop_requested = False
    last_log = 0.0
    last_csv_time = None
    solver_log_state = {"pos": 0}
    solver_metrics = {}

    try:
        check_interval = max(0.2, float(config["checkInterval"]))
        log_interval = max(1.0, float(config.get("logInterval", 30.0)))
        while proc.poll() is None:
            u_path = _find_probe_dat("U")
            a_path = _find_probe_dat("alpha.water")

            u_ok = False
            a_ok = False
            latest_t = None
            latest_u = None
            latest_da = None
            latest_alpha_values = None
            if u_path:
                u_rows = _parse_rows(u_path)
                u_samples = _series_norm_u(u_rows)
                u_ok = should_stop_metric(u_samples, config, "normU")
                if u_samples:
                    latest_t = u_samples[-1][0]
                    latest_u = u_samples[-1][1]
            if a_path:
                a_rows = _parse_rows(a_path)
                a_samples = _series_max_delta_scalar(a_rows)
                a_ok = should_stop_metric(a_samples, config, "maxDeltaAlpha")
                if a_rows:
                    latest_alpha_values = a_rows[-1][1]
                    if latest_t is None:
                        latest_t = a_rows[-1][0]
                if a_samples:
                    latest_da = a_samples[-1][1]
                    if latest_t is None:
                        latest_t = a_samples[-1][0]

            interface_hmin = None
            interface_hmax = None
            if latest_alpha_values is not None and probe_points:
                interface_hmin, interface_hmax = _estimate_interface_minmax(
                    latest_alpha_values, probe_points
                )
                if still_h_ref is None and interface_hmin is not None and interface_hmax is not None:
                    # Fallback baseline when H is unavailable: first observed midpoint.
                    still_h_ref = 0.5 * (interface_hmin + interface_hmax)

            interface_hmin_rel = None
            interface_hmax_rel = None
            if still_h_ref is not None:
                if interface_hmin is not None:
                    interface_hmin_rel = interface_hmin - still_h_ref
                if interface_hmax is not None:
                    interface_hmax_rel = interface_hmax - still_h_ref

            if csv_enabled and latest_t is not None:
                if last_csv_time is None or abs(latest_t - last_csv_time) > 1e-12:
                    _append_metrics_row(
                        csv_path,
                        {
                            "timestamp": latest_t,
                            "normU": latest_u,
                            "maxDeltaAlpha": latest_da,
                            "interfaceHeightMaxRel": interface_hmax_rel,
                            "interfaceHeightMinRel": interface_hmin_rel,
                        },
                    )
                    last_csv_time = latest_t

            if not stop_requested and u_ok and a_ok:
                print("Adaptive stop: steady state detected, requesting stop at next write.")
                update_control_dict(stop_at="writeNow")
                stop_requested = True

            if quiet_solver_output:
                _scan_solver_log(stdout_path, solver_log_state, solver_metrics)

            now = time.time()
            if latest_u is not None and (now - last_log) >= log_interval:
                if latest_da is None:
                    da_str = "n/a"
                else:
                    da_str = f"{latest_da:.3g}"
                if interface_hmin_rel is None or interface_hmax_rel is None:
                    h_str = "hRel=[n/a, n/a]"
                else:
                    h_str = f"hRel=[{interface_hmin_rel:.4g}, {interface_hmax_rel:.4g}]"
                t_str = f"{latest_t:.6g}" if latest_t is not None else "n/a"
                dt_str = (
                    f"{solver_metrics['delta_t']:.3g}" if "delta_t" in solver_metrics else "n/a"
                )
                co_max_str = (
                    f"{solver_metrics['co_max']:.3g}" if "co_max" in solver_metrics else "n/a"
                )
                ico_max_str = (
                    f"{solver_metrics['ico_max']:.3g}" if "ico_max" in solver_metrics else "n/a"
                )
                print(
                    f"[adaptive_stop] t={t_str} dt={dt_str} Co_max={co_max_str} iCo_max={ico_max_str} "
                    f"normU={latest_u:.3g} maxDeltaAlpha={da_str} {h_str}",
                    flush=True,
                )
                last_log = now
            time.sleep(check_interval)
    except KeyboardInterrupt:
        proc.send_signal(signal.SIGINT)
    finally:
        return_code = proc.wait()
        if solver_stdout is not None:
            solver_stdout.close()
        if solver_stderr is not None:
            solver_stderr.close()

    if stop_requested:
        final_time = find_latest_time()
        final_time_str = format_time(final_time)
        if final_time_str:
            update_control_dict(stop_at="endTime", end_time=final_time_str)
        else:
            update_control_dict(stop_at="endTime")

    if csv_enabled:
        fig_dir = os.path.dirname(csv_path) or "."
        fig_path = os.path.join(fig_dir, "convergence_summary.png")
        if _write_convergence_summary_figure(csv_path, fig_path):
            print(f"Adaptive stop: convergence figure -> {fig_path}")

    return return_code


if __name__ == "__main__":
    sys.exit(main())
