#!/usr/bin/env python3
import argparse
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
}

FLOAT_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")


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
        with open(path, "r") as handle:
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
        with open(path, "r") as handle:
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
    with open(CONTROL_DICT, "r") as handle:
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
    proc = subprocess.Popen(cmd)
    stop_requested = False
    last_log = 0.0

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
                if a_samples:
                    latest_da = a_samples[-1][1]
                    if latest_t is None:
                        latest_t = a_samples[-1][0]

            if not stop_requested and u_ok and a_ok:
                print("Adaptive stop: steady state detected, requesting stop at next write.")
                update_control_dict(stop_at="writeNow")
                stop_requested = True

            now = time.time()
            if latest_u is not None and (now - last_log) >= log_interval:
                if latest_da is None:
                    da_str = "n/a"
                else:
                    da_str = f"{latest_da:.3g}"
                t_str = f"{latest_t:.6g}" if latest_t is not None else "n/a"
                print(f"[adaptive_stop] t={t_str} normU={latest_u:.3g} maxDeltaAlpha={da_str}", flush=True)
                last_log = now
            time.sleep(check_interval)
    except KeyboardInterrupt:
        proc.send_signal(signal.SIGINT)
    finally:
        return_code = proc.wait()

    if stop_requested:
        final_time = find_latest_time()
        final_time_str = format_time(final_time)
        if final_time_str:
            update_control_dict(stop_at="endTime", end_time=final_time_str)
        else:
            update_control_dict(stop_at="endTime")

    return return_code


if __name__ == "__main__":
    sys.exit(main())
