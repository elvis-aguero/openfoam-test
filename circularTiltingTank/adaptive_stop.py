#!/usr/bin/env python3
import argparse
import os
import re
import signal
import subprocess
import sys
import time

FUNC_NAME = "maxU"
CONTROL_DICT = os.path.join("system", "controlDict")
ADAPTIVE_DICT = os.path.join("system", "adaptiveStopDict")

DEFAULT_CONFIG = {
    "enabled": True,
    "maxU": 1e-3,
    "minTime": 1.0,
    "window": 1.0,
    "minSamples": 5,
    "checkInterval": 2.0,
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


def parse_samples(path):
    samples = []
    try:
        with open(path, "r") as handle:
            for line in handle:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                values = [float(v) for v in FLOAT_RE.findall(line)]
                if len(values) < 2:
                    continue
                samples.append((values[0], values[-1]))
    except FileNotFoundError:
        return samples
    except Exception:
        return samples
    return samples


def should_stop(samples, config):
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
    for _, max_u in window_samples:
        if max_u > float(config["maxU"]):
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


def find_max_u_dat():
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
                if name.endswith(".dat"):
                    candidates.append(os.path.join(root, name))
        if not candidates:
            continue
        for path in candidates:
            if "fieldMinMax" in os.path.basename(path):
                return path
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

    print("Adaptive stop enabled: watching max(U) for steady state.")
    proc = subprocess.Popen(cmd)
    stop_requested = False

    try:
        check_interval = max(0.2, float(config["checkInterval"]))
        while proc.poll() is None:
            data_path = find_max_u_dat()
            if data_path:
                samples = parse_samples(data_path)
                if not stop_requested and should_stop(samples, config):
                    print(
                        "Adaptive stop: steady state detected, requesting stop at next write."
                    )
                    update_control_dict(stop_at="writeNow")
                    stop_requested = True
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
