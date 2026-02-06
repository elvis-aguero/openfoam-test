#!/usr/bin/env python3
import argparse
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

def _coerce_value(raw):
    try:
        if "." in raw or "e" in raw.lower():
            return float(raw)
        return int(raw)
    except Exception:
        return raw

def _match_filter(params, key, raw_value):
    if params is None:
        return False
    if key not in params:
        return False
    desired = _coerce_value(raw_value)
    actual = params.get(key)
    if isinstance(desired, (int, float)):
        try:
            return float(actual) == float(desired)
        except Exception:
            return False
    return str(actual) == str(desired)

def _apply_filters(params, filters):
    for key, raw_value in filters:
        if not _match_filter(params, key, raw_value):
            return False
    return True

def _build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Plot L2 RMS vs a chosen x-axis parameter for selected simulations.\n\n"
            "Selection (Option A): use --filter KEY=VALUE repeatedly to include only\n"
            "cases whose case_params.json matches all filters.\n"
            "Examples:\n"
            "  --filter H=0.0083 --filter mesher=snappy --filter n_cpus=8\n\n"
            "X-axis: use --x KEY where KEY is any top-level key in case_params.json,\n"
            "or 'tilt_deg' (also inferred from case directory name if missing).\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--filter",
        action="append",
        default=[],
        help="Filter cases by case_params.json (repeatable), e.g. --filter H=0.0083",
    )
    parser.add_argument(
        "--x",
        required=True,
        help="X-axis key from case_params.json (or 'tilt_deg')",
    )
    parser.add_argument(
        "--metric",
        default="l2_rms",
        help="Metric key from comparison_summary.csv (default: l2_rms)",
    )
    parser.add_argument(
        "--summary",
        default=os.path.join("postProcessing", "interface_compare", "comparison_summary.csv"),
        help="Relative path to comparison_summary.csv inside each case dir",
    )
    parser.add_argument(
        "--out-prefix",
        default=os.path.join("postProcessing", "l2_rms_vs_x"),
        help="Output prefix for CSV/PNG (default: postProcessing/l2_rms_vs_x)",
    )
    return parser

def _parse_filters(raw_filters):
    filters = []
    for item in raw_filters:
        if "=" not in item:
            raise ValueError(f"Invalid --filter '{item}', expected KEY=VALUE")
        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            raise ValueError(f"Invalid --filter '{item}', empty key")
        filters.append((key, value))
    return filters

def _get_x_value(params, case_dir, x_key):
    if params and x_key in params:
        return params.get(x_key)
    if x_key == "tilt_deg":
        return _extract_tilt(case_dir)
    return None

def main():
    parser = _build_parser()
    args = parser.parse_args()
    try:
        filters = _parse_filters(args.filter)
    except ValueError as e:
        print(e)
        return 2

    cases = sorted([d for d in glob.glob("case_*") if os.path.isdir(d)])
    rows = []
    for case in cases:
        params = _load_case_params(case)
        if not _apply_filters(params, filters):
            continue
        x_val = _get_x_value(params, case, args.x)
        if x_val is None:
            continue
        summary = os.path.join(case, args.summary)
        l2 = _read_metric(summary, args.metric)
        if l2 is None:
            continue
        rows.append((x_val, l2, case))

    if not rows:
        print("No matching cases with metric found.")
        return 1

    # If multiple cases share the same tilt, keep all but sort by tilt then name.
    rows.sort(key=lambda r: (r[0], r[2]))

    out_csv = f"{args.out_prefix}.csv"
    os.makedirs("postProcessing", exist_ok=True)
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([args.x, args.metric, "case_dir"])
        for x_val, l2, case in rows:
            w.writerow([x_val, l2, case])

    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"Wrote {out_csv}. Matplotlib not available: {e}")
        return 0

    xs = [r[0] for r in rows]
    ys = [r[1] for r in rows]

    plt.figure(figsize=(6, 4))
    plt.plot(xs, ys, marker="o", linestyle="-")
    plt.xlabel(args.x)
    plt.ylabel(args.metric)
    plt.title(f"{args.metric} vs {args.x}")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    out_png = f"{args.out_prefix}.png"
    plt.savefig(out_png, dpi=200)
    print(f"Wrote {out_csv} and {out_png}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
