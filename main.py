#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
import argparse
import threading

def _patch_vtk_for_pyvista():
    try:
        import vtkmodules.vtkFiltersSources as vfs
    except Exception:
        return
    if hasattr(vfs, "vtkCapsuleSource"):
        return
    try:
        class vtkCapsuleSource(vfs.vtkSphereSource):
            pass
        vfs.vtkCapsuleSource = vtkCapsuleSource
    except Exception:
        pass

def _import_pyvista():
    try:
        import pyvista as pv
        return pv
    except Exception:
        _patch_vtk_for_pyvista()
        import pyvista as pv
        return pv

# --- Dependency Management ---
def ensure_dependencies():
    """Check and install required Python packages with robust venv detection."""
    # Building/running cases only needs the standard library. Allow skipping any
    # venv/pip work (useful on HPC/login nodes) unless post-processing is used.
    if os.environ.get("SLOSHING_SKIP_DEPS", "0").strip().lower() in ("1", "yes", "true", "on"):
        return

    base_dir = os.path.dirname(os.path.abspath(__file__))
    venv_path = os.path.join(base_dir, "sloshing")
    restarted = os.environ.get("SLOSHING_ENV_RESTARTED") == "1"

    # Allow users to opt out of venv auto-switching entirely.
    if os.environ.get("SLOSHING_NO_VENV", "0").strip().lower() in ("1", "yes", "true", "on"):
        return

    # Robust detection of whether we are running in the 'sloshing' venv
    in_venv = False
    active_venv = os.environ.get("VIRTUAL_ENV")
    if active_venv and os.path.exists(active_venv) and os.path.exists(venv_path):
        try:
            if os.path.samefile(active_venv, venv_path):
                in_venv = True
        except: pass
    
    if not in_venv:
        try:
            if os.path.exists(venv_path) and os.path.samefile(sys.prefix, venv_path):
                in_venv = True
        except: pass

    # Get venv python/pip paths
    if sys.platform == "win32":
        pip_path = os.path.join(venv_path, "Scripts", "pip")
        python_path = os.path.join(venv_path, "Scripts", "python")
    else:
        pip_path = os.path.join(venv_path, "bin", "pip")
        python_path = os.path.join(venv_path, "bin", "python")
        if not os.path.exists(python_path):
            python_path = os.path.join(venv_path, "bin", "python3")

    # If a venv already exists, always jump into it immediately.
    # This avoids the noisy "Missing numpy" path when the user runs `python main.py`
    # with a system python that lacks the deps.
    if not in_venv and os.path.exists(venv_path) and not restarted:
        os.environ["SLOSHING_ENV_RESTARTED"] = "1"
        os.execv(python_path, [python_path] + sys.argv)

    # Ensure venv exists (if we're not in it yet, we will create then re-exec).
    if not os.path.exists(venv_path):
        print(f"Creating virtual environment: {venv_path}")
        subprocess.run([sys.executable, "-m", "venv", venv_path], check=True)
        os.environ["SLOSHING_ENV_RESTARTED"] = "1"
        os.execv(python_path, [python_path] + sys.argv)

    # We are inside the venv: verify packages without importing them (fast/quiet).
    try:
        import importlib.util

        required = ["numpy", "scipy", "matplotlib", "imageio", "imageio_ffmpeg", "h5py"]
        missing = [m for m in required if importlib.util.find_spec(m) is None]
        if missing:
            print(f"Installing python deps into {venv_path} (missing: {', '.join(missing)})...")
            req_file = os.path.join(base_dir, "requirements.txt")
            subprocess.run([pip_path, "install", "--upgrade", "pip", "-q"], check=False, capture_output=True)
            res = subprocess.run([pip_path, "install", "-r", req_file, "-q"], capture_output=True)
            if res.returncode != 0:
                stderr = (res.stderr or b"").decode("utf-8", errors="ignore").strip().splitlines()
                tail = "\n".join(stderr[-30:]) if stderr else "pip failed (no stderr)"
                raise RuntimeError(tail)

        # Optional: PyVista (may be unavailable on some systems)
        if importlib.util.find_spec("pyvista") is None:
            print("Note: PyVista not available; OpenFOAM interface extraction will be skipped.")
        return
    except Exception as e:
        if in_venv or restarted:
            print(f"\n❌ Python dependency setup failed: {e}")
            print(f"   Executable: {sys.executable}")
            print(f"   Venv Path:  {venv_path}")
            print("\n   If this is a broken venv, delete it and rerun:")
            print(f"   rm -rf {venv_path}")
            sys.exit(1)
        raise
    except Exception as e:
        print(f"\n❌ Unexpected error during dependency check: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

# Run dependency check
ensure_dependencies()

import math
import itertools
import re
import json
import tempfile

# --- Constants & Defaults ---
TEMPLATE_DIR = "circularTiltingTank"
VIDEO_FPS = 30
DEFAULTS = {
    "H": 0.01,
    "D": 0.0083,
    "mesh": 0.0005,
    "mesher": "gmsh",
    "geo": "flat",
    "tilt_deg": 5.0,
    "duration": 50.0,
    "dt": 0.1,
    "n_cpus": 1,
}

# --- Utility Functions ---

def parse_range(s):
    """
    Parses a MATLAB-style range (start:step:end) or comma-separated list.
    Returns a list of floats.
    """
    s = s.strip()
    if ':' in s:
        parts = s.split(':')
        if len(parts) == 2:
            start, end = float(parts[0]), float(parts[1])
            step = 1.0
        elif len(parts) == 3:
            start, step, end = float(parts[0]), float(parts[1]), float(parts[2])
        else:
            raise ValueError(f"Invalid range format: {s}")
        # Generate range
        vals = []
        v = start
        while v <= end + 1e-9:  # Tolerance for floating point
            vals.append(round(v, 6))
            v += step
        return vals
    else:
        # Comma-separated
        return [float(x.strip()) for x in s.split(',')]

def parse_indices(s, max_idx):
    """
    Parses comma-separated indices and ranges (e.g., "1, 3-5, 7").
    Returns a list of 0-indexed integers.
    """
    indices = set()
    for part in s.split(','):
        part = part.strip()
        if '-' in part:
            start, end = part.split('-')
            for i in range(int(start), int(end) + 1):
                if 1 <= i <= max_idx:
                    indices.add(i - 1)
        else:
            i = int(part)
            if 1 <= i <= max_idx:
                indices.add(i - 1)
    return sorted(list(indices))

def format_time(hours):
    """Formats hours as an HH:MM:SS Slurm time string."""
    if hours is None:
        return "00:00:00"
    total_minutes = int(math.ceil(max(hours, 0.0) * 60.0))
    h = total_minutes // 60
    m = total_minutes % 60
    return f"{h:02d}:{m:02d}:00"

def _patch_alpha_water_bc(case_dir):
    path = os.path.join(case_dir, "0", "alpha.water")
    if not os.path.exists(path):
        return
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    content = "".join(lines)
    if "AlphaContactAngle" not in content and "constantAlphaContactAngle" not in content:
        return
    out = []
    in_walls = False
    for line in lines:
        if re.match(r"\s*walls\s*\{", line):
            in_walls = True
            out.append(line)
            continue
        if in_walls:
            if re.match(r"\s*\}", line):
                in_walls = False
                out.append(line)
                continue
            if re.match(r"\s*type\s+", line):
                prefix = re.match(r"^(\s*)", line).group(1)
                out.append(f"{prefix}type            contactAngle;\n")
                continue
        out.append(line)
    with open(path, "w") as f:
        f.write("".join(out))

def _write_functions_dict(case_dir, params):
    """
    Writes a minimal, portable functionObjects file.
    We only rely on `probes` (widely available) to avoid per-version syntax issues.
    """
    H = float(params.get("H", DEFAULTS["H"]))
    D = float(params.get("D", DEFAULTS["D"]))
    R = 0.5 * D

    # Probe points inside the cylinder.
    # Bias sampling around the expected interface height (~H/2) so we can detect
    # interface stillness (alpha stops changing), not just small velocities.
    thetas = [0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi]
    r_ring = 0.45 * R

    points = []
    z_levels = [0.25 * H, 0.45 * H, 0.5 * H, 0.55 * H, 0.75 * H]
    for z in z_levels:
        zc = min(max(z, 0.05 * H), 0.95 * H)
        points.append((0.0, 0.0, zc))
        for th in thetas:
            points.append((r_ring * math.cos(th), r_ring * math.sin(th), zc))

    functions_path = os.path.join(case_dir, "system", "functions")
    content = [
        "/*--------------------------------*- C++ -*----------------------------------*\\",
        "  =========                 |",
        "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox",
        "   \\\\    /   O peration     | Website:  https://openfoam.org",
        "    \\\\  /    A nd           | Version:  13",
        "     \\\\/     M anipulation  |",
        "\\*---------------------------------------------------------------------------*/",
        "FoamFile",
        "{",
        "    format      ascii;",
        "    class       dictionary;",
        "    location    \"system\";",
        "    object      functions;",
        "}",
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",
        "",
        "probesU",
        "{",
        "    type            probes;",
        "    libs            (\"libsampling.so\");",
        "    writeControl    timeStep;",
        "    writeInterval   5;",
        "    fixedLocations  true;",
        "    fields",
        "    (",
        "        U",
        "        alpha.water",
        "    );",
        "    probeLocations",
        "    (",
    ]
    for x, y, z in points:
        content.append(f"        ({x:.8g} {y:.8g} {z:.8g})")
    content += [
        "    );",
        "}",
        "",
        "// ************************************************************************* //",
        "",
    ]

    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    with open(functions_path, "w") as f:
        f.write("\n".join(content))

def _write_case_params(case_dir, params):
    path = os.path.join(case_dir, "case_params.json")
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(params, f, indent=2, sort_keys=True)
    except Exception as e:
        print(f"  ⚠️  Failed to write case_params.json: {e}")

def _load_case_params(case_dir):
    """
    Read-only access to case parameters.
    Preference order:
      1) case_params.json (authoritative build-time params)
      2) parse from folder name (fallback for legacy cases)
    """
    path = os.path.join(case_dir, "case_params.json")
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                return json.load(f)
        except Exception as e:
            print(f"  ⚠️  Failed to read case_params.json, falling back to name: {e}")
    return parse_case_params(os.path.basename(case_dir))

def _ensure_functions_dict(case_dir):
    params = _load_case_params(case_dir)
    _write_functions_dict(case_dir, params)

def _patch_control_dict_for_speed(case_dir, params):
    control_path = os.path.join(case_dir, "system", "controlDict")
    if not os.path.exists(control_path):
        return
    with open(control_path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    # Moderate, stability-first time stepping for VOF.
    content = re.sub(r'(^\s*maxCo\s+)[^;]+;', r'\g<1>1.0;', content, flags=re.M)
    content = re.sub(r'(^\s*maxAlphaCo\s+)[^;]+;', r'\g<1>0.5;', content, flags=re.M)

    # `dt` is treated as the maximum allowed timestep (maxDeltaT).
    max_dt = float(params.get("dt", DEFAULTS["dt"]))
    content = re.sub(r'(^\s*maxDeltaT\s+)[^;]+;', r'\g<1>' + f"{max_dt:g}" + ';', content, flags=re.M)

    # Start small to avoid first-step blow-ups; adjustable time step will ramp up.
    dt0 = min(max_dt, 1e-4)
    content = re.sub(r'(^\s*deltaT\s+)[^;]+;', r'\g<1>' + f"{dt0:g}" + ';', content, flags=re.M)

    # Ensure we always write something even if we stop early (and keep I/O low).
    # `stopAt writeNow` will still force a final write at steady state.
    duration = float(params.get("duration", DEFAULTS["duration"]))
    write_interval = 0.25
    content = re.sub(r'(^\s*endTime\s+)[^;]+;', r'\g<1>' + f"{duration:g}" + ';', content, flags=re.M)
    content = re.sub(r'(^\s*writeControl\s+)[^;]+;', r'\g<1>adjustableRunTime;', content, flags=re.M)
    content = re.sub(r'(^\s*writeInterval\s+)[^;]+;', r'\g<1>' + f"{write_interval:g}" + ';', content, flags=re.M)
    content = re.sub(r'(^\s*purgeWrite\s+)[^;]+;', r'\g<1>0;', content, flags=re.M)
    with open(control_path, "w") as f:
        f.write(content)

def _patch_fvsolution_prefpoint(case_dir, params):
    """
    Ensure pRefPoint is inside the domain.

    The template previously used z=0.15 (from a larger tank). For small H,
    this can place pRefPoint outside the mesh and destabilize the pressure
    solution (can lead to SIGFPE inside GAMG/PCG).
    """
    fv_path = os.path.join(case_dir, "system", "fvSolution")
    if not os.path.exists(fv_path):
        return
    H = float(params.get("H", DEFAULTS["H"]))
    z = 0.5 * H
    with open(fv_path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    content2, _ = re.subn(
        r'(^\s*pRefPoint\s*)\([^)]*\)\s*;',
        r'\g<1>' + f"(0 0 {z:.6g});",
        content,
        flags=re.M,
    )
    # Match p_rgh atmosphere p0=0 (gauge) to avoid artificial driving pressure.
    content2 = re.sub(
        r'(^\s*pRefValue\s+)[^;]+;',
        r'\g<1>0;',
        content2,
        flags=re.M,
    )
    if content2 != content:
        with open(fv_path, "w") as f:
            f.write(content2)

def _patch_fvsolution_for_stability(case_dir):
    """
    Conservative defaults that tend to keep VOF bounded.
    We keep this minimal to avoid version-specific dictionary pitfalls.
    """
    fv_path = os.path.join(case_dir, "system", "fvSolution")
    if not os.path.exists(fv_path):
        return
    with open(fv_path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    # alpha.water block: more correction + subcycling.
    # We only patch if the keys exist in the file to avoid injecting unsupported syntax.
    content = re.sub(r'(^\s*nCorrectors\s+)[^;]+;', r'\g<1>2;', content, flags=re.M)
    content = re.sub(r'(^\s*nSubCycles\s+)[^;]+;', r'\g<1>2;', content, flags=re.M)

    # PIMPLE block: at least 2 correctors improves robustness.
    content = re.sub(r'(^\s*nCorrectors\s+)[^;]+;', r'\g<1>2;', content, flags=re.M)

    with open(fv_path, "w") as f:
        f.write(content)

def _check_mesh_quality_gmsh(case_dir, msh_path, target_lc):
    try:
        from mesh_quality import analyze_msh2, write_summary
    except Exception as e:
        print(f"  ⚠️  Mesh quality check skipped (cannot import mesh_quality): {e}")
        return {"ok": True, "summary": None}
    if not os.path.exists(msh_path):
        return {"ok": True, "summary": None}
    summary = analyze_msh2(msh_path)
    out_path = os.path.join(case_dir, "postProcessing", "mesh_quality.json")
    write_summary(summary, out_path)

    # Warn aggressively if tiny elements exist; they force tiny deltaT and huge runtime.
    min_edge = summary.min_edge
    if min_edge is None:
        return {"ok": True, "summary": summary}
    ratio = (min_edge / target_lc) if target_lc > 0 else 1.0
    ok = True
    if ratio < 0.3:
        ok = False
        print(
            f"  ⚠️  Mesh warning: min edge {min_edge:.3g}m is {ratio:.2f}x target lc={target_lc:g}m; "
            "expect very small deltaT and very slow runs."
        )
    if summary.max_aspect_ratio is not None and summary.max_aspect_ratio > 20:
        ok = False
        print(
            f"  ⚠️  Mesh warning: max aspect ratio ~{summary.max_aspect_ratio:.1f} (high); "
            "may hurt stability/time step."
        )
    # Also print element count to expose accidental over-refinement.
    if summary.n_tets > 0:
        print(f"  Mesh: {summary.n_tets:,} tetrahedra (nodes: {summary.n_nodes:,}).")
    return {"ok": ok, "summary": summary}

def _preflight_mesh_quality(params):
    """Build a temporary Gmsh mesh and warn if it produces tiny elements."""
    if params.get("mesher", "gmsh") != "gmsh":
        return {"ok": True, "summary": None}
    gmsh_path = shutil.which("gmsh")
    if not gmsh_path:
        print("Mesh preflight: gmsh not found; skipping mesh-quality check.")
        return {"ok": True, "summary": None}
    mesh_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), TEMPLATE_DIR, "generate_mesh.py")
    if not os.path.exists(mesh_script):
        return {"ok": True, "summary": None}
    with tempfile.TemporaryDirectory(prefix="openfoam_meshcheck_") as tmpdir:
        try:
            subprocess.run(
                [
                    sys.executable,
                    mesh_script,
                    str(params["H"]),
                    str(params["D"]),
                    str(params["mesh"]),
                    params["geo"],
                    params.get("mesher", "gmsh"),
                ],
                cwd=tmpdir,
                check=True,
                capture_output=True,
            )
            subprocess.run(
                ["gmsh", "-3", "cylinder.geo", "-format", "msh2", "-o", "cylinder.msh"],
                cwd=tmpdir,
                check=True,
                capture_output=True,
            )
            # Write any JSON into the tempdir so we don't clutter the repo.
            return _check_mesh_quality_gmsh(tmpdir, os.path.join(tmpdir, "cylinder.msh"), float(params["mesh"]))
        except subprocess.CalledProcessError as e:
            msg = (e.stderr or b"").decode("utf-8", errors="ignore").strip()
            print(f"Mesh preflight: failed ({msg[:200]})")
        except Exception as e:
            print(f"Mesh preflight: failed ({e})")
    return {"ok": True, "summary": None}

def _start_mesh_preflight_async(params):
    """
    Fire-and-forget mesh preflight so the UI doesn't block.
    Prints a completion message when done.
    """
    def _worker():
        result = _preflight_mesh_quality(params)
        status = "ok" if result.get("ok", True) else "warnings"
        print(f"Mesh preflight complete ({status}).")

    t = threading.Thread(target=_worker, name="mesh_preflight", daemon=True)
    t.start()
    return t

def get_case_name(params):
    """Generates a unique case folder name from parameters."""
    return (
        f"case_H{params['H']}_D{params['D']}_{params['geo']}_tilt_"
        f"T{params['tilt_deg']}_m{params['mesh']}"
    )

def _list_time_folders(path):
    try:
        return sorted(
            [d for d in os.listdir(path) if d.replace(".", "", 1).isdigit()],
            key=lambda x: float(x),
        )
    except FileNotFoundError:
        return []

def _latest_time_in_dir(path):
    times = [float(d) for d in _list_time_folders(path)]
    return max(times) if times else None

def _get_latest_time(case_dir):
    # Prefer serial latest time if present (excluding 0)
    times = [t for t in _list_time_folders(case_dir) if t != "0"]
    if times:
        return max(float(t) for t in times)
    # Fall back to parallel latest time across processors
    latest = None
    for pd in [d for d in os.listdir(case_dir) if d.startswith("processor")]:
        ppath = os.path.join(case_dir, pd)
        t = _latest_time_in_dir(ppath)
        if t is None or t == 0:
            continue
        latest = t if latest is None else max(latest, t)
    return latest

def is_case_done(case_dir, duration):
    """Checks if the simulation for this case is complete."""
    latest = _get_latest_time(case_dir)
    if latest is None:
        return False
    return latest >= float(duration) - 1e-9

def _get_case_duration(case_dir):
    params = _load_case_params(case_dir)
    return float(params.get("duration", DEFAULTS["duration"]))

def _get_case_status(case_dir):
    duration = _get_case_duration(case_dir)
    if is_case_done(case_dir, duration):
        return "FINISHED"
    if has_case_progress(case_dir):
        return "UNFINISHED"
    return "NEW"

def has_case_progress(case_dir):
    """Checks if the case has any progress (output folders or processor dirs)."""
    # Check for serial time folders (excluding '0')
    time_folders = [d for d in os.listdir(case_dir) if d.replace('.','',1).isdigit() and d != '0']
    if time_folders:
        return True
    # Check for parallel time folders (excluding '0')
    proc_dirs = [d for d in os.listdir(case_dir) if d.startswith("processor")]
    for pd in proc_dirs:
        ppath = os.path.join(case_dir, pd)
        try:
            times = [d for d in os.listdir(ppath) if d.replace('.','',1).isdigit() and d != '0']
            if times:
                return True
        except FileNotFoundError:
            continue
    return False

def parse_case_params(case_name):
    """Extracts parameters from a case folder name."""
    # Format: case_H{H}_D{D}_{geo}_tilt_T{tilt}[_d{duration}]_m{mesh}
    # We now skip duration in the name, but support parsing it for backward compatibility.
    match = re.match(
        r'case_H([\d.]+)_D([\d.]+)_(\w+)_tilt_T([\d.]+)(?:_d([\d.]+))?_m([\d.]+)',
        case_name
    )
    if not match:
        return DEFAULTS.copy()

    # If duration group (5) is None, use default.
    duration = float(match.group(5)) if match.group(5) else DEFAULTS["duration"]

    return {
        "H": float(match.group(1)),
        "D": float(match.group(2)),
        "geo": match.group(3),
        "tilt_deg": float(match.group(4)),
        "duration": duration,
        "mesh": float(match.group(6)),
        "dt": DEFAULTS["dt"],
        "mesher": DEFAULTS["mesher"],
    }

def _load_mesh_quality_summary(case_dir: str):
    path = os.path.join(case_dir, "postProcessing", "mesh_quality.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            return json.load(f)
    except Exception:
        return None

def _estimate_effective_dt(params, dx=None, max_co=1.0, max_alpha_co=0.5):
    """
    Estimate the timestep the solver is likely to run at (order-of-magnitude).
    Used for Oscar walltime sizing, since runtime scales with number of steps.
    """
    H = float(params["H"])
    dt_max = float(params.get("dt", DEFAULTS["dt"]))
    
    # Gravity-wave velocity scale.
    u_gravity = math.sqrt(9.81 * max(H, 1e-9))
    
    # For small scales, surface tension and numerical interface noise often
    # drive higher local velocities than gravity alone. We add a conservative 
    # floor for u_est (0.5 m/s) to avoid over-optimistic timesteps.
    u_est = max(0.5, u_gravity)
    
    dx_est = float(dx) if dx else float(params["mesh"])

    dt_co = max_co * dx_est / u_est
    dt_alpha = max_alpha_co * dx_est / u_est
    dt_eff = min(dt_max, dt_co, dt_alpha)
    # Avoid nonsense from bad inputs; this is only for estimation.
    return max(dt_eff, 1e-7)

def _read_control_dict_values(case_dir: str):
    path = os.path.join(case_dir, "system", "controlDict")
    if not os.path.exists(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
    except Exception:
        return {}
    values = {}
    for key in ("maxCo", "maxAlphaCo", "maxDeltaT", "deltaT"):
        m = re.search(rf"^\s*{key}\s+([^;]+);", content, flags=re.M)
        if m:
            try:
                values[key] = float(m.group(1))
            except ValueError:
                pass
    return values

def _read_phase_sigma(case_dir: str):
    path = os.path.join(case_dir, "constant", "phaseProperties")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
    except Exception:
        return None
    m = re.search(r"^\s*sigma\s+([^;]+);", content, flags=re.M)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None

def _read_rho_nu(case_dir: str, phase: str):
    path = os.path.join(case_dir, "constant", f"physicalProperties.{phase}")
    if not os.path.exists(path):
        return None, None
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
    except Exception:
        return None, None
    def _read_value(key):
        m = re.search(rf"^\s*{key}\s+([^;]+);", content, flags=re.M)
        if not m:
            return None
        try:
            return float(m.group(1))
        except ValueError:
            return None
    return _read_value("rho"), _read_value("nu")

def _estimate_capillary_dt(dx, rho1, rho2, sigma, c_sigma=0.063):
    if not dx or not rho1 or not rho2 or not sigma:
        return None
    return c_sigma * math.sqrt(((rho1 + rho2) * (dx ** 3)) / sigma)

def _estimate_viscous_dt(dx, nu):
    if not dx or not nu:
        return None
    return (dx ** 2) / (4.0 * nu)

def _read_n_cpus_from_decompose(case_dir: str, default: int = 1) -> int:
    """Read numberOfSubdomains from decomposeParDict, falling back to default."""
    path = os.path.join(case_dir, "system", "decomposeParDict")
    if not os.path.exists(path):
        return default
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
    except Exception:
        return default
    match = re.search(r'numberOfSubdomains\s+([\d.]+)\s*;', content)
    if not match:
        return default
    try:
        return max(1, int(float(match.group(1))))
    except ValueError:
        return default

def _write_mesh_tool(case_dir: str, tool: str) -> None:
    path = os.path.join(case_dir, "system", "meshTool")
    try:
        with open(path, "w", encoding="utf-8") as f:
            f.write(f"{tool}\n")
    except Exception:
        pass

def _read_mesh_tool(case_dir: str, default: str = "gmsh") -> str:
    path = os.path.join(case_dir, "system", "meshTool")
    if not os.path.exists(path):
        snappy_dict = os.path.join(case_dir, "system", "snappyHexMeshDict")
        if os.path.exists(snappy_dict):
            return "snappy"
        return default
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            val = f.read().strip().lower()
        return val or default
    except Exception:
        return default

def estimate_resources(params, case_dir=None, mesh_summary=None):
    """
    Estimates required CPUs, memory, and wall-clock time.
    Model is based on number of cells * number of timesteps (dt matters).
    """
    h, d, mesh_size = float(params["H"]), float(params["D"]), float(params["mesh"])
    duration = float(params["duration"])

    # Prefer measured mesh info when available (post-build), otherwise estimate.
    ms = mesh_summary
    if ms is None and case_dir:
        ms = _load_mesh_quality_summary(case_dir)

    if ms and hasattr(ms, "n_tets"):
        n_cells = float(getattr(ms, "n_tets") or 0.0)
    elif isinstance(ms, dict) and ms.get("n_tets"):
        n_cells = float(ms["n_tets"])
    else:
        vol = math.pi * ((d / 2.0) ** 2) * h
        cell_vol = mesh_size**3
        n_cells = vol / max(cell_vol, 1e-30)
    n_cells = max(n_cells, 1.0)

    if ms and hasattr(ms, "min_edge"):
        dx = getattr(ms, "min_edge")
    elif isinstance(ms, dict):
        dx = ms.get("min_edge")
    else:
        dx = None
    dx = dx or mesh_size
    control_vals = _read_control_dict_values(case_dir) if case_dir else {}
    max_co = float(control_vals.get("maxCo", 1.0))
    max_alpha_co = float(control_vals.get("maxAlphaCo", 0.5))
    dt_max = float(control_vals.get("maxDeltaT", params.get("dt", DEFAULTS["dt"])))
    params_eff = params.copy()
    params_eff["dt"] = dt_max
    dt_eff = _estimate_effective_dt(params_eff, dx=dx, max_co=max_co, max_alpha_co=max_alpha_co)

    # Add capillary + viscous stability limits when phase properties are available.
    if case_dir:
        sigma = _read_phase_sigma(case_dir)
        rho_w, nu_w = _read_rho_nu(case_dir, "water")
        rho_a, nu_a = _read_rho_nu(case_dir, "air")
    else:
        sigma = _read_phase_sigma(TEMPLATE_DIR)
        rho_w, nu_w = _read_rho_nu(TEMPLATE_DIR, "water")
        rho_a, nu_a = _read_rho_nu(TEMPLATE_DIR, "air")

    if sigma and rho_w and rho_a:
        dt_sigma = _estimate_capillary_dt(dx, rho_w, rho_a, sigma)
        if dt_sigma:
            dt_eff = min(dt_eff, dt_sigma)

    nu_max = None
    if nu_w and nu_a:
        nu_max = max(nu_w, nu_a)
    elif nu_w:
        nu_max = nu_w
    elif nu_a:
        nu_max = nu_a

    if nu_max:
        dt_nu = _estimate_viscous_dt(dx, nu_max)
        if dt_nu:
            dt_eff = min(dt_eff, dt_nu)
    n_steps = max(1.0, duration / dt_eff)

    # Calibrated from observed Oscar runs: ~0.0016 CPU-hr per (Mcell-step) in this repo.
    cpu_hr_per_mcell_step = 0.0016
    total_cpu_hours = cpu_hr_per_mcell_step * (n_cells / 1e6) * n_steps

    # Buffers for variability, I/O, and aggressive VOF sub-cycling.
    # Increased from 3.0 to 6.0 based on observed small-scale VOF performance.
    total_cpu_hours *= 6.0

    # Suggest CPUs to target ~2-4 hours wall time (but avoid over-parallelization).
    target_wall_h = 3.0
    suggested_cpus = max(1, int(math.ceil(total_cpu_hours / target_wall_h)))

    # Efficiency guard: keep >=15k cells/core if possible.
    suggested_cpus = min(suggested_cpus, max(1, int(n_cells / 15000)))
    suggested_cpus = min(suggested_cpus, 32)
    if suggested_cpus > 1:
        suggested_cpus = 2 ** math.floor(math.log2(suggested_cpus))

    wall_clock_hours = total_cpu_hours / suggested_cpus
    safe_hours = wall_clock_hours * 2.0 + 1.0
    safe_hours = max(safe_hours, 1.0)
    safe_hours = min(safe_hours, 24.0)

    time_limit = format_time(safe_hours)

    # Memory: conservative scaling.
    mem_gb = (n_cells / 100000.0) * 0.2 + 2.0
    mem_gb = max(4.0, math.ceil(mem_gb))
    return f"{int(mem_gb)}G", time_limit, n_cells, suggested_cpus

# --- Core Actions ---

def setup_case(params):
    """Creates the case directory and runs setup scripts."""
    case_name = get_case_name(params)
    
    if os.path.exists(case_name):
        print(f"  ⚠️  {case_name} already exists. Skipping.")
        return case_name
    
    print(f"  📂 Creating: {case_name}")
    shutil.copytree(TEMPLATE_DIR, case_name)
    
    # Ensure writable
    for root, dirs, files in os.walk(case_name):
        for d in dirs:
            os.chmod(os.path.join(root, d), 0o777)
        for f in files:
            os.chmod(os.path.join(root, f), 0o666)

    _patch_alpha_water_bc(case_name)
    _ensure_functions_dict(case_name)
    _patch_fvsolution_prefpoint(case_name, params)
    _patch_fvsolution_for_stability(case_name)

    cwd = os.path.join(os.getcwd(), case_name)
    
    # Static tilt (rotated gravity + zero motion)
    subprocess.run([
        sys.executable, "generate_tilt.py",
        str(params["tilt_deg"]), str(params["duration"]), str(params["dt"])
    ], cwd=cwd, check=True, capture_output=True)
    
    # Fields
    subprocess.run([sys.executable, "update_setFields.py", str(params['H'])], 
                   cwd=cwd, check=True, capture_output=True)
    
    # Mesh Geometry
    subprocess.run([
        sys.executable, "generate_mesh.py", 
        str(params['H']), str(params['D']), str(params['mesh']), params['geo'], params.get("mesher", "gmsh")
    ], cwd=cwd, check=True, capture_output=True)

    _write_mesh_tool(case_name, params.get("mesher", "gmsh"))
    
    # Run Gmsh
    gmsh_path = shutil.which("gmsh")
    if gmsh_path and params.get("mesher", "gmsh") == "gmsh":
        subprocess.run([
            "gmsh", "-3", "cylinder.geo", "-format", "msh2", "-o", "cylinder.msh"
        ], cwd=cwd, check=True, capture_output=True)
        _check_mesh_quality_gmsh(case_name, os.path.join(cwd, "cylinder.msh"), float(params["mesh"]))
    elif params.get("mesher", "gmsh") == "gmsh":
        print("  ❌ gmsh not found in PATH. Cannot generate mesh.")

    # Parallel Setup (Inject numberOfSubdomains)
    if params.get('n_cpus', 1) > 1:
        decomp_path = os.path.join(cwd, "system", "decomposeParDict")
        if os.path.exists(decomp_path):
            with open(decomp_path, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
            content = re.sub(
                r'numberOfSubdomains\s+[\d.]+\s*;',
                f'numberOfSubdomains {int(params["n_cpus"])};',
                content,
            )
            with open(decomp_path, 'w') as f:
                f.write(content)

    # Update controlDict endTime
    control_path = os.path.join(cwd, "system", "controlDict")
    if os.path.exists(control_path):
        with open(control_path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
        content = re.sub(r'endTime\s+[\d.]+;', f'endTime {params["duration"]};', content)
        # Use dt as the maximum dt target; start smaller to keep the first step stable.
        dt_max = float(params["dt"])
        dt0 = min(dt_max, 1e-4)
        content = re.sub(r'deltaT\s+[\d.]+;', f'deltaT {dt0};', content)
        with open(control_path, 'w') as f:
            f.write(content)
    _patch_control_dict_for_speed(case_name, params)

    # Persist build-time parameters for read-only usage during run.
    _write_case_params(case_name, params)
        
    return case_name

def run_case_local(case_name, n_cpus=1):
    """Runs simulation locally."""
    _patch_alpha_water_bc(case_name)
    _ensure_functions_dict(case_name)
    params = _load_case_params(case_name)
    mesh_tool = _read_mesh_tool(case_name, default=params.get("mesher", "gmsh"))
    _patch_fvsolution_prefpoint(case_name, params)
    _patch_fvsolution_for_stability(case_name)
    shutil.copy2(os.path.join(TEMPLATE_DIR, "adaptive_stop.py"), os.path.join(case_name, "adaptive_stop.py"))
    # Check for existing progress (non-zero time folders only)
    has_progress = has_case_progress(case_name)
            
    if has_progress:
        print(f"  🏃 Resuming {case_name} (CPUs={n_cpus})...")
        subprocess.run(
            ["make", "-C", case_name, "resume", f"N_CPUS={n_cpus}", f"MESH_TOOL={mesh_tool}", "ADAPTIVE_STOP=1"],
            check=True,
        )
    else:
        print(f"  🏃 Running {case_name} (CPUs={n_cpus})...")
        subprocess.run(
            ["make", "-C", case_name, "run", f"N_CPUS={n_cpus}", f"MESH_TOOL={mesh_tool}", "ADAPTIVE_STOP=1"],
            check=True,
        )

def run_case_oscar(case_name, params, is_oscar):
    """Submits job to Slurm on Oscar."""
    _patch_alpha_water_bc(case_name)
    _ensure_functions_dict(case_name)
    _patch_fvsolution_prefpoint(case_name, params)
    _patch_fvsolution_for_stability(case_name)
    shutil.copy2(os.path.join(TEMPLATE_DIR, "adaptive_stop.py"), os.path.join(case_name, "adaptive_stop.py"))
    mem, time_limit, n_cells, _ = estimate_resources(params, case_dir=case_name)
    mesh_tool = _read_mesh_tool(case_name, default=params.get("mesher", "gmsh"))
    
    # Read the ACTUAL number of subdomains from the case folder
    # This is the single source of truth for parallel runs
    n_cpus = 1
    decomp_path = os.path.join(case_name, "system", "decomposeParDict")
    if os.path.exists(decomp_path):
        with open(decomp_path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
            match = re.search(r'numberOfSubdomains\s+([\d.]+)\s*;', content)
            if match:
                n_cpus = int(float(match.group(1)))

    slurm_dir = os.path.join(case_name, "slurm")
    os.makedirs(slurm_dir, exist_ok=True)
    script_path = os.path.join(slurm_dir, "run_simulation.slurm")
    
    header = [
        "#!/usr/bin/env bash",
        f"#SBATCH -J {case_name}",
        "#SBATCH -p batch",
        "#SBATCH -N 1",
        f"#SBATCH -n {n_cpus}",
        f"#SBATCH --time={time_limit}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH -o {slurm_dir}/slurm.%j.out",
        f"#SBATCH -e {slurm_dir}/slurm.%j.err",
        "#SBATCH --mail-type=END",
        "#SBATCH --mail-user=elvis_vera@brown.edu",
        "",
        "set -euo pipefail",
        "export OMP_NUM_THREADS=1",
        "",
        f"echo 'Case: {case_name}'",
        "# Check if we are resuming (parallel processors or serial time folders existed)",
        f"cd {case_name}",
        "has_progress=0",
        "if (ls -d [0-9]* 2>/dev/null | grep -v '^0$' | grep -q .); then",
        "    has_progress=1",
        "else",
        "    for p in processor*; do",
        "        [ -d \"$p\" ] || continue",
        "        if ls -d \"$p\"/[0-9]* 2>/dev/null | grep -v '/0$' | grep -q .; then",
        "            has_progress=1",
        "            break",
        "        fi",
        "    done",
        "fi",
        "if [ \"$has_progress\" -eq 1 ]; then",
        "    echo 'Found existing progress. Resuming simulation...'",
        f"    make resume OSCAR=1 N_CPUS={n_cpus} MESH_TOOL={mesh_tool} ADAPTIVE_STOP=1",
        "else",
        "    echo 'Starting fresh simulation...'",
        f"    make run OSCAR=1 N_CPUS={n_cpus} MESH_TOOL={mesh_tool} ADAPTIVE_STOP=1",
        "fi",
        "echo 'End: $(date)'"
    ]
    
    with open(script_path, "w") as f:
        f.write("\n".join(header))
    
    print(f"  🚀 Submitting {case_name} ({n_cpus} CPUs, {mem}, {time_limit})...")
    subprocess.run(["sbatch", script_path], check=True)

# --- Menu System ---

# Human-readable labels for parameters
PARAM_LABELS = {
    "H": "Height (m)",
    "D": "Diameter (m)",
    "mesh": "Mesh Size (m)",
    "mesher": "Mesher",
    "geo": "Geometry",
    "tilt_deg": "Tilt Angle (deg)",
    "duration": "Duration (s)",
    "dt": "Time Step (s)",
    "n_cpus": "Parallel CPUs (1=serial)",
}

GEO_OPTIONS = ["flat", "cap"]
MESHER_OPTIONS = ["gmsh", "snappy"]

def display_config(current_values, sweeps):
    """Displays the current configuration with any overrides."""
    print("\nCurrent Configuration:")
    param_keys = list(DEFAULTS.keys())
    for i, k in enumerate(param_keys):
        label = PARAM_LABELS.get(k, k)
        if k in sweeps:
            val_str = str(sweeps[k])
            print(f"  {i+1}) {label:25}: {val_str} (SWEEP)")
        else:
            print(f"  {i+1}) {label:25}: {current_values[k]}")

def menu_build_cases(is_oscar):
    """Submenu 1: Build Case Setups"""
    print("\n--- Build Case Setups ---")
    
    current_values = DEFAULTS.copy()
    sweeps = {}
    param_keys = list(DEFAULTS.keys())
    
    while True:
        display_config(current_values, sweeps)
        print("\nOptions: Enter number to edit, 'done' to build, 'cancel' to abort.")
        
        user_input = input("Select: ").strip()
        
        if user_input.lower() == 'cancel':
            print("Cancelled.")
            return
        
        if user_input.lower() == 'done':
            break
        
        # Parse selection
        param = None
        if user_input.isdigit():
            idx = int(user_input) - 1
            if 0 <= idx < len(param_keys):
                param = param_keys[idx]
        else:
            match = [k for k in DEFAULTS if k.lower() == user_input.lower()]
            if match:
                param = match[0]
        
        if not param:
            print(f"  Invalid selection: {user_input}")
            continue
        
        # Special handling for 'geo' (categorical)
        if param == 'geo':
            print(f"\n  Select geometry:")
            for i, opt in enumerate(GEO_OPTIONS):
                print(f"    {i+1}) {opt}")
            geo_input = input("  Choice (or comma-separated for sweep, e.g., '1,2'): ").strip()
            try:
                if ',' in geo_input:
                    indices = [int(x.strip()) - 1 for x in geo_input.split(',')]
                    sweeps[param] = [GEO_OPTIONS[i] for i in indices]
                else:
                    idx = int(geo_input) - 1
                    current_values[param] = GEO_OPTIONS[idx]
                    if param in sweeps:
                        del sweeps[param]
            except (ValueError, IndexError):
                print("  Invalid choice.")
            continue

        if param == 'mesher':
            print(f"\n  Select mesher:")
            for i, opt in enumerate(MESHER_OPTIONS):
                print(f"    {i+1}) {opt}")
            mesher_input = input("  Choice (or comma-separated for sweep, e.g., '1,2'): ").strip()
            try:
                if ',' in mesher_input:
                    indices = [int(x.strip()) - 1 for x in mesher_input.split(',')]
                    sweeps[param] = [MESHER_OPTIONS[i] for i in indices]
                else:
                    idx = int(mesher_input) - 1
                    current_values[param] = MESHER_OPTIONS[idx]
                    if param in sweeps:
                        del sweeps[param]
            except (ValueError, IndexError):
                print("  Invalid choice.")
            continue
        
        # Numeric parameters
        label = PARAM_LABELS.get(param, param)
        val_str = input(f"  Enter value(s) for '{label}' (single or sweep, e.g., 0.1 or 0.1:0.05:0.2): ").strip()
        try:
            vals = parse_range(val_str)
            if len(vals) == 1:
                current_values[param] = vals[0]
                if param in sweeps:
                    del sweeps[param]
            else:
                sweeps[param] = vals
        except ValueError as e:
            print(f"  ❌ Error: {e}")
    
    # Confirmation
    display_config(current_values, sweeps)
    
    # Build param_sets
    if not sweeps:
        param_sets = [current_values.copy()]
    else:
        lengths = [len(v) for v in sweeps.values()]
        
        if len(set(lengths)) == 1:
            print(f"\n✅ All sweep lists are length {lengths[0]}. Using ZIP mode.")
            keys = list(sweeps.keys())
            param_sets = []
            for i in range(lengths[0]):
                p = current_values.copy()
                for k in keys:
                    p[k] = sweeps[k][i]
                param_sets.append(p)
        else:
            total = 1
            for l in lengths:
                total *= l
            confirm = input(f"\n⚠️  Sweep lists have different lengths. This will generate {total} cases (Cartesian Product). Continue? (y/n): ").strip().lower()
            if confirm != 'y':
                print("Cancelled.")
                return
            
            keys = list(sweeps.keys())
            combos = list(itertools.product(*[sweeps[k] for k in keys]))
            param_sets = []
            for combo in combos:
                p = current_values.copy()
                for i, k in enumerate(keys):
                    p[k] = combo[i]
                param_sets.append(p)
    
    # Final Case Review & Resource Estimation
    print("\n" + "="*40)
    print("   Final Review & Resource Estimation")
    print("="*40)
    
    # Calculate for the first case in param_sets to show representative estimate.
    sample_params = param_sets[0]
    # Preflight mesh quality in the background so we don't block the UI.
    print("Mesh preflight started in background (results will print when ready).")
    _start_mesh_preflight_async(sample_params)
    mem, time_limit, n_cells, suggested_cpus = estimate_resources(
        sample_params, mesh_summary=None
    )
    
    print(f"Total Cases to Build: {len(param_sets)}")
    print(f"Estimated Cells per Case: {int(n_cells):,}")
    print(f"Suggested Wall-Clock Time: {time_limit}")
    print(f"Suggested Parallelization: {suggested_cpus} CPUs")

    if suggested_cpus > 1 and current_values['n_cpus'] == 1:
        print(f"\n💡 [RECOMMENDED] Multi-processing is highly recommended for this cell count.")
        use_multi = input(f"   Enable parallel execution with {suggested_cpus} CPUs? (y/n): ").strip().lower()
        if use_multi == 'y':
            for p in param_sets:
                p['n_cpus'] = suggested_cpus
    
    # Final confirmation
    confirm = input(f"\nConfirm building {len(param_sets)} case(s)? (y/n): ").strip().lower()
    if confirm != 'y':
        print("Cancelled.")
        return
    
    print(f"\nGenerating {len(param_sets)} case(s)...")
    for params in param_sets:
        setup_case(params)
    print("✅ Done building cases.")

def menu_run_cases(is_oscar):
    """Submenu 2: Run Cases"""
    print("\n--- Run Cases ---")
    
    cases = sorted([d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('case_')])
    if not cases:
        print("No cases found. Use 'Build Case Setups' first.")
        return
    
    # Display cases with status
    print("Available Cases:")
    for i, c in enumerate(cases):
        status = _get_case_status(c)
        status_str = f"({status})" if status != "NEW" else ""
        print(f"  {i+1}) {c} {status_str}")
    
    idx_str = input("\nEnter case indices to run (e.g., 1, 3-5, all): ").strip().lower()
    if idx_str == 'all':
        indices = list(range(len(cases)))
    else:
        indices = parse_indices(idx_str, len(cases))
    
    if not indices:
        print("No valid indices selected.")
        return
    
    print(f"\nRunning {len(indices)} case(s)...")
    
    has_openfoam = shutil.which("foamRun") is not None
    
    for i in indices:
        case_name = cases[i]
        params = _load_case_params(case_name)
        
        if is_oscar:
            run_case_oscar(case_name, params, is_oscar)
        elif has_openfoam:
            # Honor the case's decomposeParDict (source of truth for parallel runs).
            n_cpus = _read_n_cpus_from_decompose(case_name, default=1)
            run_case_local(case_name, n_cpus=n_cpus)
        else:
            print(f"  ❌ OpenFOAM not installed. Cannot run {case_name} locally.")

def generate_video(case_dir):
    """Generates a video from OpenFOAM results using VTK (PyVista-free)."""
    import imageio
    import numpy as np

    from vtkmodules.vtkIOGeometry import vtkOpenFOAMReader
    from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline
    from vtkmodules.vtkFiltersCore import vtkContourFilter, vtkCellDataToPointData
    from vtkmodules.vtkFiltersModeling import vtkOutlineFilter
    import vtkmodules.vtkRenderingOpenGL2
    from vtkmodules.vtkFiltersGeometry import vtkDataSetSurfaceFilter
    from vtkmodules.vtkCommonTransforms import vtkTransform
    from vtkmodules.vtkFiltersGeneral import vtkTransformFilter
    from vtkmodules.vtkRenderingCore import vtkTextActor
    from vtkmodules.vtkRenderingCore import (
        vtkActor,
        vtkPolyDataMapper,
        vtkRenderer,
        vtkRenderWindow,
        vtkWindowToImageFilter,
    )
    from vtkmodules.vtkCommonColor import vtkNamedColors
    from vtkmodules.util import numpy_support

    print(f"  🎬 Generating video for {case_dir} using PyVista...")
    
    foam_file = os.path.join(case_dir, "case.foam")
    if not os.path.exists(foam_file):
        # Create empty .foam file if it doesn't exist (PyVista needs it)
        with open(foam_file, 'w') as f:
            pass
            
    try:
        reader = vtkOpenFOAMReader()
        reader.SetFileName(foam_file)
        reader.UpdateInformation()
    except Exception as e:
        print(f"  ❌ Error loading OpenFOAM case: {e}")
        return False

    # Get time values from pipeline information
    info = reader.GetOutputInformation(0)
    if info and info.Has(vtkStreamingDemandDrivenPipeline.TIME_STEPS()):
        time_values = list(info.Get(vtkStreamingDemandDrivenPipeline.TIME_STEPS()))
    else:
        time_values = [0.0]
    print(f"  Found {len(time_values)} timesteps.")
    duration = max(time_values[-1] - time_values[0], 0.0) if time_values else 0.0
    if duration <= 0.0:
        duration = DEFAULTS["duration"]
    n_frames = int(round(duration * VIDEO_FPS)) + 1
    frame_times = [time_values[0] + i * (duration / max(n_frames - 1, 1)) for i in range(n_frames)]
    
    # Setup Output in CASE folder
    results_dir = os.path.join(case_dir, "postProcessing")
    os.makedirs(results_dir, exist_ok=True)

    # 1. Generate 3D Moving Mesh Video
    print("    - Generating 3D perspective video...")
    colors = vtkNamedColors()
    renderer = vtkRenderer()
    renderer.SetBackground(colors.GetColor3d("White"))
    render_window = vtkRenderWindow()
    offscreen = os.environ.get("SLOSHING_OFFSCREEN", "0") == "1"
    render_window.SetOffScreenRendering(1 if offscreen else 0)
    render_window.SetSize(1280, 720)
    render_window.AddRenderer(renderer)

    # Determine bounds from the first timestep for a stable camera
    internal_mesh0 = None
    if time_values:
        info.Set(vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(), time_values[0])
        reader.Update()
        output0 = reader.GetOutput()
        if output0 and output0.GetNumberOfBlocks() > 0:
            internal_mesh0 = output0.GetBlock(0)

    if internal_mesh0 is not None:
        xmin, xmax, ymin, ymax, zmin, zmax = internal_mesh0.GetBounds()
        cx = 0.5 * (xmin + xmax)
        cy = 0.5 * (ymin + ymax)
        cz = 0.5 * (zmin + zmax)
        span = max(xmax - xmin, ymax - ymin, zmax - zmin, 1e-6)
        camera = renderer.GetActiveCamera()
        camera.SetPosition(cx, cy - 2.5 * span, cz + 1.2 * span)
        camera.SetFocalPoint(cx, cy, cz)
        camera.SetViewUp(0.0, 0.0, 1.0)
        renderer.ResetCameraClippingRange()

    # Lab-frame motion parameters from case name
    import re
    match = re.search(r'_R([\d.]+)_f([\d.]+)', os.path.basename(case_dir))
    orbital_radius = float(match.group(1)) if match else 0.0
    omega = 2.0 * math.pi * float(match.group(2)) if match else 0.0

    # Time ticker
    time_actor = vtkTextActor()
    time_actor.GetTextProperty().SetFontSize(20)
    time_actor.GetTextProperty().SetColor(0.0, 0.0, 0.0)
    time_actor.SetDisplayPosition(20, 20)
    renderer.AddActor2D(time_actor)
    
    # Use a distinct name to avoid confusion with old runs
    video_filename = "video_3d_render.mp4"
    video_path_3d = os.path.join(results_dir, video_filename)
    print(f"    - Target video path: {os.path.abspath(video_path_3d)}")
    
    try:
        with imageio.get_writer(video_path_3d, fps=VIDEO_FPS, macro_block_size=None) as writer:
            for i, t in enumerate(frame_times):
                info.Set(vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(), t)
                reader.Update()
                output = reader.GetOutput()

                if output and output.GetNumberOfBlocks() > 0:
                    internal_mesh = output.GetBlock(0)
                    if internal_mesh is None:
                        continue

                    renderer.RemoveAllViewProps()
                    renderer.AddActor2D(time_actor)

                    # Lab-frame translation of tank and fluid
                    tx = orbital_radius * math.cos(omega * t)
                    ty = orbital_radius * math.sin(omega * t)
                    transform = vtkTransform()
                    transform.Translate(tx, ty, 0.0)
                    tf = vtkTransformFilter()
                    tf.SetInputData(internal_mesh)
                    tf.SetTransform(transform)
                    tf.Update()
                    internal_mesh_tf = tf.GetOutput()

                    # Convert cell data to point data for contouring
                    cell_to_point = vtkCellDataToPointData()
                    cell_to_point.SetInputData(internal_mesh_tf)
                    cell_to_point.Update()
                    mesh_point = cell_to_point.GetOutput()

                    # Water surface isosurface (alpha.water = 0.5)
                    if mesh_point.GetPointData().HasArray("alpha.water"):
                        mesh_point.GetPointData().SetActiveScalars("alpha.water")
                        contour = vtkContourFilter()
                        contour.SetInputData(mesh_point)
                        contour.SetValue(0, 0.5)
                        contour.Update()

                        contour_mapper = vtkPolyDataMapper()
                        contour_mapper.SetInputConnection(contour.GetOutputPort())
                        contour_actor = vtkActor()
                        contour_actor.SetMapper(contour_mapper)
                        contour_actor.GetProperty().SetColor(colors.GetColor3d("DeepSkyBlue"))
                        contour_actor.GetProperty().SetSpecular(0.5)
                        renderer.AddActor(contour_actor)

                    # Tank walls (surface)
                    surface = vtkDataSetSurfaceFilter()
                    surface.SetInputData(internal_mesh_tf)
                    surface.Update()
                    surface_mapper = vtkPolyDataMapper()
                    surface_mapper.SetInputConnection(surface.GetOutputPort())
                    surface_actor = vtkActor()
                    surface_actor.SetMapper(surface_mapper)
                    surface_actor.GetProperty().SetColor(colors.GetColor3d("Black"))
                    surface_actor.GetProperty().SetOpacity(0.15)
                    surface_actor.GetProperty().SetRepresentationToWireframe()
                    renderer.AddActor(surface_actor)

                    # Time ticker
                    time_actor.SetInput(f"t = {t:.2f} s")

                    render_window.Render()
                    w2i = vtkWindowToImageFilter()
                    w2i.SetInput(render_window)
                    w2i.Update()
                    vtk_image = w2i.GetOutput()
                    width, height, _ = vtk_image.GetDimensions()
                    arr = numpy_support.vtk_to_numpy(vtk_image.GetPointData().GetScalars())
                    img = arr.reshape(height, width, -1)
                    img = np.flipud(img)
                    writer.append_data(img)
                    
                if (i+1) % 20 == 0:
                     print(f"      Rendered 3D frame {i+1}/{len(frame_times)}")
        print(f"      ✅ Saved: {video_filename}")
    except Exception as e:
        print(f"      ❌ Error saving 3D video: {e}")

    # 2. Generate Dashboard Video
    # We use a helper from potential_flow to avoid duplicate code
    sys.path.insert(0, 'utils')
    try:
        from potential_flow import generate_dashboard_animation
        csv_path = os.path.join(results_dir, "interface", "wall_elevation.csv")
        if os.path.exists(csv_path):
            print("    - Generating dashboard analysis video...")
            # Detect R for plotting
            import re
            match = re.search(r'_D([\d.]+)_', os.path.basename(case_dir))
            R_val = float(match.group(1))/2.0 if match else 0.1
            
            generate_dashboard_animation(csv_path, case_dir, R_val, duration=duration, fps=VIDEO_FPS)
            # Find and rename the file generated by the helper
            dash_src = os.path.join(case_dir, "postProcessing", "potential_flow", "potential_flow_dashboard.mp4")
            dash_dst = os.path.join(results_dir, "animation_dashboard_openfoam.mp4")
            # Potential flow helper saves to potential_flow subfolder, we move it up
            if os.path.exists(dash_src):
                if os.path.exists(dash_dst): os.remove(dash_dst)
                os.rename(dash_src, dash_dst)
                print(f"      ✅ Saved: animation_dashboard_openfoam.mp4")
    except Exception as e:
        print(f"      ⚠️  Could not generate dashboard: {e}")

    return True

def generate_lateral_video(case_dir):
    """Generates a lateral view video of the simulation using PyVista."""
    pv = _import_pyvista()
    import imageio
    import numpy as np

    print(f"  🎬 Generating lateral video for {case_dir}...")
    
    foam_file = os.path.join(case_dir, "case.foam")
    if not os.path.exists(foam_file):
        with open(foam_file, 'w') as f:
            pass
            
    try:
        reader = pv.POpenFOAMReader(foam_file)
    except Exception as e:
        print(f"  ❌ Error loading OpenFOAM case: {e}")
        return False

    time_values = reader.time_values
    if not time_values:
        print("  ⚠️  No timesteps found.")
        return False
        
    # Output path
    results_dir = os.path.join(case_dir, "postProcessing")
    os.makedirs(results_dir, exist_ok=True)
    video_path = os.path.join(results_dir, "video_lateral.mp4")
    
    # Setup Plotter
    # Use off_screen=True to support headless environments (requires xvfb if no display)
    off_screen = True 
    plotter = pv.Plotter(off_screen=off_screen, window_size=(1024, 768))
    plotter.set_background('white')
    
    # Determine bounds from the last timestep (usually fully developed)
    reader.set_active_time_value(time_values[-1])
    mesh = reader.read()
    if mesh.n_blocks > 0:
        internal_mesh = mesh[0]
    else:
        internal_mesh = mesh
        
    if internal_mesh is None:
        print("  ❌ Could not read mesh.")
        return False
        
    bounds = internal_mesh.bounds # (xmin, xmax, ymin, ymax, zmin, zmax)
    center = [(bounds[0]+bounds[1])/2, (bounds[2]+bounds[3])/2, (bounds[4]+bounds[5])/2]
    
    # Camera setup: Side view (Look from -Y towards +Y)
    # Position: x=center, y=negative large, z=center
    span = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
    camera_pos = [center[0], bounds[2] - span*2.5, center[2]]
    
    plotter.camera.position = camera_pos
    plotter.camera.focal_point = center
    plotter.camera.up = (0, 0, 1)
    
    print(f"    - Target video path: {os.path.abspath(video_path)}")
    print(f"    - Processing {len(time_values)} frames...")
    
    try:
        with imageio.get_writer(video_path, fps=VIDEO_FPS) as writer:
            for i, t in enumerate(time_values):
                plotter.clear()
                # Re-add light and background? PyVista clears actors but keeps renderer settings usually.
                # Background is property of renderer.
                
                reader.set_active_time_value(t)
                mesh = reader.read()
                if mesh.n_blocks > 0:
                    internal_mesh = mesh[0]
                else:
                    internal_mesh = mesh
                
                # 1. Interface (alpha.water = 0.5)
                if 'alpha.water' in internal_mesh.cell_data:
                    mesh_point = internal_mesh.cell_data_to_point_data()
                    try:
                        iso = mesh_point.contour(isosurfaces=[0.5], scalars='alpha.water')
                        plotter.add_mesh(iso, color='blue', opacity=0.9, smooth_shading=True, lighting=True)
                    except:
                        pass # Isosurface might fail if alpha is all 0 or 1
                
                # 2. Tank outline (wireframe)
                plotter.add_mesh(internal_mesh.outline(), color='black')
                
                # 3. Time Text
                plotter.add_text(f"t = {t:.2f} s", position='upper_left', color='black', font_size=12)
                
                # Render
                img = plotter.screenshot(return_img=True)
                writer.append_data(img)
                
                if (i+1) % 20 == 0:
                    print(f"      Rendered {i+1}/{len(time_values)}")
        print(f"      ✅ Saved: {os.path.basename(video_path)}")
    except Exception as e:
        print(f"      ❌ Error generating video: {e}") 
        return False
        
    return True
        
def extract_interface(case_dir):
    """Extracts the water-air interface (alpha.water=0.5) using PyVista."""
    pv = _import_pyvista()
    import numpy as np
    
    print(f"  📊 Extracting interface for {case_dir} using PyVista...")
    
    foam_file = os.path.join(case_dir, "case.foam")
    if not os.path.exists(foam_file):
        with open(foam_file, 'w') as f:
            pass
            
    try:
        reader = pv.POpenFOAMReader(foam_file)
    except Exception as e:
         print(f"  ❌ Error loading OpenFOAM case: {e}")
         return False

    time_values = reader.time_values
    
    # Setup Output in CASE folder
    results_dir = os.path.join(case_dir, "postProcessing", "interface")
    os.makedirs(results_dir, exist_ok=True)
    
    csv_summary = ["time,max_z,min_z,mean_z,num_points"]
    csv_wall = ["time,theta,zeta_wall"] # For dashboard
    
    # Parse R from case name
    import re
    match = re.search(r'_D([\d.]+)_', os.path.basename(case_dir))
    R_target = float(match.group(1))/2.0 if match else 0.1

    print(f"  Processing {len(time_values)} timesteps (R={R_target})...")
    
    for i, t in enumerate(time_values):
        reader.set_active_time_value(t)
        mesh = reader.read()
        
        if mesh.n_blocks > 0:
            internal_mesh = mesh[0]
            if 'alpha.water' in internal_mesh.cell_data:
                mesh_point = internal_mesh.cell_data_to_point_data()
                try:
                    isosurface = mesh_point.contour(isosurfaces=[0.5], scalars='alpha.water')
                    
                    # Save VTP
                    vtp_file = os.path.join(results_dir, f'interface_t{t:.6f}.vtp')
                    isosurface.save(vtp_file)
                    
                    if isosurface.n_points > 0:
                        pts = isosurface.points
                        z_coords = pts[:, 2]
                        # Aggregate Stats
                        csv_summary.append(f"{t},{np.max(z_coords)},{np.min(z_coords)},{np.mean(z_coords)},{len(pts)}")
                        
                        # Extract Wall elevation profile for dashboard
                        # We project points to (r, theta) and pick points near r=R
                        r = np.sqrt(pts[:,0]**2 + pts[:,1]**2)
                        # Find points near the wall (within 2% margin)
                        wall_mask = r > (R_target * 0.98)
                        if np.any(wall_mask):
                            wall_pts = pts[wall_mask]
                            wall_thetas = np.arctan2(wall_pts[:,1], wall_pts[:,0])
                            # Bin by theta to get a clean profile
                            n_bins = 64
                            bins = np.linspace(-np.pi, np.pi, n_bins+1)
                            for b in range(n_bins):
                                bin_mask = (wall_thetas >= bins[b]) & (wall_thetas < bins[b+1])
                                if np.any(bin_mask):
                                    z_bin = np.mean(wall_pts[bin_mask, 2])
                                    theta_bin = (bins[b] + bins[b+1])/2.0
                                    csv_wall.append(f"{t},{theta_bin},{z_bin}")
                    else:
                        csv_summary.append(f"{t},0,0,0,0")
                except:
                    csv_summary.append(f"{t},0,0,0,0")
            else:
                csv_summary.append(f"{t},0,0,0,0")
        else:
            csv_summary.append(f"{t},0,0,0,0")
            
        if (i+1) % 20 == 0:
            print(f"    Processed {i+1}/{len(time_values)}")
            
    # Save CSVs
    with open(os.path.join(results_dir, 'interface_summary.csv'), 'w') as f:
        f.write('\n'.join(csv_summary))
    with open(os.path.join(results_dir, 'wall_elevation.csv'), 'w') as f:
        f.write('\n'.join(csv_wall))
        
    print(f"  ✅ Extraction complete.")
    return True

def _read_scalar_value(path, key, default=None):
    if not os.path.exists(path):
        return default
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    match = re.search(rf'^\s*{re.escape(key)}\s+([-+0-9.eE]+);', content, re.M)
    if not match:
        return default
    try:
        return float(match.group(1))
    except ValueError:
        return default

def _read_g_vector(path):
    if not os.path.exists(path):
        return (0.0, 0.0, -9.81)
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    match = re.search(r'^\s*value\s+\(([^)]+)\);', content, re.M)
    if not match:
        return (0.0, 0.0, -9.81)
    parts = match.group(1).split()
    if len(parts) != 3:
        return (0.0, 0.0, -9.81)
    try:
        return (float(parts[0]), float(parts[1]), float(parts[2]))
    except ValueError:
        return (0.0, 0.0, -9.81)

def _read_contact_angle(path):
    return _read_scalar_value(path, "theta0", 90.0)

def _save_points_csv(path, points):
    header = "x,y,z"
    lines = [header]
    for x, y, z in points:
        lines.append(f"{x},{y},{z}")
    with open(path, "w") as f:
        f.write("\n".join(lines))

def _extract_openfoam_interface_latest(case_dir, results_dir):
    pv = _import_pyvista()
    import numpy as np

    foam_file = os.path.join(case_dir, "case.foam")
    if not os.path.exists(foam_file):
        with open(foam_file, "w") as f:
            pass

    try:
        reader = pv.POpenFOAMReader(foam_file)
    except Exception as e:
        print(f"  ❌ Error loading OpenFOAM case: {e}")
        return None, None

    time_values = reader.time_values
    if not time_values:
        print("  ⚠️  No OpenFOAM time values found.")
        return None, None

    t = max(time_values)
    reader.set_active_time_value(t)
    mesh = reader.read()
    if hasattr(mesh, "n_blocks") and mesh.n_blocks > 0:
        internal_mesh = mesh[0]
    else:
        internal_mesh = mesh
    if internal_mesh is None:
        print("  ⚠️  OpenFOAM mesh missing.")
        return None, t

    if "alpha.water" not in internal_mesh.cell_data:
        print("  ⚠️  OpenFOAM alpha.water not found.")
        return None, t

    mesh_point = internal_mesh.cell_data_to_point_data()
    try:
        isosurface = mesh_point.contour(isosurfaces=[0.5], scalars="alpha.water")
    except Exception as e:
        print(f"  ❌ Error extracting iso-surface: {e}")
        return None, t

    if isosurface.n_points == 0:
        print("  ⚠️  OpenFOAM iso-surface has no points.")
        return None, t

    vtp_file = os.path.join(results_dir, f"openfoam_interface_t{t:.6f}.vtp")
    isosurface.save(vtp_file)
    csv_file = os.path.join(results_dir, f"openfoam_interface_t{t:.6f}.csv")
    _save_points_csv(csv_file, isosurface.points)
    return isosurface.points, t

def _extract_analytical_interface(case_dir, results_dir):
    import numpy as np
    from yl_nonlin import yl_nonlin

    params = parse_case_params(os.path.basename(case_dir))
    H = params.get("H", DEFAULTS["H"])
    R = params.get("D", DEFAULTS["D"]) / 2.0

    rho = _read_scalar_value(os.path.join(case_dir, "constant", "physicalProperties.water"), "rho", 1000.0)
    sigma = _read_scalar_value(os.path.join(case_dir, "constant", "phaseProperties"), "sigma", 0.072)
    gx, gy, gz = _read_g_vector(os.path.join(case_dir, "constant", "g"))
    g_horizontal = math.hypot(gx, gy)
    g_vertical = abs(gz) if abs(gz) > 1e-12 else math.sqrt(gx * gx + gy * gy + gz * gz)
    if g_vertical <= 0:
        g_vertical = 9.81
    F = g_horizontal / g_vertical if g_horizontal > 0 else 0.0

    thetac_deg = _read_contact_angle(os.path.join(case_dir, "0", "alpha.water"))

    a_orbit = 1.0
    if F > 0:
        omega = math.sqrt(F * g_vertical / a_orbit)
        omega_rpm = omega * 60.0 / (2.0 * math.pi)
    else:
        omega_rpm = 0.0

    area, hL, pts = yl_nonlin(
        rho,
        sigma,
        g_vertical,
        omega_rpm,
        a_orbit,
        R,
        thetac_deg,
        hmax=0.02 * R,
    )

    if g_horizontal > 0:
        phi = math.atan2(gy, gx)
        x = pts[:, 0]
        y = pts[:, 1]
        x_rot = x * math.cos(phi) - y * math.sin(phi)
        y_rot = x * math.sin(phi) + y * math.cos(phi)
        pts = np.column_stack((x_rot, y_rot, pts[:, 2]))

    if H:
        z_mean = float(np.mean(pts[:, 2]))
        pts[:, 2] += (0.5 * H - z_mean)

    csv_file = os.path.join(results_dir, "analytical_interface.csv")
    _save_points_csv(csv_file, pts)
    try:
        pv = _import_pyvista()
        poly = pv.PolyData(pts)
        poly.save(os.path.join(results_dir, "analytical_interface.vtp"))
    except Exception:
        pass

    return pts, area, hL

def _compute_l2_between_interfaces(analytic_pts, openfoam_pts):
    import numpy as np
    from scipy.interpolate import griddata

    if analytic_pts is None or openfoam_pts is None:
        return None, 0
    if len(analytic_pts) == 0 or len(openfoam_pts) == 0:
        return None, 0

    xy_sim = openfoam_pts[:, :2]
    z_sim = openfoam_pts[:, 2]
    xy_ref = analytic_pts[:, :2]
    z_ref = analytic_pts[:, 2]

    z_interp = griddata(xy_sim, z_sim, xy_ref, method="linear")
    if np.any(np.isnan(z_interp)):
        z_near = griddata(xy_sim, z_sim, xy_ref, method="nearest")
        z_interp = np.where(np.isnan(z_interp), z_near, z_interp)

    valid = np.isfinite(z_interp)
    if not np.any(valid):
        return None, 0

    diff = z_interp[valid] - z_ref[valid]
    l2_rms = float(np.sqrt(np.mean(diff * diff)))
    l2_sum = float(np.sqrt(np.sum(diff * diff)))
    return {"l2_rms": l2_rms, "l2_sum": l2_sum}, int(np.sum(valid))

def compare_interfaces(case_dir):
    print(f"  📊 Comparing analytical and OpenFOAM interfaces for {case_dir}...")
    results_dir = os.path.join(case_dir, "postProcessing", "interface_compare")
    os.makedirs(results_dir, exist_ok=True)

    analytic_pts, area, hL = _extract_analytical_interface(case_dir, results_dir)
    openfoam_pts, t = _extract_openfoam_interface_latest(case_dir, results_dir)

    l2_info, n_samples = _compute_l2_between_interfaces(analytic_pts, openfoam_pts)
    l2_rms_cap = ""
    if l2_info and l2_info.get("l2_rms") is not None:
        rho = _read_scalar_value(os.path.join(case_dir, "constant", "physicalProperties.water"), "rho", 1000.0)
        sigma = _read_scalar_value(os.path.join(case_dir, "constant", "phaseProperties"), "sigma", 0.072)
        gx, gy, gz = _read_g_vector(os.path.join(case_dir, "constant", "g"))
        g_vertical = abs(gz) if abs(gz) > 1e-12 else math.sqrt(gx * gx + gy * gy + gz * gz)
        if g_vertical <= 0:
            g_vertical = 9.81
        lc = math.sqrt(sigma / (rho * g_vertical)) if rho > 0 else 0.0
        if lc > 0:
            l2_rms_cap = l2_info["l2_rms"] / lc

    summary = [
        "metric,value",
        f"openfoam_time,{t if t is not None else ''}",
        f"analytical_area,{area}",
        f"analytical_hL,{hL}",
        f"l2_rms,{l2_info['l2_rms'] if l2_info else ''}",
        f"l2_rms_capillary,{l2_rms_cap}",
        f"l2_sum,{l2_info['l2_sum'] if l2_info else ''}",
        f"num_samples,{n_samples}",
    ]
    with open(os.path.join(results_dir, "comparison_summary.csv"), "w") as f:
        f.write("\n".join(summary))

    if openfoam_pts is None:
        print("  ⚠️  OpenFOAM interface not found; comparison summary saved anyway.")
    else:
        print("  ✅ Comparison complete.")
    return True

def menu_postprocess(is_oscar):
    """Submenu 3: Postprocess"""
    print("\n" + "="*60)
    print("  POSTPROCESS MENU")
    print("="*60)
    
    cases = sorted([d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('case_')])
    if not cases:
        print("No cases found.")
        return
    
    # Display cases
    print("\nAvailable Cases:")
    for i, c in enumerate(cases):
        status = "(DONE)" if is_case_done(c, DEFAULTS['duration']) else ""
        print(f"  {i+1}) {c} {status}")
    
    print("\n" + "-"*60)
    print("Select Action:")
    print("  1) Compare Interfaces (Analytical vs OpenFOAM)")
    print("  2) Generate Lateral View Video")
    print("  Q) Back to Main Menu")
    print("-"*60)
    
    choice = input("\nAction: ").strip().lower()
    
    if choice == '1':
        print("\n→ Interface Comparison (Analytical vs OpenFOAM)")
        idx_str = input("  Enter case numbers (e.g., 1, 3-5, all): ").strip().lower()
        if idx_str == 'all':
            indices = list(range(len(cases)))
        else:
            indices = parse_indices(idx_str, len(cases))
        
        if not indices:
            print("No valid indices selected.")
            return
        
        print(f"\nComparing interfaces for {len(indices)} case(s)...")
        for i in indices:
            if is_oscar:
                if i == indices[0]:
                    submit = input("\n⚠️  Post-processing detected. Submit as Slurm job? (y/n): ").strip().lower()
                    if submit == 'y':
                        for idx in indices:
                            run_postprocess_oscar(cases[idx], "compare")
                        return
            compare_interfaces(cases[i])
    
    elif choice == '2':
        print("\n→ Generate Lateral View Video")
        idx_str = input("  Enter case numbers (e.g., 1, 3-5, all): ").strip().lower()
        if idx_str == 'all':
            indices = list(range(len(cases)))
        else:
            indices = parse_indices(idx_str, len(cases))
        
        if not indices:
            print("No valid indices selected.")
            return
            
        print(f"\nGenerating videos for {len(indices)} case(s)...")
        for i in indices:
            if is_oscar:
                if i == indices[0]:
                    submit = input("\n⚠️  Post-processing detected. Submit as Slurm job? (y/n): ").strip().lower()
                    if submit == 'y':
                        for idx in indices:
                            run_postprocess_oscar(cases[idx], "video")
                        return
            generate_lateral_video(cases[i])
    elif choice == 'q':
        return

def run_postprocess_oscar(case_name, action):
    """Submits a post-processing job to Slurm."""
    slurm_dir = os.path.join(case_name, "slurm")
    os.makedirs(slurm_dir, exist_ok=True)
    script_path = os.path.join(slurm_dir, f"postprocess_{action}.slurm")
    
    header = [
        "#!/usr/bin/env bash",
        f"#SBATCH -J post_{action}_{case_name}",
        "#SBATCH -p batch",
        "#SBATCH -N 1",
        "#SBATCH -n 1",
        "#SBATCH --time=01:00:00",
        "#SBATCH --mem=8G",
        f"#SBATCH -o {slurm_dir}/slurm_postprocessing.log",
        "#SBATCH --open-mode=append",
        "",
        "set -euo pipefail",
        "",
        "# --- Load Consistent Python Module ---",
        "module load python/3.13",
        "",
        "# --- Activate Shared Environment ---",
        "VENV_DIR=\"sloshing\"",
        "if [ ! -d \"$VENV_DIR\" ]; then",
        "  echo \"📦 Venv not found. Creating $VENV_DIR on compute node...\"",
        "  python3 -m venv $VENV_DIR",
        "  source $VENV_DIR/bin/activate",
        "  pip install --upgrade pip",
        "  pip install -r requirements.txt",
        "else",
        "  source $VENV_DIR/bin/activate",
        "fi",
        "# -----------------------------------",
        "",
        "echo '------------------------------------------------------------'",
        f"echo 'Action: {action} | Case: {case_name}'",
        f"echo 'Date: $(date)'",
        f"echo 'Python: $(which python)'",  # Debug print
        "export SLOSHING_OFFSCREEN=1",
        "export VTK_DEFAULT_RENDER_WINDOW_OFFSCREEN=1",
        f"cd {case_name}",
        "if command -v xvfb-run >/dev/null 2>&1; then",
        f"  xvfb-run -s \"-screen 0 1280x720x24\" python ../main.py --headless --case . --action {action}",
        "else",
        f"  python ../main.py --headless --case . --action {action}",
        "fi",
        "echo 'End: $(date)'",
        "echo '------------------------------------------------------------'",
        ""
    ]
    
    os.makedirs(slurm_dir, exist_ok=True)
    
    with open(script_path, "w") as f:
        f.write("\n".join(header))
    
    print(f"  🚀 Submitting post-processing job for {case_name} ({action})...")
    subprocess.run(["sbatch", script_path], check=True)

def main_menu():
    """Main entry point."""
    print("\n" + "="*40)
    print("   Tilting Tank Manager")
    print("="*40)
    
    oscar_input = input("Are you on Oscar? (y/n): ").strip().lower()
    is_oscar = oscar_input == 'y'
    
    while True:
        print("\n--- Main Menu ---")
        print("1) Build Case Setups")
        print("2) Run Cases")
        print("3) Postprocess Cases")
        print("Q) Quit")
        
        choice = input("\nSelect an option: ").strip().lower()
        
        if choice == '1':
            menu_build_cases(is_oscar)
        elif choice == '2':
            menu_run_cases(is_oscar)
        elif choice == '3':
            menu_postprocess(is_oscar)
        elif choice == 'q':
            print("Goodbye!")
            break
        else:
            print("Invalid option.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--headless", action="store_true", help="Run without menu")
    parser.add_argument("--case", type=str, help="Case directory for headless mode")
    parser.add_argument("--action", type=str, choices=["compare", "video"], help="Action for headless mode")
    
    args = parser.parse_args()
    
    if args.headless:
        if not args.case or not args.action:
            print("Error: --case and --action are required in headless mode.")
            sys.exit(1)
        
        if args.action == "compare":
            compare_interfaces(args.case)
        elif args.action == "video":
            generate_lateral_video(args.case)
    else:
        main_menu()
