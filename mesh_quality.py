#!/usr/bin/env python3
"""
Lightweight Gmsh MSH v2.2 quality checks.

Goal: catch tiny elements (min edge << target lc) which force tiny deltaT and
exploding runtime. This avoids relying on OpenFOAM `checkMesh` availability.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, asdict


@dataclass
class MeshQualitySummary:
    msh_path: str
    n_nodes: int
    n_elements: int
    n_tets: int
    min_edge: float | None
    max_edge: float | None
    mean_edge: float | None
    min_volume: float | None
    max_aspect_ratio: float | None


def _tet_volume(a, b, c, d) -> float:
    ax, ay, az = a
    bx, by, bz = b
    cx, cy, cz = c
    dx, dy, dz = d
    ab = (bx - ax, by - ay, bz - az)
    ac = (cx - ax, cy - ay, cz - az)
    ad = (dx - ax, dy - ay, dz - az)
    det = (
        ab[0] * (ac[1] * ad[2] - ac[2] * ad[1])
        - ab[1] * (ac[0] * ad[2] - ac[2] * ad[0])
        + ab[2] * (ac[0] * ad[1] - ac[1] * ad[0])
    )
    return abs(det) / 6.0


def _dist(p, q) -> float:
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    dz = p[2] - q[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def analyze_msh2(msh_path: str) -> MeshQualitySummary:
    nodes: dict[int, tuple[float, float, float]] = {}
    tets: list[tuple[int, int, int, int]] = []
    n_elements = 0

    with open(msh_path, "r") as f:
        lines = iter(f)
        for raw in lines:
            line = raw.strip()
            if line == "$Nodes":
                n_nodes = int(next(lines).strip())
                for _ in range(n_nodes):
                    parts = next(lines).split()
                    idx = int(parts[0])
                    nodes[idx] = (float(parts[1]), float(parts[2]), float(parts[3]))
                # consume end
                while True:
                    if next(lines).strip() == "$EndNodes":
                        break
            elif line == "$Elements":
                n_elements = int(next(lines).strip())
                for _ in range(n_elements):
                    parts = next(lines).split()
                    # id, type, ntags, tags..., nodeIds...
                    etype = int(parts[1])
                    ntags = int(parts[2])
                    node_ids = [int(x) for x in parts[3 + ntags :]]
                    # 4-node tetrahedron in msh2 is type 4
                    if etype == 4 and len(node_ids) == 4:
                        tets.append((node_ids[0], node_ids[1], node_ids[2], node_ids[3]))
                while True:
                    if next(lines).strip() == "$EndElements":
                        break

    min_edge = None
    max_edge = None
    sum_edge = 0.0
    count_edge = 0
    min_volume = None
    max_ar = None

    for n1, n2, n3, n4 in tets:
        p1 = nodes.get(n1)
        p2 = nodes.get(n2)
        p3 = nodes.get(n3)
        p4 = nodes.get(n4)
        if not (p1 and p2 and p3 and p4):
            continue
        edges = (
            _dist(p1, p2),
            _dist(p1, p3),
            _dist(p1, p4),
            _dist(p2, p3),
            _dist(p2, p4),
            _dist(p3, p4),
        )
        e_min = min(edges)
        e_max = max(edges)
        vol = _tet_volume(p1, p2, p3, p4)
        ar = (e_max / e_min) if e_min > 0 else float("inf")

        if min_edge is None or e_min < min_edge:
            min_edge = e_min
        if max_edge is None or e_max > max_edge:
            max_edge = e_max
        if min_volume is None or vol < min_volume:
            min_volume = vol
        if max_ar is None or ar > max_ar:
            max_ar = ar
        sum_edge += sum(edges)
        count_edge += 6

    mean_edge = (sum_edge / count_edge) if count_edge else None

    return MeshQualitySummary(
        msh_path=os.path.abspath(msh_path),
        n_nodes=len(nodes),
        n_elements=n_elements,
        n_tets=len(tets),
        min_edge=min_edge,
        max_edge=max_edge,
        mean_edge=mean_edge,
        min_volume=min_volume,
        max_aspect_ratio=max_ar,
    )


def write_summary(summary: MeshQualitySummary, out_path: str) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(asdict(summary), f, indent=2, sort_keys=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("msh", help="Path to Gmsh .msh (v2.2) file")
    parser.add_argument("--out", help="Write JSON summary to this path")
    args = parser.parse_args()

    s = analyze_msh2(args.msh)
    print(json.dumps(asdict(s), indent=2, sort_keys=True))
    if args.out:
        write_summary(s, args.out)
