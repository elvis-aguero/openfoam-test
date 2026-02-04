
import sys
import os
import math

def generate_geo_flat(H, D, mesh_size):
    R = D / 2.0
    
    geo_content = f"""
// Inputs (Flat Bottom)
H = {H};
R = {R};
lc = {mesh_size};

// Geometry using OpenCASCADE kernel
SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;

// Cylinder
v1 = newv;
Cylinder(v1) = {{0, 0, 0, 0, 0, H, R}};

// Mesh Quality
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;
// Hard-disable automatic size variation that can create tiny elements.
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

// Physical Groups
eps = 1e-3;

// Find Top Surface (at z=H)
Surface_Top[] = Surface In BoundingBox {{ -R-eps, -R-eps, H-eps, R+eps, R+eps, H+eps }};
Physical Surface("atmosphere") = Surface_Top[];

// All other surfaces are "walls"
All_Surfs[] = Surface "*";
Wall_Surfs[] = {{}};

For i In {{0 : #All_Surfs[]-1}}
    id = All_Surfs[i];
    is_top = 0;
    For j In {{0 : #Surface_Top[]-1}}
        If (id == Surface_Top[j])
            is_top = 1;
        EndIf
    EndFor
    
    If (is_top == 0)
        Wall_Surfs[] += {{id}};
    EndIf
EndFor

Physical Surface("walls") = Wall_Surfs[];
Physical Volume("internalMesh") = v1;
"""
    return geo_content

def generate_geo_cap(H, D, mesh_size):
    R = D / 2.0
    
    geo_content = f"""
// Inputs (Spherical Cap)
H = {H};
R = {R};
lc = {mesh_size};

SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;

// Cylinder Body from z=0 to z=H
v1 = newv;
Cylinder(v1) = {{0, 0, 0, 0, 0, H, R}};

// Spherical Cap at z=0 (radius R)
v2 = newv;
Sphere(v2) = {{0, 0, 0, R}};

// Fuse them
v_fused[] = BooleanUnion{{ Volume{{v1}}; Delete; }}{{ Volume{{v2}}; Delete; }};

// Mesh Quality
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;
// Hard-disable automatic size variation that can create tiny elements.
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

eps = 1e-3;

// Find Top Surface (at z=H)
Surface_Top[] = Surface In BoundingBox {{ -R-eps, -R-eps, H-eps, R+eps, R+eps, H+eps }};
Physical Surface("atmosphere") = Surface_Top[];

// All other surfaces are "walls"
All_Surfs[] = Surface "*";
Wall_Surfs[] = {{}};

For i In {{0 : #All_Surfs[]-1}}
    id = All_Surfs[i];
    is_top = 0;
    For j In {{0 : #Surface_Top[]-1}}
        If (id == Surface_Top[j])
            is_top = 1;
        EndIf
    EndFor
    
    If (is_top == 0)
        Wall_Surfs[] += {{id}};
    EndIf
EndFor

Physical Surface("walls") = Wall_Surfs[];
Physical Volume("internalMesh") = v_fused[0];
"""
    return geo_content

def _write_ascii_stl(path, solids):
    with open(path, "w", encoding="utf-8") as f:
        for name, triangles in solids:
            f.write(f"solid {name}\n")
            for tri in triangles:
                f.write("  facet normal 0 0 0\n")
                f.write("    outer loop\n")
                for vx, vy, vz in tri:
                    f.write(f"      vertex {vx:.8e} {vy:.8e} {vz:.8e}\n")
                f.write("    endloop\n")
                f.write("  endfacet\n")
            f.write(f"endsolid {name}\n")

def _generate_cylinder_stl(H, D, mesh_size, out_path, geo_type="flat"):
    R = D / 2.0
    circumference = 2.0 * math.pi * R
    n_theta = max(32, int(round(circumference / max(mesh_size, 1e-9))))

    bottom = []
    top = []
    for i in range(n_theta):
        th = 2.0 * math.pi * i / n_theta
        x = R * math.cos(th)
        y = R * math.sin(th)
        bottom.append((x, y, 0.0))
        top.append((x, y, H))

    wall_tris = []
    for i in range(n_theta):
        j = (i + 1) % n_theta
        p0 = bottom[i]
        p1 = bottom[j]
        p2 = top[j]
        p3 = top[i]
        wall_tris.append((p0, p1, p2))
        wall_tris.append((p0, p2, p3))

    bottom_tris = []
    if geo_type == "flat":
        center_bottom = (0.0, 0.0, 0.0)
        for i in range(n_theta):
            j = (i + 1) % n_theta
            bottom_tris.append((center_bottom, bottom[j], bottom[i]))
    else:
        # Spherical cap with sphere centered at origin, radius R (lower hemisphere).
        n_phi = max(16, n_theta // 2)
        phi_vals = [math.pi / 2.0 + (math.pi / 2.0) * i / n_phi for i in range(n_phi + 1)]
        for i in range(n_theta):
            j = (i + 1) % n_theta
            th0 = 2.0 * math.pi * i / n_theta
            th1 = 2.0 * math.pi * j / n_theta
            for k in range(n_phi):
                ph0 = phi_vals[k]
                ph1 = phi_vals[k + 1]
                p00 = (R * math.cos(th0) * math.sin(ph0), R * math.sin(th0) * math.sin(ph0), R * math.cos(ph0))
                p10 = (R * math.cos(th1) * math.sin(ph0), R * math.sin(th1) * math.sin(ph0), R * math.cos(ph0))
                p01 = (R * math.cos(th0) * math.sin(ph1), R * math.sin(th0) * math.sin(ph1), R * math.cos(ph1))
                p11 = (R * math.cos(th1) * math.sin(ph1), R * math.sin(th1) * math.sin(ph1), R * math.cos(ph1))
                bottom_tris.append((p00, p10, p11))
                bottom_tris.append((p00, p11, p01))

    center_top = (0.0, 0.0, H)
    top_tris = []
    for i in range(n_theta):
        j = (i + 1) % n_theta
        top_tris.append((center_top, top[i], top[j]))

    solids = [
        ("walls", wall_tris + bottom_tris),
        ("atmosphere", top_tris),
    ]
    _write_ascii_stl(out_path, solids)

def _write_block_mesh_dict(path, H, D, mesh_size, geo_type="flat"):
    R = D / 2.0
    pad = max(2.0 * mesh_size, 0.01 * R)
    zmin_base = 0.0 if geo_type == "flat" else -R
    xmin, xmax = -R - pad, R + pad
    ymin, ymax = -R - pad, R + pad
    zmin, zmax = zmin_base - pad, H + pad

    # Base mesh is intentionally coarser: one refinement level halves cell size.
    base_dx = max(mesh_size * 2.0, 1e-9)
    nx = max(1, int(round((xmax - xmin) / base_dx)))
    ny = max(1, int(round((ymax - ymin) / base_dx)))
    nz = max(1, int(round((zmax - zmin) / base_dx)))

    content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    ({xmin} {ymin} {zmin})
    ({xmax} {ymin} {zmin})
    ({xmax} {ymax} {zmin})
    ({xmin} {ymax} {zmin})
    ({xmin} {ymin} {zmax})
    ({xmax} {ymin} {zmax})
    ({xmax} {ymax} {zmax})
    ({xmin} {ymax} {zmax})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    background
    {{
        type patch;
        faces
        (
            (0 1 2 3)
            (4 5 6 7)
            (0 1 5 4)
            (1 2 6 5)
            (2 3 7 6)
            (3 0 4 7)
        );
    }}
);

mergePatchPairs
(
);

// ************************************************************************* //
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)

def _write_surface_features_dict(path):
    content = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    object      surfaceFeaturesDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

surfaces ("cylinder.stl");
includedAngle   150;
writeFeatureEdgeMesh yes;

// ************************************************************************* //
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)

def _write_snappy_hex_mesh_dict(path, H, geo_type="flat"):
    zmin_base = 0.0 if geo_type == "flat" else -0.5 * H
    location_in_mesh = (0.0, 0.0, max(zmin_base + 0.1 * H, 0.25 * H))
    content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

castellatedMesh true;
snap            true;
addLayers       false;
mergeTolerance  1e-6;

geometry
{{
    cylinder
    {{
        type triSurfaceMesh;
        file "cylinder.stl";
        name cylinder;
        regions
        {{
            atmosphere {{ name atmosphere; }}
            walls       {{ name walls; }}
        }}
    }}
}};

castellatedMeshControls
{{
    maxLocalCells       2000000;
    maxGlobalCells      20000000;
    minRefinementCells  0;
    nCellsBetweenLevels 2;
    allowFreeStandingZoneFaces true;

    features
    (
        {{
            file "cylinder.eMesh";
            level 1;
        }}
    );

    refinementSurfaces
    {{
        cylinder
        {{
            level (1 1);
            regions
            {{
                walls
                {{
                    level (1 1);
                    patchInfo {{ type wall; }}
                }}
                atmosphere
                {{
                    level (1 1);
                    patchInfo {{ type patch; }}
                }}
            }}
        }}
    }}

    resolveFeatureAngle 30;
    locationInMesh ({location_in_mesh[0]} {location_in_mesh[1]} {location_in_mesh[2]});
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance   2.0;
    nSolveIter  30;
    nRelaxIter  5;
}}

addLayersControls
{{
    relativeSizes       true;
    layers              {{ }};
}}

meshQualityControls
{{
    maxNonOrtho         65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave          80;
    minVol              1e-13;
    minArea             1e-13;
    minTetQuality       1e-9;
    minTwist            0.02;
    minDeterminant      0.001;
    minFaceWeight       0.02;
    minVolRatio         0.01;
    minTriangleTwist    -1;
    nSmoothScale        4;
    errorReduction      0.75;
}}

// ************************************************************************* //
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)

def _write_topo_set_dict(path, H, D, mesh_size, geo_type="flat"):
    R = D / 2.0
    pad = max(2.0 * mesh_size, 0.01 * R)
    zmin_base = 0.0 if geo_type == "flat" else -R
    xmin, xmax = -R - pad, R + pad
    ymin, ymax = -R - pad, R + pad
    zmin, zmax = zmin_base - pad, H + pad
    content = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

actions
(
    {{
        name    allCells;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     ({xmin} {ymin} {zmin}) ({xmax} {ymax} {zmax});
    }}
    {{
        name    internalMesh;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        set     allCells;
    }}
);

// ************************************************************************* //
""".format(xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 generate_mesh.py <Height> <Diameter> <MeshSize> [GeometryType: flat/cap] [Mesher: gmsh/snappy]")
        sys.exit(1)

    H = float(sys.argv[1])
    D = float(sys.argv[2])
    lc = float(sys.argv[3])

    geo_type = "flat"
    if len(sys.argv) >= 5:
        geo_type = sys.argv[4]

    mesher = "gmsh"
    if len(sys.argv) >= 6:
        mesher = sys.argv[5].strip().lower()

    if mesher == "snappy":
        tri_dir = os.path.join("constant", "triSurface")
        os.makedirs(tri_dir, exist_ok=True)
        stl_path = os.path.join(tri_dir, "cylinder.stl")
        _generate_cylinder_stl(H, D, lc, stl_path, geo_type=geo_type)
        _write_block_mesh_dict(os.path.join("system", "blockMeshDict"), H, D, lc, geo_type=geo_type)
        _write_surface_features_dict(os.path.join("system", "surfaceFeaturesDict"))
        _write_snappy_hex_mesh_dict(os.path.join("system", "snappyHexMeshDict"), H, geo_type=geo_type)
        _write_topo_set_dict(os.path.join("system", "topoSetDict"), H, D, lc, geo_type=geo_type)
        print(f"Generated snappyHexMesh inputs at {stl_path}.")
        sys.exit(0)

    if geo_type == "cap":
        geo_content = generate_geo_cap(H, D, lc)
        print(f"Generating Spherical Cap Geometry (D={D}, H={H})")
    else:
        geo_content = generate_geo_flat(H, D, lc)
        print(f"Generating Flat Bottom Geometry (D={D}, H={H})")

    with open("cylinder.geo", "w") as f:
        f.write(geo_content)

    print("cylinder.geo created.")
