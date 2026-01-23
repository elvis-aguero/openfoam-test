
import sys

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

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 generate_mesh.py <Height> <Diameter> <MeshSize> [GeometryType: flat/cap]")
        sys.exit(1)
        
    H = float(sys.argv[1])
    D = float(sys.argv[2])
    lc = float(sys.argv[3])
    
    # Default to flat
    geo_type = "flat"
    if len(sys.argv) >= 5:
        geo_type = sys.argv[4]
    
    if geo_type == "cap":
        geo_content = generate_geo_cap(H, D, lc)
        print(f"Generating Spherical Cap Geometry (D={D}, H={H})")
    else:
        geo_content = generate_geo_flat(H, D, lc)
        print(f"Generating Flat Bottom Geometry (D={D}, H={H})")
    
    with open("cylinder.geo", "w") as f:
        f.write(geo_content)
    
    print("cylinder.geo created.")
