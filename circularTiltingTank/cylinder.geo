
// Inputs (Spherical Cap)
H = 1.0;
R = 0.25;
lc = 0.05;

SetFactory("OpenCASCADE");

// Cylinder Body from z=0 to z=H
v1 = newv;
Cylinder(v1) = {0, 0, 0, 0, 0, H, R};

// Spherical Cap at z=0 (radius R)
// This sphere will fuse with the cylinder.
// The top half (z > 0) is inside the cylinder.
// The bottom half (z < 0) forms the hemispherical cap.
v2 = newv;
Sphere(v2) = {0, 0, 0, R};

// Fuse them
v_fused[] = BooleanUnion{ Volume{v1}; Delete; }{ Volume{v2}; Delete; };

// Mesh Quality
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

eps = 1e-3;

// Physical Surfaces

// Top Surface (at z=H)
Surface_Top[] = Surface In BoundingBox { -R-eps, -R-eps, H-eps, R+eps, R+eps, H+eps };
Physical Surface("atmosphere") = Surface_Top[];

// Walls (Everything else)
All_Surfaces[] = Surface "*";
Wall_Surfaces[] = {};

For i In {0 : #All_Surfaces[]-1}
    id = All_Surfaces[i];
    is_top = 0;
    
    For j In {0 : #Surface_Top[]-1}
        If (id == Surface_Top[j])
            is_top = 1;
        EndIf
    EndFor
    
    If (is_top == 0)
        Wall_Surfaces[] += {id};
    EndIf
EndFor

Physical Surface("walls") = Wall_Surfaces[];
Physical Volume("internalMesh") = v_fused[];
