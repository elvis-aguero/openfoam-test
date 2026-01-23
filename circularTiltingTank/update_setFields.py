
import sys

def update_setfields(h):
    content = fr"""/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      setFieldsDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defaultValues
{{
    alpha.water 0;
}}

zones
{{
    water
    {{
        type box;
        box (-100 -100 -1.0) (100 100 {h/2});
        values
        {{
            alpha.water 1;
        }}
    }}
}}

// ************************************************************************* //
"""
    with open("system/setFieldsDict", "w") as f:
        f.write(content)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 update_setFields.py <H>")
        sys.exit(1)
        
    h = float(sys.argv[1])
    update_setfields(h)
    print(f"system/setFieldsDict updated with water level {h/2}")
