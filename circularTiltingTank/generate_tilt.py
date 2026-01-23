import math
import sys

G = 9.81


def write_gravity(theta_deg, output_file="constant/g"):
    theta = math.radians(theta_deg)
    gy = G * math.sin(theta)
    gz = -G * math.cos(theta)
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
    class       uniformDimensionedVectorField;
    location    "constant";
    object      g;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -2 0 0 0 0];
value           (0 {gy:.6g} {gz:.6g});


// ************************************************************************* //
"""
    with open(output_file, "w") as f:
        f.write(content)


def write_static_motion(duration, dt, output_file="constant/6DoF.dat"):
    n_steps = int(duration / dt) + 1
    with open(output_file, "w") as file:
        file.write(f"{n_steps}\n(\n")
        for i in range(n_steps):
            ti = i * dt
            file.write(f"({ti:.6g} (0 0 0) (0 0 0))\n")
        file.write(")\n")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 generate_tilt.py <theta_deg> <duration> <dt>")
        sys.exit(1)

    theta_deg = float(sys.argv[1])
    duration = float(sys.argv[2])
    dt = float(sys.argv[3])

    write_gravity(theta_deg)
    write_static_motion(duration, dt)
    print(f"Tilted gravity written to constant/g with theta={theta_deg} deg")
    print(f"Static motion file written to constant/6DoF.dat for duration={duration}s, dt={dt}s")
