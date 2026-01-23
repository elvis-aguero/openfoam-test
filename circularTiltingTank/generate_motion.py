
import sys
import math

def smootherstep(t):
    """Smootherstep function: 6t^5 - 15t^4 + 10t^3"""
    return t * t * t * (t * (t * 6 - 15) + 10)

def generate_motion(r_max, f, duration, dt, ramp_duration, output_file):
    # Number of steps
    n_steps = int(duration / dt) + 1
    
    with open(output_file, 'w') as file:
        file.write(f"{n_steps}\n(\n")
        
        for i in range(n_steps):
            ti = i * dt
            
            # Calculate current radius with soft start
            if ti < ramp_duration:
                # Normalized time for ramp
                tau = ti / ramp_duration
                # Apply smootherstep to radius
                r = r_max * smootherstep(tau)
            else:
                r = r_max

            # Circular motion
            theta = 2 * math.pi * f * ti
            x = r * math.cos(theta)
            y = r * math.sin(theta)
            z = 0.0
            
            # Rotation (none)
            rx = 0.0
            ry = 0.0
            rz = 0.0
            
            # Format: (time (tx ty tz) (rx ry rz))
            file.write(f"({ti:.6g} ({x:.6g} {y:.6g} {z:.6g}) ({rx:.6g} {ry:.6g} {rz:.6g}))\n")
            
        file.write(")\n")

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python3 generate_motion.py <radius> <frequency> <duration> <dt> <ramp_duration>")
        sys.exit(1)
        
    r = float(sys.argv[1])
    f = float(sys.argv[2])
    duration = float(sys.argv[3])
    dt = float(sys.argv[4])
    ramp_duration = float(sys.argv[5])
    
    output_file = "constant/6DoF.dat"
    
    # Handle default ramp duration (if passed as negative, use 10% of duration)
    ramp_arg = float(sys.argv[5])
    if ramp_arg < 0:
        ramp_duration = duration * 0.1
        print(f"Using default ramp duration: 10% of duration = {ramp_duration}s")
    else:
        ramp_duration = ramp_arg

    generate_motion(r, f, duration, dt, ramp_duration, output_file)
    print(f"Motion file written to {output_file} with ramp duration {ramp_duration}s")
