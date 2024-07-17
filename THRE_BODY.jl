using DifferentialEquations
using Plots

# Define the three-body problem function
function three_body!(du, u, p, t)
    G = p[1]  # Gravitational constant
    m1 = p[2] # Mass of body 1
    m2 = p[3] # Mass of body 2
    m3 = p[4] # Mass of body 3
    
    # Unpack positions and velocities
    x1, y1, vx1, vy1 = u[1:4]
    x2, y2, vx2, vy2 = u[5:8]
    x3, y3, vx3, vy3 = u[9:12]
    
    # Distances between bodies
    r12 = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    r13 = sqrt((x1 - x3)^2 + (y1 - y3)^2)
    r23 = sqrt((x2 - x3)^2 + (y2 - y3)^2)
    
    # Gravitational forces
    f12 = G * m1 * m2 / r12^3
    f13 = G * m1 * m3 / r13^3
    f23 = G * m2 * m3 / r23^3
    
    # Accelerations
    ax1 = -f12 * (x1 - x2) - f13 * (x1 - x3)
    ay1 = -f12 * (y1 - y2) - f13 * (y1 - y3)
    
    ax2 = f12 * (x1 - x2) - f23 * (x2 - x3)
    ay2 = f12 * (y1 - y2) - f23 * (y2 - y3)
    
    ax3 = f13 * (x1 - x3) + f23 * (x2 - x3)
    ay3 = f13 * (y1 - y3) + f23 * (y2 - y3)
    
    # Update derivatives
    du[1] = vx1
    du[2] = vy1
    du[3] = ax1
    du[4] = ay1
    
    du[5] = vx2
    du[6] = vy2
    du[7] = ax2
    du[8] = ay2
    
    du[9] = vx3
    du[10] = vy3
    du[11] = ax3
    du[12] = ay3
end

# Initial conditions and parameters
G = 1.0      # Gravitational constant
m1 = 1.0     # Mass of body 1
m2 = 1.0     # Mass of body 2
m3 = 1.0     # Mass of body 3

# Initial positions and velocities (adjust as needed)
u0 = [
    1.0, 0.0, 0.0, 0.5,   # x1, y1, vx1, vy1
    -0.5, 0.866, -0.25, 0.0,  # x2, y2, vx2, vy2
    -0.5, -0.866, -0.25, 0.0  # x3, y3, vx3, vy3
]
