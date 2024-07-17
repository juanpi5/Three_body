using Plots
using LinearAlgebra

#constants 
m1 = 10.0
m2 = 10.0
m3 = 1e-17

function forces(pi, pj, pk, mi, mj, mk)

    #r, distance between two planets 
    #p, position of each planet 
    #m, mass of each planet
    r_ij = norm(pj-pi)
    r_ik = norm(pk-pi)
    F = ((mi*mj)*(pj-pi)/r_ij^3 + (mi*mk)*(pk-pi)/r_ik^3)

end

function derivatives(t, u)

    #initial conditions:

    pos1, pos2, pos3, vel1, vel2, vel3 = u[1:6]

    #Forces:
    F1, F2, F3 = forces(pos1, pos2, pos3, m1, m2, m3), forces(pos2, pos1, pos3, m2, m1, m3), forces(pos3, pos1, pos2, m3, m1, m2)

    #Accelerations:
    acc1, acc2, acc3 = F1/m1, F2/m2, F3/m3

    #derivatives with respect to u
    du = Vector([vel1, vel2, vel3, acc1, acc2, acc3])

    return du
end


#Runge-Kutta
function Runge_kutta(u, t, dt) 

    k1 = derivatives(t, u);
    k2 = derivatives(t + dt/2, [u[1] + dt/2*k1[1], u[2] + dt/2*k1[2], u[3] + dt/2*k1[3], u[4] + dt/2*k1[4], u[5] + dt/2*k1[5], u[6] + dt/2*k1[6]]);
    k3 = derivatives(t + dt/2, [u[1] + dt/2*k2[1], u[2] + dt/2*k2[2], u[3] + dt/2*k2[3], u[4] + dt/2*k2[4], u[5] + dt/2*k2[5], u[6] + dt/2*k2[6]]);
    k4 = derivatives(t + dt, [u[1] + dt*k3[1], u[2] + dt*k3[2], u[3] + dt*k3[3], u[4] + dt*k3[4], u[5] + dt*k3[5], u[6] + dt*k3[6]]);

    for i in 1:6
        u[i] += dt/6 * (k1[i]+2*k2[i]+2*k3[i]+k4[i])
    end

end



function updated_positions() 

    #parameters

    p1 = [3.0, 0.0]
    p2 = [-3.0, 0.0]
    p3 = [0.0, 0.0] 

    v1 = [0.0, 1.0]
    v2 = [0.0, -1.0]
    v3 = [0.0, 0.0]

    u0 = Vector([p1, p2, p3, v1, v2, v3])
    t0 = 0
    dt = 0.1
    t_max = 100
    steps = convert(Int64, t_max/dt)

    poss_1 = Vector{Array{Float64}}(undef, steps)
    poss_2 = Vector{Array{Float64}}(undef, steps)
    poss_3 = Vector{Array{Float64}}(undef, steps)

    for i in 1:steps

        poss_1[i] = copy(u0[1])
        poss_2[i] = copy(u0[2])
        poss_3[i] = copy(u0[3])

        Runge_kutta(u0, t0, dt)

    end

# Running the simulation

    @gif for i in 1:steps
        plot(map(p -> p[1], poss_1[1:i]), map(p -> p[2], poss_1[1:i]), label="m1", xlims=(-20, 20), ylims=(-20, 20), linecolor=:red)
        plot!(map(p -> p[1], poss_2[1:i]), map(p -> p[2], poss_2[1:i]), label="m2", linecolor=:blue)
        plot!(map(p -> p[1], poss_3[1:i]), map(p -> p[2], poss_3[1:i]), label="m3", linecolor=:green)
        scatter!([poss_1[i][1]], [poss_1[i][2]], label="", markercolor=:red)
        scatter!([poss_2[i][1]], [poss_2[i][2]], label="", markercolor=:blue)
        scatter!([poss_3[i][1]], [poss_3[i][2]], label="", markercolor=:green)
    end

end

anim = updated_positions()
