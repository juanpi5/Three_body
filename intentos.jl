using Plots
using LinearAlgebra
#constants 
G = 6.6743e-11 
m1 = 5.972e24
m2 = 7.347e22
m3 = 1.988e30

function distance(pi, pj)
    r_ij = sqrt((pj[1]-pi[1])^2 + (pj[2]-pi[2])^2)
    return r_ij
end

function forces(pi, pj, pk, mi, mj, mk)
    #r, distance between two planets 
    #p, position of each planet 
    #m, mass of each planet
    r_ij = norm(pj-pi)
    r_ik = norm(pk-pi)
    F = G*((mi*mj)*(pj-pi)/r_ij^3 + (mi*mk)*(pk-pi)/r_ik^3)

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

    p1 = [1.496e11, 0]
    p2 = [0, 0]
    p3 = [1.5e11, 3.84e8] 

    v1 = [3e4, 3e4]
    v2 = [-1.1e3, -1.1e3]
    v3 = [1.4e5, 2.6e3]

    u0 = Vector([p1, p2, p3, v1, v2, v3])
    t0 = 0
    dt = 0.1
    t_max = 50

    poss_1 = Vector{Array{Float64}}(undef, convert(Int64, t_max/dt))
    poss_2 = Vector{Array{Float64}}(undef, convert(Int64, t_max/dt))
    poss_3 = Vector{Array{Float64}}(undef, convert(Int64, t_max/dt))

    for dt in 1:t_max

        poss_1[dt] = copy(u0[1])
        poss_2[dt] = copy(u0[2])
        poss_3[dt] = copy(u0[3])
        Runge_kutta(u0, t0, dt)

    end

# Running the simulation
    tp = 0
    
    @gif for dt in 1:t_max

        tp += 1 

            plot(poss_1[1:tp,1], poss_1[1:tp,2], label="m1", xlims=(-5,5), ylims=(-5,5),  linecolor=:red)
            plot!(poss_2[1:tp,1], poss_2[1:tp,2], label="m2", xlims=(-5,5), ylims=(-5,5), linecolor=:blue)
            plot!(poss_3[1:tp,1], poss_3[1:tp,2], label="m3", xlims=(-5,5), ylims=(-5,5), linecolor=:green)
            scatter!([poss_1[tp,1]], [poss_1[tp,2]], label="", xlims=(-5,5), ylims=(-5,5), markercolor=:red)
            scatter!([poss_2[tp,1]], [poss_2[tp,2]], label="", xlims=(-5,5), ylims=(-5,5), markercolor=:blue)
            scatter!([poss_3[tp,1]], [poss_3[tp,2]], label="", xlims=(-5,5), ylims=(-5,5), markercolor=:green)

    end

end
