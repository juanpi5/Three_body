using Javis

#constants 
G = 6.6743e-11 
m1 = 5.972e24
m2 = 7.347e22
m3 = 1.988e30

function distance(pi, pj)
    r_ij = sqrt((pj[1]-pi[1])^2 + (pj[2]-pi[2])^2)
    return r_ij
end

function forces(r_ij, r_ik, pi, pj, pk, mi, mj, mk)
    #r, distance between two planets 
    #p, position of each planet 
    #m, mass of each planet
    
    F = ([
        G*((mi*mj)*(pj[1]-pi[1])/r_ij^3 + (mi*mk)*(pk[1]-pi[1])/r_ik),
        G*((mi*mj)*(pj[2]-pi[2])/r_ij^3 + (mi*mk)*(pk[2]-pi[2])/r_ik)
    ])

end

function derivatives(t, u)

    #initial conditions

    pos1, pos2, pos3, vel1, vel2, vel3 = u[1:6]

    #distance between centers
    r12, r13, r23 = distance(pos1, pos2), distance(pos1, pos3), distance(pos2, pos3)

    #Forces:
    F1, F2, F3 = forces(r12, r13, pos1, pos2, pos3, m1, m2, m3), forces(r12, r23, pos2, pos1, pos3, m2, m1, m3), forces(r13, r23, pos3, pos1, pos2, m3, m1, m2)

    #Accelerations
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


    p1 = [1.496e11, 0]
    p2 = [0, 0]
    p3 = [1.5e11, 3.84e8] 

    v1 = [3e4, 3e4]
    v2 = [-1.1e3, -1.1e3]
    v3 = [1.4e5, 2.6e3]

    u0 = Vector([p1, p2, p3, v1, v2, v3])
    t0 = 0
    dt = 0.1
    
    

    for t in 1:60
        Runge_kutta(u0, t0, dt)
        println(u0)
        println("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")
    end

