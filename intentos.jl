using Javis

#constants 
G = 6.6743e-11 
m1 = 5.972e24
m2 = 7.347e22
m3 = 1.988e30


function derivatives(t, u)

    #initial conditions

    pos1 = u[1]
    pos2 = u[2]
    pos3 = u[3]
    vel1 = u[4]
    vel2 = u[5]
    vel3 = u[6]

    #distance between centers
    r12 = sqrt((pos2[1]-pos1[1])^2 + (pos2[2]-pos1[2])^2)
    r13 = sqrt((pos3[1]-pos1[1])^2 + (pos3[2]-pos1[2])^2)
    r23 = sqrt((pos3[1]-pos2[1])^2 + (pos3[2]-pos2[2])^2)


    #Forces:
    f12 = [
        G*(m1*m2)*(pos2[1]-pos1[1])/r12^3,
        G*(m1*m2)*(pos2[2]-pos1[2])/r12^3
    ]

    f21 = f12 * -1

    f13 = [
        G*(m1*m3)*(pos3[1]-pos1[1])/r13^3,
        G*(m1*m3)*(pos3[2]-pos1[2])/r13^3
    ]

    f31 = f13 * -1

    f23 = [
        G*(m2*m3)*(pos3[1]-pos2[1])/r23^3,
        G*(m2*m3)*(pos3[2]-pos2[2])/r23^3
    ]


    f32 = f23 * -1

    #total forces
    F1 = f12 + f13
    F2 = f21 + f23
    F3 = f31 + f32

    #Accelerations
    acc1 = F1/m1
    acc2 = F2/m2
    acc3 = F3/m3

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
    end

