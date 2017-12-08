import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

def frw(eta, w, p):
    L, k, Omega_Lambda, Omega_m, H_0 = p

    a, r, rdot, t, tdot, phi = w
    # r, rdot, t, tdot, phi, E, R = w

    # spacetime quantities
    R = r*a
    E = -k*r**2
    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)
    R_t = r*a_t
    R_r = a
    R_rt = a_t
    R_rr = 0
    E_r = 2*k*r
    
    # coordinates
    phidot = L / R**2
    tddot = -R*R_t*phidot**2 - rdot**2/(1+E)*R_rt*R_r
    rddot = -2*R_rt/R_r*tdot*rdot + (E_r/(1+E)-R_rr/R_r)*rdot**2 + (1+E)*R/R_r*phidot**2

    # print(r, rdot, rddot)
    return [
        a_t *tdot,
        rdot,
        rddot,
        tdot,
        tddot,
        phidot,
    ]

def frw2(w, eta, p):
    L, k, Omega_Lambda, Omega_m, H_0 = p

    a, r, rdot, t, tdot, phi = w
    # r, rdot, t, tdot, phi, E, R = w

    # spacetime quantities
    R = r*a
    E = -k*r**2
    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)
    R_t = r*a_t
    R_r = a
    R_rt = a_t
    R_rr = 0
    E_r = 2*k*r
    
    # coordinates
    phidot = L / R**2
    tddot = -R*R_t*phidot**2 - rdot**2/(1+E)*R_rt*R_r
    rddot = -2*R_rt/R_r*tdot*rdot + (E_r/(1+E)-R_rr/R_r)*rdot**2 + (1+E)*R/R_r*phidot**2

    # print(r, rdot, rddot)
    return [
        a_t *tdot,
        rdot,
        rddot,
        tdot,
        tddot,
        phidot,
    ]

count = 0
def kottler(eta, w, p):
    E, L, M, Omega_Lambda, Omega_m, H_0 = p
    r_h, t, r, rdot, phi = w

    Lambda = 3*Omega_Lambda*H_0**2
    f = 1 - 2*M/r - Lambda/3*r**2
    rddot = L**2 * (r - 3*M) / r**4
    tdot = E / f
    phidot = L/r**2
    r_h_t = (1 - 2*M/r_h - Lambda/3*r_h**2) * np.sqrt(2*M/r_h + Lambda/3*r_h**2)

    global count
    if count < 20:
        print(r_h, r, rdot, rddot, L, phi, phidot)
        count += 1

    return [
        r_h_t*tdot,
        tdot,
        rdot,
        rddot,
        phidot,
    ]

INTEGRATOR = 'vode'
def solve():

    k = 0
    M = 1.
    # H_0 = 7.33e-27
    # H_0 = 1e-10
    H_0 = 1e-3
    Omega_Lambda = 0
    Omega_m = 1 - Omega_Lambda

    a0 = 1
    initial_a = a0
    # initial_r = 10e17
    initial_r = 150
    initial_phi = 0.
    initial_t = 0.

    angle_to_horizontal = 0.15
    initial_rdot = -initial_r/10
    initial_phidot = -np.tan(angle_to_horizontal) * initial_rdot / initial_r
    print("initial velocities:", initial_rdot, initial_phidot, initial_tdot)
    initial_R = initial_a*initial_r
    initial_tdot = -np.sqrt(initial_a**2/(1-k*initial_r**2)*initial_rdot**2 + initial_R**2*initial_phidot**2)


    rho = Omega_m*3*H_0**2/(8*np.pi)
    r_h = 1/initial_a*(3*M/(4*np.pi*rho))**(1./3)
    print('r_h:', r_h)

    Lambda = 3*Omega_Lambda*H_0**2
    L_frw = (initial_a*initial_r)**2*initial_phidot

    solver_frw = spi.ode(frw).set_integrator(INTEGRATOR)

    p_frw = [L_frw, k, Omega_Lambda, Omega_m, H_0]
    initial = [initial_a, initial_r, initial_rdot, initial_t, initial_tdot, initial_phi]

    eta = np.arange(0, 100, 0.1)
    # sol_frw_only = spi.odeint(frw2, initial, eta, args=(p_frw,), printmessg=1, mxstep=500000)

    solver_frw.set_initial_value(initial, 0).set_f_params(p_frw)
    sol = []
    dt = 0.1
    while solver_frw.successful():
        solver_frw.integrate(solver_frw.t + dt)
        sol.append(list(solver_frw.y))
        if solver_frw.y[1] <= r_h:
            last = solver_frw.y
            break
    
    solver_kottler = spi.ode(kottler).set_integrator(INTEGRATOR)
    r_out = last[0] * last[1]
    initial_t = 0
    initial_r = r_out
    initial_phi = last[5]
    initial_rh = initial_r
    f = 1 - 2*M/r_out + Lambda/3*r_out**2
    etadot = last[4] / last[0] # conformal time

    # print('kottler rdot start:', last[2])
    # a, r, rdot, t, tdot, phi
    initial_rdot = last[0] * (np.sqrt(1 - f)*etadot + last[2])
    initial_tdot = last[0]/f*(etadot + np.sqrt(1-f)*last[2])
    initial_phidot = L_frw / last[1]**2
    L_kottler = initial_phidot *initial_r**2

    initial_kottler = [initial_rh, initial_t, initial_r, initial_rdot, initial_phi]
    E = f*initial_tdot
    p_kottler = [E, L_kottler, M, Omega_Lambda, Omega_m, H_0]

    solver_kottler.set_initial_value(initial_kottler, 0).set_f_params(p_kottler)
    sol_kottler = []
    while solver_kottler.successful():
        # dt = 0.001
        dt = 100
        solver_kottler.integrate(solver_kottler.t + dt, step=True)
        sol_kottler.append(list(solver_kottler.y))
        # print("step: ", solver_kottler.y)
        if solver_kottler.y[2] > solver_kottler.y[0]:
            last = solver_kottler.y
            break

    # r_h, t, r, rdot, phi = w
    initial_phi = last[4]
    initial_r = r_h
    initial_a = last[0] / r_h
    print("scale factor on exit:", initial_a)
    initial_phidot = L_kottler / last[2]**2
    f = 1-2*M/last[2]-Lambda/3*last[2]**2
    last_tdot = E/f
    initial_rdot = 1/initial_a*(1/f*last[3] - np.sqrt(1-f)*last_tdot)
    initial_etadot = 1/initial_a*(last_tdot - np.sqrt(1-f)/f*last[3])
    initial_tdot = initial_etadot * initial_a

    # print(p_frw)
    p_frw[0] = initial_a**2 * initial_r**2*initial_phidot # change the L_frw
    # print(p_frw)
    # p_frw = [L_frw, k, Omega_Lambda, Omega_m, H_0]
    # r_h, t, r, rdot, phi = w

    if Omega_Lambda:
        initial_t = 2/(3*H_0*np.sqrt(Omega_Lambda))*np.arcsin(np.sqrt(Omega_Lambda/(1-Omega_Lambda))*(initial_r/a0/r_h)**(3/2))
    else:
        initial_t = 2/(3*H_0)*(initial_r/a0/r_h)**(3/2)
    # a, t, r, rdot, phi 

    # check if its going to cross the axis
    initial_ydot = initial_rdot*np.sin(initial_phi) + initial_r*np.cos(initial_phi)*initial_phidot
    if initial_ydot > 0:
        print("light ray is not going to cross the axis, shoot it closer to the black hole")

    initial_frw2 = [initial_a, initial_r, initial_rdot, initial_t, initial_tdot, initial_phi]
    solver_frw2 = spi.ode(frw).set_integrator(INTEGRATOR)
    solver_frw2.set_initial_value(initial_frw2, 0).set_f_params(p_frw)

    t_end = 100
    print("=====")
    while solver_frw2.successful() and solver_frw2.t < t_end:
        dt = 0.1
        solver_frw2.integrate(solver_frw2.t + dt)
        sol.append(list(solver_frw2.y))
        if solver_frw2.y[1] * np.sin(solver_frw2.y[5]) < 0:  # stop when it crosses the axis
            break


    sol = np.array(sol)
    r = sol[:,1]
    phi = sol[:,5]

    x = r * np.cos(phi)
    y = r * np.sin(phi)
    plt.plot(x, y, 'bo')

    # the frw only line
    # r_frw = sol_frw_only[:,1]
    # phi_frw = sol_frw_only[:,5]
    # plt.plot(r_frw*np.cos(phi_frw), r_frw*np.sin(phi_frw), 'r-')

    swiss_cheese_hole = plt.Circle((0., 0.), r_h, color='grey', fill=False, zorder=10)
    axes = plt.gca()
    axes.add_artist(swiss_cheese_hole)
    axes.set_aspect('equal', adjustable='box')
    # plot the center
    plt.plot(0, 0, 'ro')

    lim = 200
    axes.set_xlim([-lim, lim])
    axes.set_ylim([-lim, lim])



    plt.figure()
    sol_kottler = np.array(sol_kottler)
    r_kottler = sol_kottler[:,2]
    phi_kottler = sol_kottler[:,4]
    x_kottler = r_kottler * np.cos(phi_kottler)
    y_kottler = r_kottler * np.sin(phi_kottler)
    axes = plt.gca()
    swiss_cheese_hole = plt.Circle((0., 0.), r_h, color='grey', fill=False, zorder=10)
    axes.add_artist(swiss_cheese_hole)
    axes.set_aspect('equal', adjustable='box')
    plt.plot(x_kottler, y_kottler, 'g-')
    lim = 200
    axes.set_xlim([-lim, lim])
    axes.set_ylim([-lim, lim])

solve()
plt.show()