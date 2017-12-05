# current version with all the printouts

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

def get_R(M, Lambda, R_h, t, r, print_msg=False):
    if Lambda == 0:
        return (9/2*M*(t-r+R_h)**2)**(1/3)
    else:
        if print_msg:
            print(t, r)
        return (6*M/Lambda*(np.sinh(np.sqrt(3*Lambda/4)*(t-r+R_h)))**2)**(1/3)

kottler = False
L_kottler = None

def ltb(w, eta, p):
    global kottler
    global L_kottler
    L_frw, k, Omega_Lambda, Omega_m, H_0, R_h, M = p
    r, rdot, t, tdot, phi = w

    # H = H_0 * np.sqrt(Omega_Lambda + Omega_m/a**3)
    # a_t = H*a

    rho_0 = 3*H_0**2/(8*np.pi)
    a_0 = 1
    a = (9/4)**(1/3)*(8*np.pi/3*rho_0*a_0**3)**(1/3)*t**(2/3)
    a_t = (2/3)*(9/4)**(1/3)*(8*np.pi/3*rho_0*a_0**3)**(1/3)*t**(-1/3)


    if np.abs(r) < R_h:
        if not kottler:
            print("enter")
            L_kottler = L_frw
        #     Lambda = 3*Omega_Lambda*H**2
        #     R = get_R(M, Lambda, R_h, t, r)
        #     print(rdot, phi, t, r)
        kottler = True
    else:
        if kottler:
            print("exit")
        kottler = False

    if kottler:
        Lambda = 3*Omega_Lambda*H_0**2
        R = np.abs((9/2*M*(t-r+R_h)**2)**(1/3))


        R_r = -np.sqrt(2*M/R + Lambda/3*R**2)
        # R_r = (6*M/Lambda)**(1/3)*(2/3)*(np.sinh(np.sqrt(3*Lambda/4)*(t-r)))**(-1/3)*np.cosh(np.sqrt(3*Lambda/4)*(t-r))*(-np.sqrt(3*Lambda/4))
        R_t = -R_r
        # R_rr = Lambda/3*R-2*M/R**2
        R_rr = 1/2*(-2*M/R**2+2*Lambda/3*R)
        # R_rr = -1/2*(2*M/R+Lambda/3*R**2)**(-1/2)*(-2*M/R**2+2*Lambda/3*R)*R_r
        R_rt = R_t/2/R_r*(-2*M/R**2+2*Lambda*R/3)
        E = 0
        E_r = 0

        phidot = L_kottler/R**2
        tddot = -R_r*R_rt*rdot**2 - R*R_t*phidot**2
        rddot = -2*R_rt/R_r*tdot*rdot + (E_r/2/(1+E)-R_rr/R_r)*rdot**2 + (1+E)*R/R_r*phidot**2

        # tddot = -0.5*(-2*M*R_t/R**2+2*Lambda/3*R*R_t)*rdot**2 - R*R_t*phidot**2
        # rddot = 1/(2*M/R+Lambda/3*R**2)*(-2*M*rdot/R**2*(-rdot*R_r-tdot*R_t) - Lambda/3*rdot*(2*rdot*R_r+2*tdot*R_t) + rdot**2/2*(-2*M*R_r/R**2+2*Lambda/3*R*R_r))

        print(eta, R, t, M, r, rdot, rddot)
        print((1+E)*R/R_r*phidot**2)
        # print(r, rdot, rddot, R, R_rt)
        # print(t, r, rdot, rddot, R, R_rr)

        return [
            rdot,
            rddot,
            tdot,
            tddot,
            phidot,
        ]

    # tddot = -R*R_t*phidot**2 - rdot**2/(1+E)*R_rt*R_r
    # rddot = -2*R_rt/R_r*tdot*rdot + (E_r/2/(1+E)-R_rr/R_r)*rdot**2 + (1+E)*R/R_r*phidot**2


    else:
        R = a*r
        R_r = a
        R_t = r*a_t
        R_rr = 0
        R_rt = a_t
        E = 0
        E_r = 0

        phidot = L_frw/R**2
        tddot = -R*R_t*phidot**2 - rdot**2/(1+E)*R_rt*R_r
        rddot = -2*R_rt/R_r*tdot*rdot + (E_r/2/(1+E)-R_rr/R_r)*rdot**2 + (1+E)*R/R_r*phidot**2

        # print(r, rdot, rddot)
        print(eta, R, t, M, r, rdot, rddot, a)

        return [
            rdot,
            rddot,
            tdot,
            tddot,
            phidot,
        ]


def solve():
    eta = np.arange(0, 100, 0.005)
    
    k = 0
    M = 1
    H_0 = 7.33e-27
    Omega_Lambda = 0.
    Omega_m = 1 - Omega_Lambda
    
    initial_a = 1
    initial_r = 10e17
    initial_phi = 0.
    rho = Omega_m*3*H_0**2/(8*np.pi)
    initial_t = (6*np.pi*rho)**(-1/2)

    initial_rdot = -1e17
    initial_phidot = 0.03
    
    initial_R = initial_a*initial_r
    initial_tdot = np.sqrt(initial_a**2/(1-k*initial_r**2)*initial_rdot**2 + initial_R**2*initial_phidot**2)
    initial = [initial_r, initial_rdot, initial_t, initial_tdot, initial_phi]
    print('initial_tdot:', initial_tdot)

    rho = Omega_m*3*H_0**2/(8*np.pi)
    R_h = 1/initial_a*(3*M/(4*np.pi*rho))**(1./3)
    print('r_h:', R_h)

    Lambda = 3*Omega_Lambda*H_0**2
    L_frw = (initial_a*initial_r)**2*initial_phidot
    # L_kottler = get_R(M, Lambda, R_h, initial_t, initial_r)**2*initial_phidot

    p = [L_frw, k, Omega_Lambda, Omega_m, H_0, R_h, M]

    # solver = spi.ode(ltb)
    # solver.set_integrator('vode')
    # solver.set_initial_value(initial, 1).set_f_params(p)
    # t_end = 50
    # sol = []
    # count = 0
    # while solver.successful() and solver.t < t_end:
    #     dt = 0.5
    #     solver.integrate(dt)
    #     if count < 2:
    #         print(solver.y)
    #     sol.append(list(solver.y))
    #     count += 1

    # sol = np.array(sol)
    # print(sol)

    sol, info = spi.odeint(ltb, initial, eta, args=(p,), printmessg=1, full_output=1, mxstep=500000)

    r = sol[:,0]
    # r = np.trim_zeros(r, 'b')
    phi = sol[:,4]

    print(r)

    print(r[-1])
    print(phi[-1])

    x = r * np.cos(phi)
    y = r * np.sin(phi)
    # plt.plot(x, y, 'b-')
    plt.plot(x, y, marker='.', color='b', linestyle='-')

    # black_hole = plt.Circle((0., 0.), 2*M, color='black', fill=True, zorder=10)
    swiss_cheese_hole = plt.Circle((0., 0.), R_h, color='grey', fill=False, zorder=10)
    axes = plt.gca()
    # axes.add_artist(black_hole)
    axes.add_artist(swiss_cheese_hole)
    axes.set_aspect('equal', adjustable='box')
    # plot the center
    plt.plot(0, 0, 'ro')

    lim = 12e17
    axes.set_xlim([-lim, lim])
    axes.set_ylim([-lim, lim])


solve()
plt.show()