# formulas according to fleury paper that seems to be working, but don't think it's right

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt


kottler = False
np.seterr(all='warn', over='raise')

# def ode_solver(f, initial, t, args):
#     current = np.array(initial)
#     results = []
#     for i in range(len(t)-1):
#         dt = t[i+1] - t[i]
#         derivatives = np.array(f(current, t[i], *args))
#         current = derivatives*dt + current
#         results.append(current)

#     return np.array(results)

def ltb(w, eta, p):
    global kottler
    L, k, Omega_Lambda, Omega_m, H_0, R_h, M = p
    r, rdot, t, tdot, phi, a, a_tilde = w
    rho_tilde = 0

    H = H_0*np.sqrt(Omega_m/a**3 + Omega_Lambda)
    a_t = H*a
    a_tilde_t = a_t
    a_tilde_dot = None

    if r < R_h:
        # enter the Kottler hole
        if not kottler:
            print("entering======================")
        kottler = True
        rho_tilde = 3*M/(4*np.pi*(a*R_h)**3)


    if r > R_h:
        if kottler:
            print("exiting======================")
        # exit the Kottler hole
        kottler = False

    if kottler:
        Lambda = Omega_Lambda*3*H_0**2
        rho_tilde = 3*M/(4*np.pi*(r)**3)
        H_tilde = np.sqrt(8*np.pi*rho_tilde/3*(a/a_tilde)**3 + Lambda/3)
        R = a_tilde*r
        # R_r = a_tilde*H_tilde*r
        R_r = a_tilde
        a_tilde_t = a_tilde*H_tilde
        R_t = r*a_tilde_t
        a_tilde_r = (a_tilde*H_tilde*r-a_tilde)/r
        a_tilde_tt = (-4*np.pi/3*rho_tilde + Lambda/3)*a_tilde
        R_rt = r*a_tilde_tt
        R_rr = a_tilde_t
        # R_rt = -1./2*(2*M/R+Lambda/3.*R**2)*(-2*M/R**2+2*Lambda/3*R)*R_t
        # R_rr = -1./2*(2*M/R+Lambda/3.*R**2)*(-2*M/R**2+2*Lambda/3*R)*R_r
        E = 0
        E_r = 0
        a_tilde_dot = a_tilde_r*rdot + a_tilde_t*tdot

    else:
        R = a*r
        R_t = a_t*r
        R_r = a
        R_rt = a_t
        R_rr = 0
        E = -k*r**2
        E_r = -2*k*r

    phidot = L/R**2
    tddot = -R*R_t*phidot**2 - rdot**2/(1+E)*R_rt*R_r
    rddot = -2*R_rt/R_r*tdot*rdot + (E_r/2/(1+E) - R_rr/R_r)*rdot**2 + (1+E)*R/R_r*phidot**2
    adot = a_t*tdot
    if not a_tilde_dot:
        a_tilde_dot = adot
    

    # print(r, rdot, rddot, R_r)

    return [
        rdot,
        rddot,
        tdot,
        tddot,
        phidot,
        adot,
        a_tilde_dot,
    ]



def frw(w, eta, p):
    global kottler
    L, k, Omega_Lambda, Omega_m, H_0, R_h, M = p

    # differential equations
    r, rdot, t, tdot, phi, a, a_tilde = w

    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)

    R = a*r
    R_t = r*a_t
    R_r = a
    R_rt = a_t
    R_rr = 0
    E_r = 2*k*r
    E = k*r**2
    
    # coordinates
    phidot = L / R**2
    Edot = E_r * rdot
    tddot = -R*R_t*phidot**2 - rdot**2/(1+E)*R_rt*R_r
    rddot = -2*R_rt/R_r*tdot*rdot + (E_r/2/(1+E)-R_rr/R_r)*rdot**2 + (1+E)*R/R_r*phidot**2
    Rdot = R_t*tdot + R_r*rdot


    return [
        rdot,
        rddot,
        tdot,
        tddot,
        phidot,
        a_t*tdot,
        a_t*tdot,
    ]


def solve():
    # eta = np.linspace(100, 200000, 5000)
    # eta = np.arange(1, 10000, 0.01)
    # eta = np.linspace(1, 100000, 10000)
    # eta = (1e7 - np.logspace(1, 7, 10000))[::-1]
    eta = np.arange(1, 100, 0.1)

    k = 0

    initial_a = 1
    initial_r = 2e+18
    initial_phi = 0.
    initial_rdot = -1e17
    initial_phidot = 0.01
    initial_tdot = np.sqrt(initial_a**2*initial_rdot**2 + initial_a**2*initial_r**2*initial_phidot**2)
    print(initial_tdot)
    initial_t = 1.

    # r, rdot, t, tdot, phi, E, R
    initial = [initial_r, initial_rdot, initial_t, initial_tdot, initial_phi, initial_a, initial_a]

    L = initial_r**2*initial_phidot
    Omega_Lambda = 0.7
    Omega_m = 1 - Omega_Lambda
    H_0 = 7.33e-27
    M = 1
    rho_c = 3*H_0**2/8/np.pi
    R_h = 1/initial_a*(3*M/(4*Omega_m*rho_c*np.pi))**(1/3)
    print(R_h)
    p = [L, k, Omega_Lambda, Omega_m, H_0, R_h, M]
    sol = spi.odeint(ltb, initial, eta, args=(p,))
    r = sol[:,0]
    phi = sol[:,4]
    a = sol[:,5]
    # print(r)
    # print(phi)
    # print(rdot)


    # sol_frw = spi.odeint(frw, initial, eta, args=(p,), mxstep=500000)
    # r_frw = sol_frw[:,0]
    # phi_frw = sol_frw[:,4]
    # x_frw = r_frw * np.cos(phi_frw)
    # y_frw = r_frw * np.sin(phi_frw)
    # plt.plot(x_frw, y_frw, "r-")


    x = r * np.cos(phi)
    y = r * np.sin(phi)

    plt.plot(x, y, "b-")

    # plot the center
    plt.plot(0, 0, 'ro')

    # black_hole = plt.Circle((0., 0.), 2*M, color='black', fill=True, zorder=10)
    swiss_cheese_hole = plt.Circle((0., 0.), R_h, color='grey', fill=False, zorder=10)
    axes = plt.gca()
    # axes.add_artist(black_hole)
    axes.add_artist(swiss_cheese_hole)
    axes.set_aspect('equal', adjustable='box')

    lim = 2.5e+18
    axes.set_xlim([-lim, lim])
    axes.set_ylim([-lim, lim])

solve()
plt.show()