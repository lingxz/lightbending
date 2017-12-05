import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

def ltb(w, eta, p):
    L, k, Omega_Lambda, Omega_m, H_0 = p

    # differential equations
    r, rdot, t, tdot, phi, E, R = w

    # spacetime quantities
    a = R / r
    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)
    R_t = r*a_t
    R_r = a
    R_rt = a_t
    R_rr = 0
    E_r = 2*k*r
    
    # coordinates
    phidot = L / R**2
    Edot = E_r * rdot
    tddot = -R*R_t*phidot**2 - rdot**2/(1+E)*R_rt*R_r
    rddot = -2*R_rt/R_r*tdot*rdot + (E_r/(1+E)-R_rr/R_r)*rdot**2 + (1+E)*R/R_r*phidot**2
    Rdot = R_t*tdot + R_r*rdot

    print(r, rdot, rddot)
    return [
        rdot,
        rddot,
        tdot,
        tddot,
        phidot,
        Edot,
        Rdot,
    ]


def solve():
    eta = np.arange(1, 100, 0.1)

    k = 0

    initial_a = 1
    initial_r = 150
    initial_phi = 0.
    initial_rdot = -20.
    initial_phidot = 0.05
    initial_tdot = np.sqrt(initial_a**2*initial_rdot**2 + initial_a**2*initial_r**2*initial_phidot**2)
    initial_R = initial_r*initial_a
    initial_E = k*initial_r**2

    # r, rdot, t, tdot, phi, E, R
    initial = [initial_r, initial_rdot, 0, initial_tdot, initial_phi, initial_E, initial_R]
    print(initial)

    L = initial_r**2*initial_phidot
    Omega_Lambda, Omega_m = 0.1, .9
    H_0 = 1e-3
    p = [L, k, Omega_Lambda, Omega_m, H_0]
    sol = spi.odeint(ltb, initial, eta, args=(p,), mxstep=500000)
    r = sol[:,0]
    phi = sol[:,4]

    print(r)
    print(phi)


    x = r * np.cos(phi)
    y = r * np.sin(phi)
    plt.plot(x, y, "b-")

    # plot the center
    plt.plot(0, 0, 'ro')
    axes = plt.gca()

    lim = 200
    axes.set_xlim([-lim, lim])
    axes.set_ylim([-lim, lim])

solve()
plt.show()