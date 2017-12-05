import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

def frw(w, eta, p):
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

# def frw(w, eta, p):
#     L, k, Omega_Lambda, Omega_m, H_0 = p

#     a, r, rdot, t, tdot, phi = w

#     a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)

#     phidot = L / r**2
#     tddot = -a*r**2*phidot**2*a_t - a*rdot**2*a_t/ (1 - k*r**2)
#     rddot = r*phidot**2 - k*rdot**2 - 2*a_t/a*rdot*tdot

#     print(r, rdot, rddot)
#     return [
#         a_t*tdot,
#         rdot,
#         rddot,
#         tdot,
#         tddot,
#         phidot
#     ]


def solve():
    eta = np.arange(1, 100, 0.1)

    k = 0

    initial_a = 1
    initial_r = 150
    initial_phi = 0.
    initial_rdot = -20
    initial_phidot = 0.05
    initial_tdot = np.sqrt(initial_a**2*initial_rdot**2 + initial_a**2*initial_r**2*initial_phidot**2)

    # a, r, rdot, t, tdot, phi
    initial = [initial_a, initial_r, initial_rdot, 0, initial_tdot, initial_phi]
    print(initial)

    L = initial_r**2*initial_phidot
    Omega_Lambda, Omega_m = 0.1, 0.9
    H_0 = 1e-3
    p = [L, k, Omega_Lambda, Omega_m, H_0]
    sol = spi.odeint(frw, initial, eta, args=(p,), mxstep=500000)
    r = sol[:,1]
    phi = sol[:,5]

    # print(r)
    # print(phi)


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
