import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

def ltb(w, eta, p):
    L, k, Omega_Lambda, Omega_m, H_0 = p

    # differential equations
    r, rdot, t, tdot, phi, a = w
    
    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)
    
    # coordinates
    phidot = L / (r*a)**2
    tddot = -a*r**2*a_t*phidot**2 - rdot**2*a*a_t
    rddot = -2*a_t/a*tdot*rdot + r*phidot**2

    # print(r, rdot, rddot)
    return [
        rdot,
        rddot,
        tdot,
        tddot,
        phidot,
        a_t*tdot,
    ]


def solve():
    k = 0
    eta = np.arange(1, 100, 0.1)

    initial_a = 1
    initial_r = 150
    initial_phi = 0.
    initial_rdot = -20.
    initial_phidot = 0.05
    initial_tdot = np.sqrt(initial_a**2*initial_rdot**2 + initial_a**2*initial_r**2*initial_phidot**2)


    # r, rdot, t, tdot, phi, E, R
    initial = [initial_r, initial_rdot, 0, initial_tdot, initial_phi, initial_a]
    print(initial)

    L = initial_r**2*initial_phidot
    Omega_Lambda, Omega_m = 0.1, .9
    H_0 = 1e-3
    p = [L, k, Omega_Lambda, Omega_m, H_0]
    sol = spi.odeint(ltb, initial, eta, args=(p,), mxstep=500000)
    r = sol[:,0]
    phi = sol[:,4]

    a = sol[:,5]
    rdot = sol[:,1]
    phidot = L/(a*r)**2
    angles = np.arctan((rdot*np.sin(phi) + r*np.cos(phi)*phidot)/(rdot*np.cos(phi) - r*np.sin(phi)*phidot))

    print("angles:")
    print(angles)

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