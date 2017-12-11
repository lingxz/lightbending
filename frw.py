import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

def frw(eta, w, p):
    L, k, Omega_Lambda, Omega_m, H_0 = p

    a, r, rdot, t, phi = w

    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)

    phidot = L / (a*r)**2
    tdot = np.sqrt(a**2*rdot**2+a**2*r**2*phidot**2)
    rddot = r*phidot**2 - k*rdot**2 - 2*a_t/a*rdot*tdot

    print(phidot, r, r*np.cos(phi), get_angle(r, phi, rdot, phidot))

    return [
        a_t*tdot,
        rdot,
        rddot,
        tdot,
        phidot,
    ]

def get_angle(r, phi, rdot, phidot):
    res = np.arctan((rdot*np.sin(phi)+r*np.cos(phi)*phidot)/(rdot*np.cos(phi)-r*np.sin(phi)*phidot))
    return res

def solve():
    eta = np.arange(0, 30, 0.1)

    k = 0

    # initial_a = 1
    # initial_r = 150
    # initial_phi = 0.
    # initial_rdot = -20
    # initial_phidot = 0.05
    # initial_tdot = np.sqrt(initial_a**2*initial_rdot**2 + initial_a**2*initial_r**2*initial_phidot**2)

    initial_a = 1
    initial_r = 398.9463488840856
    initial_phi = np.pi
    initial_rdot = -39.89463488840856
    initial_phidot = -8.00000000017e-07
    initial_tdot = -39.8946348897

    # a, r, rdot, t, tdot, phi
    initial = [initial_a, initial_r, initial_rdot, 0, initial_phi]
    print(initial)

    L = initial_a**2*initial_r**2*initial_phidot
    Omega_Lambda, Omega_m = 0, 1.
    H_0 = 1e-3
    p = [L, k, Omega_Lambda, Omega_m, H_0]

    r = spi.ode(frw)
    r.set_initial_value(initial).set_f_params(p).set_integrator('vode', atol=1e-100, rtol=1e-15)
    results = []
    dt = 0.1
    r_h = 1.0599032239320836
    while r.successful():
        res = r.integrate(r.t+dt)
        results.append(list(res))
        if r.y[1] < r_h:
            # print(r.t)
            break
        # if r.t > 30:
        #     break

    results = np.array(results)
    sol = results


    # sol = spi.odeint(lambda x, t, *args: frw(t, x, *args), initial, eta, args=(p,))
    # print("========")
    # print(sol)
    r = sol[:,1]
    phi = sol[:,4]

    # print(r)
    # print(phi)
    angles = get_angle(r, phi, sol[:,2], L/(sol[:,0]*r)**2)
    print(angles)
    plt.plot(r * np.cos(phi), angles)
    plt.show()


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
# plt.show()
