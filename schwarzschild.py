import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

G = 6.67e-11
M_sun = 1.989e30
R_sun = 695700e3
c = 299792458
plt.style.use('ggplot')

def schw_null_geodesics(w, t, p):
    r, rdot, phi = w
    M, L = p
    phidot = L / r**2
    return [
        rdot,
        L**2 * (r - 3*M) / r**4,
        phidot
    ]

def schw_null_geodesics2(t, w, p):
    r, rdot, phi = w
    M, L = p
    phidot = L / r**2
    return [
        rdot,
        L**2 * (r - 3*M) / r**4,
        phidot
    ]

def expected_schw_light_bending(r, M):
    return 4. * M / r

def solve(max_t, initial, p):
    t = np.arange(0, max_t, 0.01)

    r = spi.ode(schw_null_geodesics2)
    r.set_initial_value(initial).set_f_params(p).set_integrator('vode')
    results = []
    dt = 0.01
    # res = r.integrate(max_t)
    # while r.successful() and r.t < max_t:
    #     res = r.integrate(r.t+dt)
    #     results.append(res)
    # results = np.array(results)
    # sol = results

    sol = spi.odeint(schw_null_geodesics, initial, t, args=(p,),)

    # plot it
    r = sol[:,0]
    phi = sol[:,2]
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    plt.plot(x, y, color="b")

    # get angle to horizontal
    num_entries = len(t)
    gradient = (y[num_entries*4//5] - y[num_entries-1]) / (x[num_entries*4//5] - x[num_entries-1])
    deflection = np.arctan(-gradient)

    # # unperturbed orbit
    # unperturbed_phi = np.linspace(np.pi * 1/5,  np.pi * 4/5)
    # unperturbed_r = initial[0] / np.sin(unperturbed_phi)
    # unperturbed_x = unperturbed_r * np.cos(unperturbed_phi)
    # unperturbed_y = unperturbed_r * np.sin(unperturbed_phi)
    # # plt.plot(unperturbed_x, unperturbed_y, linestyle='--', color='gray')

    return deflection

def main():
    M = 1.
    plt.figure(dpi=300)
    # initial_x = -150
    max_t = 1000.
    # bs = np.arange(20., 100., 1.)
    bs = np.array([4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40])
    # bs = np.array([1000.])
    # E = -1.
    initial_xs = bs*0 - 10
    # initial_xs = bs*0 - 50
    deflections = []
    for b in bs:
        initial_x = -b
        initial_r = np.sqrt(b**2 + initial_x**2)
        initial_phi = np.arccos(initial_x / initial_r)

        # L = b/E
        # f = 1. - 2*M/initial_r
        # initial_tdot = -E/f
        # initial_rdot = initial_tdot * np.cos(initial_phi)
        # initial_phidot = initial_tdot / initial_r * np.sqrt(1 - (initial_rdot/initial_tdot)**2)

        initial_rdot = np.cos(initial_phi)
        initial_phidot = -np.sqrt((1 - initial_rdot**2) / initial_r**2)
        L = initial_r**2 * initial_phidot

        d = solve(max_t, [initial_r, initial_rdot, initial_phi], [M, L])
        deflections.append(d)

    axes = plt.gca()
    lim = max(-initial_xs)
    # lim = 50
    axes.set_xlim([-lim,lim])
    axes.set_ylim([-lim,lim])
    # plot the blackhole
    circle = plt.Circle((0., 0.), 2*M, color='black', fill=True, zorder=10)
    plt.legend(bs)
    axes.add_artist(circle)
    axes.set_aspect('equal', adjustable='box')


    # plot deflection angles
    expected = expected_schw_light_bending(bs, M)
    # plt.figure(dpi=300)

    # plt.plot(bs, expected, 'r+')
    # plt.plot(bs, deflections, 'b+')
    # plt.errorbar(bs, deflections, yerr=np.array(deflections)*0.05, fmt="none", ecolor='b', elinewidth='1.')
    # plt.ylabel('Deflection angle')
    # plt.xlabel('Distance of closest approach')
    # plt.legend([r'Theoretical deflection ($\frac{4M}{R}$)', 'Numerical result'])


    # plt.savefig('images/deflections_trans.png')
    percentage_errors = np.absolute(expected - deflections) / deflections
    print(percentage_errors)
    # plt.gca().set_aspect('equal', adjustable='box')


# main()
# plt.show()