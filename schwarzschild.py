# code for generating the graphs in https://theconfused.me/blog/numerical-integration-of-light-paths-in-a-schwarzschild-metric/
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi


def schw_null_geodesics(w, t, p):
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


def solve(t, initial, p):
    sol = spi.odeint(schw_null_geodesics, initial, t, args=(p,),)
    r = sol[:,0]
    phi = sol[:,2]
    return r, phi


def calculate_deflection(x, y, t):
    num_entries = len(t)
    gradient = (y[num_entries*4//5] - y[num_entries-1]) / (x[num_entries*4//5] - x[num_entries-1])
    return np.arctan(-gradient)


def calculate_initial_conditions(b, initial_x, M):
    # calculate initial conditions for each impact parameter
    initial_r = np.sqrt(b**2 + initial_x**2)
    initial_phi = np.arccos(initial_x / initial_r)

    f = 1 - 2*M/initial_r
    new_rdot = -np.sqrt(f**2 / (1 + f * (np.tan(initial_phi))**2))
    new_phidot = -new_rdot * np.tan(initial_phi) / initial_r
    # initial_rdot = new_rdot
    # initial_phidot = new_phidot

    initial_rdot = np.cos(initial_phi)
    initial_phidot = -np.sqrt((1 - initial_rdot**2) / initial_r**2)

    print("=====")
    print(initial_rdot, initial_phidot)
    print(new_rdot, new_phidot)

    L = initial_r**2 * initial_phidot
    return [initial_r, initial_rdot, initial_phi], L


def plot_light_rays():
    plt.figure()
    M = 1.  # natural units
    initial_x = -50
    max_t = 10000.

    # impact parameters to plot
    bs = np.array([4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40])
    # time points to evaluate
    t = np.arange(0, max_t, 0.01)

    for b in bs:
        initial, L = calculate_initial_conditions(b, initial_x, M)
        r, phi = solve(t, initial, [M, L])
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        plt.plot(x, y)

    axes = plt.gca()
    lim = -initial_x
    axes.set_xlim([-lim, lim])
    axes.set_ylim([-lim, lim])

    # plot the black hole, schwarzschild radius is 2M
    circle = plt.Circle((0., 0.), 2*M, color='black', fill=True, zorder=10)
    plt.legend(bs)
    axes.add_artist(circle)
    axes.set_aspect('equal', adjustable='box')

def plot_deflection_angles():
    plt.figure()
    M = 1.  # natural units
    initial_x = -50
    max_t = 10000.

    # impact parameters to plot
    bs = np.linspace(20, 100, 40)
    # time points to evaluate
    t = np.arange(0, max_t, 0.01)

    deflections = []
    for b in bs:
        initial, L = calculate_initial_conditions(b, initial_x, M)
        r, phi = solve(t, initial, [M, L])
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        deflections.append(calculate_deflection(x, y, t))

    # plot deflections angles against analytical predicted deflections angles
    expected = expected_schw_light_bending(bs, M)
    plt.plot(bs, expected, 'ro')
    plt.plot(bs, deflections, 'b+')
    plt.ylabel('Deflection angle')
    plt.xlabel('Distance of closest approach')
    plt.legend([r'Theoretical deflection ($\frac{4M}{R}$)', 'Numerical result'])


if __name__ == '__main__':
    plt.style.use('ggplot')
    plot_light_rays()
    plot_deflection_angles()
    plt.show()