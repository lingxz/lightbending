import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

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

def solve(initial, p, r_h):

    r = spi.ode(schw_null_geodesics2)
    r.set_initial_value(initial).set_f_params(p).set_integrator('vode')
    results = []
    dt = 1e-5
    # res = r.integrate(max_t)
    while r.successful():
        res = r.integrate(r.t+dt)
        results.append(res)
        if r.y[0] > r_h and r.y[2] < np.pi/2:
            break
    results = np.array(results)
    sol = results

    # sol = spi.odeint(schw_null_geodesics, initial, t, args=(p,),)

    # plot it
    r = sol[:,0]
    rdot = sol[:,1]
    phi = sol[:,2]
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    plt.plot(x, y, color="b")

    # get angle to horizontal
    exiting_angle = r[-1]*p[1]/r[-1]**2/rdot[-1]
    print("exiting_angle: ", exiting_angle)
    # num_entries = len(t)
    # gradient = (y[num_entries*4//5] - y[num_entries-1]) / (x[num_entries*4//5] - x[num_entries-1])
    # deflection = np.arctan(-gradient)

    # # unperturbed orbit
    # unperturbed_phi = np.linspace(np.pi * 1/5,  np.pi * 4/5)
    # unperturbed_r = initial[0] / np.sin(unperturbed_phi)
    # unperturbed_x = unperturbed_r * np.cos(unperturbed_phi)
    # unperturbed_y = unperturbed_r * np.sin(unperturbed_phi)
    # # plt.plot(unperturbed_x, unperturbed_y, linestyle='--', color='gray')

    # return deflection

# r, rdot, phi
initial = [0.00096342312817610052, -0.043884038009890776, 3.1396295238764198]
p = [3.240440699935191e-11, -6.5801956132742141e-08]
r_h = 0.000962923396595
solve(initial, p, r_h)