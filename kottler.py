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
    M, L, Lambda = p
    phidot = L / r**2
    return [
        rdot,
        L**2 * (r - 3*M) / r**4,
        phidot
    ]

def kottler(w, t, p):
    E, L, M, Omega_Lambda, Omega_m, H_0 = p
    t, r, rdot, phi = w

    Lambda = 3*Omega_Lambda*H_0**2
    f = 1 - 2*M/r - Lambda/3*r**2
    rddot = L**2 * (r - 3*M) / r**4
    tdot = E / f
    phidot = L/r**2
    # r_h_t = (1 - 2*M/r_h - Lambda/3*r_h**2) * np.sqrt(2*M/r_h + Lambda/3*r_h**2)

    return [
        # r_h_t*tdot,
        tdot,
        rdot,
        rddot,
        phidot,
    ]

def expected_schw_light_bending(r, M):
    return 4. * M / r

def get_angle(r, phi, rdot, phidot):
    res = np.arctan((rdot*np.sin(phi)+r*np.cos(phi)*phidot)/(rdot*np.cos(phi)-r*np.sin(phi)*phidot))
    return res

def solve(max_t, initial, p, phi_plane):
    t = np.arange(0, max_t, 0.05)

    sol = spi.odeint(schw_null_geodesics, initial, t, args=(p,),)

    # plot it
    r = sol[:,0]
    phi = sol[:,2]
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    # plt.plot(x, y, color="b")

    rdot = sol[:,1]
    phidot = p[1]/r**2
    deflection = get_angle(r, phi, rdot, phidot)
    print("deflection:", deflection[-1])

    # # get angle to horizontal
    # isclose = np.isclose(phi, np.full(phi.shape, phi_plane), atol=1e-5)
    # index = np.argwhere(isclose == True)[0,0]
    # dr = r[index+1] - r[index]
    # dphi = phi[index+1] - phi[index]
    # A = dr/dphi

    # M, L, Lambda = p
    # f = 1 - 2*M/r[index] - Lambda / 3. * r[index]**2
    # deflection = r[index] * dphi / dr * np.sqrt(f)

    # # cosphi = (1./f *) +  
    # # deflection = np.arctan((f**0.5 * r[index]) / np.abs(A)) - phi_plane
    # # deflection = np.arccos(np.abs(A) / (A**2 + f*r[index]**2)**0.5) - phi_plane

    # # num_entries = len(t)
    # # gradient = (y[num_entries*4//5] - y[num_entries-1]) / (x[num_entries*4//5] - x[num_entries-1])
    # # deflection = np.arctan(-gradient)
    # return -deflection

M = 1
b = 1000
Lambda = 10**(-40)
max_t = 5000
initial_x = -b
initial_r = np.sqrt(b**2 + initial_x**2)
initial_phi = np.arccos(initial_x / initial_r)
initial_rdot = np.cos(initial_phi)
initial_phidot = -np.sqrt((1 - initial_rdot**2) / initial_r**2)
initial_tdot = 1
f = 1-2*M/initial_r
L = initial_r**2 * initial_phidot
E = f*initial_tdot
p_original = [M, L, Lambda]
# p = [E, L, M, 0, 1, 1]
# initial = [0, initial_r, initial_rdot, initial_phi]
initial_orig = [initial_r, initial_rdot, initial_phi]
phi_plane = 0.
solve(max_t, initial_orig, p_original, 0)
print("expected:", expected_schw_light_bending(b, M))
# plt.plot(0, 0, 'ro')
# lim = 10
# axes = plt.gca()
# axes.set_xlim([-lim, lim])
# axes.set_ylim([-lim, lim])
# plt.show()


def main():
    M = 1.
    plt.figure(dpi=300)
    # initial_x = -150
    max_t = 90000.
    bs = np.arange(500., 550., 1.)
    # bs = np.array([4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40])
    # bs = np.array([40.])
    # bs = np.array([1000.])
    # E = -1.
    # initial_xs = bs*0 - 10
    deflections = []
    phi_plane = 0.
    Lambda = 10**(-40)
    for b in bs:
        initial_x = -b
        initial_r = np.sqrt(b**2 + initial_x**2)
        initial_phi = np.arccos(initial_x / initial_r)
        initial_rdot = np.cos(initial_phi)
        initial_phidot = -np.sqrt((1 - initial_rdot**2) / initial_r**2)
        L = initial_r**2 * initial_phidot
        d = solve(max_t, [initial_r, initial_rdot, initial_phi], [M, L, Lambda], phi_plane)
        deflections.append(d)

    Lambda = 50**(-40)
    deflections2 = []
    for b in bs:
        initial_x = -b
        initial_r = np.sqrt(b**2 + initial_x**2)
        initial_phi = np.arccos(initial_x / initial_r)
        initial_rdot = np.cos(initial_phi)
        initial_phidot = -np.sqrt((1 - initial_rdot**2) / initial_r**2)
        L = initial_r**2 * initial_phidot
        d = solve(max_t, [initial_r, initial_rdot, initial_phi], [M, L, Lambda], phi_plane)
        deflections2.append(d)


    # plot deflection angles
    expected = expected_schw_light_bending(bs, M)
    # print(deflections)
    # print(expected)
    # plt.figure(dpi=300)

    plt.plot(bs, expected, 'r+')
    plt.plot(bs, deflections, 'b+')
    plt.plot(bs, deflections2, 'g+')
    # plt.errorbar(bs, deflections, yerr=np.array(deflections)*0.05, fmt="none", ecolor='b', elinewidth='1.')
    plt.ylabel('Deflection angle')
    plt.xlabel('Distance of closest approach')
    plt.legend([r'Theoretical deflection ($\frac{4M}{R}$)', 'Numerical result'])

    plt.savefig('images/kottler_deflections.png')
    percentage_errors = np.absolute(expected - deflections) / deflections
    print(percentage_errors)

    # plt.gca().set_aspect('equal', adjustable='box')


