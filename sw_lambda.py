'''
origin at lens
'''

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import time


length_scale = 3.086e23 # mega parsec
INTEGRATOR = 'vode'
INTEGRATOR_PARAMS = {
    'atol': 1e-110, 
    # 'atol': 0,
    # 'rtol': 0,
    'rtol': 1e-15,
    'nsteps': 100000000,
    # 'method': 'bdf',
}
# H_0 = 7.33e-27
H_0 = 7.56e-27 * length_scale
# Omega_Lambda = 0
# Omega_m = 1 - Omega_Lambda
M = 0.5e15 / length_scale
print("M: ", M)

def frw(eta, w, p):
    L, k, Omega_Lambda, Omega_m, H_0 = p

    a, r, rdot, t, phi = w

    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)

    phidot = L / (a*r)**2
    # tddot = -a*r**2*phidot**2*a_t - a*rdot**2*a_t/ (1 - k*r**2)
    tdot = -np.sqrt(a**2*rdot**2+a**2*r**2*phidot**2)
    rddot = r*phidot**2 - k*rdot**2 - 2*a_t/a*rdot*tdot

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

def kottler(eta, w, p):
    E, L, M, Omega_Lambda, Omega_m, H_0 = p
    r_h, t, r, rdot, phi = w

    Lambda = 3*Omega_Lambda*H_0**2
    f = 1 - 2*M/r - Lambda/3*r**2
    rddot = L**2 * (r - 3*M) / r**4
    tdot = E / f
    phidot = L/r**2
    r_h_t = (1 - 2*M/r_h - Lambda/3*r_h**2) * np.sqrt(2*M/r_h + Lambda/3*r_h**2)

    return [
        r_h_t*tdot,
        tdot,
        rdot,
        rddot,
        phidot,
    ]

def solve(angle_to_horizontal, comoving_lens=1e25, plot=True, Omega_Lambda=0):
    k = 0
    a0 = 1
    initial_a = a0
    # initial_r = 10e17
    initial_r = comoving_lens
    initial_phi = np.pi
    initial_t = 0.
    Omega_m = 1 - Omega_Lambda

    # initial_tdot = -1.
    # initial_rdot = -np.sqrt(initial_tdot**2/initial_a**2/(1+(np.tan(angle_to_horizontal))**2))
    # initial_phidot = initial_rdot * np.tan(angle_to_horizontal) / initial_r

    initial_rdot = -initial_r
    initial_phidot = np.tan(angle_to_horizontal) * initial_rdot / initial_r
    initial_R = initial_a*initial_r
    initial_tdot = -np.sqrt(initial_a**2/(1-k*initial_r**2)*initial_rdot**2 + initial_R**2*initial_phidot**2)
    
    if plot:
        print("initial velocities:", initial_rdot, initial_phidot, initial_tdot)

    rho = Omega_m*3*H_0**2/(8*np.pi)
    r_h = 1/initial_a*(3*M/(4*np.pi*rho))**(1./3)
    if plot:
        print('r_h:', r_h, "\n")

    Lambda = 3*Omega_Lambda*H_0**2
    L_frw = (initial_a*initial_r)**2*initial_phidot

    solver_frw = spi.ode(frw).set_integrator(INTEGRATOR, **INTEGRATOR_PARAMS)

    p_frw = [L_frw, k, Omega_Lambda, Omega_m, H_0]
    initial = [initial_a, initial_r, initial_rdot, initial_t, initial_phi]
    
    # save for later frw straight line propagation
    initial0 = initial

    # eta = np.arange(0, 100, 0.1)
    # sol_frw_only = spi.odeint(frw2, initial, eta, args=(p_frw,), printmessg=1, mxstep=500000)

    solver_frw.set_initial_value(initial, 0).set_f_params(p_frw)
    sol = []
    dt = 5e-7
    while solver_frw.successful():
        solver_frw.integrate(solver_frw.t + dt, step=False)
        sol.append(list(solver_frw.y))
        if solver_frw.y[1] <= r_h:
            last = solver_frw.y
            break

    sol = sol
    last = sol[-1]
    # print("FRW coordinates, before entering kottler:")
    sol1 = np.array(sol)
    # angle = get_angle(sol1[:,1], sol1[:,4], sol1[:,2], L_frw/(sol1[:,0]*sol1[:,1])**2)
    # plt.plot(sol1[:,1]*np.cos(sol1[:,4]), angle)
    # plt.show()

    if plot:
        frw_angle_before_entering = get_angle(last[1], last[4], last[2], L_frw/(last[0]*last[1])**2)
        print("Angle in FRW::")
        print(frw_angle_before_entering)
        print(last)
        print("\n")
    solver_kottler = spi.ode(kottler).set_integrator(INTEGRATOR, **INTEGRATOR_PARAMS)
    # a, r, rdot, t, phi = w
    r_out = last[0] * last[1]
    tdot_out = -np.sqrt(last[0]**2*last[2]**2+last[0]**2*last[1]**2*L_frw/(last[0]*last[1])**2)
    initial_t = 0
    initial_r = r_out
    initial_phi = last[4]
    initial_rh = initial_r
    f = 1 - 2*M/r_out - Lambda/3*r_out**2
    etadot = tdot_out / last[0] # conformal time

    # a, r, rdot, t, phi
    initial_rdot = last[0] * (np.sqrt(1 - f)*etadot + last[2])
    initial_tdot = last[0]/f*(etadot + np.sqrt(1-f)*last[2])
    initial_phidot = L_frw / (last[0]*last[1])**2
    L_kottler = initial_phidot *initial_r**2

    initial_kottler = [initial_rh, initial_t, initial_r, initial_rdot, initial_phi]
    E = f*initial_tdot
    p_kottler = [E, L_kottler, M, Omega_Lambda, Omega_m, H_0]

    solver_kottler.set_initial_value(initial_kottler, 0).set_f_params(p_kottler)

    if plot:
        print("kottler initial:")
        print("initial_rh, initial_t, initial_r, initial_rdot, initial_phi")
        print(initial_kottler)

    angle_before_entering = get_angle(initial_r, initial_phi, initial_rdot, initial_phidot)
    if plot:
        print("light ray angle before entering kottler hole: ", angle_before_entering)
        print("initial conditions on entering kottler hole: ", initial_kottler)
        print("initial params on entering kottler hole: ", p_kottler)
    sol_kottler = []
    first_time = True
    prev_r = np.inf
    while solver_kottler.successful():
        # dt = 1e-6
        solver_kottler.integrate(solver_kottler.t + 5e-7, step=False)
        sol_kottler.append(list(solver_kottler.y))
        # if solver_kottler.y[2] > prev_r and first_time:
        #     if plot:
        #         print("turning point in kottler metric:", solver_kottler.y[2])
        #     first_time = False
        # else:
        #     prev_r = solver_kottler.y[2]
        if solver_kottler.y[4] < np.pi/2 and first_time:
            if plot:
                print("turning point in kottler metric:", solver_kottler.y[2])
            first_time = False
        if solver_kottler.y[2] > solver_kottler.y[0]:
            last = solver_kottler.y
            break

    angle_after_exiting = get_angle(last[2], last[4], last[3], L_kottler/last[2]**2)
    if plot:
        print("light ray angle after exiting kottler hole: ", angle_after_exiting)
        print("bending angle in kottler: ", angle_before_entering - angle_after_exiting)
        print("\n")

    # r_h, t, r, rdot, phi = w
    initial_phi = last[4]
    initial_r = r_h
    initial_a = last[2] / r_h
    if plot:
        print("scale factor on exit:", initial_a)
    initial_phidot = L_kottler / last[2]**2
    f = 1-2*M/last[2]-Lambda/3*last[2]**2
    last_tdot = E/f
    initial_rdot = 1/initial_a*(1/f*last[3] - np.sqrt(1-f)*last_tdot)
    initial_etadot = 1/initial_a*(last_tdot - np.sqrt(1-f)/f*last[3])
    initial_tdot = initial_etadot * initial_a


    p_frw[0] = initial_a**2 * initial_r**2*initial_phidot # change the L_frw
    # p_frw = [L_frw, k, Omega_Lambda, Omega_m, H_0]
    # r_h, t, r, rdot, phi = w

    if Omega_Lambda:
        initial_t = 0
        # initial_t = 2/(3*H_0*np.sqrt(Omega_Lambda))*np.arcsin(np.sqrt(Omega_Lambda/(1-Omega_Lambda))*(last[2]/a0/r_h)**(3/2))
    else:
        initial_t = 2/(3*H_0)*(last[2]/a0/r_h)**(3/2)
    # a, t, r, rdot, phi 

    frw_angle_after_exiting = get_angle(initial_r, initial_phi, initial_rdot, initial_phidot)
    if plot:
        print("Angle after exiting FRW:", frw_angle_after_exiting)
        print("bending angle in FRW: ", frw_angle_before_entering - frw_angle_after_exiting)
        print("\n")

    # check if its going to cross the axis
    initial_ydot = initial_rdot*np.sin(initial_phi) + initial_r*np.cos(initial_phi)*initial_phidot
    if initial_ydot > 0:
        print("light ray is not going to cross the axis, decrease angle_to_horizontal")
        print("initial angle to horizontal: ", angle_to_horizontal)
        print("----")

    if initial_r*np.sin(initial_phi) < 0:
        print("light ray is bent too much by the hole, increase angle_to_horizontal")
        print("initial angle to horizontal: ", angle_to_horizontal)
        print("----")

    initial_frw2 = [initial_a, initial_r, initial_rdot, initial_t, initial_phi]
    if plot:
        print("FRW2 initial: ", initial_frw2)
        print(initial_r*np.cos(initial_phi), initial_r*np.sin(initial_phi))
    solver_frw2 = spi.ode(frw).set_integrator(INTEGRATOR, **INTEGRATOR_PARAMS)
    solver_frw2.set_initial_value(initial_frw2, 0).set_f_params(p_frw)

    while solver_frw2.successful():
        # dt = 1e-6
        solver_frw2.integrate(solver_frw2.t + dt, step=False)
        sol.append(list(solver_frw2.y))
        if solver_frw2.y[1] * np.sin(solver_frw2.y[4]) < 0:  # stop when it crosses the axis
            break

    sol = np.array(sol)
    r = sol[:,1]
    phi = sol[:,4]
    a = sol[:,0]

    # s = spi.ode(frw)
    # p_s = [0, k, Omega_Lambda, Omega_m, H_0]
    # source_a = None
    # # initial0[1] = 1e-100
    # initial_s = [1, 1e-8, comoving_lens, 0, 0]
    # # initial = [initial_a, initial_r, initial_rdot, initial_t, initial_phi]
    # s.set_f_params(p_s).set_initial_value(initial_s).set_integrator(INTEGRATOR, **INTEGRATOR_PARAMS)
    # while s.successful():
    #     s.integrate(s.t + dt)
    #     if s.y[1] > (r[-1] + comoving_lens):
    #         source_a = s.y[0]
    #         break

    if plot:
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        plt.plot(x, y, 'bo')
        
        print("axis crossing points", r[-1], phi[-1])
        swiss_cheese_hole = plt.Circle((0., 0.), r_h, color='grey', fill=False, zorder=10)
        axes = plt.gca()
        axes.add_artist(swiss_cheese_hole)
        axes.set_aspect('equal', adjustable='box')
        # plot the center
        plt.plot(0, 0, 'ro')

        lim = 0.5
        axes.set_xlim([-lim, lim])
        axes.set_ylim([-lim, lim])

        plt.figure()
        sol_kottler = np.array(sol_kottler)
        r_kottler = sol_kottler[:,2]
        phi_kottler = sol_kottler[:,4]
        x_kottler = r_kottler * np.cos(phi_kottler)
        y_kottler = r_kottler * np.sin(phi_kottler)
        axes = plt.gca()
        swiss_cheese_hole = plt.Circle((0., 0.), r_h, color='grey', fill=False, zorder=10)
        axes.add_artist(swiss_cheese_hole)
        axes.set_aspect('equal', adjustable='box')
        plt.plot(x_kottler, y_kottler, 'g-')
        lim = 200
        axes.set_xlim([-lim, lim])
        axes.set_ylim([-lim, lim])

    # print("diff between as", source_a, a[-1])
    return r[-1], a[-1]
    # return r[-1], source_a

# [ 0.18075975  0.05122652  0.23824122  0.31020329  0.20010044  0.39768539
#   0.43731914  0.46608477  0.37572935  0.39980157]

def get_distances(z, Omega_Lambda=0):
    Omega_m = 1 - Omega_Lambda
    def integrand(z):
        return 1/np.sqrt(Omega_m*(1+z)**3 + Omega_Lambda)
    integral, error = spi.quad(integrand, 0, z)
    comoving = integral/H_0
    dang = comoving/(1+z)
    return comoving, dang

def calc_theta(D_LS, D_L, D_S):
    return np.sqrt(4*M*D_LS/D_L/D_S/1.01)

from tqdm import tqdm
def main():
    start = time.time()
    # thetas = np.linspace(25, 35, 10)*10**(-6)
    theta = 25e-6
    # print("thetas", thetas)
    # thetas = np.array([15e-6])
    om_lambdas = np.linspace(0, 0.9, 10)
    # om = 0
    z_lens = 0.01
    a_lens = 1/(z_lens+1)
    ds = []
    dls = []
    rs = []
    dl = []
    # for theta in tqdm(thetas):
    for om in tqdm(om_lambdas):
        # print("omega_lambda:", om)
        comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om)
        dl.append(dang_lens)
        # print("lens distances: ", comoving_lens, dang_lens)
        r, a = solve(theta, plot=False, comoving_lens=comoving_lens, Omega_Lambda=om)
        d_s = a*(r + comoving_lens)
        d_ls = a/a_lens * r
        ds.append(d_s)
        dls.append(d_ls)
        # rs.append(r)
    ds = np.array(ds)
    dls = np.array(dls)
    dl = np.array(dl)
    numerical_thetas = calc_theta(dls, dl, ds)
    # alpha = 4*M/dang_lens/thetas
    # theoretical_rs = thetas/alpha*comoving_lens/(1-thetas/alpha/(1+z_lens))

    # print("=====")
    # rs = np.array(rs) / 1.01
    # print(rs)
    # print(theoretical_rs)
    # print((rs - theoretical_rs)/theoretical_rs * 100)
    print(numerical_thetas)
    percentage_errors = np.abs(numerical_thetas - theta)/theta
    # percentage_errors = np.abs(th - ds)/ds
    print(percentage_errors*100)
    print("Time taken: {}".format(time.time() - start))

main()