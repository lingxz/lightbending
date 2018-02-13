'''
origin at lens
'''

# diff_lambdas_small2: latest, diff z_lens, 100 per Omega_Lambda, but wrong DL
# diff_lambdas_small: same z_lens, diff z_source, 100 per Omega_Lambda
# diff_lambdas_small3: same as diff_lambdas_small2, but with correct DL

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import time
import pandas as pd
from scipy.special import erf
from scipy import misc


length_scale = 3.086e22 # mega parsec
INTEGRATOR = 'vode'
INTEGRATOR_PARAMS = {
    'atol': 1e-120, 
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
# M = 0.5e15 / length_scale
M = 1474e12 / length_scale
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

def ltb_Efun(r):
    rho_outside_now = ?? 
    return lambda r: (rho_outside_now*8*np.pi/3)**2*r**2*(1-3*M/(4*np.pi*r**3*rho_outside_now))


def mass_inside_hole(r, rh, sigma=None):
    if not sigma:
        sigma = rh/10.
    integral, error = spi.quad(lambda x: rho_m(0, x,rh,sigma), 0, r)
    return integral*4*np.pi


def normalization_factor(sigma, rh):
    c = 1/(2*sigma**2)
    B = np.sqrt(np.sqrt(np.pi))*erf(np.sqrt(c)*rh)/(4*c**(3/2)) - rh*np.exp(-c*rh**2)/2*c
    B = np.sqrt(2*np.pi/2*sigma**3)
    return M/(4*np.pi*B)

def rho_m(t, r, rh, sigma=None):
    if not sigma:
        sigma = rh/10.
    A = normalization_factor(sigma, rh)
    return A*np.exp(-r**2/(2*sigma**2))
    # om_m = 1
    # return om_m*3*H_0**2/(8*np.pi)


def ltb(eta, w, p):
    L, Omega_Lambda, Omega_m, H_0, rh, a_at_border = p  # M is mass within a comoving sphere at a given r
    Lambda = 3*Omega_Lambda*H_0**2

    a, t, r, rdot, phi, R, E = w

    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)

    ##############
    # a_at_border_t = a_at_border * H_0 * np.sqrt(Omega_m/a_at_border**3 + Omega_Lambda)
    # H0_at_border = a_at_border_t/a_at_border
    # tau = 3*Lambda/4*t
    # R0 = lambda x: x
    # R = R0(r)*(np.cosh(tau) + np.sqrt(3/Lambda)*ltb_hubble(r)*np.sinh(tau))
    # R_t = 
    #################

    # spacetime quantities
    mass = mass_inside_hole(r, rh)
    mprime = misc.derivative(mass_inside_hole, r, n=1, args=(rh,))

    E = 0
    E_r = 0
    # E = 2*(1/2*H0_at_border**2*r**2 - 1/(6*np.pi)*mass/r)
    # E_r = 2*H0_at_border**2*r - 1/(3*np.pi)*(r*mprime - mass)/r**2

    # print(2*mass/R, 1/3*Lambda*R**2, E)
    R_t = np.sqrt(2*mass/R + 1/3*Lambda*R**2 + E)

    R_r = mprime/R**2/(4*np.pi*rho_m())

    R_r = -E_r*R/(8*np.pi*rho_m(t, r, rh)*R**2 + Lambda*r**2 + E + 3*R_t**2)
    R_rt = -2*R_t*R_r/R
    R_rr = 0


    phidot = L/R**2
    tdot = -np.sqrt(R_r**2/(1+E)*rdot**2 + R**2*phidot**2)
    # print(R_rt, tdot, rdot, R_r)
    # print(-2*R_rt*tdot*rdot/R_r, (E_r/2/(1+E) - R_rr/R_r)*rdot**2, (1+E)*R/R_r*phidot**2)
    rddot = -2*R_rt*tdot*rdot/R_r + (E_r/2/(1+E) - R_rr/R_r)*rdot**2 + (1+E)*R/R_r*phidot**2
    
    # r_h, t, r, rdot, phi, R, E
    return [
        a_t*tdot,
        tdot,
        rdot,
        rddot,
        phidot,
        R_r*rdot + R_t*tdot,
        E_r*rdot,
    ]


def solve(angle_to_horizontal, comoving_lens=1e25, plot=True, Omega_Lambda=0, dt=5e-7):
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
    if r_h > initial_r:
        print("Starting point is inside the hole! Make the hole smaller or bring the lens further away.")
        return
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
    # sol = []
    while solver_frw.successful():
        solver_frw.integrate(solver_frw.t + dt, step=False)
        # sol.append(list(solver_frw.y))
        if solver_frw.y[1] <= r_h:
            last = solver_frw.y
            break

    if plot:
        frw_angle_before_entering = get_angle(last[1], last[4], last[2], L_frw/(last[0]*last[1])**2)
        print("Angle in FRW::")
        print(frw_angle_before_entering)
        print(last)
        print("\n")


    ############## preparing numbers to enter hole ##############################

    solver_ltb = spi.ode(ltb).set_integrator(INTEGRATOR, **INTEGRATOR_PARAMS)
    # frw: a, r, rdot, t, phi = w
    initial_a = last[0]
    initial_r = last[1]
    initial_rdot = last[2]
    initial_t = last[3]
    initial_phi = last[4]
    initial_R = last[0]*last[1]
    initial_phidot = L_frw/last[1]**2
    L_ltb = initial_R**2*initial_phidot
    initial_E = 0
    # initial_a_t = initial_a * H_0 * np.sqrt(Omega_m/initial_a**3 + Omega_Lambda)
    # initial_E = 2*(1/2*(initial_a_t/initial_a)**2*r_h**2 - 1/(6*np.pi)*M/r_h)
    # a, t, r, rdot, phi, R, E = w
    initial_ltb = [initial_a, initial_t, initial_r, initial_rdot, initial_phi, initial_R, initial_E]
    p_ltb = [L_ltb, Omega_Lambda, Omega_m, H_0, r_h, last[0]]

    solver_ltb.set_initial_value(initial_ltb, 0).set_f_params(p_ltb)


    if plot:
        print("ltb initial:")
        print("a, t, r, rdot, phi, R, E")
        print(initial_ltb)

    # if plot:
    #     angle_before_entering = get_angle(initial_r, initial_phi, initial_rdot, initial_phidot)
    #     print("light ray angle before entering LTB hole: ", angle_before_entering)
    #     print("initial conditions on entering LTB hole: ", initial_ltb)
    #     print("initial params on entering LTB hole: ", p_ltb)



    ############## end of preparing numbers to enter hole ##############################
    #######################################


    ############## in the hole ##############################

    first_time = True
    while solver_ltb.successful():
        solver_ltb.integrate(solver_ltb.t + dt, step=False)
        if solver_ltb.y[2] > r_h:
            last = solver_ltb.y
            break

    ############## exit the hole ##############################
    #########################################


    # angle_after_exiting = get_angle(last[2], last[4], last[3], L_kottler/last[2]**2)
    # if plot:
    #     print("light ray angle after exiting kottler hole: ", angle_after_exiting)
    #     print("bending angle in kottler: ", angle_before_entering - angle_after_exiting)
    #     print("\n")

    ############## converting numbers to go back into cheese ##############################

    
    # frw: a, r, rdot, t, phi = w
    # LTB: a, t, r, rdot, phi, R, E = w
    initial_a = last[0]
    initial_r = last[2]
    initial_rdot = last[3]
    initial_t = last[1]
    initial_phi = last[4]

    initial_phidot = L_ltb / last[5]**2
    initial_rdot = last[3]
    initial_tdot = -np.sqrt(initial_a**2/(1+last[6])*initial_rdot**2 + last[5]**2*initial_phidot**2)

    p_frw[0] = initial_a**2 * initial_r**2*initial_phidot # change the L_frw

    ############## finished converting numbers to go back into cheese ##############################
    ##################################################################

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

    sol = []
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
    return np.sqrt(4*M*D_LS/D_L/D_S)

def theta2rls_flat(theta, z_lens):
    rl, dl = get_distances(z_lens, Omega_Lambda=0)
    return 4*M*rl/(4*M-dl*theta**2)

def rs2theta(rs, rl, dl):
    return np.sqrt(4*M*(rs-rl)/rs/dl)

def rs2redshift_flat(rs):
    return 1/(1-H_0*rs/2)**2 -1

from tqdm import tqdm

def main():
    start = time.time()
    om_lambdas = np.linspace(0, 0.99, 1)
    # z_lens_all = np.linspace(0.05, 0.2, 10)
    z_lens_all = np.linspace(0.05, 0.2, 1)
    # z_lens = 0.1
    # a_lens = 1/(z_lens+1)
    # start_thetas = np.linspace(0.7e-5, 2e-5, 100)
    start_thetas = np.array([1e-5]*1)
    source_rs = np.array([theta2rls_flat(th1, z1) for th1, z1 in zip(start_thetas, z_lens_all)])
    # source_rs = theta2rls_flat(start_thetas, z_lens)
    source_zs = rs2redshift_flat(source_rs)
    print("source_zs: ", source_zs)
    # step_size = 6.12244897959e-07
    # step_size = 9.63265306122e-07
    step_size = 4.306122e-07
    step_size = 5e-7
    # step_size = 6.45454545455e-07
    first = True
    for source_z, z_lens in tqdm(list(zip(source_zs, z_lens_all))):
        rs = []
        thetas = []
        source_rs_array = []
        numerical_thetas = []
        dl = []
        dls = []
        ds = []
        for om in om_lambdas:
            comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om)
            source_r, dang_r = get_distances(source_z, Omega_Lambda=om)
            theta = rs2theta(source_r, comoving_lens, dang_lens)
            # dls, dl, ds
            print("calculated theta", source_r, dang_r)
            r, a = solve(theta, plot=False, comoving_lens=comoving_lens, Omega_Lambda=om, dt=step_size)
            thetas.append(theta)
            rs.append(r+comoving_lens)
            source_rs_array.append(source_r)
            num_theta = rs2theta(r+comoving_lens, comoving_lens, dang_lens)
            # num_theta2 = calc_theta(a*r, dang_lens, a*(r+comoving_lens))
            numerical_thetas.append(num_theta)
            ds.append((r+comoving_lens)*a)
            dls.append(r*a)
            dl.append(dang_lens)
        rs = np.array(rs)
        thetas = np.array(thetas)
        source_rs_array = np.array(source_rs_array)
        numerical_thetas = np.array(numerical_thetas)
        ds = np.array(ds)
        dls = np.array(dls)
        dl = np.array(dl)
        df = pd.DataFrame({'rs': rs, 'DL': dl, 'DLS': dls, 'DS': ds,'theta': thetas, 'rs_initial': source_rs_array, 'om_lambdas': om_lambdas, 'numerical_thetas': numerical_thetas, 'step': [step_size]*len(thetas)})
        filename = 'data/dopri5_diff_lambdas2.csv'
        # if first:
        #     df.to_csv(filename, index=False)
        #     first = False
        # else:
        #     df.to_csv(filename, index=False, header=False, mode='a')
    print("Time taken: {}".format(time.time() - start))

def main2():
    start = time.time()
    # thetas = np.linspace(25, 35, 10)*10**(-6)
    theta = 5e-6
    # print("thetas", thetas)
    # thetas = np.array([15e-6])
    # om_lambdas = np.linspace(0, 0.9, 1)
    om_lambdas = np.array([0])
    # om = 0
    z_lens = 0.05
    a_lens = 1/(z_lens+1)
    ds = []
    dls = []
    rs = []
    dl = []
    step_size = 1e-5
    # for theta in tqdm(thetas):
    for om in tqdm(om_lambdas):
        # print("omega_lambda:", om)
        comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om)
        print("lens r", comoving_lens)
        dl.append(dang_lens)
        # print("lens distances: ", comoving_lens, dang_lens)
        r, a = solve(theta, plot=True, comoving_lens=comoving_lens, Omega_Lambda=om, dt=step_size)
        d_s = a*(r + comoving_lens)
        d_ls = a * r
        ds.append(d_s)
        dls.append(d_ls)
        # rs.append(r)
    ds = np.array(ds)
    dls = np.array(dls)
    dl = np.array(dl)
    numerical_thetas = calc_theta(dls, dl, ds)
    percentage_errors = np.abs(numerical_thetas - theta)/theta*100
    # percentage_errors = np.abs(th - ds)/ds
    print(percentage_errors)
    print("Time taken: {}".format(time.time() - start))
    plt.plot(om_lambdas, percentage_errors, 'ro')

main2()