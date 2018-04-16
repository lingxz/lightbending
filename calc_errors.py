'''
sw_lambda, but generalized for curved universe
dimensions are in k
origin at lens
'''

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import time
import pandas as pd
from astropy.cosmology import LambdaCDM

length_scale = 3.086e22 # mega parsec
INTEGRATOR = 'vode'
INTEGRATOR_PARAMS = {
    'atol': 1e-110, 
    # 'atol': 0,
    # 'rtol': 0,
    'rtol': 1e-15,
    'nsteps': 100000000,
    'method': 'bdf',
}

tols = {
    'atol': 1e-120,
    'rtol': 1e-12,
    'method': 'RK45',
}

# H_0 = 7.33e-27
# H_0 = 7.56e-27 * length_scale
H_0 = 70*1000/(3.086e22)/299792458 * length_scale # 70km/s/Mpc
# Omega_Lambda = 0
# Omega_m = 1 - Omega_Lambda
# M = 0.5e15 / length_scale
M = 1474e13 / length_scale
print("M: ", M)

def frw(eta, w, p):
    L, Omega_k, Omega_Lambda, Omega_m, H_0 = p

    a, r, rdot, t, phi = w

    # Omega_k = 1 - Omega_Lambda - Omega_m
    k = omk2k(Omega_k, H_0)
    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda + Omega_k/a**2)

    phidot = L / (a*r)**2
    # tddot = -a*r**2*phidot**2*a_t - a*rdot**2*a_t/ (1 - k*r**2)
    # tdot = -np.sqrt(a**2*rdot**2+a**2*r**2*phidot**2)
    tdot = -np.sqrt(a**2*rdot**2/(1-k*r**2)+a**2*r**2*phidot**2)

    # rddot = r*phidot**2 - k*rdot**2 - 2*a_t/a*rdot*tdot
    rddot = (1-k*r**2)*r*phidot**2 - k*r*rdot**2/(1-k*r**2) - 2*a_t/a*rdot*tdot

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

def get_angle_frw_curved(r, phi, rdot, phidot, k):
    if phi > np.pi/2:
        angle = np.pi - phi
    else:
        angle = phi
    # res = np.arctan((rdot*np.sin(phi)+r*np.cos(phi)*phidot)/(rdot*np.cos(phi)-r*np.sin(phi)*phidot))
    res = np.arctan(np.abs(r*phidot/rdot)*np.sqrt(1-k*r**2)) - angle
    # res = np.arctan(np.abs(r*phidot/rdot)) - angle
    return res

def r2chi(k, r):
    if k == 0:
        return r
    if k > 0:
        return np.arcsin(np.sqrt(k)*r)/np.sqrt(k)
    if k < 0:
        return np.arcsinh(np.sqrt(-k)*r)/np.sqrt(-k)

def chi2r(k, chi):
    if k == 0:
        return chi
    if k > 0:
        return np.sin(chi*np.sqrt(k))/np.sqrt(k)
    if k < 0:
        return np.sinh(chi*np.sqrt(-k))/np.sqrt(-k)

def omk2k(om_k, H0):
    return -H0**2*om_k

def kottler(eta, w, p):
    E, L, M, Omega_Lambda, Omega_m, H_0, k, rh = p
    r_h, t, r, rdot, phi = w

    Lambda = 3*Omega_Lambda*H_0**2
    f = 1 - 2*M/r - Lambda/3*r**2
    rddot = L**2 * (r - 3*M) / r**4
    tdot = E / f
    phidot = L/r**2

    Ah = 1 - 2*M/r_h - Lambda/3*r_h**2
    r_h_t = Ah*np.sqrt(1-Ah/(1-k*rh**2))
    # r_h_t = (1 - 2*M/r_h - Lambda/3*r_h**2) * np.sqrt(2*M/r_h + Lambda/3*r_h**2)

    return [
        r_h_t*tdot,
        tdot,
        rdot,
        rddot,
        phidot,
    ]


def solve(angle_to_horizontal, comoving_lens=None, plot=True, Omega_Lambda=None, Omega_m=None, dt=None):
    a0 = 1
    initial_a = a0
    # initial_r = 10e17
    initial_r = comoving_lens
    initial_phi = np.pi
    initial_t = 0.
    Omega_k = 1 - Omega_m - Omega_Lambda
    k = omk2k(Omega_k, H_0)
    # initial_tdot = -1.
    # initial_rdot = -np.sqrt(initial_tdot**2/initial_a**2/(1+(np.tan(angle_to_horizontal))**2))
    # initial_phidot = initial_rdot * np.tan(angle_to_horizontal) / initial_r

    initial_rdot = -initial_r
    # print("1-kr^2", np.sqrt(1-k*initial_r**2))
    if plot:
        print("om_k", Omega_k, k)
        print("1-kr^2", np.sqrt(1-k*initial_r**2))
    initial_phidot = np.tan(angle_to_horizontal)*initial_rdot/initial_r/np.sqrt(1-k*initial_r**2)
    # initial_phidot = np.tan(angle_to_horizontal) * initial_rdot / initial_r
    initial_tdot = -np.sqrt(initial_a**2*initial_rdot**2/(1-k*initial_r**2) + initial_a**2*initial_r**2*initial_phidot**2)

    if plot:
        print("initial velocities:", initial_rdot, initial_phidot, initial_tdot)

    rho = Omega_m*3*H_0**2/(8*np.pi)
    r_h = 1/initial_a*(3*M/(4*np.pi*rho))**(1./3)
    if plot:
        print('r_h:', r_h, "\n")
    chi_h = r2chi(k, r_h)

    Lambda = 3*Omega_Lambda*H_0**2
    L_frw = (initial_a*initial_r)**2*initial_phidot

    p_frw = [L_frw, Omega_k, Omega_Lambda, Omega_m, H_0]
    initial = [initial_a, initial_r, initial_rdot, initial_t, initial_phi]
    
    # save for later frw straight line propagation
    initial0 = initial

    def reach_hole(t, y): return r2chi(k, y[1]) - chi_h
    reach_hole.terminal = True
    reach_hole.direction = -1

    sol = spi.solve_ivp(lambda t, y: frw(t, y, p_frw), [0, 10], initial, events=reach_hole, **tols)
    last = sol.y[:,-1]


    frw_angle_before_entering = get_angle_frw_curved(last[1], last[4], last[2], L_frw/(last[0]*last[1])**2, k)
    # frw_angle_before_entering = get_angle(last[1], last[4], last[2], L_frw/(last[0]*last[1])**2)

    if plot:
        print("Angle in FRW::")
        print(frw_angle_before_entering)
        print(last)
        print("\n")

    solver_kottler = spi.ode(kottler).set_integrator(INTEGRATOR, **INTEGRATOR_PARAMS)
    # a, r, rdot, t, phi = w
    r_out = last[0] * last[1]
    exit_rh = r_out
    tdot_out = -np.sqrt(last[0]**2*last[2]**2+last[0]**2*last[1]**2*L_frw/(last[0]*last[1])**2)
    initial_t = 0
    initial_r = r_out
    initial_phi = last[4]
    initial_rh = initial_r
    f = 1 - 2*M/r_out - Lambda/3*r_out**2
    etadot = tdot_out / last[0] # conformal time

    # a, r, rdot, t, phi
    frw_r = last[1]
    frw_a = last[0]
    if plot:
        print("sqrt", 1-k*frw_r**2, 2*M/r_out+Lambda/3*r_out**2-k*frw_r**2)
    jacobian = np.matrix([[1/f*np.sqrt(1-k*frw_r**2), frw_a/f/np.sqrt(1-k*frw_r**2)*np.sqrt(1-f-k*frw_r**2)], [np.sqrt(1-f-k*frw_r**2), frw_a]])
    vels = np.matmul(jacobian, np.array([tdot_out, last[2]]))
    initial_rdot, initial_tdot = vels[0, 1], vels[0, 0]

    # initial_rdot = last[0] * (np.sqrt(1 - f)*etadot + last[2]/np.sqrt(1-k*initial_r**2))
    # initial_tdot = last[0]/f*(etadot + np.sqrt(1-f)*last[2]/np.sqrt(1-k*initial_r**2))

    initial_phidot = L_frw / (last[0]*last[1])**2
    L_kottler = initial_phidot *initial_r**2

    initial_kottler = [initial_rh, initial_t, initial_r, initial_rdot, initial_phi]
    E = f*initial_tdot
    p_kottler = [E, L_kottler, M, Omega_Lambda, Omega_m, H_0, k, r_h]

    solver_kottler.set_initial_value(initial_kottler, 0).set_f_params(p_kottler)

    if plot:
        print("kottler initial:")
        print("initial_rh, initial_t, initial_r, initial_rdot, initial_phi")
        print(initial_kottler)
        print("p_kottler")
        print(p_kottler)

    angle_before_entering = get_angle(initial_r, initial_phi, initial_rdot, initial_phidot)
    if plot:
        print("light ray angle before entering kottler hole: ", angle_before_entering)
        print("initial conditions on entering kottler hole: ", initial_kottler)
        print("initial params on entering kottler hole: ", p_kottler)
    

    initial_kottler[2] -= 1e-10
    def exit_hole(t, y): return y[2] - y[0]
    def turning_point(t, y): return y[4] - np.pi/2 
    exit_hole.terminal = True
    exit_hole.direction = 1
    turning_point.terminal = False
    turning_point.direction = -1
    sol_kottler = spi.solve_ivp(lambda t, y: kottler(t, y, p_kottler), [0, 10], initial_kottler, dense_output=True, events=[turning_point, exit_hole], **tols)
    last = sol_kottler.y[:,-1]
    if plot:
        print("sol kottler", sol_kottler.t_events)
        print("turning point in kottler", sol_kottler.sol(sol_kottler.t_events[0])[2])
        print()

    # first_time = True
    # prev_r = np.inf
    # while solver_kottler.successful():
    #     solver_kottler.integrate(solver_kottler.t + dt, step=False)
    #     if solver_kottler.y[4] < np.pi/2 and first_time:
    #         if plot:
    #             print("turning point in kottler metric:", solver_kottler.y[2])
    #         first_time = False
    #     if solver_kottler.y[2] > solver_kottler.y[0]:
    #         last = solver_kottler.y
    #         break




    if plot:
        print("last in Kottler", last)
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

    # a, r, rdot, t, phi
    frw_r = r_h
    frw_a = initial_a
    jacobian = np.matrix([[1/f*np.sqrt(1-k*frw_r**2), frw_a/f/np.sqrt(1-k*frw_r**2)*np.sqrt(1-f-k*frw_r**2)], [np.sqrt(1-k*frw_r**2-f), frw_a]])
    inv_jac = np.linalg.inv(jacobian)
    vels = np.matmul(inv_jac, np.array([last_tdot, last[3]]))
    initial_rdot, initial_tdot = vels[0, 1], vels[0, 0]

    # initial_rdot = 1/initial_a*(1/f*last[3] - np.sqrt(1-f)*last_tdot) * np.sqrt(1-k*initial_r**2)
    # initial_etadot = 1/initial_a*(last_tdot - np.sqrt(1-f)/f*last[3])
    # initial_tdot = initial_etadot * initial_a


    p_frw[0] = initial_a**2 * initial_r**2*initial_phidot # change the L_frw
    # p_frw = [L_frw, Omega_k, Omega_Lambda, Omega_m, H_0]
    # r_h, t, r, rdot, phi = w

    initial_t = 0
    # a, t, r, rdot, phi 

    # frw_angle_after_exiting = get_angle(initial_r, initial_phi, initial_rdot, initial_phidot)
    frw_angle_after_exiting = get_angle_frw_curved(initial_r, initial_phi, initial_rdot, initial_phidot, k)

    if plot:
        print("Angle after exiting FRW:", frw_angle_after_exiting)
        print("bending angle in FRW: ", frw_angle_before_entering + frw_angle_after_exiting)
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

    enter_phi = initial_phi
    alpha = frw_angle_before_entering + frw_angle_after_exiting

    def reach_axis(t, y): return y[4] - 0
    reach_axis.terminal = True
    reach_axis.direction = -1
    sol2 = spi.solve_ivp(lambda t, y: frw(t, y, p_frw), [0, 10], initial_frw2, events=reach_axis, **tols)
    last = sol2.y[:,-1]

    if plot:
        print("last one frw", last)

    # frw_angle_after_exiting = get_angle_frw_curved(initial_r, initial_phi, initial_rdot, initial_phidot, k)

    # frw_angle_after_exiting = get_angle_frw_curved(last[1], last[4], last[2], p_frw[0]/(last[0]*last[1])**2, k)
    # alpha = angle_to_horizontal + frw_angle_after_exiting
    alpha  = get_angle_frw_curved(last[1], last[4], last[2], p_frw[0]/(last[0]*last[1])**2, k) + angle_to_horizontal
    return alpha, exit_rh, enter_phi, last[1], last[0]


def get_distances(z, Omega_Lambda=None, Omega_m=None):
    Omega_k = 1 - Omega_Lambda - Omega_m
    k = omk2k(Omega_k, H_0)
    # print(Omega_k, Omega_m, Omega_Lambda)
    def integrand(z):
        return 1/np.sqrt(Omega_m*(1+z)**3 + Omega_Lambda + Omega_k*(1+z)**2)
    integral, error = spi.quad(integrand, 0, z)
    chi = integral/H_0
    # comoving = chi
    comoving = chi2r(k, chi)
    dang = comoving/(1+z)
    return comoving, dang


from tqdm import tqdm
def main():
    start = time.time()
    om_lambdas = np.linspace(0., 0.99, 1)
    z_lens_all = np.linspace(.5, 1.2, 1)
    om_m = 1
    om = 0
    om_k = 0

    start_thetas = np.array([1/3600/180*np.pi]*50) # 1 arcsec
    step_size = 1e-9
    first = True
    filename = 'data/errors.csv'
    to_file = True
    print(filename)
    thetas = []
    dl = []
    ms = []
    com_lens = []
    exit_rhs = []
    enter_phis = []
    alphas = []
    om_ks = []

    raw_rs = []
    z_lens = 0.5
    theta = 1/3600/180*np.pi
    tolerances = np.linspace(1e-14, 1e-12, 1000)
    for t in tqdm(tolerances):
        global tols
        tols['rtol'] = t
        k = omk2k(om_k, H_0)
        ms.append(M)

        # comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om, Omega_m=om_m)
        cosmo = LambdaCDM(H0=70, Om0=om_m, Ode0=om)
        dang_lens = cosmo.angular_diameter_distance(z_lens).value
        comoving_lens = cosmo.comoving_transverse_distance(z_lens).value
        # print("comoving_lens", comoving_lens)

        dl.append(dang_lens)
        com_lens.append(comoving_lens)
        thetas.append(theta)
        alpha, exit_rh, enter_phi, raw_r, source_a = solve(theta, plot=False, comoving_lens=comoving_lens, Omega_Lambda=om, Omega_m=om_m, dt=step_size)
        exit_rhs.append(exit_rh)
        enter_phis.append(enter_phi)
        alphas.append(alpha)
        om_ks.append(om_k)
        raw_rs.append(raw_r)

        # R = dang_lens*theta
        # dls = cosmo.angular_diameter_distance_z1z2(z_lens, 1/source_a-1)
        # ds = cosmo.angular_diameter_distance_z1z2(0, 1/source_a-1)
        # print(dls, ds)
        # print("dang_s", get_distances(1./source_a-1, Omega_Lambda=om, Omega_m=om_m))
        # dang_s = source_a*chi2r(k, r2chi(k, raw_r) + r2chi(k, comoving_lens))
        # print("dang_s from numerical", dang_s, source_a*raw_r)
        # print("raw_rs", raw_r)

        # r0 = 1/(1/R + M/R**2 + 3/16*M**2/R**3)
        # print("r0", r0)
        # A_frw = 4*M/R + 15*np.pi*M**2/4/R**2 + 401/12*M**3/R**3
        # frw = comoving_lens/(A_frw/theta -1)
        # print("R", dang_lens*theta)
        # # alpha2 = (comoving_lens + raw_r)/raw_r*theta
        # # alpha2 = np.arctan(np.tan(theta)*chi2r(k, r2chi(k, raw_r)+r2chi(k, comoving_lens))/raw_r) + theta
        # alpha2 = chi2r(k, r2chi(k, raw_r)+r2chi(k, comoving_lens))/raw_r*theta
        # print("A_frw", A_frw, alpha)
        # print("compare", alpha/A_frw-1)
        # print(dang_s*theta/raw_r/6.62068423e-01)

        # chi_s = r2chi(k, r) + r2chi(k, comoving_lens)
        # d_s = a*chi2r(k, chi_s)
        # d_ls = a * r
        # ds.append(d_s)
        # dls.append(d_ls)

    thetas = np.array(thetas)
    dl = np.array(dl)
    ms = np.array(ms)
    com_lens = np.array(com_lens)
    exit_rhs = np.array(exit_rhs)
    enter_phis = np.array(enter_phis)
    alphas = np.array(alphas)
    om_ks = np.array(om_ks)
    raw_rs = np.array(raw_rs)

    # print("lengths", len(thetas), len(dl), len(ms), len(com_lens), len(exit_rhs), len(enter_phis), len(alphas))
    df = pd.DataFrame({
        'raw_rs': raw_rs,
        'tols': tolerances
    })
    if to_file:
        if first:
            df.to_csv(filename, index=False)
            first = False
        else:
            df.to_csv(filename, index=False, header=False, mode='a')

    print("Time taken: {}".format(time.time() - start))

main()