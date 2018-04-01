'''
origin at lens
'''

# diff_lambdas_small2: latest, diff z_lens, 100 per Omega_Lambda, but wrong DL
# diff_lambdas_small: same z_lens, diff z_source, 100 per Omega_Lambda
# diff_lambdas_small3: same as diff_lambdas_small2, but with correct DL
# diff_lambdas_small4: same as previous, but the step size in ltb solver is changed to agree with dt. Previously was fixed at 5e-7, oops!!!
# diff_lambdas_bigger_redshifts: bigger redshifts, 0.2 to 1 instead of 0.05 to 0.2.
# diff_lambdas_bigger_redshifts2： same as above, just more points

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import time
import pandas as pd
import random

length_scale = 3.086e22 # mega parsec
INTEGRATOR = 'vode'
INTEGRATOR_PARAMS = {
    # 'atol': 1e-110, 
    'atol': 1e-30,
    # 'rtol': 1e-3,
    'rtol': 3e-16,
    'nsteps': 100000000,
    # 'method': 'bdf',
}
# H_0 = 7.33e-27
H_0 = 7.56e-27 * length_scale
# Omega_Lambda = 0
# Omega_m = 1 - Omega_Lambda
# M = 0.5e15 / length_scale
M = 1474e12 / length_scale
rho_frw_initial = None
print("M: ", M)

def cart2spherical(x):
    return np.array([
        np.sqrt(x[0]**2 + x[1]**2),
        np.arctan(x[1]/x[0])
    ])

def conformal_time(a, Omega_m, Omega_Lambda):
    def tmp(a1):
        return 1/a1**2/np.sqrt(Omega_m/a1**3 + Omega_Lambda)
    result, err = spi.quad(tmp, 1, a, epsrel=1e-8, epsabs=1e-8)
    return result/H_0

def binary_search(start, end, answer, Omega_m, Omega_Lambda):
    mid = (end+start)/2
    res = conformal_time(mid, Omega_m, Omega_Lambda)
    # print(mid, res, answer)
    if np.isclose(res, answer, rtol=1e-10, atol=0):
        return mid
    if res < answer:
        return binary_search(mid, end, answer, Omega_m, Omega_Lambda)
    else:
        return binary_search(start, mid, answer, Omega_m, Omega_Lambda)


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

# def mass(r):
#     initial_rh = (3*M/(4*np.pi*rho_frw_initial))**(1./3)
#     rlimit = initial_rh/1200
#     rho_peak = 3*M/(np.pi*rlimit**3)
#     # integrand = lambda r: rho(r)*
#     # integral, error = spi.quad(lambda r1:rho(r1)*r1**2, 0, r)
#     if r > rlimit:
#         return M
#     else:
#         return 4*np.pi*(-(rho_peak/rlimit)*r**4/4 + rho_peak*r**3/3)
#     # if r > rlimit:
#     #     print(r/rlimit, M, integral)
#     # return 4*np.pi*integral

def mass(r):
    initial_rh = (3*M/(4*np.pi*rho_frw_initial))**(1./3)
    rlimit = initial_rh*0.5
    if r > rlimit:
        return M
    else:
        c = 10
        Rvir = rlimit/100
        Rs = Rvir/c
        rho0 = M/(4*np.pi*Rs**3*(np.log((Rs + rlimit)/Rs) - rlimit/(Rs + rlimit)))
        return 4*np.pi*rho0*Rs**3*(np.log((Rs+r)/Rs) - r/(Rs+r))

        # integral, error = spi.quad(lambda r1: rho(r1)*r1**2, 0, r)
        # return 4*np.pi*integral

def rho(r):
    initial_rh = (3*M/(4*np.pi*rho_frw_initial))**(1./3)
    rlimit = initial_rh*0.5
    if r > rlimit:
        return 0
    else:
        c = 10
        Rvir = rlimit/100
        Rs = Rvir/c
        rho0 = M/(4*np.pi*Rs**3*(np.log((Rs + rlimit)/Rs) - rlimit/(Rs + rlimit)))
        return rho0/(r/Rs)/(1 + r/Rs)**2

# def rho(r):
#     initial_rh = (3*M/(4*np.pi*rho_frw_initial))**(1./3)
#     rlimit = initial_rh/1200
#     # c = 10
#     # Rs = rlimit / c
#     # return rho0/(r/Rs)/(1+r/Rs)**2
#     rho_peak = 3*M/(np.pi*rlimit**3)
#     if r > rlimit:
#         return 0
#     else:
#         return rho_peak - (rho_peak/rlimit)*r

def ltb(eta, w, p):
    h, L, M, Omega_Lambda, Omega_m, H_0 = p
    r_h, t, r, rdot, phi, alpha, P = w

    current_m = mass(r)
    current_rho = rho(r)

    Lambda = 3*Omega_Lambda*H_0**2
    E = -2*current_m/r - Lambda*r**2/3
    tdot = h/alpha**2
    phidot = L/r**2

    P_r = (current_rho + P)/2/(1+E)/r*(Lambda*r**2 - 8*np.pi*P*r**2 + E)
    alpha_r = -alpha/2/(1+E)/r*(Lambda*r**2 - 8*np.pi*P*r**2 + E)
    
    E_r = 2*current_m/r**2 - 8*np.pi*current_rho*r - 2*Lambda*r/3
    # rddot = -alpha*alpha_r*(1+E)*tdot**2 + r*(1+E)*phidot**2 + E_r/2/(1+E)*r**2
    rddot = h**2/2*(E_r/alpha**2-2*(1+E)*alpha_r/alpha**3) - L**2/2*(E_r/r**2-2*(1+E)/r**3)

    r_h_t = (1 - 2*M/r_h - Lambda/3*r_h**2) * np.sqrt(2*M/r_h + Lambda/3*r_h**2)

    return [
        r_h_t*tdot,
        tdot,
        rdot,
        rddot,
        phidot,
        alpha_r*rdot,
        P_r*rdot,
    ]

def solve(angle_to_horizontal, comoving_lens=None, plot=True, Omega_Lambda=0, dt=None):
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

    rho_frw = Omega_m*3*H_0**2/(8*np.pi)
    r_h = 1/initial_a*(3*M/(4*np.pi*rho_frw))**(1./3)
    global rho_frw_initial
    rho_frw_initial = rho_frw
    if r_h >= initial_r:
        print("Starting point is inside the hole! Make the hole smaller or bring the lens further away.")
        return
    if plot:
        print("rho_frw_initial", rho_frw_initial)
        print('r_h:', r_h, "\n")

    Lambda = 3*Omega_Lambda*H_0**2
    L_frw = (initial_a*initial_r)**2*initial_phidot

    p_frw = [L_frw, k, Omega_Lambda, Omega_m, H_0]
    initial = [initial_a, initial_r, initial_rdot, initial_t, initial_phi]
    
    # save for later frw straight line propagation
    initial0 = initial

    x_initial = np.array([0, 0])
    velocity = np.sqrt((initial_rdot*np.cos(initial_phi)-initial_r*np.sin(initial_phi)*initial_phidot)**2 + (initial_rdot*np.sin(initial_phi)+initial_r*np.cos(initial_phi)*initial_phidot)**2)
    direction = np.array([
        initial_rdot*np.cos(initial_phi)-initial_r*np.sin(initial_phi)*initial_phidot, 
        initial_rdot*np.sin(initial_phi)+initial_r*np.cos(initial_phi)*initial_phidot,
    ])
    direction = direction / velocity

    # direction = np.array([
    #     velocity*np.cos(angle_to_horizontal),
    #     velocity*np.sin(angle_to_horizontal),
    # ])

    delta_eta = np.roots([direction[0]**2+direction[1]**2, -2*direction[0]*comoving_lens, comoving_lens**2-r_h**2])

    if plot:
        print("delta_eta", delta_eta)
    delta_eta = sorted(delta_eta[delta_eta > 0])[0]
    eta_initial = 0
    eta_out = eta_initial - delta_eta
    if plot:
        print("eta_out", eta_out, delta_eta)

    x_final = x_initial + (eta_initial - eta_out)*direction
    x_final[0] -= comoving_lens
    r_final, phi_final = cart2spherical(x_final)
    phi_final = np.pi + phi_final  # change from 4th to second quadrant
    # print("conformal_time", conformal_time(0.9997899999999995, Omega_m, Omega_Lambda))
    a_final = binary_search(0.5, 1, eta_out, Omega_m, Omega_Lambda)
    if plot:
        print("a_final result", conformal_time(a_final, Omega_m, Omega_Lambda))
    phidot_final = L_frw/(a_final*r_final)**2
    rdot_final = (r_final*np.cos(phi_final)*phidot_final+r_final*np.sin(phi_final)*phidot_final*np.tan(angle_to_horizontal))/(np.cos(phi_final)*np.tan(angle_to_horizontal)-np.sin(phi_final))
    # rdot_final = r_final*phidot_final/np.tan(angle_to_horizontal)
    tdot_final = -np.sqrt(a_final**2*rdot_final**2 + a_final**2*r_final**2*phidot_final**2)
    t_final = 0 # not important

    if plot:
        print("phi_final", phi_final)
        print("a_final", a_final)

    last = [a_final, r_final, rdot_final, t_final, phi_final]

    if plot:
        frw_angle_before_entering = get_angle(last[1], last[4], last[2], L_frw/(last[0]*last[1])**2)
        print("Angle in FRW::")
        print(frw_angle_before_entering)
        print(last)
        print("\n")


    solver_ltb = spi.ode(ltb).set_integrator(INTEGRATOR, **INTEGRATOR_PARAMS)
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
    initial_rdot = last[0] * (np.sqrt(1 - f)*etadot + last[2])
    initial_tdot = last[0]/f*(etadot + np.sqrt(1-f)*last[2])
    initial_phidot = L_frw / (last[0]*last[1])**2
    L_ltb = initial_phidot *initial_r**2

    initial_alpha = np.sqrt(f)
    initial_P = 0
    # r_h, t, r, rdot, phi, alpha, P
    initial_ltb = [initial_rh, initial_t, initial_r, initial_rdot, initial_phi, initial_alpha, initial_P]
    h = f*initial_tdot
    p_ltb = [h, L_ltb, M, Omega_Lambda, Omega_m, H_0]

    solver_ltb.set_initial_value(initial_ltb, 0).set_f_params(p_ltb)

    if plot:
        print("LTB initial:")
        print("initial_rh, initial_t, initial_r, initial_rdot, initial_phi, initial_alpha, initial_P")
        print(initial_ltb)

    angle_before_entering = get_angle(initial_r, initial_phi, initial_rdot, initial_phidot)
    if plot:
        print("light ray angle before entering ltb hole: ", angle_before_entering)
        print("initial conditions on entering ltb hole: ", initial_ltb)
        print("initial params on entering ltb hole: ", p_ltb)

    initial_ltb[2] -= 1e-10
    def exit_hole(t, y): return y[2] - y[0]
    def turning_point(t, y): return y[4] - np.pi/2 
    exit_hole.terminal = True
    exit_hole.direction = 1
    turning_point.terminal = False
    turning_point.direction = -1
    sol_ltb = spi.solve_ivp(lambda t, y: ltb(t, y, p_ltb), [0, 100], initial_ltb, dense_output=True, method='RK45', events=[turning_point, exit_hole], rtol=1e-20, atol=1e-30)
    last = sol_ltb.y[:,-1]
    if plot:
        print("sol ltb", sol_ltb.t_events)
        print(sol_ltb)
        print("turning point in ltb", sol_ltb.sol(sol_ltb.t_events[0])[2])
        print()


    # # sol_ltb = []
    # first_time = True
    # prev_r = np.inf
    # while solver_ltb.successful():
    #     solver_ltb.integrate(solver_ltb.t + dt, step=False)
    #     if solver_ltb.y[4] < np.pi/2 and first_time:
    #         # dt = 5e-9
    #         if plot:
    #             print("turning point in ltb metric:", solver_ltb.y[2])
    #         first_time = False
    #     if solver_ltb.y[2] > solver_ltb.y[0]:
    #         last = solver_ltb.y
    #         break

    angle_after_exiting = get_angle(last[2], last[4], last[3], L_ltb/last[2]**2)
    if plot:
        print("last ltb", last)
        print("light ray angle after exiting ltb hole: ", angle_after_exiting)
        print("bending angle in ltb: ", angle_before_entering - angle_after_exiting)
        print("\n")

    # r_h, t, r, rdot, phi = w
    initial_phi = last[4]
    initial_r = r_h
    initial_a = last[2] / r_h
    if plot:
        print("time spent in hole: ", last[1])
        print("scale factor on exit:", initial_a)
        print("exit r in FRW:", initial_r)

    initial_phidot = L_ltb / last[2]**2
    f = 1-2*M/last[2]-Lambda/3*last[2]**2
    last_tdot = h/f
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
    initial_t = 0
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
        print("initial_ydot", initial_ydot)
        print("----")

    if initial_r*np.sin(initial_phi) < 0:
        print("light ray is bent too much by the hole, increase angle_to_horizontal")
        print("initial angle to horizontal: ", angle_to_horizontal)
        print("----")

    ans = initial_r*np.sin(initial_phi)/np.tan(np.abs(frw_angle_after_exiting)) + initial_r*np.cos(initial_phi)
    return ans, exit_rh


    return ans, exit_rh
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
    om_lambdas = np.linspace(0., 0.9, 1)
    # om_lambdas = np.array([0.99])
    z_lens_all = np.linspace(0.1, 0.5, 1)
    # z_lens = 0.1
    # a_lens = 1/(z_lens+1)
    # start_thetas = np.linspace(0.7e-5, 2e-5, 100)
    start_thetas = np.array([8e-7]*50)
    step_size = 1e-8
    first = True
    filename = 'data/ltb_cheese_throwaway.csv'
    print(filename)
    for theta, z_lens in tqdm(list(zip(start_thetas, z_lens_all))):
        rs = []
        thetas = []
        dl = []
        ms = []
        raw_rs = []
        com_lens = []
        exit_rhs = []
        for om in tqdm(om_lambdas):
            # global M
            # rho = (1-om)*3*H_0**2/(8*np.pi)
            # M = 4/3*np.pi*rho*r_h**3
            ms.append(M)
            comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om)
            # print("lens r", comoving_lens, theta)
            com_lens.append(comoving_lens)
            r, exit_rh = solve(theta, plot=True, comoving_lens=comoving_lens, Omega_Lambda=om, dt=step_size)
            exit_rhs.append(exit_rh)
            # print("result", r)

            # A_frw = 4*M/(dang_lens*theta) + 15*np.pi*M**2/4/(dang_lens*theta)**2 + 401/12*M**3/(dang_lens*theta)**3
            # frw = comoving_lens/(A_frw/theta -1)
            # print("R", dang_lens*theta)
            # print("compare", r, frw, r/frw-1)
            # print("A_frw", A_frw)

            thetas.append(theta)
            raw_rs.append(r)
            rs.append(r+comoving_lens)
            dl.append(dang_lens)
        rs = np.array(rs)
        thetas = np.array(thetas)
        dl = np.array(dl)
        ms = np.array(ms)
        raw_rs = np.array(raw_rs)
        com_lens = np.array(com_lens)
        exit_rhs = np.array(exit_rhs)

        # DL,DLS,DS,M,numerical_thetas,om_lambdas,rs,rs_initial,step,theta,z_lens,comoving_lens,raw_rs
        df = pd.DataFrame({
            'DL': dl,
            'M': ms, 
            'numerical_thetas': thetas,
            'om_lambdas': om_lambdas, 
            'rs': rs, 
            # 'rs_initial': [2]*len(thetas), 
            'step': [step_size]*len(thetas),
            'theta': thetas, 
            'z_lens': [z_lens]*len(thetas),
            'comoving_lens': com_lens,
            'raw_rs': raw_rs,
            'exit_rhs': exit_rhs,
        })
        # df = df[['DL','M','comoving_lens','numerical_thetas','om_lambdas','raw_rs','rs','rs_initial','step','theta','z_lens']]
        if first:
            df.to_csv(filename, index=False)
            first = False
        else:
            df.to_csv(filename, index=False, header=False, mode='a')
    print("Time taken: {}".format(time.time() - start))

main()