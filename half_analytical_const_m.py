'''
origin at lens
'''

# diff_lambdas_small2: latest, diff z_lens, 100 per Omega_Lambda, but wrong DL
# diff_lambdas_small: same z_lens, diff z_source, 100 per Omega_Lambda
# diff_lambdas_small3: same as diff_lambdas_small2, but with correct DL
# diff_lambdas_small4: same as previous, but the step size in kottler solver is changed to agree with dt. Previously was fixed at 5e-7, oops!!!
# diff_lambdas_bigger_redshifts: bigger redshifts, 0.2 to 1 instead of 0.05 to 0.2.
# diff_lambdas_bigger_redshifts2： same as above, just more points

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import time
import pandas as pd


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
# H_0 = 7.33e-27
H_0 = 7.56e-27 * length_scale
# Omega_Lambda = 0
# Omega_m = 1 - Omega_Lambda
# M = 0.5e15 / length_scale
M = 1474e12 / length_scale
tdot_to_steps_ratio = None

print("M: ", M)

# def frw(eta, w, p):
#     L, k, Omega_Lambda, Omega_m, H_0 = p

#     a, r, rdot, t, phi = w

#     a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda)

#     phidot = L / (a*r)**2
#     # tddot = -a*r**2*phidot**2*a_t - a*rdot**2*a_t/ (1 - k*r**2)
#     tdot = -np.sqrt(a**2*rdot**2+a**2*r**2*phidot**2)
#     rddot = r*phidot**2 - k*rdot**2 - 2*a_t/a*rdot*tdot

#     return [
#         a_t*tdot,
#         rdot,
#         rddot,
#         tdot,
#         phidot,
#     ]

def get_angle(r, phi, rdot, phidot):
    res = np.arctan((rdot*np.sin(phi)+r*np.cos(phi)*phidot)/(rdot*np.cos(phi)-r*np.sin(phi)*phidot))
    return res

def kottler(eta, w, p):
    E, L, M, Omega_Lambda, Omega_m, H_0 = p
    r_h, t, r, rdot, phi = w

    Lambda = 3*Omega_Lambda*H_0**2
    # Lambda = 0
    f = 1 - 2*M/r - Lambda/3*r**2
    rddot = L**2 * (r - 3*M) / r**4
    tdot = E / f
    phidot = L/r**2
    r_h_t = (1 - 2*M/r_h - Lambda/3*r_h**2) * np.sqrt(2*M/r_h + Lambda/3*r_h**2)
    # r_h_t = (1 - 2*M/r_h) * np.sqrt(2*M/r_h)

    return [
        r_h_t*tdot,
        tdot,
        rdot,
        rddot,
        phidot,
    ]

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

def solve(angle_to_horizontal, comoving_lens=None, plot=True, Omega_Lambda=0, dt=None):
    k = 0
    a0 = 1
    initial_a = a0
    # initial_r = 10e17
    initial_r = comoving_lens
    initial_phi = np.pi
    initial_t = 0.
    Omega_m = 1 - Omega_Lambda



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
    phidot_final = L_frw/(a_final*r_final)**2
    rdot_final = (r_final*np.cos(phi_final)*phidot_final+r_final*np.sin(phi_final)*phidot_final*np.tan(angle_to_horizontal))/(np.cos(phi_final)*np.tan(angle_to_horizontal)-np.sin(phi_final))
    # rdot_final = r_final*phidot_final/np.tan(angle_to_horizontal)
    tdot_final = -np.sqrt(a_final**2*rdot_final**2 + a_final**2*r_final**2*phidot_final**2)
    t_final = 0 # not important

    if plot:
        print("phi_final", phi_final)
        print("a_final", a_final)

    last = [a_final, r_final, rdot_final, t_final, phi_final]

    # res = np.arctan((rdot*np.sin(phi)+r*np.cos(phi)*phidot)/(rdot*np.cos(phi)-r*np.sin(phi)*phidot))


    if plot:
        print("angle_to_horizontal", angle_to_horizontal)
        frw_angle_before_entering = get_angle(last[1], last[4], last[2], L_frw/(last[0]*last[1])**2)
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
    initial_rdot = last[0] * (np.sqrt(1 - f)*etadot + last[2])
    initial_tdot = last[0]/f*(etadot + np.sqrt(1-f)*last[2])
    initial_phidot = L_frw / (last[0]*last[1])**2
    L_kottler = initial_phidot *initial_r**2

    initial_kottler = [initial_rh, initial_t, initial_r, initial_rdot, initial_phi]
    E = f*initial_tdot
    p_kottler = [E, L_kottler, M, Omega_Lambda, Omega_m, H_0]

    # global tdot_to_steps_ratio
    # if not tdot_to_steps_ratio:
    #     tdot_to_steps_ratio = initial_tdot/dt
    # else:
    #     dt = initial_tdot / tdot_to_steps_ratio
    # print("tdot_to_steps_ratio", tdot_to_steps_ratio, dt)
    # print()

    solver_kottler.set_initial_value(initial_kottler, 0).set_f_params(p_kottler)

    if plot:
        print("kottler initial:")
        print("initial_rh, initial_t, initial_r, initial_rdot, initial_phi")
        print(initial_kottler)
        print("initial_tdot", initial_tdot)

    angle_before_entering = get_angle(initial_r, initial_phi, initial_rdot, initial_phidot)
    if plot:
        print("light ray angle before entering kottler hole: ", angle_before_entering)
        print("initial conditions on entering kottler hole: ", initial_kottler)
        print("initial params on entering kottler hole: ", p_kottler)

    first_time = True
    prev_r = np.inf
    while solver_kottler.successful():
        solver_kottler.integrate(solver_kottler.t + dt, step=False)
        # sol_kottler.append(list(solver_kottler.y))

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
        # if solver_kottler.y[4] < np.pi/2 and np.isclose(solver_kottler.y[2], solver_kottler.y[0], rtol=1e-5):
        #     last = solver_kottler.y
        #     print("here??")
        #     break
        if solver_kottler.y[2] > solver_kottler.y[0]:
            # print("exit hole at", solver_kottler.y[2], solver_kottler.y[0], solver_kottler.y[2]-solver_kottler.y[0])
            last = solver_kottler.y
            break

    angle_after_exiting = get_angle(last[2], last[4], last[3], L_kottler/last[2]**2)
    if plot:
        print("light ray angle after exiting kottler hole: ", angle_after_exiting)
        print("bending angle in kottler: ", angle_before_entering - angle_after_exiting)
        print("\n")

    if plot:
        print("time in hole: ", last[1])

    # r_h, t, r, rdot, phi = w
    initial_phi = last[4]
    initial_r = r_h
    initial_a = last[2] / r_h
    if plot:
        print("scale factor on exit:", initial_a)
        print("exit r in FRW:", initial_r)
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
        print("----")

    if initial_r*np.sin(initial_phi) < 0:
        print("light ray is bent too much by the hole, increase angle_to_horizontal")
        print("initial angle to horizontal: ", angle_to_horizontal)
        print("----")

    ans = initial_r*np.sin(initial_phi)/np.tan(np.abs(frw_angle_after_exiting)) + initial_r*np.cos(initial_phi)
    return ans, exit_rh
    # return r[-1], source_a

# [ 0.18075975  0.05122652  0.23824122  0.31020329  0.20010044  0.39768539
#   0.43731914  0.46608477  0.37572935  0.39980157]

def get_distances(z, Omega_Lambda=0):
    Omega_m = 1 - Omega_Lambda
    def integrand(z):
        return 1/np.sqrt(Omega_m*(1+z)**3 + Omega_Lambda)
    integral, error = spi.quad(integrand, 0, z, epsrel=1e-8, epsabs=1e-8)
    comoving = integral/H_0
    dang = comoving/(1+z)
    return comoving, dang

def calc_theta(D_LS, D_L, D_S):
    return np.sqrt(4*M*D_LS/D_L/D_S)

def theta2rls_flat(theta, z_lens):
    rl, dl = get_distances(z_lens, Omega_Lambda=0)
    return 4*calculate_mass(0)*rl/(4*calculate_mass(0)-dl*theta**2)

def rs2theta(rs, rl, dl):
    return np.sqrt(4*M*(rs-rl)/rs/dl)

def rs2redshift_flat(rs):
    return 1/(1-H_0*rs/2)**2 -1

def calculate_mass(om):
    rho = (1-om)*3*H_0**2/(8*np.pi)
    return 4/3*np.pi*rho*r_h**3

from tqdm import tqdm

def main():
    start = time.time()
    om_lambdas = np.linspace(0, 0.995, 50)
    # om_lambdas = np.array([0.99])
    z_lens_all = np.linspace(0.5, 1., 50)
    # z_lens = 0.1
    # a_lens = 1/(z_lens+1)
    # start_thetas = np.linspace(0.7e-5, 2e-5, 100)
    start_thetas = np.array([7e-7]*50)
    # source_rs = np.array([theta2rls_flat(th1, z1) for th1, z1 in zip(start_thetas, z_lens_all)])
    # source_rs = theta2rls_flat(start_thetas, z_lens)
    
    # source_zs = rs2redshift_flat(source_rs)
    # print("source_zs: ", source_zs)

    # step_size = 4.67346938776e-07
    step_size = 1e-9
    first = True
    filename = 'data/half_analytical_const_m3.csv'
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
            ms.append(M)
            comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om)
            # source_r, dang_r = get_distances(source_z, Omega_Lambda=om)
            # theta = rs2theta(source_r, comoving_lens, dang_lens)
            # print("lens r", comoving_lens, theta)
            com_lens.append(comoving_lens)
            r, exit_rh = solve(theta, plot=False, comoving_lens=comoving_lens, Omega_Lambda=om, dt=step_size)
            exit_rhs.append(exit_rh)
            # A_frw = 4*M/(dang_lens*theta) + 15*np.pi*M**2/4/(dang_lens*theta)**2 + 401/12*M**3/(dang_lens*theta)**3
            # frw = comoving_lens/(A_frw/theta -1)
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
            # 'rs_initial': source_rs_array, 
            'step': [step_size]*len(thetas),
            'theta': thetas, 
            'z_lens': [z_lens]*len(thetas),
            'comoving_lens': com_lens,
            'raw_rs': raw_rs,
            'exit_rhs': exit_rhs,
        })
        # df = df[['DL','DLS','DS','M','numerical_thetas','om_lambdas','rs','rs_initial','step','theta','z_lens','comoving_lens','raw_rs']]
        if first:
            df.to_csv(filename, index=False)
            first = False
        else:
            df.to_csv(filename, index=False, header=False, mode='a')
    print("Time taken: {}".format(time.time() - start))

main()
# main2()
# plt.show()
# plt.savefig('images/lambda.png')