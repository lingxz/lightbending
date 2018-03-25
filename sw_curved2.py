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
    rddot = (1-k*r**2)*r*phidot**2 - k*rdot**2/(1-k*r**2) - 2*a_t/a*rdot*tdot

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
        return np.sin(np.sqrt(k)*chi)/np.sqrt(k)
    if k < 0:
        return np.sinh(np.sqrt(-k)*chi)/np.sqrt(-k)

def omk2k(om_k, H0):
    return H0**2*om_k

def kottler2(eta, w, p):
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

def omega_lambda2lambda(Omega_Lambda):
    return 3*Omega_Lambda*H_0**2

def kantowski_alpha(R, phi, Omega_Lambda):
    r0 = 1/(1/R + M/R**2)
    Lambda = omega_lambda2lambda(Omega_Lambda)
    rs = 2*M
    first_term = (rs/2/r0)*np.cos(phi)*(-4*(np.cos(phi))**2 - 12*np.cos(phi)*np.sin(phi)*np.sqrt(Lambda*r0**2/3+rs/r0*(np.sin(phi))**3) + Lambda*r0**2*(8/3-20/3*(np.sin(phi))**2))
    second_term = (rs/2/r0)**2*(15/4*(2*phi-np.pi) + np.cos(phi)*(4+33/2*np.sin(phi)-4*(np.sin(phi))**2+19*(np.sin(phi))**3-64*(np.sin(phi))**5) - 12*np.log(np.tan(phi/2))*(np.sin(phi))**3)
    return first_term + second_term

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
    a0 = 1
    initial_a = a0
    # initial_r = 10e17
    initial_r = comoving_lens
    initial_phi = np.pi
    initial_t = 0.
    Omega_m = 1 - Omega_Lambda
    Omega_k = 0
    k = omk2k(Omega_k, H_0)


    initial_rdot = -initial_r
    # initial_phidot = np.tan(angle_to_horizontal) * initial_rdot / initial_r
    initial_phidot = np.tan(angle_to_horizontal) * initial_rdot / initial_r * np.sqrt(1-k*initial_r**2)
    initial_R = initial_a*initial_r
    initial_tdot = -np.sqrt(initial_a**2/(1-k*initial_r**2)*initial_rdot**2 + initial_R**2*initial_phidot**2)

    if plot:
        print("initial velocities:", initial_rdot, initial_phidot, initial_tdot)

    rho = Omega_m*3*H_0**2/(8*np.pi)
    r_h = 1/initial_a*(3*M/(4*np.pi*rho))**(1./3)
    chi_h = r2chi(k, r_h)
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

    def reach_hole(t, y): return r2chi(k, y[1]) - chi_h
    reach_hole.terminal = True
    reach_hole.direction = -1

    sol = spi.solve_ivp(lambda t, y: frw(t, y, p_frw), [0, 10], initial, method='RK45', events=reach_hole, rtol=1e-15, atol=1e-120)
    last = sol.y[:,-1]


    # res = np.arctan((rdot*np.sin(phi)+r*np.cos(phi)*phidot)/(rdot*np.cos(phi)-r*np.sin(phi)*phidot))

    frw_angle_before_entering = get_angle_frw_curved(last[1], last[4], last[2], L_frw/(last[0]*last[1])**2, k)
    if plot:
        print("angle_to_horizontal", angle_to_horizontal)
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
    jacobian = np.matrix([[1/f*np.sqrt(1-k*frw_r**2), frw_a/f/np.sqrt(1-k*frw_r**2)*np.sqrt(1-f-k*frw_r**2)], [np.sqrt(1-k*frw_r**2-f), frw_a]])
    vels = np.matmul(jacobian, np.array([tdot_out, last[2]]))
    initial_rdot, initial_tdot = vels[0, 1], vels[0, 0]


    # initial_rdot = last[0] * (np.sqrt(1 - f)*etadot + last[2])
    # initial_tdot = last[0]/f*(etadot + np.sqrt(1-f)*last[2])


    initial_phidot = L_frw / (last[0]*last[1])**2
    L_kottler = initial_phidot *initial_r**2

    initial_kottler = [initial_rh, initial_t, initial_r, initial_rdot, initial_phi]
    E = f*initial_tdot
    p_kottler = [E, L_kottler, M, Omega_Lambda, Omega_m, H_0, k, r_h]

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
        print("p_kottler")
        print(p_kottler)

    angle_before_entering = get_angle_frw_curved(initial_r, initial_phi, initial_rdot, initial_phidot, k)
    if plot:
        print("light ray angle before entering kottler hole: ", angle_before_entering)
        print("initial conditions on entering kottler hole: ", initial_kottler)
        print("initial params on entering kottler hole: ", p_kottler)

    first_time = True
    prev_r = np.inf
    while solver_kottler.successful():
        solver_kottler.integrate(solver_kottler.t + dt, step=False)
        if solver_kottler.y[4] < np.pi/2 and first_time:
            if plot:
                print("turning point in kottler metric:", solver_kottler.y[2])
            first_time = False
        if solver_kottler.y[2] > solver_kottler.y[0]:
            last = solver_kottler.y
            break
    print("last in Kottler", last)
    angle_after_exiting = get_angle(last[2], last[4], last[3], L_kottler/last[2]**2)
    if plot:
        print("light ray angle after exiting kottler hole: ", angle_after_exiting)
        print("bending angle in kottler: ", angle_before_entering - angle_after_exiting)
        print("\n")

    if plot:
        print("time in hole: ", last[1])

    exit_phi = last[4]
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

    # a, r, rdot, t, phi
    frw_r = r_h
    frw_a = initial_a
    jacobian = np.matrix([[1/f*np.sqrt(1-k*frw_r**2), frw_a/f/np.sqrt(1-k*frw_r**2)*np.sqrt(1-f-k*frw_r**2)], [np.sqrt(1-k*frw_r**2-f), frw_a]])
    inv_jac = np.linalg.inv(jacobian)
    vels = np.matmul(inv_jac, np.array([last_tdot, last[3]]))
    initial_rdot, initial_tdot = vels[0, 1], vels[0, 0]

    # initial_rdot = 1/initial_a*(1/f*last[3] - np.sqrt(1-f)*last_tdot)
    # initial_etadot = 1/initial_a*(last_tdot - np.sqrt(1-f)/f*last[3])
    # initial_tdot = initial_etadot * initial_a


    p_frw[0] = initial_a**2 * initial_r**2*initial_phidot # change the L_frw
    # p_frw = [L_frw, k, Omega_Lambda, Omega_m, H_0]
    # r_h, t, r, rdot, phi = w

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

    alpha = frw_angle_before_entering - frw_angle_after_exiting
    ans2 = comoving_lens*np.tan(angle_to_horizontal)/np.tan(np.abs(frw_angle_after_exiting))
    ans3 = comoving_lens*angle_to_horizontal/(np.abs(frw_angle_after_exiting))
    ans = initial_r*np.sin(initial_phi)/np.tan(np.abs(frw_angle_after_exiting)) + initial_r*np.cos(initial_phi)
    # print("anss", ans, ans3, ans2)
    return ans, exit_rh, exit_phi, alpha
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
    om_lambdas = np.linspace(0.5, 0.99, 1)
    # om_lambdas = np.array([0.])
    z_lens_all = np.linspace(1., 1.2, 1)
    # z_lens = 0.1
    # a_lens = 1/(z_lens+1)
    # start_thetas = np.linspace(0.7e-5, 2e-5, 100)
    start_thetas = np.array([1e-5]*50)
    # source_rs = np.array([theta2rls_flat(th1, z1) for th1, z1 in zip(start_thetas, z_lens_all)])
    # source_rs = theta2rls_flat(start_thetas, z_lens)
    
    # source_zs = rs2redshift_flat(source_rs)
    # print("source_zs: ", source_zs)

    # step_size = 4.67346938776e-07
    step_size = 1e-9
    first = True
    filename = 'data/half_analytical_throwaway.csv'
    print(filename)
    for theta, z_lens in tqdm(list(zip(start_thetas, z_lens_all))):
        rs = []
        thetas = []
        dl = []
        ms = []
        raw_rs = []
        com_lens = []
        exit_rhs = []
        exit_phis = []
        alphas = []
        for om in tqdm(om_lambdas):
            ms.append(M)
            comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om)
            # source_r, dang_r = get_distances(source_z, Omega_Lambda=om)
            # theta = rs2theta(source_r, comoving_lens, dang_lens)
            # print("lens r", comoving_lens, theta)
            com_lens.append(comoving_lens)
            r, exit_rh, exit_phi, alpha = solve(theta, plot=True, comoving_lens=comoving_lens, Omega_Lambda=om, dt=step_size)
            exit_rhs.append(exit_rh)
            exit_phis.append(exit_phi)
            alphas.append(alpha)

            A_frw = 4*M/(dang_lens*theta) + 15*np.pi*M**2/4/(dang_lens*theta)**2 + 401/12*M**3/(dang_lens*theta)**3
            print("A_frw", A_frw, alpha, alpha/A_frw-1)
            k_alpha = -kantowski_alpha(dang_lens*theta, exit_phi, om)
            print("k_alpha", k_alpha)

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
        exit_phis = np.array(exit_phis)
        alphas = np.array(alphas)

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
            'enter_phis': exit_phis,
            'alphas': alphas,
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