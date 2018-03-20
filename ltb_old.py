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


length_scale = 3.086e22 # mega parsec
INTEGRATOR = 'vode'
INTEGRATOR_PARAMS = {
    # 'atol': 1e-110, 
    'atol': 1e-20,
    # 'rtol': 1e-3,
    'rtol': 1e-15,
    'nsteps': 100000000,
    # 'method': 'adams',
}
# H_0 = 7.33e-27
H_0 = 7.56e-27 * length_scale
# Omega_Lambda = 0
# Omega_m = 1 - Omega_Lambda
# M = 0.5e15 / length_scale
M = 1474e12 / length_scale
rho_frw_initial = None
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
    rlimit = initial_rh*0.7
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
    rlimit = initial_rh*0.7
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
        if solver_frw.y[4] < np.pi/2:
            print("Light ray is not going to enter the hole, reduce theta!")
        if solver_frw.y[1] <= r_h:
            last = solver_frw.y
            break

    # sol = sol
    # last = sol[-1]
    # print("FRW coordinates, before entering ltb:")
    
    # sol1 = np.array(sol)

    # angle = get_angle(sol1[:,1], sol1[:,4], sol1[:,2], L_frw/(sol1[:,0]*sol1[:,1])**2)
    # plt.plot(sol1[:,1]*np.cos(sol1[:,4]), angle)
    # plt.show()

    if plot:
        frw_angle_before_entering = get_angle(last[1], last[4], last[2], L_frw/(last[0]*last[1])**2)
        print("Angle in FRW::")
        print(frw_angle_before_entering)
        print(last)
        print("\n")


    solver_ltb = spi.ode(ltb).set_integrator(INTEGRATOR, **INTEGRATOR_PARAMS)
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
    # sol_ltb = []
    first_time = True
    prev_r = np.inf
    while solver_ltb.successful():
        # dt = 1e-6
        solver_ltb.integrate(solver_ltb.t + dt, step=False)
        # sol_ltb.append(list(solver_ltb.y))

        # if solver_ltb.y[2] > prev_r and first_time:
        #     if plot:
        #         print("turning point in ltb metric:", solver_ltb.y[2])
        #     first_time = False
        # else:
        #     prev_r = solver_ltb.y[2]
        if solver_ltb.y[4] < np.pi/2 and first_time:
            if plot:
                print("turning point in ltb metric:", solver_ltb.y[2])
            first_time = False
        if solver_ltb.y[2] > solver_ltb.y[0]:
            last = solver_ltb.y
            break

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

    # if plot:
    #     x = r * np.cos(phi)
    #     y = r * np.sin(phi)
    #     plt.plot(x, y, 'bo')
        
    #     print("axis crossing points", r[-1], phi[-1])
    #     swiss_cheese_hole = plt.Circle((0., 0.), r_h, color='grey', fill=False, zorder=10)
    #     axes = plt.gca()
    #     axes.add_artist(swiss_cheese_hole)
    #     axes.set_aspect('equal', adjustable='box')
    #     # plot the center
    #     plt.plot(0, 0, 'ro')

    #     lim = 0.5
    #     axes.set_xlim([-lim, lim])
    #     axes.set_ylim([-lim, lim])

    #     plt.figure()
    #     sol_ltb = np.array(sol_ltb)
    #     r_ltb = sol_ltb[:,2]
    #     phi_ltb = sol_ltb[:,4]
    #     x_ltb = r_ltb * np.cos(phi_ltb)
    #     y_ltb = r_ltb * np.sin(phi_ltb)
    #     axes = plt.gca()
    #     swiss_cheese_hole = plt.Circle((0., 0.), r_h, color='grey', fill=False, zorder=10)
    #     axes.add_artist(swiss_cheese_hole)
    #     axes.set_aspect('equal', adjustable='box')
    #     plt.plot(x_ltb, y_ltb, 'g-')
    #     lim = 200
    #     axes.set_xlim([-lim, lim])
    #     axes.set_ylim([-lim, lim])

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
    om_lambdas = np.linspace(0., 0.99, 50)[:-1]
    # z_lens_all = np.linspace(0.05, 0.2, 10)
    z_lens_all = np.linspace(0.2, 1, 20)
    # z_lens = 0.1
    # a_lens = 1/(z_lens+1)
    # start_thetas = np.linspace(0.7e-5, 2e-5, 100)
    start_thetas = np.array([1e-6]*50)
    source_rs = np.array([theta2rls_flat(th1, z1) for th1, z1 in zip(start_thetas, z_lens_all)])
    # source_rs = theta2rls_flat(start_thetas, z_lens)
    source_zs = rs2redshift_flat(source_rs)
    print("source_zs: ", source_zs)
    # step_size = 6.12244897959e-07
    # step_size = 9.63265306122e-07
    step_size = 1e-07
    # step_size = 5e-7
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
        filename = 'data/diff_lambdas_ltb_cheese.csv'
        if first:
            df.to_csv(filename, index=False)
            first = False
        else:
            df.to_csv(filename, index=False, header=False, mode='a')
    print("Time taken: {}".format(time.time() - start))

def main2():
    start = time.time()
    om_lambdas = np.linspace(0, 0.99, 1)
    # z_lens_all = np.linspace(0.05, 0.2, 10)
    z_lens_all = np.linspace(0.2, 1, 1)
    start_thetas = np.array([8e-7]*50)
    source_rs = np.array([theta2rls_flat(th1, z1) for th1, z1 in zip(start_thetas, z_lens_all)])
    # source_rs = theta2rls_flat(start_thetas, z_lens)
    source_zs = rs2redshift_flat(source_rs)
    print("source_zs: ", source_zs)
    step_size = 1e-7
    first = True
    for source_z, z_lens in list(zip(source_zs, z_lens_all)):
        rs = []
        thetas = []
        source_rs_array = []
        numerical_thetas = []
        dl = []
        dls = []
        ds = []
        for om in om_lambdas:
            comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om)
            print("DL", dang_lens)
            source_r, dang_r = get_distances(source_z, Omega_Lambda=om)
            theta = rs2theta(source_r, comoving_lens, dang_lens)
            # dls, dl, ds
            r, a = solve(theta, plot=True, comoving_lens=comoving_lens, Omega_Lambda=om, dt=step_size)
            thetas.append(theta)
            rs.append(r+comoving_lens)
            print("result", r)
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
        print("ds: ", ds)
        print("dls: ", dls)
    print("Time taken: {}".format(time.time() - start))

main2()