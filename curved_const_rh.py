'''
sw_lambda, but generalized for curved universe
dimensions are in k
origin at lens
keeps constant rh when varying Lambda
'''

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import time
import pandas as pd
from astropy.cosmology import LambdaCDM

length_scale = 3.086e22 # mega parsec

tols = {
    'atol': 1e-120,
    'rtol': 1e-14,
    'method': 'RK45',
}

H_0 = 70*1000/(3.086e22)/299792458 * length_scale # 70km/s/Mpc

M_initial = 1474e13 / length_scale
rho_initial = 3*H_0**2/(8*np.pi)
r_h = (3*M_initial/(4*np.pi*rho_initial))**(1./3)
M = None
print("r_h", r_h)
print("M: ", M)

def frw(eta, w, p):
    L, Omega_k, Omega_Lambda, Omega_m, H_0 = p

    a, r, rdot, t, phi = w

    k = omk2k(Omega_k, H_0)
    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda + Omega_k/a**2)

    phidot = L / (a*r)**2
    tdot = -np.sqrt(a**2*rdot**2/(1-k*r**2)+a**2*r**2*phidot**2)
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

# didn't check this function properly
# it is not used in generating any results
def get_angle_frw_curved(r, phi, rdot, phidot, k):
    if phi > np.pi/2:
        angle = np.pi - phi
    else:
        angle = phi
    # res = np.arctan((rdot*np.sin(phi)+r*np.cos(phi)*phidot)/(rdot*np.cos(phi)-r*np.sin(phi)*phidot))
    res = np.arctan(np.abs(r*phidot/rdot)*np.sqrt(1-k*r**2)) - angle
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

    return [
        r_h_t*tdot,
        tdot,
        rdot,
        rddot,
        phidot,
    ]


def solve(angle_to_horizontal, comoving_lens=None, verbose=True, Omega_Lambda=None, Omega_m=None):
    a0 = 1
    initial_a = a0
    # initial_r = 10e17
    initial_r = comoving_lens
    initial_phi = np.pi
    initial_t = 0.
    Omega_k = 1 - Omega_m - Omega_Lambda
    k = omk2k(Omega_k, H_0)


    ######################### frw1 initial conditions ##################################
    initial_rdot = -initial_r
    if verbose:
        print("om_k", Omega_k, k)
        print("1-kr^2", np.sqrt(1-k*initial_r**2))
    initial_phidot = np.tan(angle_to_horizontal)*initial_rdot/initial_r/np.sqrt(1-k*initial_r**2)
    # initial_phidot = np.tan(angle_to_horizontal) * initial_rdot / initial_r
    initial_tdot = -np.sqrt(initial_a**2*initial_rdot**2/(1-k*initial_r**2) + initial_a**2*initial_r**2*initial_phidot**2)

    if verbose:
        print("initial velocities:", initial_rdot, initial_phidot, initial_tdot)

    # rho = Omega_m*3*H_0**2/(8*np.pi)
    # r_h = 1/initial_a*(3*M/(4*np.pi*rho))**(1./3)
    if verbose:
        print('r_h:', r_h, "\n")
    chi_h = r2chi(k, r_h)

    Lambda = 3*Omega_Lambda*H_0**2
    L_frw = (initial_a*initial_r)**2*initial_phidot

    p_frw = [L_frw, Omega_k, Omega_Lambda, Omega_m, H_0]
    initial = [initial_a, initial_r, initial_rdot, initial_t, initial_phi]

    ######################### end frw1 initial conditions ##################################


    ######################### frw1 integration ##################################
    def reach_hole(t, y): return r2chi(k, y[1]) - chi_h
    reach_hole.terminal = True
    reach_hole.direction = -1

    sol = spi.solve_ivp(lambda t, y: frw(t, y, p_frw), [0, 10], initial, events=reach_hole, **tols)
    last = sol.y[:,-1]

    ######################### end frw1 integration ##################################

    frw_angle_before_entering = get_angle_frw_curved(last[1], last[4], last[2], L_frw/(last[0]*last[1])**2, k)

    if verbose:
        print("Angle in FRW::")
        print(frw_angle_before_entering)
        print(last)
        print("\n")


    ######################### conversion to kottler ##################################
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

    jacobian = np.matrix([[1/f*np.sqrt(1-k*frw_r**2), frw_a/f/np.sqrt(1-k*frw_r**2)*np.sqrt(1-f-k*frw_r**2)], [np.sqrt(1-f-k*frw_r**2), frw_a]])
    vels = np.matmul(jacobian, np.array([tdot_out, last[2]]))
    initial_rdot, initial_tdot = vels[0, 1], vels[0, 0]

    initial_phidot = L_frw / (last[0]*last[1])**2
    L_kottler = initial_phidot *initial_r**2

    initial_kottler = [initial_rh, initial_t, initial_r, initial_rdot, initial_phi]
    E = f*initial_tdot
    p_kottler = [E, L_kottler, M, Omega_Lambda, Omega_m, H_0, k, r_h]

    if verbose:
        print("kottler initial:")
        print("initial_rh, initial_t, initial_r, initial_rdot, initial_phi")
        print(initial_kottler)
        print("p_kottler")
        print(p_kottler)

    angle_before_entering = get_angle(initial_r, initial_phi, initial_rdot, initial_phidot)
    if verbose:
        print("light ray angle before entering kottler hole: ", angle_before_entering)
        print("initial conditions on entering kottler hole: ", initial_kottler)
        print("initial params on entering kottler hole: ", p_kottler)
    


    ######################### kottler integration ##################################
    
    # this is so that it doesn't hit the boundary condition initially
    initial_kottler[2] -= 1e-10

    def exit_hole(t, y): return y[2] - y[0]
    def turning_point(t, y): return y[4] - np.pi/2 
    exit_hole.terminal = True
    exit_hole.direction = 1
    turning_point.terminal = False
    turning_point.direction = -1
    sol_kottler = spi.solve_ivp(lambda t, y: kottler(t, y, p_kottler), [0, 10], initial_kottler, dense_output=True, events=[turning_point, exit_hole], **tols)
    last = sol_kottler.y[:,-1]
    if verbose:
        print("sol kottler", sol_kottler.t_events)
        print("turning point in kottler", sol_kottler.sol(sol_kottler.t_events[0])[2])
        print()


    if verbose:
        print("last in Kottler", last)
    angle_after_exiting = get_angle(last[2], last[4], last[3], L_kottler/last[2]**2)



    ######################### conversion back to FRW ##################################
    # r_h, t, r, rdot, phi = w
    initial_phi = last[4]
    initial_r = r_h
    initial_a = last[2] / r_h

    if verbose:
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

    # change the L_frw, I don't think it is needed, but doing it just in case
    p_frw[0] = initial_a**2 * initial_r**2*initial_phidot

    initial_t = 0
    # a, t, r, rdot, phi 

    # this angle is not used for results, just to get a feel of what it is when debugging
    frw_angle_after_exiting = get_angle_frw_curved(initial_r, initial_phi, initial_rdot, initial_phidot, k)

    if verbose:
        print("Angle after exiting FRW:", frw_angle_after_exiting)
        print("bending angle in FRW: ", frw_angle_before_entering + frw_angle_after_exiting)
        print("\n")

    # check if it is going to cross the axis, inform if otherwise
    # doesn't throw an error if it isn't
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
    if verbose:
        print("FRW2 initial: ", initial_frw2)

    enter_phi = initial_phi
    alpha = frw_angle_before_entering + frw_angle_after_exiting


    ######################### FRW integration ##################################
    def reach_axis(t, y): return y[4] - 0
    reach_axis.terminal = True
    reach_axis.direction = -1
    sol2 = spi.solve_ivp(lambda t, y: frw(t, y, p_frw), [0, 10], initial_frw2, events=reach_axis, **tols)
    last = sol2.y[:,-1]

    if verbose:
        print("last one frw", last)

    ######################### end FRW integration, return everything ##################################
    # alpha is not used, I also dunno why I'm returning it
    alpha  = get_angle_frw_curved(last[1], last[4], last[2], p_frw[0]/(last[0]*last[1])**2, k) + angle_to_horizontal
    return alpha, exit_rh, enter_phi, last[1], last[0]


def get_distances(z, Omega_Lambda=None, Omega_m=None):
    Omega_k = 1 - Omega_Lambda - Omega_m
    k = omk2k(Omega_k, H_0)
    def integrand(z):
        return 1/np.sqrt(Omega_m*(1+z)**3 + Omega_Lambda + Omega_k*(1+z)**2)
    integral, error = spi.quad(integrand, 0, z)
    chi = integral/H_0
    # comoving = chi
    comoving = chi2r(k, chi)
    dang = comoving/(1+z)
    return comoving, dang


def binary_search(start, end, answer, Omega_m, Omega_Lambda):
    mid = (end+start)/2
    res = get_distances(mid, Omega_Lambda, Omega_m)[1]
    # print(mid, res, answer)
    if np.isclose(res, answer, rtol=1e-14, atol=1e-14):
        return mid
    if res < answer:
        return binary_search(mid, end, answer, Omega_m, Omega_Lambda)
    else:
        return binary_search(start, mid, answer, Omega_m, Omega_Lambda)


# just a library to wrap around iterables to show a progress bar 
# it's a great library!
from tqdm import tqdm

def main(om_k = 0., to_file=True):
    start = time.time()

    # change this to edit the number of data points you want
    om_lambdas = np.linspace(0., 0.8, 50)

    # z of the lens
    z_lens_initial = 0.5
    # om_k = 0.1

    one_arcsec = 1/3600/180*np.pi
    # start_thetas = np.array([1/3600/180*np.pi]*50) # 1 arcsec
    filename = 'data/curved_const_rh_output.csv'
    first = True
    to_file = True
    print(filename)

    ## block A, see below
    # cosmo = LambdaCDM(H0=70, Om0=1, Ode0=0)
    # dang_lens = cosmo.angular_diameter_distance(z_lens_initial).value
    # comoving_lens = cosmo.comoving_transverse_distance(z_lens_initial).value
    # theta = one_arcsec

    z_lens = z_lens_initial
    theta = one_arcsec

    thetas = []
    dl = []
    ms = []
    com_lens = []
    exit_rhs = []
    enter_phis = []
    alphas = []
    om_ks = []
    raw_rs = []
    for om in tqdm(om_lambdas):
        om_m = 1 - om - om_k

        global M
        rho = (1-om)*3*H_0**2/(8*np.pi)
        M = 4/3*np.pi*rho*r_h**3
        ms.append(M)
        k = omk2k(om_k, H_0)

        # comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om, Omega_m=om_m)
        
        # the block below achieves the same effect as the above line, 
        # but using astropy function. I checked they are the same,
        # but you can use the below code if you have more faith in a tested library function (:
        cosmo = LambdaCDM(H0=70, Om0=om_m, Ode0=om)
        dang_lens = cosmo.angular_diameter_distance(z_lens).value
        comoving_lens = cosmo.comoving_transverse_distance(z_lens).value

        ######################################################

        ## you can also fix lens angular diameter distance (dang_lens)
        ## and calculate z_lens for different Lambdas
        ## it finds the z_lens through a binary search
        ## you have to uncomment the stuff below and block A from above
        # z_lens = binary_search(0, z_lens_initial, dang_lens, om_m, om)
        # comoving_lens = dang_lens * (1+z_lens)

        dl.append(dang_lens)
        com_lens.append(comoving_lens)
        thetas.append(theta)
        alpha, exit_rh, enter_phi, raw_r, source_a = solve(theta, verbose=False, comoving_lens=comoving_lens, Omega_Lambda=om, Omega_m=om_m)
        exit_rhs.append(exit_rh)
        enter_phis.append(enter_phi)
        alphas.append(alpha)
        om_ks.append(om_k)
        raw_rs.append(raw_r)

    thetas = np.array(thetas)
    dl = np.array(dl)
    ms = np.array(ms)
    com_lens = np.array(com_lens)
    exit_rhs = np.array(exit_rhs)
    enter_phis = np.array(enter_phis)
    alphas = np.array(alphas)
    om_ks = np.array(om_ks)
    raw_rs = np.array(raw_rs)

    df = pd.DataFrame({
        'DL': dl,
        'M': ms, 
        'om_lambdas': om_lambdas, 
        'theta': thetas, 
        'comoving_lens': com_lens,
        'exit_rhs': exit_rhs,    # for calculating Rindler & Ishak predictions
        'enter_phis': enter_phis,  # for calculating kantowski's predictions
        # 'alphas': alphas,
        'om_ks': om_ks,
        'raw_rs': raw_rs,
    })
    if to_file:
        # this conditional is not really needed,
        # was left over from a time when I had another outer loop
        if first:
            df.to_csv(filename, index=False)
            first = False
        else:
            df.to_csv(filename, index=False, header=False, mode='a')

    print("Time taken: {}".format(time.time() - start))
    return df

def main_multiple_omk():
    current = None
    filename = "curvedpyoutput_fixedr_withk.csv"
    # omks = np.linspace(0., 0.001, 10)
    omks = [0, 0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.05]
    for om_k in omks:
        df = main(om_k=om_k, to_file=False)
        df['om_k'] = om_k
        if current is None:
            df.to_csv(filename, index=False)
        else:
            df.to_csv(filename, index=False, header=False, mode='a')
        current = df


# main()
main_multiple_omk()