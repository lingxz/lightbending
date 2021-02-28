'''
curved universe
dimensions are in k
origin at lens
'''

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import time
import pandas as pd
from astropy.cosmology import LambdaCDM

length_scale = 3.086e22 # 1 mega parsec

tols = {
    'atol': 1e-120,
    'rtol': 1e-14,
    'method': 'RK45',
}

H_0 = 70*1000/(3.086e22)/299792458 * length_scale # 70km/s/Mpc
M = 1474e13 / length_scale

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
    R_h, T, R, Rdot, phi = w

    Lambda = 3*Omega_Lambda*H_0**2
    f = 1 - 2*M/R - Lambda/3*R**2
    Rddot = L**2 * (R - 3*M) / R**4
    Tdot = E / f
    phidot = L/R**2

    fh = 1 - 2*M/R_h - Lambda/3*R_h**2
    R_h_t = fh*np.sqrt(1-fh/(1-k*rh**2))

    return [
        R_h_t*Tdot,
        Tdot,
        Rdot,
        Rddot,
        phidot,
    ]

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

def kottler_with_frw(eta, w, p):
    E, L, M, Omega_Lambda, Omega_m, H_0, k, rh, L_frw, Omega_k = p
    R_h, T, R, Rdot, phi, a, r, rdot, t, phi_frw = w

    Lambda = 3*Omega_Lambda*H_0**2
    f = 1 - 2*M/R - Lambda/3*R**2
    Rddot = L**2 * (R - 3*M) / R**4
    Tdot = E / f
    phidot = L/R**2

    fh = 1 - 2*M/R_h - Lambda/3*R_h**2
    R_h_t = fh*np.sqrt(1-fh/(1-k*rh**2))

    # frw
    a_t = a * H_0 * np.sqrt(Omega_m/a**3 + Omega_Lambda + Omega_k/a**2)
    phidot_frw = L_frw / (a*r)**2
    tdot = -np.sqrt(a**2*rdot**2/(1-k*r**2)+a**2*r**2*phidot_frw**2)
    rddot = (1-k*r**2)*r*phidot_frw**2 - k*r*rdot**2/(1-k*r**2) - 2*a_t/a*rdot*tdot

    return [
        R_h_t*Tdot,
        Tdot,
        Rdot,
        Rddot,
        phidot,
        
        # frw
        a_t*tdot,
        rdot,
        rddot,
        tdot,
        phidot_frw,
    ]

class SwissCheese:
    # class to encapsulate the solution for light propagation in a swiss cheese model.
    def __init__(self, M, Omega_Lambda, Omega_m, comoving_lens, angle_to_horizontal):
        self.Omega_Lambda = Omega_Lambda
        self.Omega_m = Omega_m
        self.Omega_k = 1 - Omega_Lambda - Omega_m
        self.Lambda = 3*Omega_Lambda*H_0**2
        self.comoving_lens = comoving_lens # r coordinate of the lens
        self.angle_to_horizontal = angle_to_horizontal
        self.k = omk2k(self.Omega_k, H_0)
        
        initial_a = 1
        rho_0 = Omega_m*3*H_0**2/(8*np.pi) # critical density
        self.M = M # mass of the hole
        self.r_h = 1/initial_a*(3*M/(4*np.pi*rho_0))**(1./ 3) # r coordinate of the hole, based on the mass
        self.chi_h = r2chi(self.k, self.r_h)

        # params that will be set during the run
        self.frw_parameters = None
        self.frw_initial = None
        self.frw_final = None
        self.kottler_parameters = None
        self.kottler_initial = None
        self.kottler_final = None
        self.frw_initial_right = None  # for the FRW integration to the right of the hole.
        self.frw_final_right = None
    
    def run(self):
        self.frw_initial_conditions_and_parameters()
        self.integrate_frw()
        self.convert_frw_to_kottler_coordinates()
        self.integrate_kottler()
        self.convert_kottler_to_frw_coordinates()
        self.integrate_frw_right()

    
    def frw_initial_conditions_and_parameters(self):
        ## Initial conditions to start the integration
        initial_r = self.comoving_lens
        initial_rdot = -initial_r
        initial_phi = np.pi
        initial_t = 0.
        initial_a = 1
        initial_phidot = np.tan(self.angle_to_horizontal)*initial_rdot/initial_r/np.sqrt(1-self.k*initial_r**2)
        initial_tdot = -np.sqrt(initial_a**2*initial_rdot**2/(1-self.k*initial_r**2) + initial_a**2*initial_r**2*initial_phidot**2)
        initial = [initial_a, initial_r, initial_rdot, initial_t, initial_phi]

        ## Parameters needed for FRW integration
        L = (initial_a*initial_r)**2*initial_phidot
        parameters = [L, self.Omega_k, self.Omega_Lambda, self.Omega_m, H_0]
        self.frw_initial = initial
        self.frw_parameters = parameters
        return initial, parameters
    
    def integrate_frw(self):
        initial = self.frw_initial
        parameters = self.frw_parameters
        def reach_hole(t, y): return r2chi(self.k, y[1]) - self.chi_h
        reach_hole.terminal = True
        reach_hole.direction = -1
        sol = spi.solve_ivp(lambda t, y: frw(t, y, parameters), [0, 10], initial, events=reach_hole, **tols)
        self.frw_final = sol.y_events[0][0]

    def convert_frw_to_kottler_coordinates(self):
        a, r, rdot, t, phi = self.frw_final
        L_frw = self.frw_parameters[0]
        R_out = a * r
        exit_Rh = R_out
        Tdot_out = -np.sqrt(a**2*rdot**2+a**2*r**2*L_frw/(a*r)**2)
        initial_T = 0
        initial_R = R_out
        initial_phi = phi
        initial_Rh = initial_R
        f = 1 - 2*self.M/R_out - self.Lambda/3*R_out**2
        etadot = Tdot_out / a # conformal time

        jacobian = np.matrix([[1/f*np.sqrt(1-self.k*r**2), a/f/np.sqrt(1-self.k*r**2)*np.sqrt(1-f-self.k*r**2)], [np.sqrt(1-f-self.k*r**2), a]])
        vels = np.matmul(jacobian, np.array([Tdot_out, rdot]))
        initial_Rdot, initial_Tdot = vels[0, 1], vels[0, 0]

        initial_phidot = L_frw / (a*r)**2
        L_kottler = initial_phidot *initial_R**2

        self.kottler_initial = [initial_Rh, initial_T, initial_R, initial_Rdot, initial_phi]
        E = f*initial_Tdot
        self.kottler_parameters = [E, L_kottler, self.M, self.Omega_Lambda, self.Omega_m, H_0, self.k, self.r_h]

    def integrate_kottler(self):
        # this is so that it doesn't hit the boundary condition initially
        initial = self.kottler_initial[:]
        parameters = self.kottler_parameters
        initial[2] -= 1e-10

        def exit_hole(t, y): return y[2] - y[0]
        def turning_point(t, y): return y[4] - np.pi/2 
        exit_hole.terminal = True
        exit_hole.direction = 1  # when value goes from negative to positive
        turning_point.terminal = False
        turning_point.direction = -1
        sol_kottler = spi.solve_ivp(lambda t, y: kottler(t, y, parameters), [0, 10], initial, dense_output=True, events=[turning_point, exit_hole], **tols)
        self.kottler_final = sol_kottler.y_events[1][0]
    
    def convert_kottler_to_frw_coordinates(self):
        # r_h, t, r, rdot, phi = w
        Rh, T, R, Rdot, phi = self.kottler_final
        L_kottler = self.kottler_parameters[1]
        E_kottler = self.kottler_parameters[0]
        initial_phi = phi
        initial_r = self.r_h
        initial_a = R / self.r_h

        initial_phidot = L_kottler / R**2
        f = 1-2*M/R-self.Lambda/3*R**2
        Tdot = E_kottler/f

        # a, r, rdot, t, phi
        frw_r = self.r_h
        frw_a = initial_a
        jacobian = np.matrix([[1/f*np.sqrt(1-self.k*frw_r**2), frw_a/f/np.sqrt(1-self.k*frw_r**2)*np.sqrt(1-f-self.k*frw_r**2)], [np.sqrt(1-self.k*frw_r**2-f), frw_a]])
        inv_jac = np.linalg.inv(jacobian)
        vels = np.matmul(inv_jac, np.array([Tdot, Rdot]))
        initial_rdot, initial_tdot = vels[0, 1], vels[0, 0]

        # change the L_frw, I don't think it is needed, but doing it just in case
        self.frw_parameters[0] = initial_a**2 * initial_r**2*initial_phidot

        initial_t = 0
        # a, t, r, rdot, phi 

        # check if it is going to cross the axis, inform if otherwise
        # doesn't throw an error if it isn't
        initial_ydot = initial_rdot*np.sin(initial_phi) + initial_r*np.cos(initial_phi)*initial_phidot
        if initial_ydot > 0:
            print("light ray is not going to cross the axis, decrease angle_to_horizontal")
            # print("initial angle to horizontal: ", angle_to_horizontal)
            print("----")

        if initial_r*np.sin(initial_phi) < 0:
            print("light ray is bent too much by the hole, increase angle_to_horizontal")
            # print("initial angle to horizontal: ", angle_to_horizontal)
            print("----")

        self.frw_initial_right = [initial_a, initial_r, initial_rdot, initial_t, initial_phi]

        # enter_phi = initial_phi
        # alpha = frw_angle_before_entering + frw_angle_after_exiting
    
    def integrate_frw_right(self):
        ######################### FRW integration ##################################
        def reach_axis(t, y): return y[4] - 0
        reach_axis.terminal = True
        reach_axis.direction = -1
        sol = spi.solve_ivp(lambda t, y: frw(t, y, self.frw_parameters), [0, 10], self.frw_initial_right, events=reach_axis, **tols)
        self.frw_final_right = sol.y_events[0][0]
    
    def debug_logs(self):
        print("===== Initial parameters =====")
        self._print_with_names([self.k, self.Omega_k, self.r_h, self.angle_to_horizontal, self.comoving_lens], ["k", "Omega_k", "r_h", "angle_to_horizontal", "comoving_lens"])
        print("===== Initial FRW coordinates =====")
        self._print_with_names(self.frw_initial, ["a", "r", "rdot", "t", "phi"])
        self._print_with_names(self.frw_parameters, ["L", "Omega_k", "Omega_Lambda", "Omega_m", "H_0"])

        print("===== Final FRW coordinates =====")
        self._print_with_names(self.frw_final, ["a", "r", "rdot", "t", "phi"])
        # TODO clean this up
        frw_angle_before_entering_hole = get_angle_frw_curved(self.frw_final[1], self.frw_final[4], self.frw_final[2], self.frw_parameters[0]/(self.frw_final[0]*self.frw_final[1])**2, self.k)
        print("FRW angle before entering hole:", frw_angle_before_entering_hole)
        print()

        print("===== Initial Kottler coordinates =====")
        self._print_with_names(self.kottler_initial, ["R_h", "T", "R", "Rdot", "phi"])
        self._print_with_names(self.kottler_parameters, ["E", "L", "M", "Omega_Lambda", "Omega_m", "H_0", "k", "r_h"])
        
        print("===== Final Kottler coordinates =====")
        self._print_with_names(self.kottler_final, ["R_h", "T", "R", "Rdot", "phi"])
        
        print("===== Initial FRW coordinates (right) =====")
        self._print_with_names(self.frw_initial_right, ["a", "r", "rdot", "t", "phi"])
        print()
        
        print("===== Final FRW coordinates (right) =====")
        self._print_with_names(self.frw_final_right, ["a", "r", "rdot", "t", "phi"])
    
    def get(self, attribute, param):
        arr = getattr(self, attribute)
        frw_coords_to_index = ["a", "r", "rdot", "t", "phi"]
        frw_params_to_index = ["L", "Omega_k", "Omega_Lambda", "Omega_m", "H_0"]
        kottler_coords_to_index = ["R_h", "T", "R", "Rdot", "phi"]
        kottler_params_to_index = ["E", "L", "M", "Omega_Lambda", "Omega_m", "H_0", "k", "r_h"]
        
        if attribute in ["frw_initial", "frw_initial_right", "frw_final", "frw_final_right"]:
            index = frw_coords_to_index.index(param)
        elif attribute in ["kottler_initial", "kottler_final"]:
            index = kottler_coords_to_index.index(param)
        elif attribute in ["frw_parameters"]:
            index = frw_params_to_index.index(param)
        elif attribute in ["kottler_parameters"]:
            index = kottler_params_to_index.index(param)
        return arr[index]


    @staticmethod    
    def _print_with_names(values, names):
        for name, value in zip(names, values):
            print("{}: {}".format(name, value))


class NoHoleFRWIntegration(SwissCheese):
    def integrate_frw(self):
        initial = self.frw_initial
        parameters = self.frw_parameters
        def reach_hole(t, y): return y[4] - np.pi/2
        reach_hole.terminal = True
        reach_hole.direction = -1
        sol = spi.solve_ivp(lambda t, y: frw(t, y, parameters), [0, 10], initial, events=reach_hole, **tols)
        self.frw_final = sol.y_events[0][0]


class SwissCheeseKottlerModification(SwissCheese):
    def integrate_kottler(self):
        initial = self.kottler_initial[:]
        initial = list(initial) + list(self.frw_final)
        parameters = self.kottler_parameters
        parameters = list(parameters) + [self.get("frw_parameters", "L"), self.Omega_k]
        def exit_hole(t, y): return y[2] - y[0]
        exit_hole.terminal = False
        exit_hole.direction = 0 # when value goes from negative to positive
        
        # integrate past the hole before terminating so that the other boundary change of exiting hole is captured
        def way_past_hole(t, y): return y[2] - y[0] - initial[2]/ 2
        way_past_hole.terminal = True
        way_past_hole.direction = 1

        sol_kottler = spi.solve_ivp(lambda t, y: kottler_with_frw(t, y, parameters), [0, 10], initial, dense_output=True, events=[exit_hole, way_past_hole], **tols)
        self.integrated_a = sol_kottler.y_events[0][-1][5]
        self.kottler_final = list(sol_kottler.y_events[0][-1])[:5]


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


from tqdm import tqdm

def main(om_k = 0., filename=None, model=None):
    if model is None:
        model = SwissCheese
    start = time.time()

    # change this to edit the number of data points you want
    # om_lambdas = np.linspace(0., 0.99, 50)
    om_lambdas = np.linspace(0., 0.49, 50)
    # om_lambdas = [0]

    # z of the lens
    z_lens_initial = 0.5

    one_arcsec = 1/3600/180*np.pi
    # start_thetas = np.array([1/3600/180*np.pi]*50) # 1 arcsec
    # filename = 'data/curvedpyoutput_k{}.csv'.format(om_k)
    print(filename)

    ## block A, see below
    # cosmo = LambdaCDM(H0=70, Om0=1, Ode0=0)
    # dang_lens = cosmo.angular_diameter_distance(z_lens_initial).value
    # comoving_lens = cosmo.comoving_transverse_distance(z_lens_initial).value
    # theta = one_arcsec

    z_lens = z_lens_initial
    theta = one_arcsec

    thetas = [] # fixed
    dl = [] # from the z_lens, so fixed
    ms = []
    com_lens = []
    exit_rhs = []
    enter_phis = []
    om_ks = []
    raw_rs = []
    for om in tqdm(om_lambdas):
        om_m = 1 - om - om_k
        assert om_m > 0
        assert om_m <= 1
        k = omk2k(om_k, H_0)
        ms.append(M)


        # comoving_lens, dang_lens = get_distances(z_lens, Omega_Lambda=om, Omega_m=om_m)

        # astropy function to do the same thing
        cosmo = LambdaCDM(H0=70, Om0=om_m, Ode0=om)
        dang_lens = cosmo.angular_diameter_distance(z_lens).value
        comoving_lens = cosmo.comoving_transverse_distance(z_lens).value


        # numerical
        solution = model(
            M = M,
            Omega_Lambda = om,
            Omega_m = om_m,
            comoving_lens = comoving_lens,
            angle_to_horizontal = theta)
        solution.run()

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

        exit_rh, enter_phi, raw_r, source_a = solution.kottler_initial[2], solution.frw_initial_right[4], solution.frw_final_right[1], solution.frw_final_right[0]

        # exit_rh, enter_phi, raw_r, source_a = solve(theta, verbose=False, comoving_lens=comoving_lens, Omega_Lambda=om, Omega_m=om_m)
        exit_rhs.append(exit_rh)
        enter_phis.append(enter_phi)
        om_ks.append(om_k)
        raw_rs.append(raw_r)

    thetas = np.array(thetas)
    dl = np.array(dl)
    ms = np.array(ms)
    com_lens = np.array(com_lens)
    exit_rhs = np.array(exit_rhs)
    enter_phis = np.array(enter_phis)
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
        'om_ks': om_ks,
        'raw_rs': raw_rs,
    })
    first = True
    if filename is not None:
        if first:
            df.to_csv(filename, index=False)
            first = False
        else:
            df.to_csv(filename, index=False, header=False, mode='a')

    print("Time taken: {}".format(time.time() - start))
    return df


def main_multiple_omk(filename, model=None):
    current = None
    # omks = np.linspace(0., 0.001, 10)
    # omks = [0, 0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.05]
    omks = [0.1, 0.3, 0.5, 0.8, 0.9]
    omks = np.linspace(0., 0.5, 10)
    for om_k in omks:
        df = main(om_k=om_k, filename=None, model=model)
        df['om_k'] = om_k
        if current is None:
            df.to_csv(filename, index=False)
        else:
            df.to_csv(filename, index=False, header=False, mode='a')
        current = df


def compare_with_analytical():
    # analytical
    z_lens = 0.5
    one_arcsec = 1/3600/180*np.pi
    theta = one_arcsec
    om_m = 0.1
    om_lambda = 0.5
    cosmo = LambdaCDM(H0=70, Om0=om_m, Ode0=om_lambda)
    dang_lens = cosmo.angular_diameter_distance(z_lens).value
    comoving_lens = cosmo.comoving_transverse_distance(z_lens).value

    # numerical
    solution = NoHoleFRWIntegration(
        M = M,
        Omega_Lambda = om_lambda,
        Omega_m = om_m,
        comoving_lens = comoving_lens,
        angle_to_horizontal = theta)
    solution.frw_initial_conditions_and_parameters()
    solution.integrate_frw()

    # compare
    # print("numerical", solution.frw_final[0] *solution.frw_final[1] / theta)  # a * r / theta
    print("numerical", solution.get("frw_final", "a") * solution.get("frw_final", "r") / theta)
    print("dang_lens", dang_lens)


def compare_kottler_hole_size_with_frw_size():
    # analytical
    z_lens = 0.5
    one_arcsec = 1/3600/180*np.pi
    theta = one_arcsec
    om_m = 0.1
    om_lambda = 0.1
    cosmo = LambdaCDM(H0=70, Om0=om_m, Ode0=om_lambda)
    dang_lens = cosmo.angular_diameter_distance(z_lens).value
    comoving_lens = cosmo.comoving_transverse_distance(z_lens).value

    # numerical
    solution = SwissCheeseKottlerModification(
        M = M,
        Omega_Lambda = om_lambda,
        Omega_m = om_m,
        comoving_lens = comoving_lens,
        angle_to_horizontal = theta)
    solution.run()

    print(solution.integrated_a * solution.r_h, solution.get("kottler_final", "R_h"), solution.integrated_a, solution.get("frw_initial_right", "a"), (1 - solution.integrated_a/solution.get("frw_initial_right", "a")) * 100)
    # print(solution.integrated_a, solution.get("kottler_final", "R_h"), solution.get("kottler_parameters", "r_h"), solution.r_h)
    

if __name__ == '__main__':
    # main(filename="data/curvedpyoutput_new.csv")
    # compare_with_analytical()
    compare_kottler_hole_size_with_frw_size()
    # main_multiple_omk(filename = "data/curvedpyoutput_withk4.csv", model=SwissCheeseKottlerModification)
    # main_multiple_omk(filename = "data/curvedpyoutput_withk5.csv", model=SwissCheese) # normal
    