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
from utils import r2chi, chi2r, binary_search
from astropy.cosmology import LambdaCDM

length_scale = 3.086e22 # 1 mega parsec

tols = {
    'atol': 1e-120,
    'rtol': 1e-14,
    'method': 'RK45',
}

# mass of sun in kg * G / c**2 / length_scale (for convenience)
M_sun = 1.98847e30 * (6.67408e-11) / 299792458**2 / length_scale
one_arcsec = 1/3600/180*np.pi
H_0 = 70*1000/(3.086e22)/299792458 * length_scale # 70km/s/Mpc
M = M_sun * 1e13

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
    # fh = 1 - 2*M/R_h
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
        try:
            self.integrate_frw_right()
        except Exception as e:
            print("Did not run last FRW integration")

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
        # def reach_hole(t, y): return y[1] - self.r_h
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
        self.kottler_final = sol_kottler.y_events[1][-1]

        # ###### distance travelled calculation
        # rs = sol_kottler.y[2]
        # phis = sol_kottler.y[4]
        # dr = np.diff(rs)
        # dphi = np.diff(phis)
        # rs = rs[:-1]
        # distance_travelled = np.sum(np.sqrt(dr**2 + (rs*dphi)**2))
        # print("Distance travelled in Kottler hole:", distance_travelled)

        # extra
        self.kottler_turning_point = sol_kottler.y_events[0][-1]
        self.kottler_closest_approach = self.kottler_turning_point[2]

    def convert_kottler_to_frw_coordinates(self):
        # r_h, t, r, rdot, phi = w
        Rh, T, R, Rdot, phi = self.kottler_final
        L_kottler = self.kottler_parameters[1]
        E_kottler = self.kottler_parameters[0]
        initial_phi = phi
        initial_r = self.r_h
        initial_a = R / self.r_h

        initial_phidot = L_kottler / R**2
        f = 1-2*self.M/R-self.Lambda/3*R**2
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

        self.frw_initial_right = [initial_a, initial_r, initial_rdot, initial_t, initial_phi]

        # check if it is going to cross the axis, inform if otherwise
        # doesn't throw an error if it isn't
        initial_ydot = initial_rdot*np.sin(initial_phi) + initial_r*np.cos(initial_phi)*initial_phidot
        initial_xdot = initial_rdot*np.cos(initial_phi) - initial_r*np.sin(initial_phi)*initial_phidot
        self.frw_right_initial_ydot = initial_ydot
        self.frw_right_initial_xdot = initial_xdot
        if initial_ydot > 0:
            print("light ray is not going to cross the axis, decrease angle_to_horizontal")
            # print("initial angle to horizontal: ", angle_to_horizontal)
            print("----")

        if initial_r*np.sin(initial_phi) < 0:
            print("light ray is bent too much by the hole, increase angle_to_horizontal")
            # print("initial angle to horizontal: ", angle_to_horizontal)
            print("----")


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


class SwissCheeseNoExpansion(SwissCheese):
    def integrate_kottler(self):
        # this is so that it doesn't hit the boundary condition initially
        initial = self.kottler_initial[:]
        parameters = self.kottler_parameters
        initial[2] -= 1e-10

        def exit_hole(t, y): return y[2] - self.kottler_initial[0]
        def turning_point(t, y): return y[4] - np.pi/2 
        exit_hole.terminal = True
        exit_hole.direction = 1  # when value goes from negative to positive
        turning_point.terminal = False
        turning_point.direction = -1
        sol_kottler = spi.solve_ivp(lambda t, y: kottler(t, y, parameters), [0, 10], initial, dense_output=True, events=[turning_point, exit_hole], **tols)
        self.kottler_final = sol_kottler.y_events[1][-1]
        
        # extra
        self.kottler_turning_point = sol_kottler.y_events[0][-1]
        self.kottler_closest_approach = self.kottler_turning_point[2]

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


# def binary_search(start, end, answer, Omega_m, Omega_Lambda):
#     mid = (end+start)/2
#     res = get_distances(mid, Omega_Lambda, Omega_m)[1]
#     # print(mid, res, answer)
#     if np.isclose(res, answer, rtol=1e-14, atol=1e-14):
#         return mid
#     if res < answer:
#         return binary_search(mid, end, answer, Omega_m, Omega_Lambda)
#     else:
#         return binary_search(start, mid, answer, Omega_m, Omega_Lambda)


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

        rs = chi2r(k, r2chi(k, raw_r) + r2chi(k, solution.comoving_lens))
        R = theta * dang_lens
        alpha_numerical = theta * rs / raw_r
        alpha_schwarzschild = 4*solution.M/R + 15*np.pi*solution.M**2/4/R**2 + 401/12*solution.M**3/R**3
        print("====")
        print(1 - alpha_numerical / alpha_schwarzschild)
        print(theta, comoving_lens, raw_r)

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

    print("comoving!!!", comoving_lens)
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


# schucker case is flat, integrate two lines, see where they meet
# repro of table 4 of Schucker 2009
def schucker(alpha, alpha_prime, z_lens, om_lambda=0.77, mass=M_sun*1e13):
    om_m = 1 - om_lambda
    cosmo = LambdaCDM(H0=70, Om0=om_m, Ode0=om_lambda)

    dang_lens = cosmo.angular_diameter_distance(z_lens).value
    comoving_lens = cosmo.comoving_transverse_distance(z_lens).value
    print(dang_lens)

    # keep dang_lens the same
    dang_lens = 1505.44763
    z_lens = binary_search(cosmo, 0, 2., dang_lens)
    comoving_lens = cosmo.comoving_transverse_distance(z_lens).value

    # lower ray
    solution1 = SwissCheese(
        M = mass,
        Omega_Lambda = om_lambda,
        Omega_m = om_m,
        comoving_lens=comoving_lens,
        angle_to_horizontal = alpha
    )

    # upper ray
    solution2 = SwissCheese(
        M = mass,
        Omega_Lambda = om_lambda,
        Omega_m = om_m,
        comoving_lens=comoving_lens,
        angle_to_horizontal = alpha_prime
    )
    solution1.run()
    solution2.run()


    # find intersection between the two lines
    from sympy import Point, Line
    # upper ray
    p1 = Point(solution2.get("frw_initial_right", "r") * np.cos(solution2.get("frw_initial_right", "phi")), solution2.get("frw_initial_right", "r") * np.sin(solution2.get("frw_initial_right", "phi")))
    p2 = Point(solution2.get("frw_final_right", "r") * np.cos(solution2.get("frw_final_right", "phi")), solution2.get("frw_final_right", "r") * np.sin(solution2.get("frw_final_right", "phi")))
    upper_line = Line(p1, p2)
    # lower ray
    x1 = solution1.get("frw_initial_right", "r") * np.cos(solution1.get("frw_initial_right", "phi"))
    y1 = solution1.get("frw_initial_right", "r") * np.sin(solution1.get("frw_initial_right", "phi"))
    m = solution1.frw_right_initial_ydot / solution2.frw_right_initial_xdot
    x2 = solution2.get("frw_final_right", "r")
    y2 = m * (x2 - x1) + y1
    lower_line = Line(Point(x1, -y1), Point(x2, -y2))

    intersect = tuple(upper_line.intersection(lower_line)[0].evalf())
    phi = np.arctan(float(intersect[1]) / float(intersect[0]))
    phi_in_arcsec = abs(phi) * 180 / np.pi * 3600
    bending = 4 * mass / (alpha_prime * dang_lens)

    print("============", phi_in_arcsec, solution1.r_h, intersect, solution2.get("frw_final", "phi"), comoving_lens, dang_lens)

def swiss_chess_single_run(om_lambda, mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec):
    # analytical
    om_m = 1 - om_lambda
    cosmo = LambdaCDM(H0=70, Om0=om_m, Ode0=om_lambda)
    dang_lens = cosmo.angular_diameter_distance(z_lens).value
    comoving_lens = cosmo.comoving_transverse_distance(z_lens).value

    # numerical
    solution = SwissCheese(
        M = mass,
        Omega_Lambda = om_lambda,
        Omega_m = om_m,
        comoving_lens = comoving_lens,
        angle_to_horizontal = theta)
    solution.run()

    print("kantowski rb:", solution.get("kottler_initial", "R"))

    # check distance of closest approach
    R = theta * dang_lens
    r0 = 1/(1/R + solution.M/R**2 + 3/16*solution.M**2/R**3)
    # print(solution.kottler_closest_approach, r0)
    
    phi_tilda = solution.get("kottler_final", "phi")
    print("phi_tilda in degrees", phi_tilda *  180/np.pi)
    kantowski, kantowski_error = kantowski_alpha(r0, solution.M, phi_tilda, solution.Lambda)
    kantowski = np.abs(kantowski)
    # alpha_numerical = solution.angle_to_horizontal * (solution.get("frw_final_right", "r") + solution.comoving_lens) / solution.get("frw_final_right", "r")
    alpha_schwarzschild = 4*solution.M/R + 15*np.pi*solution.M**2/4/R**2 + 401/12*solution.M**3/R**3
    # print(alpha_numerical, kantowski, alpha_schwarzschild)
    return {
        # "fractional_deviation": np.abs(1 - kantowski / alpha_numerical),
        # "higher_order_kantowski": kantowski_error / alpha_numerical,
        # "numerical_fractional_deviation_from_schwarzschild": np.abs(1 - alpha_numerical / alpha_schwarzschild),
        "kantowski_fractional_deviation_from_schwarzschild": np.abs(1 - kantowski / alpha_schwarzschild)
    }

def swiss_chess_curved_single_run(om_k, om_lambda, mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec, const_dl=True):
    k = omk2k(om_k, H_0)
    om_m = 1 - om_lambda - om_k

    cosmo = LambdaCDM(H0=70, Om0=om_m, Ode0=om_lambda)
    dang_lens = cosmo.angular_diameter_distance(z_lens).value
    # this is NOT chi, this is r, which is the comoving coordinate distance, or comoving transverse distance.
    comoving_coordinate_lens = cosmo.comoving_transverse_distance(z_lens).value
    # this is chi for clarity (but it's not needed)
    chi = cosmo.comoving_distance(z_lens).value

    #### testing out keeping dang_lens constant instead of z_lens
    if const_dl:
        dang_lens = 1130
        z_lens = binary_search(cosmo, 0, 2., dang_lens)
        comoving_coordinate_lens = cosmo.comoving_transverse_distance(z_lens).value

    # print("comoving_lens", comoving_coordinate_lens)
    # numerical
    solution = SwissCheese(
        M = mass,
        Omega_Lambda = om_lambda,
        Omega_m = om_m,
        comoving_lens = comoving_coordinate_lens,
        angle_to_horizontal = theta)
    solution.run()

    # check distance of closest approach
    R = theta * dang_lens    
    r0 = 1/(1/R + solution.M/R**2 + 3/16*solution.M**2/R**3)
    phi_tilda = solution.get("kottler_final", "phi")

    rs = chi2r(k, r2chi(k, solution.get("frw_final_right", "r")) + r2chi(k, solution.comoving_lens))
    alpha_numerical = solution.angle_to_horizontal * rs / solution.get("frw_final_right", "r")
    alpha_schwarzschild = 4*solution.M/R + 15*np.pi*solution.M**2/4/R**2 + 401/12*solution.M**3/R**3
    # Ishak&Rindler 2010 The relevance of the Cosmological constant for lensing eq 31
    alpha_rindler_2010 = 4 * solution.M / R + 15 * np.pi/4 * solution.M**2/R**2 + 305 /12*solution.M**3/R**3 - solution.Lambda*R/3 * solution.r_h
    # Rindler&Ishak 2007 The contribution of the cosmological constant to the relativistic bending of light revisited eq 17
    alpha_rindler_2007 = 4 * solution.M / R + 15 * np.pi/4 * solution.M**2/R**2 - solution.Lambda * R**3 / solution.M / 6
    kantowski, kantowski_error = kantowski_alpha(r0, solution.M, phi_tilda, solution.Lambda)
    kantowski = np.abs(kantowski)

    # adot = H_0 * np.sqrt(om_m/0.9**3 + om_lambda + om_k / 0.9**2)
    # print("fake adot", om_lambda, adot)

    return {
        "numerical_fractional_deviation_from_kantowski": alpha_numerical / kantowski  - 1,
        "higher_order_kantowski": kantowski_error / alpha_numerical,
        "numerical_fractional_deviation_from_schwarzschild": alpha_numerical / alpha_schwarzschild - 1,
        "kantowski_fractional_deviation_from_schwarzschild": kantowski / alpha_schwarzschild - 1,
        "rindler_2010_fractional_deviation_from_schwarzschild": alpha_rindler_2010 / alpha_schwarzschild - 1,
        "rindler_2007_fractional_deviation_from_schwarzschild": alpha_rindler_2007 / alpha_schwarzschild - 1,
        "phi_tilda": phi_tilda,
        "Rh_enter": solution.get("kottler_final", "R_h"),
        "Rh_exit": solution.get("kottler_initial", "R_h")
    }

# this produces figure 4
# Plots fractional deviation from schwarzschild for various om_k, om_lambda
def swiss_chess_curved_multi_run(mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec):
    om_ks = [-0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4]
    for om_k in om_ks:
        results = []
        om_lambdas = np.linspace(0., 1 - om_k - 0.1, 30)
        for om_lambda in om_lambdas:
            result = swiss_chess_curved_single_run(om_k, om_lambda, mass, z_lens, theta)
            results.append(result)
    
        fractional_deviations = [r["numerical_fractional_deviation_from_schwarzschild"] for r in results]
        plt.plot(om_lambdas, np.array(fractional_deviations) * 1e5, label=r'$\Omega_k = {}$'.format(om_k))

        plt.xlabel(r'$\Omega_{\Lambda}$')
        plt.ylabel(r'Fractional deviation of $\alpha$ / $10^{-5}$')
        plt.legend()

# plots Fig 6, but need to comment out the enter and exit Rhs.
def swiss_chess_multi_run_rh(om_m=0.5, mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec):
    results = []
    # om_k = 0
    om_lambdas = np.linspace(0., 0.5, 10)
    for om_lambda in om_lambdas:
        om_k = 1 - om_m - om_lambda
        # om_m = 1 - om_lambda - om_k
        result = swiss_chess_curved_single_run(om_k, om_lambda, mass, z_lens, theta)
        results.append(result)

    # plot enter and exit Rhs
    # enter_rhs = [r["Rh_enter"] for r in results]
    # plt.plot(om_lambdas, enter_rhs, label="enter_rhs")
    # exit_rhs = [r["Rh_exit"] for r in results]
    # plt.plot(om_lambdas, exit_rhs, label="exit_rhs")

    plt.xlabel("omega_lambda")
    plt.legend()

# constant omega_m, omega_k used to compensate for change in omega_m
def swiss_chess_curved_multi_run_schwarzschild_const_omega_m(om_m=0.5, mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec):
    results = []
    om_lambdas = np.linspace(0., 1.3, 30)
    for om_lambda in om_lambdas:
        om_k = 1 - om_m - om_lambda
        result = swiss_chess_curved_single_run(om_k, om_lambda, mass, z_lens, theta)
        results.append(result)

    #### plot normal fractional deviation
    fractional_deviations = [r["numerical_fractional_deviation_from_schwarzschild"] for r in results]
    plt.plot(om_lambdas, np.array(fractional_deviations) * 1e5, label=r"$\Omega_m$ = constant = {}".format(om_m))

    plt.xlabel(r"$\Omega_{\Lambda}$")
    plt.legend()


def swiss_chess_rindler(mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec):
    results = []
    om_m = 0.5
    om_lambdas = np.linspace(0., 0.99, 50)
    for om_lambda in om_lambdas:
        om_k = 1 - om_m - om_lambda
        # flat
        om_m = 1 - om_lambda
        om_k = 0
        result = swiss_chess_curved_single_run(om_k, om_lambda, mass, z_lens, theta)
        results.append(result)

    #### plot normal fractional deviation
    fractional_deviations = [r["numerical_fractional_deviation_from_schwarzschild"] for r in results]
    plt.plot(om_lambdas, np.array(fractional_deviations)*10**5, label="Numerical deviations".format(om_m))

    rindler_fractional_deviations = [r["rindler_2010_fractional_deviation_from_schwarzschild"] for r in results]
    plt.plot(om_lambdas, np.array(rindler_fractional_deviations)*10**5, label="Ishak & Rindler deviations".format(om_m))

    kantowski_fractional_deviations = [r["kantowski_fractional_deviation_from_schwarzschild"] for r in results]
    plt.plot(om_lambdas, np.array(kantowski_fractional_deviations)*10**5, label="Kantowski deviations".format(om_m))

    # plt.ylabel(r"Fractional deviation of $\alpha$ / $10^{-5}$")
    plt.xlabel(r"$\Omega_{\Lambda}$")
    plt.legend()

# constant omega_m, omega_k used to compensate for change in omega_m
# same as swiss_chess_curved_multi_run_schwarzschild_const_omega_m but plotted against kantowski
def swiss_chess_curved_multi_run_kantowski_const_omega_m(om_m=0.5, mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec):
    results = []
    om_lambdas = np.linspace(0., 1.4, 30)
    for om_lambda in om_lambdas:
        om_k = 1 - om_m - om_lambda
        result = swiss_chess_curved_single_run(om_k, om_lambda, mass, z_lens, theta)
        results.append(result)

    fractional_deviations = [r["numerical_fractional_deviation_from_kantowski"] for r in results]
    plt.plot(om_lambdas, np.array(fractional_deviations) * 10**5, label=r"$\Omega_m$ = constant = {}".format(om_m))

    plt.xlabel(r"$\Omega_{\Lambda}$")
    plt.legend()


def swiss_chess_curved_multi_run_kantowski(mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec, plot_errors=True):
    results = []
    om_ks = [-0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4]
    # om_ks = [-0.5, 0., 0.5]
    ax = plt.gca()
    for om_k in om_ks:
        results = []
        om_lambdas = np.linspace(0., 1 - om_k - 0.1, 30)
        for om_lambda in om_lambdas:
            result = swiss_chess_curved_single_run(om_k, om_lambda, mass, z_lens, theta)
            results.append(result)

        fractional_deviations = [r["numerical_fractional_deviation_from_kantowski"] for r in results]
        kantowski_higher_order = [r["higher_order_kantowski"] for r in results]
        color = next(ax._get_lines.prop_cycler)['color']
        plt.plot(om_lambdas, np.array(fractional_deviations) * 10**5, color=color, label=r'$\Omega_k = {}$'.format(om_k))
        if plot_errors:
            plt.plot(om_lambdas, kantowski_higher_order, "+", color=color, label='Higher order term(Omega_k = {})'.format(om_k))

        ##### plotting phi_tilda
        # phi_tildas = [r["phi_tilda"] for r in results]
        # plt.plot(om_lambdas, phi_tildas, label='phi_tilda (Omega_k = {})'.format(om_k))

        plt.xlabel(r"$\Omega_{\Lambda}$")
        plt.ylabel(r'Fractional deviation of $\alpha$ / $10^{-5}$')
        plt.legend()

def cosec(phi):
    return 1./np.cos(phi)

def kantowski_11(phi, r0):
    rs = 2*M
    rsr0 = rs/2/r0
    print("cosec", np.sin(np.pi - phi), rs/r0)
    result = cosec(phi) * (1 - rsr0*(-1 + 2*cosec(phi) - np.sin(phi)) + rsr0**2 * (-17./4 + 15./4*(phi - np.pi/2)/np.tan(phi) + 4*(cosec(phi))**2 + 1./4 * (np.sin(phi))**2))
    error = rsr0**3
    return result, error

def kantowski_alpha(r0, M, phi, Lambda):
    rs = 2*M
    first_term = (rs/2/r0)*np.cos(phi)*(-4*(np.cos(phi))**2 - 12*np.cos(phi)*np.sin(phi)*np.sqrt(Lambda*r0**2/3+rs/r0*(np.sin(phi))**3) + Lambda*r0**2*(8/3-20/3*(np.sin(phi))**2))
    second_term = (rs/2/r0)**2*(15/4*(2*phi-np.pi) + np.cos(phi)*(4+33/2*np.sin(phi)-4*(np.sin(phi))**2+19*(np.sin(phi))**3-64*(np.sin(phi))**5) - 12*np.log(np.tan(phi/2))*(np.sin(phi))**3)
    error = (rs / r0 + Lambda * r0**2)**(5./2)
    return first_term + second_term, error


if __name__ == '__main__':
    from utils import set_typography
    set_typography(latex=True)
    import matplotlib
    matplotlib.rcParams.update({'font.size': 13})
    # plt.rc('legend', fontsize=12)
    plt.figure(figsize=(6*1.2, 5*1.2))
    # main(filename="data/curvedpyoutput_new.csv")
    # compare_with_analytical()
    # compare_kottler_hole_size_with_frw_size()
    # main_multiple_omk(filename = "data/curvedpyoutput_withk4.csv", model=SwissCheeseKottlerModification)
    # main_multiple_omk(filename = "data/curvedpyoutput_withk5.csv", model=SwissCheese) # normal

    # print(swiss_chess_single_run(0.9))

    # A1689 strong
    # print(swiss_chess_single_run(0.7, mass=M_sun*8e13, z_lens=0.18, theta=45*one_arcsec))
    
    # A1689 weak
    # print(swiss_chess_single_run(0.7, mass=M_sun*1e15, z_lens=0.18, theta=600*one_arcsec))

    # RDCS1252-2927
    # print(swiss_chess_single_run(0.7, mass=M_sun*1e15, z_lens=1.24, theta=180*one_arcsec))

    # Elliptical galaxy strong
    # print(swiss_chess_single_run(0.7, mass=M_sun*3e11, z_lens=0.5, theta=2*one_arcsec))

    # Elliptical galaxy weak
    # print(swiss_chess_single_run(0.7, mass=M_sun*1e13, z_lens=0.5, theta=70*one_arcsec))

    # curved
    # print(swiss_chess_curved_single_run(0., 0.5, mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec))
    # print(swiss_chess_curved_single_run(0.5, 0., mass=M_sun*1e13*0.858, z_lens=0.5, theta=one_arcsec))
    # print(swiss_chess_curved_single_run(0., 0.5, mass=M_sun*1e13, z_lens=0.5, theta=one_arcsec))
    # swiss_chess_curved_multi_run()
    # swiss_chess_curved_multi_run_kantowski()

    # bigger masses
    # print(swiss_chess_curved_single_run(0., 0.5, mass=M_sun*1e14, z_lens=0.5, theta=60*one_arcsec))
    ### figure 4
    # swiss_chess_curved_multi_run(mass=M_sun*1e14, theta=10*one_arcsec)
    ##### 1e18 numbers
    # m1 = M_sun*1e18
    # theta1 = 20*60*one_arcsec
    # z1 = 2.
    m1 = M_sun*1e19
    theta1 = 60*60*one_arcsec

    m1 = M_sun*1e15
    theta1 = 60*one_arcsec
    # import pprint
    # swiss_chess_curved_single_run(0, 0.5)
    # swiss_chess_curved_multi_run_schwarzschild_const_omega_m(mass=m1, theta=theta1)
    # swiss_chess_curved_multi_run(mass=m1, theta=theta1)
    # swiss_chess_curved_multi_run_kantowski_const_omega_m(mass=m1, theta=theta1)
    # swiss_chess_curved_multi_run_kantowski(mass=m1, theta=theta1, plot_errors=False)
    # swiss_chess_curved_single_run(0, 0.5, mass=m1, theta=theta1)
    # pprint.pprint(swiss_chess_curved_single_run(0, 0.5, mass=m1, z_lens=z1, theta=theta1))
    # swiss_chess_curved_multi_run_kantowski_const_omega_m(om_m=0.5, mass=m1, theta=theta1, z_lens=z1)
    # swiss_chess_curved_multi_run_kantowski(mass=m1, theta=theta1, plot_errors=False, z_lens=z1)
    # swiss_chess_curved_multi_run_schwarzschild_const_omega_m(om_m=0.5, mass=M_sun*1e14, theta=10*one_arcsec)

    # swiss_chess_curved_multi_run_schwarzschild_const_omega_m(mass=m1, theta=theta1)
    # swiss_chess_curved_multi_run(mass=m1, theta=theta1)

    # swiss_chess_curved_multi_run_kantowski_const_omega_m(om_m=0.5, mass=m1, theta=theta1)
    # swiss_chess_curved_multi_run_kantowski(mass=m1, theta=theta1, plot_errors=False)
    plt.tight_layout()


    # Fig 6
    swiss_chess_multi_run_rh(mass=m1, theta=theta1)
    swiss_chess_rindler(mass=m1, theta=theta1)

    # swiss_chess_multi_run_rh()
    # swiss_chess_rindler()
    plt.show()

    # to get these two curves on the same graph
    # swiss_chess_curved_multi_run_schwarzschild_const_omega_m(om_m=0.5)
    # swiss_chess_curved_multi_run_schwarzschild_const_omega_m(om_m=0.2)
    # swiss_chess_curved_multi_run_schwarzschild_const_omega_m(om_m=0.8)
    # swiss_chess_curved_multi_run()

    # # Same as above, but compared against kantowski
    # swiss_chess_curved_multi_run_kantowski_const_omega_m(om_m=0.5)
    # swiss_chess_curved_multi_run_kantowski_const_omega_m(om_m=0.5)
    # swiss_chess_curved_multi_run_kantowski_const_omega_m(om_m=0.8)
    
    # swiss_chess_curved_multi_run_kantowski(plot_errors=False)
    # plt.ylabel("Deviation from Kantowski prediction")
    
    # swiss_chess_multi_run_rh()


    # swiss_chess_rindler()
    # plt.show()


    # swiss_chess_curved_single_run(0., 0.5)


    # ######## Schucker 2009 table 4
    # # ## col 1
    # schucker(10*one_arcsec, 5*one_arcsec, 0.68, om_lambda=0.77, mass=1.8*M_sun*1e13)
    # # ## col 2
    # schucker(10*one_arcsec, 5*one_arcsec, 0.68, om_lambda=0.92, mass=1.8*M_sun*1e13)
    # # ## col 3
    # schucker(10*one_arcsec, 5*one_arcsec, 0.68, om_lambda=0.61, mass=1.8*M_sun*1e13)