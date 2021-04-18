'''
Simulation of integration in Butcher's paper
'''
import numpy as np
import scipy.integrate as spi

length_scale = 3.086e22 # 1 mega parsec

tols = {
    'atol': 1e-120,
    'rtol': 1e-14,
    'method': 'RK45',
}

H_0 = 70*1000/(3.086e22)/299792458 * length_scale # 70km/s/Mpc

def kottler(eta, w, p):
    L, M = p
    R, Rdot, phi = w
    Rddot = L**2 * (R - 3*M) / R**4
    phidot = L/R**2

    return [
        Rdot,
        Rddot,
        phidot,
    ]

class KottlerBase:
    def __init__(self, M, start_r, Omega_m, Omega_Lambda, angle_to_horizontal):
        self.M = M
        self.start_r = start_r
        self.Omega_Lambda = Omega_Lambda
        self.Omega_m = Omega_m
        self.Lambda = 3*Omega_Lambda*H_0**2
        self.angle_to_horizontal = angle_to_horizontal

        self.kottler_initial = None
        self.kottler_parameters = None
        self.kottler_final = None
        self.final_angle = None
        self.deflection_angle = None
        
        self.theta_obs = None

    def run(self):
        self.generate_kottler_initial()
        self.integrate_kottler()
        self.calculate_deflection()

    def generate_kottler_initial(self):
        raise NotImplementedError


    def integrate_kottler(self):
        initial = self.kottler_initial[:]
        parameters = self.kottler_parameters
        
        # Estimate distance of closest approach. 
        # Closest approach is when R is minimum, which is when Rdot is zero.
        def closest_approach(t, y): return y[1]
        closest_approach.terminal = False

        # Impact parameter is when 180deg - phi = 90deg - theta
        def impact_parameter(t, y): return np.pi - y[2] - (np.pi/2 - self.angle_to_horizontal)
        impact_parameter.terminal = False

        def reach_horizon(t, y): return y[2]
        reach_horizon.terminal = True
        # reach_horizon.direction = -1 # value goes from positive to negative.
        reach_horizon.direction = 0

        events = [impact_parameter, closest_approach, reach_horizon]
        sol = spi.solve_ivp(lambda t, y: kottler(t, y, parameters), [0, 10], initial, dense_output=True, events=events, **tols)
        self.impact_parameter = sol.y_events[0][-1][0]
        self.closest_approach = sol.y_events[1][-1][0]
        self.kottler_final = sol.y_events[-1][-1]

    def get(self, attribute, param):
        arr = getattr(self, attribute)
        kottler_coords_to_index = ["R", "Rdot", "phi"]
        kottler_params_to_index = ["L", "M", "Omega_Lambda", "Omega_m"]
        if attribute in ["kottler_initial", "kottler_final"]:
            index = kottler_coords_to_index.index(param)
        elif attribute in ["kottler_parameters"]:
            index = kottler_params_to_index.index(param)
        return arr[index]
    
    def calculate_deflection(self):
        raise NotImplementedError

    def debug_logs(self):
        print("===== Initial parameters =====")
        self._print_with_names([self.Omega_m, self.Omega_Lambda, self.start_r, self.angle_to_horizontal], ["Omega_m", "Omega_Lambda", "start_r", "angle_to_horizontal"])

        print("===== Initial Kottler coordinates =====")
        self._print_with_names(self.kottler_initial, ["R", "Rdot", "phi"])
        self._print_with_names(self.kottler_parameters, ["L", "M"])
        
        print("===== Final Kottler coordinates =====")
        self._print_with_names(self.kottler_final, ["R", "Rdot", "phi"])
        self._print_with_names([self.impact_parameter, self.closest_approach], ["Impact parameter", "Closest approach"])

        print("===== Angles =====")
        self._print_with_names([self.angle_to_horizontal, self.deflection_angle], ["angle_to_horizontal", "deflection_angle"])

    @staticmethod    
    def _print_with_names(values, names):
        for name, value in zip(names, values):
            print("{}: {}".format(name, value))

class ButcherBackward(KottlerBase):
    def generate_kottler_initial(self):
        # L, M, Omega_Lambda, Omega_m = p
        # R, Rdot, phi = w

        initial_R = self.start_r
        initial_phi = np.pi
        initial_Rdot = -initial_R

        f = 1 - 2*self.M/initial_R - self.Lambda * initial_R **2 / 3
        v = np.sqrt(1-f)
        # this exact calculation seems to give less precision, so we're using eqn 20.
        # theta_stat = np.arccos((np.cos(self.angle_to_horizontal) + v) / (1 + v * np.cos(self.angle_to_horizontal)))
        theta_stat = np.sqrt((1 - v)/(1 + v)) * self.angle_to_horizontal
        initial_phidot = np.tan(theta_stat) * initial_Rdot / initial_R / np.sqrt(f)

        # eq 20        
        # initial_phidot = np.sqrt((1-v) / (1 + v)) * initial_Rdot / initial_R / np.sqrt(f) * self.angle_to_horizontal
        
        L = initial_R**2 * initial_phidot

        self.kottler_initial = [initial_R, initial_Rdot, initial_phi]
        self.kottler_parameters = [L, self.M]


    def calculate_deflection(self):
        self.theta_obs = self.angle_to_horizontal
        Rdot = self.get("kottler_final", "Rdot")
        R = self.get("kottler_final", "R")
        L = self.get("kottler_parameters", "L")
        phidot = L / R**2

        # butcher Eq 30
        rO = self.start_r
        rs = R
        D_LS = rs
        D_S = (rO + rs) / (1 + rO*np.sqrt(self.Lambda/3))
        D_L = rO / (1 + rO * np.sqrt(self.Lambda/3))

        self.impact_parameter = D_L * self.angle_to_horizontal
        self.deflection_angle  = abs(self.angle_to_horizontal) * D_S / D_LS  # eq 37
        
        # fO = 1 - 2*self.M / self.start_r - self.Lambda *  self.start_r **2 / 3
        # vO = np.sqrt(1-fO)
        # self.deflection_angle = 2 * np.sqrt(self.M * D_S / D_L/D_LS) + 15 * np.pi * self.M * D_S / 32 / D_L/D_LS  # Eq 38


class ButcherForward(KottlerBase):
    def generate_kottler_initial(self):
        # L, M, Omega_Lambda, Omega_m = p
        # R, Rdot, phi = w

        initial_R = self.start_r
        initial_phi = np.pi
        initial_Rdot = -initial_R

        theta_stat = self.angle_to_horizontal
        f = 1 - 2*self.M/initial_R - self.Lambda * initial_R **2 / 3
        initial_phidot = np.tan(theta_stat) * initial_Rdot / initial_R / np.sqrt(f)
        
        L = initial_R**2 * initial_phidot

        self.kottler_initial = [initial_R, initial_Rdot, initial_phi]
        self.kottler_parameters = [L, self.M]
    
    def calculate_deflection(self):
        Rdot = self.get("kottler_final", "Rdot")
        R = self.get("kottler_final", "R")
        f = 1 - 2*self.M/R - self.Lambda * R**2 / 3
        L = self.get("kottler_parameters", "L")
        phidot = L / R**2

        theta_stat = np.arctan(R * phidot / Rdot * np.sqrt(f))
        v = np.sqrt(1 - f)
        
        # butcher
        # theta_obs = np.arctan(np.sin(theta_stat) * np.sqrt(1 - v**2) / (np.cos(theta_stat) - v))
        # theta_rindler = np.arccos(abs(Rdot/phidot) / np.sqrt(Rdot**2/phidot**2 + f * R**2))

        # butcher
        rs = self.start_r
        rO = R
        D_LS = rs
        D_S = (rO + rs) / (1 + rO*np.sqrt(self.Lambda/3))
        theta_obs = np.sqrt((1 + v)/(1-v))* theta_stat  # eq 19
        self.deflection_angle = abs(theta_obs) * D_S / D_LS
    
        D_L = rO / (1 + rO * np.sqrt(self.Lambda/3))
        self.impact_parameter = D_L * theta_obs
        self.theta_obs = theta_obs

class SchwarzschildForward(KottlerBase):
    def generate_kottler_initial(self):
        # L, M, Omega_Lambda, Omega_m = p
        # R, Rdot, phi = w
        initial_R = self.start_r
        initial_phi = np.pi
        initial_Rdot = -initial_R
        # schwarzschild calculation: coordinate angle
        initial_phidot = np.tan(self.angle_to_horizontal) * initial_Rdot/initial_R
        
        L = initial_R**2 * initial_phidot

        self.kottler_initial = [initial_R, initial_Rdot, initial_phi]
        self.kottler_parameters = [L, self.M]

    def calculate_deflection(self):
        Rdot = self.get("kottler_final", "Rdot")
        R = self.get("kottler_final", "R")
        f = 1 - 2*self.M/R - self.Lambda * R**2 / 3
        L = self.get("kottler_parameters", "L")
        phidot = L / R**2
        angle_now = np.abs(np.arctan(R * phidot / Rdot))
        self.deflection_angle = angle_now + self.angle_to_horizontal
        

# Note: Seems like not working yet
class RindlerIshakForward(SchwarzschildForward):
    def calculate_deflection(self):
        Rdot = self.get("kottler_final", "Rdot")
        R = self.get("kottler_final", "R")
        f = 1 - 2*self.M/R - self.Lambda * R**2 / 3
        L = self.get("kottler_parameters", "L")
        phidot = L / R**2
        theta_rindler = np.arccos(abs(Rdot/phidot) / np.sqrt(Rdot**2/phidot**2 + f * R**2))
        self.deflection_angle = 2 * theta_rindler

def single_run(om_lambda):
    M = 1474e13 / length_scale
    om_m = 1 - om_lambda
    start_r = 1782
    one_arcsec = 1/3600/180*np.pi
    theta = one_arcsec

    solution = ButcherBackward(M, start_r, om_m, om_lambda, theta)

    solution.run()

    # schwarzschild_deflection = np.abs(4 * M / solution.closest_approach + (-4 + 15*np.pi/4) * (M/solution.closest_approach)**2 + (122/3 -15*np.pi/2) * M**3 / solution.closest_approach**3)
    schwarzschild_deflection = np.abs(4 * M / solution.impact_parameter + (15*np.pi/4) * (M/solution.impact_parameter)**2)
    print("third term ratio", 401 /12 * (M/solution.impact_parameter)**3 / schwarzschild_deflection)
    if solution.Lambda > 0 and solution.theta_obs is not None:
        # higher_order_butcher = solution.theta_obs **2 / solution.Lambda/solution.start_r**2
        higher_order_butcher = solution.M **3 / solution.impact_parameter**3 / solution.Lambda/solution.start_r **2
    else:
        higher_order_butcher = 0
    print("HIGHERORDER", higher_order_butcher / solution.deflection_angle)
    print((1 - schwarzschild_deflection / solution.deflection_angle))
    return {
        "fractional_deviation": 1 - schwarzschild_deflection / solution.deflection_angle,
        "higher_order_butcher": higher_order_butcher / solution.deflection_angle
    }

def multi_run():
    results = []
    om_lambdas = np.linspace(0.1, 0.9, 30)
    for om_lambda in om_lambdas:
        result = single_run(om_lambda)
        results.append(result)
    
    fractional_deviations = [r["fractional_deviation"] for r in results]
    higher_order_butcher = [r["higher_order_butcher"] for r in results]
    import matplotlib.pyplot as plt
    plt.plot(om_lambdas, higher_order_butcher, "r+", label="Butcher higher order terms (eq 41)")
    plt.plot(om_lambdas, fractional_deviations, "b+", label="Fractional deviations from predictions")
    plt.xlabel("omega_lambda")
    plt.legend()
    plt.title("Butcher backward propagation")
    plt.savefig("butcher_backward.png")


multi_run()


