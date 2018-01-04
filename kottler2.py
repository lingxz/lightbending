import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import time

INTEGRATOR = 'vode'
INTEGRATOR_PARAMS = {
    'atol': 1e-110, 
    # 'atol': 0,
    # 'rtol': 0,
    'rtol': 1e-15,
    'nsteps': 100000000,
    # 'method': 'bdf',
}

def get_angle(r, phi, rdot, phidot):
    res = np.arctan((rdot*np.sin(phi)+r*np.cos(phi)*phidot)/(rdot*np.cos(phi)-r*np.sin(phi)*phidot))
    return res

def kottler(eta, w, p):
    E, L, M, Omega_Lambda, Omega_m, H_0 = p
    t, r, rdot, phi = w

    Lambda = 3*Omega_Lambda*H_0**2
    f = 1 - 2*M/r - Lambda/3*r**2
    rddot = L**2 * (r - 3*M) / r**4
    tdot = E / f
    phidot = L/r**2
    # r_h_t = (1 - 2*M/r_h - Lambda/3*r_h**2) * np.sqrt(2*M/r_h + Lambda/3*r_h**2)

    return [
        # r_h_t*tdot,
        tdot,
        rdot,
        rddot,
        phidot,
    ]

length_scale = 3.086e23 # mega parsec
M = 0.5e15 / length_scale
b =  0.000121188611029
# initial_x = -b
# initial_tdot = 1
# initial_r = np.sqrt(b**2 + initial_x**2)
# initial_phi = np.arccos(initial_x / initial_r)
# initial_rdot = np.cos(initial_phi)
# initial_phidot = -np.sqrt((initial_tdot**2 - initial_rdot**2) / initial_r**2)
# initial_t = 0
L = -0.0005068014501392157
initial_phidot = L/0.083290272382647856**2
initial = [0, 0.083290272382647856, -4.2969555464529909, 3.1402043118034899]
start_angle = get_angle(0.083290272382647856, 3.1402043118034899, -4.2969555464529909, initial_phidot)
# initial = [initial_t, initial_r, initial_rdot, initial_phi]
# f = 1-2*M/initial_r
# E = f*initial_tdot
# L = initial_r**2*initial_phidot**2
# p = [E, L, M, 0, 1, 1]
p = [-4.2968965738964089, -0.0005068014501392157, 1.6202203499675955e-09, 0, 1, 0.0023330160000000003]

solver = spi.ode(kottler).set_integrator(INTEGRATOR, **INTEGRATOR_PARAMS)
solver.set_initial_value(initial, 0).set_f_params(p)
sol = []
first_time = True
while solver.successful():
    dt = 2e-7
    solver.integrate(solver.t + dt)
    sol.append(list(solver.y))
    if solver.y[3] < np.pi/2 and first_time:
        print("turning point in kottler metric:", solver.y[1])
        first_time = False
    if solver.y[1] * np.sin(solver.y[3]) < 0:  # stop when it crosses the axis
        break

sol = np.array(sol)
r = sol[:,1]
phi = sol[:,3]
rdot = sol[:,2]
phidot = L/r**2
deflection = get_angle(r, phi, rdot, phidot)
print("deflection:", start_angle - deflection[-1])
print(r.min(), b)
expected = 4*M/r.min()
print("expected:", expected)
