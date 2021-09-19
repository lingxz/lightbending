import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.integrate as spi
# plt.style.use('seaborn-white')

length_scale = 3.086e22
# H_0 = 7.56e-27 * length_scale
H_0 = 70*1000/(3.086e22)/299792458 * length_scale # 70km/s/Mpc


def set_typography(latex=False):
    # pass
    from matplotlib import rc, rcParams
    font = {'family' : 'sans-serif',
        # "serif": "Computer Modern",
        # "serif": "cm",
        # 'weight' : 'bold',
        'size'   : 18
    }
    # rc('font', **font)
    rc('figure', autolayout=True)
    # rcParams['mathtext.fontset'] = 'cm'


    # if latex:
    #     rc('text', usetex=True)
    # else:
    #     rc('text', usetex=False)

def omega_lambda2lambda(Omega_Lambda):
    return 3*Omega_Lambda*H_0**2

def kantowski_alpha(R, M, phi, Omega_Lambda):
    r0 = 1/(1/R + M/R**2 + 3/16*M**2/R**3)
    # r0 = 1/(1/R + M/R**2)
    # r0 = R*(1-M/R - 3/2*(M/R)**2 - 4*(M/R)**3)
    Lambda = omega_lambda2lambda(Omega_Lambda)
    rs = 2*M
    first_term = (rs/2/r0)*np.cos(phi)*(-4*(np.cos(phi))**2 - 12*np.cos(phi)*np.sin(phi)*np.sqrt(Lambda*r0**2/3+rs/r0*(np.sin(phi))**3) + Lambda*r0**2*(8/3-20/3*(np.sin(phi))**2))
    second_term = (rs/2/r0)**2*(15/4*(2*phi-np.pi) + np.cos(phi)*(4+33/2*np.sin(phi)-4*(np.sin(phi))**2+19*(np.sin(phi))**3-64*(np.sin(phi))**5) - 12*np.log(np.tan(phi/2))*(np.sin(phi))**3)
    return first_term + second_term


# def r2chi(k, r):
#     if k == 0:
#         return r
#     if k > 0:
#         return np.arcsin(np.sqrt(k)*r)/np.sqrt(k)
#     if k < 0:
#         return np.arcsinh(np.sqrt(-k)*r)/np.sqrt(-k)

# def chi2r(k, chi):
#     if k == 0:
#         return chi
#     if k > 0:
#         return np.sin(np.sqrt(k)*chi)/np.sqrt(k)
#     if k < 0:
#         return np.sinh(np.sqrt(-k)*chi)/np.sqrt(-k)

def r2chi(k, r):
    if k == 0:
        return r
    if k > 0:
        return np.arcsin(np.sqrt(k)*r)*np.sqrt(k)
    if k < 0:
        return np.arcsinh(np.sqrt(-k)*r)*np.sqrt(-k)

def chi2r(k, chi):
    if k == 0:
        return chi
    if k > 0:
        return np.sin(chi/np.sqrt(k))/np.sqrt(k)
    if k < 0:
        return np.sinh(chi/np.sqrt(-k))/np.sqrt(-k)

def omk2k(om_k, H0):
    return -H0**2*om_k


def plot_rs_multiple_ks(filename, latex=False, filenames=None):
    set_typography(latex)

    combined = pd.read_csv(filename)
    grouped = combined.groupby(combined['om_ks'])


    count = 0
    for _, df in grouped:
        if count == 7:
            break
        current_omk = df['om_ks'].values[0]
        numerical_res =[]
        preds_frw = []
        preds_ishak = []
        preds_kantowski = []
        kant_higher_order_ratio = []

        for index, row in df.iterrows():
            M = row.M
            Lambda = 3*row.om_lambdas*H_0**2

            k = omk2k(row.om_ks, H_0)

            numerical = chi2r(k, r2chi(k, row.raw_rs)+r2chi(k, row.comoving_lens))/row.raw_rs*row.theta

            R = (row.DL*row.theta)
            A_frw = 4*M/R + 15*np.pi*M**2/4/R**2 + 401/12*M**3/R**3

            A_ishak = 4*M/R + 15*np.pi*M**2/4/R**2 + 305/12*M**3/R**3 - Lambda*R*row.exit_rhs/3

            A_kantowski = -kantowski_alpha(R, M, row.enter_phis, row.om_lambdas)

            r0 = 1/(1/R + M/R**2 + 3/16*M**2/R**3)
            # for comparing with kantowski predictions
            extra_term_ratio = np.abs((2*M/r0 + omega_lambda2lambda(row.om_lambdas)*r0**2)**(5/2)/(4*M/r0*(np.cos(row.enter_phis))**3))

            preds_frw.append(A_frw/row.theta)
            preds_ishak.append(A_ishak/row.theta)
            preds_kantowski.append(A_kantowski/row.theta)
            kant_higher_order_ratio.append(extra_term_ratio)
            numerical_res.append(numerical/row.theta)

        df['preds_frw'] = preds_frw
        df['preds_ishak'] = preds_ishak
        df['preds_kantowski'] = preds_kantowski
        df['kant_higher_order_ratio'] = kant_higher_order_ratio
        df['numerical_res'] = numerical_res

        df['numerical'] = df.numerical_res/df.preds_frw - 1
        df['ishak'] = df.preds_ishak/df.preds_frw - 1
        df['kantowski'] = df.preds_kantowski/df.preds_frw - 1
        df['numerical_kantowski'] = df.numerical_res / df.preds_kantowski - 1
        df['numerical_zero'] = (df.numerical +1 ) / (df.numerical + 1).head(1).values[0] -1

        stats = df[['numerical_res', 'om_lambdas', 'numerical', 'ishak', 'kantowski', 'numerical_kantowski', 'kant_higher_order_ratio', 'numerical_zero']].groupby('om_lambdas').agg(['mean', 'std', 'count'])
        stats.columns = [' '.join(col).strip() for col in stats.columns.values]
        stats['numerical mean std'] = stats['numerical std']/np.sqrt(stats['numerical count'])

        scale = 1e-5
        plt.plot(stats.index, stats['numerical mean']/scale, '-', label='Numerical results (Omega_k = {})'.format(current_omk))

        if latex:
            plt.xlabel(r'$\Omega_{\Lambda}$')
            plt.ylabel(r'Fractional deviation of $\frac{D_S}{D_{LS}}$ / $10^{-4}$')
        else:
            plt.xlabel('Omega_Lambda')
            plt.ylabel('Mean fractional deviation/10^-5')

        # plt.plot(stats.index, [0/scale]*len(stats.index), 'r-', label='FRW predictions')

        plt.legend(bbox_to_anchor=(1.05, 1))
        count += 1
        if filenames:
            plt.savefig(filenames[0], dpi=400, transparent=True)


def calculate_frw():
    df = pd.read_csv("curvedpyoutput_withk3.csv")
    for index, row in df.iterrows():
        M = row.M
        k = omk2k(row.om_ks, H_0)

        numerical = chi2r(k, r2chi(k, row.raw_rs)+r2chi(k, row.comoving_lens))/row.raw_rs*row.theta
        R = (row.DL*row.theta)
        A_frw = 4*M/R + 15*np.pi*M**2/4/R**2 + 401/12*M**3/R**3
        print(A_frw, numerical, numerical/A_frw - 1)
    # return theta

# calculate_frw()

def plot_rs(filename, plot_ishak=True, plot_kantowski=True, latex=False, filenames=None):
    set_typography(latex)

    df = pd.read_csv(filename)

    numerical_res =[]
    preds_frw = []
    preds_ishak = []
    preds_kantowski = []
    kant_higher_order_ratio = []

    for index, row in df.iterrows():
        M = row.M
        Lambda = 3*row.om_lambdas*H_0**2

        k = omk2k(row.om_ks, H_0)

        numerical = chi2r(k, r2chi(k, row.raw_rs)+r2chi(k, row.comoving_lens))/row.raw_rs*row.theta

        R = (row.DL*row.theta)
        A_frw = 4*M/R + 15*np.pi*M**2/4/R**2 + 401/12*M**3/R**3

        A_ishak = 4*M/R + 15*np.pi*M**2/4/R**2 + 305/12*M**3/R**3 - Lambda*R*row.exit_rhs/3

        A_kantowski = -kantowski_alpha(R, M, row.enter_phis, row.om_lambdas)

        r0 = 1/(1/R + M/R**2 + 3/16*M**2/R**3)
        # for comparing with kantowski predictions
        extra_term_ratio = np.abs((2*M/r0 + omega_lambda2lambda(row.om_lambdas)*r0**2)**(5/2)/(4*M/r0*(np.cos(row.enter_phis))**3))

        preds_frw.append(A_frw/row.theta)
        preds_ishak.append(A_ishak/row.theta)
        preds_kantowski.append(A_kantowski/row.theta)
        kant_higher_order_ratio.append(extra_term_ratio)
        numerical_res.append(numerical/row.theta)

    df['preds_frw'] = preds_frw
    df['preds_ishak'] = preds_ishak
    df['preds_kantowski'] = preds_kantowski
    df['kant_higher_order_ratio'] = kant_higher_order_ratio
    df['numerical_res'] = numerical_res

    df['numerical'] = df.numerical_res/df.preds_frw - 1
    df['ishak'] = df.preds_ishak/df.preds_frw - 1
    df['kantowski'] = df.preds_kantowski/df.preds_frw - 1
    df['numerical_kantowski'] = df.numerical_res / df.preds_kantowski - 1
    df['numerical_zero'] = (df.numerical +1 ) / (df.numerical + 1).head(1).values[0] -1

    stats = df[['om_lambdas', 'numerical', 'ishak', 'kantowski', 'numerical_kantowski', 'kant_higher_order_ratio', 'numerical_zero']].groupby('om_lambdas').agg(['mean', 'std', 'count'])
    stats.columns = [' '.join(col).strip() for col in stats.columns.values]
    stats['numerical mean std'] = stats['numerical std']/np.sqrt(stats['numerical count'])
    stats['numerical_kantowski mean std'] = stats['numerical_kantowski std'] / np.sqrt(stats['numerical_kantowski count'])
    stats['ishak mean std'] = stats['ishak std']/np.sqrt(stats['ishak count'])


    scale = 1e-5
    plt.plot(stats.index, stats['numerical mean']/scale, 'k.', label='Numerical results')

    if latex:
        plt.xlabel(r'$\Omega_{\Lambda}$')
        plt.ylabel(r'Fractional deviation of $\frac{D_S}{D_{LS}}$ / $10^{-4}$')
    else:
        plt.xlabel('Omega_Lambda')
        plt.ylabel('Mean fractional deviation/10^-5')

    if plot_ishak:
        plt.plot(stats.index, stats['ishak mean']/scale, 'g-', label='Rindler and Ishak predictions')
    if plot_kantowski:
        plt.plot(stats.index, stats['kantowski mean']/scale, 'b-', label="Kantowski predictions")
    plt.plot(stats.index, [0/scale]*len(stats.index), 'r-', label='FRW predictions')

    plt.legend()

    if filenames:
        plt.savefig(filenames[0], dpi=400, transparent=True)

    plt.figure()
    scale2 = 1e-6
    plt.plot(stats.index, -stats['numerical_kantowski mean']/scale2, '.', label='Numerical deviation from Kantowski')
    plt.errorbar(stats.index, -stats['numerical_kantowski mean']/scale2, yerr=stats['numerical_kantowski mean std']/scale2, label='__nolegend__', linestyle='none')
    if plot_kantowski:
        plt.plot(stats.index, stats['kant_higher_order_ratio mean']/scale2, label='Neglected term ratio')
    
    if latex:
        plt.xlabel(r'$\Omega_{\Lambda}$')
        plt.ylabel(r'Fractional deviation of $\frac{D_S}{D_{LS}}$ / $10^{-6}$')
    
    plt.legend()
    if filenames:
        plt.savefig(filenames[1],  dpi=400, transparent=True)

    plt.figure()
    plt.title("Cannot remember what this graph was for")
    scale3 = 1e-6
    plt.plot(stats.index, stats['numerical_kantowski mean']/scale3, '.', label='Deviations from Kantowski')
    plt.plot(stats.index, stats['numerical_zero mean']/scale3, '.', label='Deviations from numerical Schwarzschild case')
    plt.legend()


def plot_rs_ltb(filename, plot_ishak=True, plot_kantowski=True, latex=False, filenames=None):
    set_typography(latex)

    df = pd.read_csv(filename)

    preds_frw = []
    preds_ishak = []
    preds_kantowski = []
    numerical_res = []
    kant_higher_order_ratio = []
    for index, row in df.iterrows():
        M_total = row.M
        Lambda = 3*row.om_lambdas*H_0**2
        enclosed_r = row.DL*row.theta
        rho_frw_initial = (1-row.om_lambdas)*3*H_0**2/(8*np.pi)
        r_h = (3*M_total/(4*np.pi*rho_frw_initial))**(1./3)
        
        def mass(r):
            initial_rh = (3*M_total/(4*np.pi*rho_frw_initial))**(1./3)
            rlimit = initial_rh*0.5
            if r > rlimit:
                return M_total
            else:
                c = 10
                Rvir = rlimit/100
                Rs = Rvir/c
                rho0 = M_total/(4*np.pi*Rs**3*(np.log((Rs + rlimit)/Rs) - rlimit/(Rs + rlimit)))
                return 4*np.pi*rho0*Rs**3*(np.log((Rs+r)/Rs) - r/(Rs+r))

                # integral, error = spi.quad(lambda r1: rho(r1)*r1**2, 0, r)
                # return 4*np.pi*integral

        def rho(r):
            initial_rh = (3*M_total/(4*np.pi*rho_frw_initial))**(1./3)
            rlimit = initial_rh*0.5
            if r > rlimit:
                return 0
            else:
                c = 10
                Rvir = rlimit/100
                Rs = Rvir/c
                rho0 = M_total/(4*np.pi*Rs**3*(np.log((Rs + rlimit)/Rs) - rlimit/(Rs + rlimit)))
                return rho0/(r/Rs)/(1 + r/Rs)**2

        def projected_mass(r):
            initial_rh = (3*M_total/(4*np.pi*rho_frw_initial))**(1./3)
            rlimit = initial_rh*0.5
            if r > rlimit:
                return M_total
            else:
                c = 10
                Rvir = rlimit/100
                Rs = Rvir/c
                
                g = 1/(np.log(1+c) - c/(1+c))
                Rtilde = r / Rvir
                if r > Rs:
                    c_inverse = np.arccos(1/c/Rtilde)
                else:
                    c_inverse = np.arccosh(1/c/Rtilde)
                Rtilde1 = rlimit/Rvir
                return 1/(np.log(1+c*Rtilde1) - c*Rtilde1/(1+c*Rtilde1))*M_total*(c_inverse/np.abs(c**2*Rtilde**2-1)**(1/2) + np.log(c*Rtilde/2))
                # return 1/(np.log(1+c*Rtilde) - c*Rtilde/(1+c*Rtilde))*mass(r)*(c_inverse/np.abs(c**2*Rtilde**2-1)**(1/2) + np.log(c*Rtilde/2))

        R = (row.DL*row.theta)
        M = projected_mass(enclosed_r)
        numerical = (row.raw_rs + row.comoving_lens)/row.raw_rs*row.theta
        
        A_frw = 4*M/R + 15*np.pi*M**2/4/R**2 + 401/12*M**3/R**3

        A_ishak = 4*M/R + 15*np.pi*M**2/4/R**2 + 305/12*M**3/R**3 - Lambda*R*row.exit_rhs/3

        A_kantowski = -kantowski_alpha(R, M, row.enter_phis, row.om_lambdas)
        
        # r0 = 1/(1/R + M/R**2)
        r0 = 1/(1/R + M/R**2 + 3/16*M**2/R**3)
        extra_term_ratio = np.abs((2*M/r0 + omega_lambda2lambda(row.om_lambdas)*r0**2)**(5/2)/(4*M/r0*(np.cos(row.enter_phis))**3))
        

        
        preds_frw.append(A_frw/row.theta)
        preds_ishak.append(A_ishak/row.theta)
        preds_kantowski.append(A_kantowski/row.theta)
        numerical_res.append(numerical/row.theta)
        kant_higher_order_ratio.append(extra_term_ratio)


    df['preds_frw'] = preds_frw
    df['preds_ishak'] = preds_ishak
    df['preds_kantowski'] = preds_kantowski
    # df['kant_higher_order_ratio'] = kant_higher_order_ratio
    df['numerical_res'] = numerical_res
    df['kant_higher_order_ratio'] = kant_higher_order_ratio

    df['numerical'] = df.numerical_res/df.preds_frw - 1
    df['ishak'] = df.preds_ishak/df.preds_frw - 1
    df['kantowski'] = df.preds_kantowski/df.preds_frw - 1
    df['numerical_kantowski'] = df.numerical_res / df.preds_kantowski - 1

    stats = df[['om_lambdas', 'numerical', 'ishak', 'kantowski', 'numerical_kantowski', 'kant_higher_order_ratio']].groupby('om_lambdas').agg(['mean', 'std', 'count'])
    # stats = df[['om_lambdas', 'om_ks', 'numerical', 'ishak', 'kantowski', 'numerical_kantowski', 'kant_higher_order_ratio']].groupby('om_ks').agg(['mean', 'std', 'count'])
    stats.columns = [' '.join(col).strip() for col in stats.columns.values]
    stats['numerical mean std'] = stats['numerical std']/np.sqrt(stats['numerical count'])
    stats['numerical_kantowski mean std'] = stats['numerical_kantowski std'] / np.sqrt(stats['numerical_kantowski count'])
    # stats['numerical first order mean std'] = stats['numerical first order std']/np.sqrt(stats['numerical first order count'])
    stats['ishak mean std'] = stats['ishak std']/np.sqrt(stats['ishak count'])


    scale = 1e-5
    plt.plot(stats.index, stats['numerical mean']/scale, 'k.', label='Numerical results')
    # plt.plot(stats.index, stats['numerical first order mean']/scale, 'b-', label='with second order corrections')
    # plt.errorbar(stats.index, stats['numerical mean']/scale, yerr=stats['numerical mean std']/scale, linestyle='none', label='__nolegend__')

    if latex:
        plt.xlabel(r'$\Omega_{\Lambda}$')
        plt.ylabel(r'Fractional deviation of $\frac{D_S}{D_{LS}}$ / $10^{-5}$')
    else:
        plt.xlabel('Omega_Lambda')
        plt.ylabel('Mean fractional deviation/10^-5')

    if plot_ishak:
        plt.plot(stats.index, stats['ishak mean']/scale, 'g-', label='Rindler and Ishak predictions')
    if plot_kantowski:
        plt.plot(stats.index, stats['kantowski mean']/scale, 'b-', label="Kantowski predictions")
    plt.plot(stats.index, [0/scale]*len(stats.index), 'r-', label='FRW predictions')
    # plt.ylim((-0.0008, 0.0008))
    if plot_ishak or plot_kantowski:
        plt.legend()

    if filenames:
        plt.savefig(filenames[0], dpi=400, transparent=True)


    plt.figure()
    scale2 = 1e-6
    plt.plot(stats.index, -stats['numerical_kantowski mean']/scale2, '.', label='Numerical deviation from Kantowski')
    plt.errorbar(stats.index, -stats['numerical_kantowski mean']/scale2, yerr=stats['numerical_kantowski mean std']/scale2, label='__nolegend__', linestyle='none')
    if plot_kantowski:
        plt.plot(stats.index, stats['kant_higher_order_ratio mean']/scale2, label='Neglected term ratio')
    
    if latex:
        plt.xlabel(r'$\Omega_{\Lambda}$')
        plt.ylabel(r'Fractional deviation of $\frac{D_S}{D_{LS}}$ / $10^{-6}$')
    
    plt.legend()
    if filenames:
        plt.savefig(filenames[1],  dpi=400, transparent=True)


class LTB:
    M = 4.776409591704472e-08
    rho_frw_initial = 6.49705928222e-09
    rlimit_ratio = 1/0.5

    def initial_rh(self):
        return (3*self.M/(4*np.pi*self.rho_frw_initial))**(1./3)

    def rs(self):
        initial_rh = (3*self.M/(4*np.pi*self.rho_frw_initial))**(1./3)
        return np.linspace(0., initial_rh, 1000)[1:]

    def mass(self, r):
        initial_rh = (3*self.M/(4*np.pi*self.rho_frw_initial))**(1./3)
        rlimit = initial_rh/self.rlimit_ratio
        if r > rlimit:
            return self.M
        else:
            c = 10
            Rvir = rlimit/100
            Rs = Rvir/c
            rho0 = self.M/(4*np.pi*Rs**3*(np.log((Rs + rlimit)/Rs) - rlimit/(Rs + rlimit)))
            return 4*np.pi*rho0*Rs**3*(np.log((Rs+r)/Rs) - r/(Rs+r))
        # else:
        #     integral, error = spi.quad(lambda r1: self.rho(r1)*r1**2, 0, r)
        #     return 4*np.pi*integral
        integral, error = spi.quad(lambda r1: self.rho(r1)*r1**2, 0, r)
        return 4*np.pi*integral

    def rho(self, r):
        initial_rh = (3*self.M/(4*np.pi*self.rho_frw_initial))**(1./3)
        rlimit = initial_rh/self.rlimit_ratio
        if r > rlimit:
            return 0
        else:
            c = 10
            Rvir = rlimit/100
            Rs = Rvir/c
            rho0 = self.M/(4*np.pi*Rs**3*(np.log((Rs + rlimit)/Rs) - rlimit/(Rs + rlimit)))
            return rho0/(r/Rs)/(1 + r/Rs)**2

def ltb_graphs():
    set_typography(True)
    solar_mass = 1474 / 3.086e22
    ltb = LTB()
    rs = ltb.rs()
    rhos = np.array([ltb.rho(r) for r in rs])
    masses = np.array([ltb.mass(r) for r in rs])
    # plt.semilogx(rs/ltb.initial_rh(), masses/solar_mass)
    plt.semilogx(rs/ltb.initial_rh(), masses/solar_mass)

    latex = True
    plt.ylabel(r'$M / M_{\odot}$')
    plt.xlabel(r'$R/ R_{h0}$')
    plt.savefig('report/images/ltb-mass.png', dpi=400, transparent=True)
    # else:
    #     plt.xlabel('mass sun')
    #     plt.ylabel('Fractional deviation')


    plt.figure()
    plt.loglog(rs/ltb.initial_rh(), rhos/ltb.rho_frw_initial)
    plt.xlabel(r'$R/ R_{h0}$')
    plt.ylabel(r'$\rho / \rho_m$')
    plt.savefig('report/images/ltb-density.png', dpi=400, transparent=True)

    def f(P, t):
        P = P[0]
        r = ltb.initial_rh() - t
        current_rho = ltb.rho(r)
        current_m = ltb.mass(r)
        Lambda = 0
        E = -2*current_m/r - Lambda*r**2/3
        return [(current_rho + P)/2/(1+E)/r*(Lambda*r**2 - 8*np.pi*P*r**2 + E)]
    from scipy.integrate import odeint, solve_ivp, ode
    # pressure_rs = np.linspace(0., 0.6*ltb.initial_rh(), num=5000, endpoint=True)[1:]
    r = ode(lambda r, p: f(p, r)).set_initial_value(0).set_integrator('vode')

    pressure_rs = []
    pressure = []
    dt = 1/5000
    while r.successful() and r.t < 0.999*ltb.initial_rh():
        pressure.append(r.integrate(r.t+dt)[0])
        pressure_rs.append(r.t+dt)

    pressure = np.array(pressure)
    pressure_rs = np.array(pressure_rs)
    pressure_rs = ltb.initial_rh() - pressure_rs
    # pressure = odeint(f, [0], pressure_rs, rtol=1e-30, atol=1e-100)[:,0]
    # print(pressure)

    plt.figure()
    plt.semilogx(pressure_rs/ltb.initial_rh(), pressure/ltb.rho_frw_initial)
    plt.xlabel(r'$R/ R_{h0}$')
    plt.ylabel(r'$P / \rho_m$')
    plt.savefig('report/images/ltb-pressure.png', dpi=400, transparent=True)

from astropy.cosmology import LambdaCDM

def binary_search(cosmo, start, end, answer):
    mid = (end+start)/2
    res = cosmo.angular_diameter_distance(mid).value
    if np.isclose(res, answer, rtol=1e-14, atol=1e-14):
        return mid
    if res < answer:
        return binary_search(cosmo, mid, end, answer)
    else:
        return binary_search(cosmo, start, mid, answer)