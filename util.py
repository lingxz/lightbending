import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.integrate as spi
# plt.style.use('seaborn-white')

length_scale = 3.086e22
H_0 = 7.56e-27 * length_scale

def omega_lambda2lambda(Omega_Lambda):
    return 3*Omega_Lambda*H_0**2

def kantowski_alpha(R, M, phi, Omega_Lambda):
    # r0 = 1/(1/R + M/R**2 + 3/16*M**2/R**3)
    r0 = 1/(1/R + M/R**2)
    # r0 = R*(1-M/R - 3/2*(M/R)**2 - 4*(M/R)**3)
    Lambda = omega_lambda2lambda(Omega_Lambda)
    rs = 2*M
    first_term = (rs/2/r0)*np.cos(phi)*(-4*(np.cos(phi))**2 - 12*np.cos(phi)*np.sin(phi)*np.sqrt(Lambda*r0**2/3+rs/r0*(np.sin(phi))**3) + Lambda*r0**2*(8/3-20/3*(np.sin(phi))**2))
    second_term = (rs/2/r0)**2*(15/4*(2*phi-np.pi) + np.cos(phi)*(4+33/2*np.sin(phi)-4*(np.sin(phi))**2+19*(np.sin(phi))**3-64*(np.sin(phi))**5) - 12*np.log(np.tan(phi/2))*(np.sin(phi))**3)
    return first_term + second_term

def plot_alphas(filename, plot_ishak=True, plot_kantowski=True):
    df = pd.read_csv(filename)

    preds_frw = []
    preds_ishak = []
    preds_kantowski = []
    kant_higher_order_ratio = []

    for index, row in df.iterrows():
        M = row.M
        Lambda = 3*row.om_lambdas*H_0**2
        rho = (1-row.om_lambdas)*3*H_0**2/(8*np.pi)
        rh = (3*M/(4*np.pi*rho))**(1./3)

        R = (row.DL*row.theta)
        A_frw = 4*M/R + 15*np.pi*M**2/4/R**2 + 401/12*M**3/R**3
        # frw = row.comoving_lens/(A_frw/row.theta -1)

        A_ishak = 4*M/R + 15*np.pi*M**2/4/R**2 + 305/12*M**3/R**3 - Lambda*R*row.exit_rhs/3
        # ishak = row.comoving_lens/(A_ishak/row.theta -1)

        A_kantowski = -kantowski_alpha(R, M, row.enter_phis, row.om_lambdas)

        r0 = 1/(1/R + M/R**2)
        extra_term_ratio = np.abs((2*M/r0 + omega_lambda2lambda(row.om_lambdas)*r0**2)**(5/2)/(4*M/r0*(np.cos(row.enter_phis))**3))
        

        preds_frw.append(A_frw)
        preds_ishak.append(A_ishak)
        preds_kantowski.append(A_kantowski)
        kant_higher_order_ratio.append(extra_term_ratio)

    df['preds_frw'] = preds_frw
    df['preds_ishak'] = preds_ishak
    df['preds_kantowski'] = preds_kantowski
    df['kant_higher_order_ratio'] = kant_higher_order_ratio

    df['numerical'] = df.alphas/df.preds_frw - 1
    df['ishak'] = df.preds_ishak/df.preds_frw - 1
    df['kantowski'] = df.preds_kantowski/df.preds_frw - 1
    df['numerical_kantowski'] = df.alphas / df.preds_kantowski - 1

    stats = df[['om_lambdas', 'numerical', 'ishak', 'kantowski', 'numerical_kantowski', 'kant_higher_order_ratio']].groupby('om_lambdas').agg(['mean', 'std', 'count'])
    stats.columns = [' '.join(col).strip() for col in stats.columns.values]
    stats['numerical mean std'] = stats['numerical std']/np.sqrt(stats['numerical count'])
    # stats['numerical first order mean std'] = stats['numerical first order std']/np.sqrt(stats['numerical first order count'])
    stats['ishak mean std'] = stats['ishak std']/np.sqrt(stats['ishak count'])

    # stats['numerical mean'] = stats['numerical mean'] - 1
    # stats['ishak mean'] = stats['ishak mean'] - 1
    # stats['kantowski mean'] = stats['kantowski mean'] - 1

    scale = 1
    plt.plot(stats.index, stats['numerical mean']/scale, '.', label='__nolegend__')
    # plt.plot(stats.index, stats['numerical first order mean']/scale, 'b-', label='with second order corrections')
    plt.errorbar(stats.index, stats['numerical mean']/scale, yerr=stats['numerical mean std']/scale, linestyle='none', label='__nolegend__')
    plt.xlabel('Omega_Lambda')
    plt.ylabel('Mean fractional deviation/10^-6')
    if plot_ishak:
        plt.plot(stats.index, stats['ishak mean']/scale, 'g-', label='Rindler and Ishak predictions')
    if plot_kantowski:
        plt.plot(stats.index, stats['kantowski mean']/scale, 'b-', label="Kantowski predictions")
    plt.plot(stats.index, [0/scale]*len(stats.index), 'r-', label='FRW predictions')
    # plt.ylim((-0.0008, 0.0008))
    if plot_ishak or plot_kantowski:
        plt.legend()

    plt.figure()
    plt.plot(stats.index, stats['numerical_kantowski mean']/scale, '.', label='numerical to kantowski')
    if plot_kantowski:
        plt.plot(stats.index, stats['kant_higher_order_ratio mean']/scale, label='neglected term')
    plt.grid()
    plt.legend()

def plot_diff_lambdas_distances(filename, plot_ishak=True):
    df = pd.read_csv(filename)

    preds_frw = []
    preds_ishak = []
    for index, row in df.iterrows():
        M = row.M
        Lambda = 3*row.om_lambdas*H_0**2
        rho = (1-row.om_lambdas)*3*H_0**2/(8*np.pi)
        rh = (3*M/(4*np.pi*rho))**(1./3)

        try:
            exit_rh = row.exit_rhs
        except:
            exit_rh = rh
        A_frw = 4*M/(row.DL*row.theta) + 15*np.pi*M**2/4/(row.DL*row.theta)**2 + 401/12*M**3/(row.DL*row.theta)**3
        frw = row.comoving_lens/(A_frw/row.theta -1)
        
        A_ishak = 4*M/(row.DL*row.theta) + 15*np.pi*M**2/4/(row.DL*row.theta)**2 + 305/12*M**3/(row.DL*row.theta)**3 - Lambda*row.DL*row.theta*exit_rh/3
        ishak = row.comoving_lens/(A_ishak/row.theta -1)
        
        preds_frw.append(frw)
        preds_ishak.append(ishak)

    df['preds_frw'] = preds_frw
    df['preds_ishak'] = preds_ishak

    df['numerical'] = df.raw_rs/df.preds_frw
    df['ishak'] = df.preds_ishak/df.preds_frw

    stats = df[['om_lambdas', 'numerical', 'ishak']].groupby('om_lambdas').agg(['mean', 'std', 'count'])
    stats.columns = [' '.join(col).strip() for col in stats.columns.values]
    stats['numerical mean std'] = stats['numerical std']/np.sqrt(stats['numerical count'])
    stats['ishak mean std'] = stats['ishak std']/np.sqrt(stats['ishak count'])
    stats.head()

    scale = 1
    plt.plot(stats.index, stats['numerical mean']/scale, '.')
    plt.errorbar(stats.index, stats['numerical mean']/scale, yerr=stats['numerical mean std']/scale, linestyle='none')
    plt.xlabel('Omega_Lambda')
    # plt.ylabel('Mean fractional deviation/10^-5')
    if plot_ishak:
        plt.plot(stats.index, stats['ishak mean']/scale, 'g-')
    plt.plot(stats.index, [1/scale]*len(stats.index), 'r-')

    # chisquared_frw = (stats['numerical mean'] - 1)**2/(stats['numerical mean std'])**2
    # chisquared_ishak = (stats['numerical mean'] - stats['ishak mean'])**2/(stats['numerical mean std'])**2
    # p_frw = np.exp(-chisquared_frw.values.sum()/2)
    # p_ishak = np.exp(-chisquared_ishak.values.sum()/2)
    # print("p_frw/p_ishak: ", p_frw/p_ishak)
    # print("p_ishak/p_frw: ", p_ishak/p_frw)
    # # print(chisquared_frw)
    # # print(chisquared_ishak)

def plot_diff_lambdas(filename, recalculate_distances=False, scale=1e-5, plot_rindler=False):
    df = pd.read_csv(filename)

    # patch
    if recalculate_distances:
        length_scale = 3.086e22
        H_0 = 7.56e-27 * length_scale
        M = 1474e12 / length_scale

        z_lens1 = np.linspace(0.05, 0.2, 100)
        z_lens2 = []
        for z in z_lens1:
            z_lens2.extend([z]*50)

        # print(len(z_lens2), len(df.index))
        df['z_lens'] = z_lens2[:len(df.index)]

        def get_distances(z, Omega_Lambda=0):
            Omega_m = 1 - Omega_Lambda
            def integrand(z):
                return 1/np.sqrt(Omega_m*(1+z)**3 + Omega_Lambda)
            integral, error = spi.quad(integrand, 0, z)
            comoving = integral/H_0
            dang = comoving/(1+z)
            return comoving, dang

        dang_lens = []
        for index, row in df.iterrows():
            com, dang = get_distances(row.z_lens, row.om_lambdas)
            dang_lens.append(dang)

        df['DL'] = dang_lens

    M = 1474e12 / length_scale


    theta_second_order = []
    theta_rindler = []
    for index, row in df.iterrows():
        Lambda = 3*row.om_lambdas*H_0**2
        
        rho = (1-row.om_lambdas)*3*H_0**2/(8*np.pi)
        r_h = (3*M/(4*np.pi*rho))**(1./3)
        
        coeff = [row.DS + Lambda*row.DL*row.DLS*r_h/3, 0, -4*M*row.DLS/row.DL, -15*np.pi*M**2/4*row.DLS/row.DL**2, -305/12*M**3*row.DLS/row.DL**3]
        roots = np.roots(coeff)
        roots = roots[roots>0 & np.isreal(roots)]
        th = np.real(roots)
        rindler = th[np.argmin(np.abs(row.theta - th))]

        # coeff2 = [row.DS, -4*M*row.DLS/row.DL, 8*M**3*row.DLS/row.DL**3]
        coeff2 = [row.DS, 0, -4*M*row.DLS/row.DL, -15*np.pi*M**2/4*row.DLS/row.DL**2, -401/12*M**3*row.DLS/row.DL**3]
        roots2 = np.roots(coeff2)
        roots2 = roots2[roots2>0 & np.isreal(roots2)]
        th2 = np.real(roots2)
        second_order = th2[np.argmin(np.abs(row.theta - th2))]
        theta_rindler.append(rindler)
        theta_second_order.append(second_order)

    df['theta_second_order'] = theta_second_order
    df['theta_rindler'] = theta_rindler


    ## removed percentage!!

    df['percentage_diff'] = (df.theta_second_order - df.theta)/df.theta
    df['rindler_preds'] = (df.theta_rindler - df.theta)/df.theta
    # df['percentage_diff'] = (df.theta_second_order - df.theta)/df.theta*100
    # df['percentage_diff'] = (df.rs - df.rs_initial)/df.rs_initial*100

    stats = df[['om_lambdas', 'percentage_diff', 'rindler_preds']].groupby('om_lambdas').agg(['mean', 'std', 'count'])
    stats.columns = [' '.join(col).strip() for col in stats.columns.values]
    stats['percentage_diff mean std'] = stats['percentage_diff std']/np.sqrt(stats['percentage_diff count'])
    stats['rindler_preds mean std'] = stats['rindler_preds std']/np.sqrt(stats['rindler_preds count'])


    plt.plot(stats.index, stats['percentage_diff mean']/scale, '.', label='Numerical results')
    if plot_rindler:
        plt.plot(stats.index, stats['rindler_preds mean']/scale, 'g-', label='Results predicted by Rindler and Ishak')
        # plt.errorbar(stats.index, stats['rindler_preds mean']/scale, yerr=stats['rindler_preds mean std']/scale, linestyle='none')
    plt.errorbar(stats.index, stats['percentage_diff mean']/scale, yerr=stats['percentage_diff mean std']/scale, linestyle='none', label='_nolegend_')
    plt.xlabel('Omega_Lambda')
    plt.ylabel('Mean deviation/10^-5')


class LTB:
    M = 4.776409591704472e-08
    rho_frw_initial = 6.49705928222e-09
    rlimit_ratio = 1/0.7

    def initial_rh(self):
        return (3*self.M/(4*np.pi*self.rho_frw_initial))**(1./3)

    def rs(self):
        initial_rh = (3*self.M/(4*np.pi*self.rho_frw_initial))**(1./3)
        return np.linspace(0., initial_rh, 1000)

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

if __name__ == "__main__":
    ltb = LTB()
    rs = ltb.rs()[1:]
    rhos = [ltb.rho(r) for r in rs]
    masses = [ltb.mass(r) for r in rs]
    plt.semilogx(rs, masses)

    # def f(P, r):
    #     current_rho = ltb.rho(r)
    #     current_m = ltb.mass(r)
    #     Lambda = 0
    #     E = -2*current_m/r - Lambda*r**2/3
    #     return (current_rho + P)/2/(1+E)/r*(Lambda*r**2 - 8*np.pi*P*r**2 + E)
    # from scipy.integrate import odeint
    # pressure_rs = np.linspace(0.1, ltb.initial_rh(), 1000)[1:][::-1]
    # pressure = odeint(f, 0, pressure_rs)

    plt.show()
