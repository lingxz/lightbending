import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.integrate as spi
# plt.style.use('seaborn-white')

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

    length_scale = 3.086e22
    H_0 = 7.56e-27 * length_scale
    M = 1474e12 / length_scale

    def calc_theta(D_LS, D_L, D_S):
        return np.sqrt(4*M*D_LS/D_L/D_S)

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
        # if r > rlimit:
        #     return self.M
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
    def f(P, r):
        current_rho = ltb.rho(r)
        current_m = ltb.mass(r)
        Lambda = 0
        E = -2*current_m/r - Lambda*r**2/3
        return (current_rho + P)/2/(1+E)/r*(Lambda*r**2 - 8*np.pi*P*r**2 + E)
    from scipy.integrate import odeint
    pressure_rs = np.linspace(0.1, ltb.initial_rh(), 1000)[1:][::-1]
    pressure = odeint(f, 0, pressure_rs)

