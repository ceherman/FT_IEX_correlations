import m_functions as m_fun
import numpy as np
from scipy import optimize

class constants:
    """Store physical constants for function access."""
    def __init__(self):
        self.e         = 1.602e-19      # [C]
        self.Na        = 6.022e23       # [mol-1]
        self.kT        = 4.11e-21       # [J]
        self.eps0      = 8.854e-12      # [C V-1 m-1]
        self.h         = 6.626e-34      # [J s]
        self.F         = 96485          # [C mol-1]
        self.kg_per_Da = 1.6605e-27     # [kg Da-1]
        self.prot_vbar = 0.73           # [ml g-1] protein partial specific volume
        self.pK        = {'K':          [10.59, 1], \
                          'R':          [12.00, 1], \
                          'H':          [ 6.55, 1], \
                          'E':          [ 4.18, 0], \
                          'D':          [ 3.53, 0], \
                          'Y':          [10.26, 0], \
                          'N_term':     [ 7.66, 1], \
                          'C_term':     [ 3.30, 0], \
                          'sulfonate':  [ 2.30, 0], \
                          'neg_ext':    [-100 , 0], \
                          'quat_amine': [20.00, 1]}
        return

class solution:
    def __init__(self, pH, ion_str, eps=80.1, m_z=1):
        self.pH      = pH
        self.ion_str = ion_str  # [M]
        self.eps     = eps      # [-]
        self.m_z     = m_z      # (Symmetrical) electorlyte valence
        self.get_ch()
        self.get_kappa()
        self.get_cap_dif()
        return

    def get_ch(self):
        # Convert pH to bulk solution hydronium ion concentration [M]
        self.ch = 10.0**(-1.0*self.pH)
        return

    def get_kappa(self):
        # Convert ionic strength [M] into the inverse Debye length [m-1]
        self.kappa = np.sqrt(2*constants().e**2*self.ion_str*constants().Na*1.0e3/(constants().kT*self.eps*constants().eps0))
        return

    def get_cap_dif(self):
        # Compute the diffuse layer capacitance [C m-2 V-1]
        self.cap_dif = self.eps*constants().eps0*self.kappa
        return


class resin:
    def __init__(self, name, ligand, surf_dens):
        self.name      = name
        self.ligand    = ligand
        self.surf_dens = surf_dens # [mol m-2]
        self.log_pK_surf_dens()
        return

    def log_pK_surf_dens(self):
        self.log = {}
        self.log[self.ligand] = constants().pK[self.ligand].copy()
        self.log[self.ligand].append(self.surf_dens)
        return

    def m_log_pKeff(self, solution):
        self.m_log_pK_eff = m_fun.m_pK_eff(self.log, solution)
        return

    def m_get_charge_dens(self, solution):
        self.m_log_pKeff(solution)
        ch = solution.ch
        self.m_charge_dens = 0.0

        # Constant charge density (original model)
        for ligand in self.m_log_pK_eff:
            pKa_eff = self.m_log_pK_eff[ligand][3]
            Ka_eff = 10**(-1.0*pKa_eff)

            if self.m_log_pK_eff[ligand][1] == 0:  # acid
                sign = -1.0
                denom = (1.0+ch/Ka_eff)
            elif self.m_log_pK_eff[ligand][1] == 1:  # base
                sign = 1.0
                denom = (1.0+Ka_eff/ch)
            self.m_charge_dens += sign*constants().F*self.surf_dens/denom       # Would need to change for multiple ligands

        # # Charge regulation
        # get_psi0(solution, self)
        # self.m_charge_dens = LHS_1([self.psi0], solution)

        return

    def m_get_charge_dens_dl(self, solution):
        self.m_get_charge_dens(solution)
        self.m_charge_dens_dl = solution.m_z*constants().e*self.m_charge_dens/(constants().kT*constants().eps0*solution.eps*solution.kappa)
        return

    def m_get_theta_r0(self, solution):
        self.m_get_charge_dens_dl(solution)
        fit = optimize.minimize(m_fun.m_theta_r0_residual, x0=[0], args=(solution, self), bounds=[(-1000, 1000)], tol=1e-12)
        self.theta_0 = fit.x[0]
        return

    def m_get_capital_theta(self, solution, dpr):
        self.m_get_theta_r0(solution)
        self.capital_theta = m_fun.m_capital_theta(dpr, solution, self)
        return


class protein:
    # I want to specify net charge and radius

    def __init__(self, name, mass, net_charge):
        self.name           = name
        self.mass           = mass          # [g mol-1]
        self.net_charge     = net_charge
        self.get_size()
        self.m_get_charge_dens()
        return

    def get_size(self):
        self.radius = (0.75/np.pi*self.mass*constants().prot_vbar/constants().Na)**(1/3)*1.0e-2 # [m]
        self.volume = 4.0/3.0 * np.pi * self.radius**3                                          # [m3]
        self.area   = 4.0 * np.pi * self.radius**2                                              # [m2]
        return

    def m_get_charge_dens(self):
        self.m_charge_dens = self.net_charge/self.area * constants().F/constants().Na   # [C m-2]
        return

    def m_get_charge_dens_dl(self, solution):
        self.m_charge_dens_dl = solution.m_z*constants().e*self.m_charge_dens*self.radius/(constants().kT*solution.eps*constants().eps0)
        return

    def m_get_theta_p0(self, solution):
        self.m_get_charge_dens_dl(solution)
        fit = optimize.minimize(m_fun.m_theta_p0_residual, x0=[0], args=(solution, self), bounds=[(-1000, 1000)], tol=1e-12)
        self.theta_0 = fit.x[0]
        # print(self.theta_0)
        return

    def m_get_capital_theta(self, solution, dpr):
        self.m_get_theta_p0(solution)
        self.capital_theta = m_fun.m_capital_theta(dpr, solution, self)
        return

    def m_fp(self, solution):
        kap_r = solution.kappa*self.radius
        num   = 1 - 1/kap_r + (1 + 1/kap_r)*np.exp(-2.0*kap_r)
        denom = 1 + 1/kap_r - (1 + 1/kap_r)*np.exp(-2.0*kap_r)     # 2016 Guelat:  1/kap_r is (+) {agrees with Hsu and Liu 2009};        2012 Guelat:  1/kap_r is (-)
        self.fp = num/denom
        return
