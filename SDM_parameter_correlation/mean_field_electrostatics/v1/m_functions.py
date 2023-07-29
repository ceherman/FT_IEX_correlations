import charged_objects as obj
import numpy as np

def m_pK_eff(log, solution):
    a_debye = 0.5114
    B       = 0.1
    ion_str = solution.ion_str
    for aa in log:
        pK = log[aa][0]
        z = log[aa][1]
        pK_eff = pK + 2*(z-1)*(a_debye*np.sqrt(ion_str)/(1+np.sqrt(ion_str))-B*ion_str)
        log[aa].append(pK_eff)
    return log

def m_theta_p0_residual(theta_p0, solution, protein):
    lhs = protein.m_charge_dens_dl*(1.0-(np.tanh(theta_p0/4.0))**2)
    rhs = theta_p0*(1.0-(np.tanh(theta_p0/4.0))**2) + 4.0*solution.kappa*protein.radius*np.tanh(theta_p0/4.0)
    residual = abs(lhs - rhs)
    return residual

def m_theta_r0_residual(theta_r0, solution, resin):
    lhs = resin.m_charge_dens_dl
    rhs = 4.0*np.tanh(theta_r0/4.0)/(1.0-(np.tanh(theta_r0/4.0))**2)
    residual = abs(lhs - rhs)
    return residual

def m_capital_theta(dpr, solution, charged_obj):
    theta_0 = charged_obj.theta_0

    capital_theta = 4.0*np.exp(solution.kappa*dpr/2.0)*np.arctanh(np.exp(-1.0*solution.kappa*dpr/2.0)*np.tanh(theta_0/4.0))
    return capital_theta

def m_delta_g_elec(dpr, solution, protein, resin):
    # Electorstatic delta G [J mol-1]
    protein.m_get_capital_theta(solution, dpr)
    protein.m_fp(solution)
    resin.m_get_capital_theta(solution, dpr)

    cap_theta_p = protein.capital_theta
    fp          = protein.fp
    cap_theta_r = resin.capital_theta
    kappa = solution.kappa

    constant = obj.constants().Na*obj.constants().eps0*solution.eps*(obj.constants().kT**2)*np.pi*protein.radius/((solution.m_z*obj.constants().e)**2)
    if fp < 0.0:
        variable_term = np.arctan(np.sqrt(abs(fp))*np.exp(-kappa*dpr))
    else:
        variable_term = np.arctanh(np.sqrt(fp)*np.exp(-kappa*dpr))
    bracket = (-1.0*(cap_theta_p**2 + fp*cap_theta_r**2)/fp)*np.log(1-fp*np.exp(-2.0*kappa*dpr)) + 4*cap_theta_p*cap_theta_r/np.sqrt(abs(fp))*variable_term
    # NB:  Negative sign should either be on the entire first term, or it should be removed completely
    # {Negative sign:  Guelat, Khalaf, and Hsu & Liu 2009; positive sign:  Hsu & Liu 1999, and possibly also Carnie & Chan 1993}

    return constant*bracket

def m_delta_g_vdw(dpr, apr, protein):
    # VdW delta G [J mol-1]; apr [J]
    radius = protein.radius
    return -1.0*apr*obj.constants().Na/6.0*(radius/dpr + radius/(2*radius+dpr) + np.log(dpr/(2*radius+dpr)))

def m_delta_g(dpr, apr, solution, protein, resin):
    # [J mol-1]
    g_elec = m_delta_g_elec(dpr, solution, protein, resin)
    g_vdw  = m_delta_g_vdw(dpr, apr, protein)
    return g_elec + g_vdw

def m_get_Keq_dim(dpr, apr, solution, protein, resin):
    g = m_delta_g(dpr, apr, solution, protein, resin)
    mass = protein.mass*obj.constants().kg_per_Da
    return obj.constants().h/np.sqrt(2.0*np.pi*mass*obj.constants().kT)*(np.exp(-1.0*g/(obj.constants().Na*obj.constants().kT))-1.0)

def m_get_Keq_nd(surf_area, porosity, dpr, apr, solution, protein, resin):
    keq_dim = m_get_Keq_dim(dpr, apr, solution, protein, resin)
    return surf_area/(1.0 - porosity) * keq_dim

def m_get_kprime(surf_area, porosity, dpr, apr, solution, protein, resin):
    keq_nd = m_get_Keq_nd(surf_area, porosity, dpr, apr, solution, protein, resin)
    return (1.0 - porosity)/porosity * keq_nd
#
