# Refers to the files in the charge_regulation_Hubbuch folder

def m_pK_eff(log, solution):
    a_debye = 0.5114
    B       = 0.1
    ion_str = solution.ion_str
    for aa in log:
        pK = log[aa][0]
        z = log[aa][1]
        pK_eff = pK + 2*(z-1)*(a_debye*numpy.sqrt(ion_str)/(1+numpy.sqrt(ion_str))-B*ion_str)
        log[aa].append(pK_eff)
    return log

def m_theta_p0_residual(theta_p0, solution, protein):
    lhs = protein.m_charge_dens_dl*(1.0-(numpy.tanh(theta_p0/4.0))**2)
    rhs = theta_p0*(1.0-(numpy.tanh(theta_p0/4.0))**2) + 4.0*solution.kappa*protein.radius*numpy.tanh(theta_p0/4.0)
    residual = abs(lhs - rhs)
    return residual

def m_theta_r0_residual(theta_r0, solution, resin):
    lhs = resin.m_charge_dens_dl
    rhs = 4.0*numpy.tanh(theta_r0/4.0)/(1.0-(numpy.tanh(theta_r0/4.0))**2)
    residual = abs(lhs - rhs)
    return residual

def m_capital_theta(dpr, solution, charged_obj):
    theta_0 = charged_obj.theta_0

    capital_theta = 4.0*numpy.exp(solution.kappa*dpr/2.0)*numpy.arctanh(numpy.exp(-1.0*solution.kappa*dpr/2.0)*numpy.tanh(theta_0/4.0))
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

    constant = Na*eps0*solution.eps*(kT**2)*numpy.pi*protein.radius/((solution.m_z*e)**2)
    if fp < 0.0:
        variable_term = numpy.arctan(numpy.sqrt(abs(fp))*numpy.exp(-kappa*dpr))
    else:
        variable_term = numpy.arctanh(numpy.sqrt(fp)*numpy.exp(-kappa*dpr))
    bracket = ((cap_theta_p**2 + fp*cap_theta_r**2)/fp)*numpy.log(1-fp*numpy.exp(-2.0*kappa*dpr)) + 4*cap_theta_p*cap_theta_r/numpy.sqrt(abs(fp))*variable_term
    # NB:  Negative sign on first term removed
    # {This contradicts Guelat, Khalaf, and Hsu & Liu 2009; appears consistent with Hsu & Liu 1999, Carnie & Chan 1993}

    # print(constant)
    # print(((cap_theta_p**2 + fp*cap_theta_r**2)/fp)*numpy.log(1-fp*numpy.exp(-2.0*kappa*dpr)))
    # print((cap_theta_p**2 + fp*cap_theta_r**2)/fp)
    # print(numpy.log(1-fp*numpy.exp(-2.0*kappa*dpr)))
    # print(4*cap_theta_p*cap_theta_r/numpy.sqrt(abs(fp))*variable_term)
    return constant*bracket

def m_delta_g_vdw(dpr, apr, protein):
    # VdW delta G [J mol-1]; apr [J]
    radius = protein.radius
    return -1.0*apr*Na/6.0*(radius/dpr + radius/(2*radius+dpr) + numpy.log(dpr/(2*radius+dpr)))

def m_delta_g(dpr, apr, solution, protein, resin):
    # [J mol-1]
    g_elec = m_delta_g_elec(dpr, solution, protein, resin)
    g_vdw  = m_delta_g_vdw(dpr, apr, protein)
    return g_elec + g_vdw

def m_get_Keq(dpr, apr, solution, protein, resin):
    g = m_delta_g(dpr, apr, solution, protein, resin)
    mass = protein.mass*kg_per_Da
    return h/numpy.sqrt(2.0*numpy.pi*mass*kT)*(numpy.exp(-1.0*g/(Na*kT))-1.0)



def m_get_res_df(in_array, df_res, apr, solution, protein, resin):
    # in_array is [protein radius, dpr]
    protein.radius = in_array[0]
    dpr = in_array[1]

    for index in range(len(df_res['Keq'])):
        solution.ion_str = df_res.at[index, 'IS(M)']
        solution.pH = df_res.at[index, 'pH']
        solution.get_ch()
        solution.get_kappa()

        df_res.at[index, 'pred_Keq'] = m_get_Keq(dpr, apr, solution, protein, resin)
        df_res.at[index, 'residual'] = (df_res.at[index, 'Keq'] - df_res.at[index, 'pred_Keq'])
        df_res.at[index, 'log_residual'] = numpy.log10(df_res.at[index, 'Keq']/df_res.at[index, 'pred_Keq'])
    return df_res

def m_get_res_vec_log(in_array, files, apr, solution, protein, resin):
    df_master = get_df_master(files)
    df_res = m_get_res_df(in_array, df_master, apr, solution, protein, resin)
    residual = []
    for element in df_res['log_residual']:
        residual.append(element)
    return residual

def m_get_sum_sq_log(in_array, files, apr, solution, protein, resin):
    residual = m_get_res_vec_log(in_array, files, apr, solution, protein, resin)
    sum_sq = 0
    for element in residual:
        sum_sq += element**2
    return sum_sq

def m2_get_res_df(in_array, df_res, apr, solution, protein, resin):
    # in_array is [dpr]
    dpr = in_array[0]

    # # in_array is [radius]
    # dpr = 1.0e-15
    # protein.radius = in_array[0]

    for index in range(len(df_res['Keq'])):
        solution.ion_str = df_res.at[index, 'IS(M)']
        solution.pH = df_res.at[index, 'pH']
        solution.get_ch()
        solution.get_kappa()

        df_res.at[index, 'pred_Keq'] = m_get_Keq(dpr, apr, solution, protein, resin)
        df_res.at[index, 'residual'] = (df_res.at[index, 'Keq'] - df_res.at[index, 'pred_Keq'])
        df_res.at[index, 'log_residual'] = numpy.log10(df_res.at[index, 'Keq']/df_res.at[index, 'pred_Keq'])
    return df_res

def m2_get_res_vec_log(in_array, files, apr, solution, protein, resin):
    df_master = get_df_master(files)
    df_res = m2_get_res_df(in_array, df_master, apr, solution, protein, resin)
    residual = []
    for element in df_res['log_residual']:
        residual.append(element)
    return residual

def m2_get_sum_sq_log(in_array, files, apr, solution, protein, resin):
    residual = m2_get_res_vec_log(in_array, files, apr, solution, protein, resin)
    sum_sq = 0
    for element in residual:
        sum_sq += element**2
    return sum_sq
