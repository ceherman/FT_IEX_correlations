def LHS_1(psi_ar, solution):
    # Charge density [C m-2] of the diffuse layer at surface potential psi
    # return solution.eps*eps0*solution.kappa*psi_ar[0]
    return 2*solution.eps*eps0*solution.kappa*kT/e*numpy.sinh(e*psi_ar[0]/(2*kT))

def RHS_1(psi_ar, solution, charged_obj):
    # Charge density [C m-2] of the inner layer at surface potential psi
    total = 0
    c = solution.ch
    for element in charged_obj.log:
        K         = 10**(-1.0*charged_obj.log[element][0])
        charge    = charged_obj.log[element][1]
        surf_dens = charged_obj.log[element][2]
        add = surf_dens*((charge-1)/(1+(c/K)*numpy.exp(-1.0*e*psi_ar[0]/kT)) + charge/(1+(K/c)*numpy.exp(e*psi_ar[0]/kT)))
        total = total + add

    extra = 0
    if type(charged_obj)==protein:
        extra = e*charged_obj.extra_charge/(4.0*numpy.pi*(charged_obj.radius**2))
    return (e*Na*total + extra)

def residual_1(psi_ar, solution, charged_obj):
    return (abs(LHS_1(psi_ar, solution) - RHS_1(psi_ar, solution, charged_obj)))

def get_psi0(solution, charged_obj):
    fit = scipy.optimize.minimize(residual_1, x0=[0], args=(sol, charged_obj), bounds=[(-1, 1)], tol=1e-12)
    # fit = scipy.optimize.minimize_scalar(residual_1, args=(sol, charged_obj), bounds=[(-1, 1)], tol=1e-8)
    # check_res = residual_1([fit.x[0]], solution, charged_obj)
    # if check_res > 1e-5:
    #     raise Exception('psi0 issue; residual = {:.2e}'.format(check_res))
    charged_obj.psi0 = fit.x[0]
    return

def get_cap_in(solution, charged_obj):
    total = 0
    c = solution.ch
    psi = charged_obj.psi0
    for element in charged_obj.log:
        K         = 10**(-1.0*charged_obj.log[element][0])
        charge    = charged_obj.log[element][1]
        surf_dens = charged_obj.log[element][2]
        add = -1*surf_dens*(((charge-1)*(c/K)*(e/kT)*numpy.exp(-1.0*e*psi/kT))/((1+(c/K)*numpy.exp(-1.0*e*psi/kT))**2)
                         - (charge*(K/c)*(e/kT)*numpy.exp(e*psi/kT))/((1+(K/c)*numpy.exp(e*psi/kT))**2))
        total = total + add
    charged_obj.cap_in = e*Na*total
    return

def get_p(solution, charged_obj):
    charged_obj.p = solution.cap_dif/(solution.cap_dif + charged_obj.cap_in)
    return

def get_system_params(solution, protein, resin):
    solution.get_ch()
    solution.get_kappa()
    solution.get_cap_dif()

    get_psi0(solution, protein)
    get_cap_in(solution, protein)
    get_p(solution, protein)

    get_psi0(solution, resin)
    get_cap_in(solution, resin)
    get_p(solution, resin)
    return

def reset_system(solution, protein, resin):
    solution.pH = 7.0
    solution.ion_str = 0.1
    protein.radius = protein.initial_radius
    protein.log_pK_surf_dens()
    get_system_params(solution, protein, resin)
    return

def get_ener(z, solution, protein, resin):
#     z is the surface separation
    kappa = solution.kappa
    a1 = protein.radius

    psi1 = protein.psi0
    psi2 = resin.psi0
    p1 = protein.p
    p2 = resin.p
    g = (1-2*p1)*(1-2*p2)

    if g < 0:
        h = numpy.arctan(numpy.sqrt(-1.0*g)*numpy.exp(-1*kappa*z))
    else:
        h = numpy.arctanh(numpy.sqrt(g)*numpy.exp(-1*kappa*z))

    if h == numpy.nan or h == numpy.inf:
        raise Exception('NaN issue. g = {.:2f}, h = {.:2f}'.format(g, h))
    elif 1-g*numpy.exp(-2*kappa*z) < 0:
        raise Exception('g = {:.2f}, kappa = {:.2e}, z = {:.2e} \n p_prot = {:.2f} p_res = {:.2f}'.format(g, kappa, z, p1, p2))
    elif g == 0.0:
        raise Exception('g = 0 issue')
    elif z < 0:
        raise Exception('z < 0 issue. z = {:.2e}'.format(z))
    else:
        ener = 2*numpy.pi*solution.eps*eps0*a1*(((1-2*p1)*psi2**2 + (1-2*p2)*psi1**2)/(2*g)
                                                *numpy.log(1-g*numpy.exp(-2*kappa*z))
                                                +2*psi1*psi2/numpy.sqrt(numpy.abs(g))*h)

#     print('Whole expression = {:.2e}'.format(((1-2*p1)*psi2**2 + (1-2*p2)*psi1**2)/(2*g)*numpy.log(1-g*numpy.exp(-2*kappa*z))))
#     print('First part = {:.2e}'.format(((1-2*p1)*psi2**2 + (1-2*p2)*psi1**2)/(2*g)))
#     print('Log = {:.2e}'.format(numpy.log(1-g*numpy.exp(-2*kappa*z))))
    return ener

def find_ener_min(solution, protein, resin, x_up_bnd=1e-8):
    # Find the energy minimum
    res = scipy.optimize.minimize_scalar(get_ener, args=(solution, protein, resin), method='bounded', bracket=(4.0e-16, x_up_bnd), \
                                         bounds=(4.0e-16, x_up_bnd), tol=5e-16)
    return res.fun

def find_ener_min_x(solution, protein, resin, x_up_bnd=1e-8):
    # Find the energy minimum
    res = scipy.optimize.minimize_scalar(get_ener, args=(solution, protein, resin), method='bounded', bracket=(4.0e-16, x_up_bnd), \
                                         bounds=(4.0e-16, x_up_bnd), tol=5e-16)
    return res.x

def get_Keq(k, u_min):
    return k*kT/u_min*(1.0-numpy.exp(-1.0*u_min/kT))

def get_Keq_obj(k, solution, protein, resin):
    u_min = find_ener_min(solution, protein, resin)
    Keq = get_Keq(k, u_min)
    return Keq

def get_integrand(z, solution, protein, resin):
    return numpy.exp(-1.0*get_ener(z, solution, protein, resin)/kT) - 1

def get_Keq_integrate(solution, protein, resin, x_up_bnd=1e-8):
    u_min_x = find_ener_min_x(solution, protein, resin)
    res = integrate.quad(get_integrand, u_min_x,  x_up_bnd, args=(solution, protein, resin))
    return res[0]



def get_df_master(files):
    df_master = pandas.DataFrame()
    if type(files) == str:
        df_exp = pandas.read_csv(files)
        df_temp = pandas.DataFrame(numpy.nan, index=range(len(df_exp['Keq'])), columns=['pred_Keq', 'residual', 'log_residual'])
        df_res = pandas.concat([df_exp, df_temp], axis=1)
        df_master = df_res
    else:
        for file in files:
            df_exp = pandas.read_csv(file)
            df_temp = pandas.DataFrame(numpy.nan, index=range(len(df_exp['Keq'])), columns=['pred_Keq', 'residual', 'log_residual'])
            df_res = pandas.concat([df_exp, df_temp], axis=1)
            if file == files[0]:
                df_master = df_res
            else:
                df_master = pandas.concat([df_master, df_res], axis=0, ignore_index=True)
    return df_master

def get_res_df(in_array, df_res, solution, protein, resin):
    # in_array is [k (final model parameter), protein radius]

    protein.radius = in_array[1]
    protein.log_pK_surf_dens() #

    for index in range(len(df_res['Keq'])):
        sol.ion_str = df_res.at[index, 'IS(M)']
        sol.pH = df_res.at[index, 'pH']
        get_system_params(solution, protein, resin)

        u_min = find_ener_min(solution, protein, resin)
        df_res.at[index, 'pred_Keq'] = get_Keq(in_array[0], u_min)
        df_res.at[index, 'residual'] = (df_res.at[index, 'Keq'] - df_res.at[index, 'pred_Keq'])
        df_res.at[index, 'log_residual'] = numpy.log10(df_res.at[index, 'Keq']/df_res.at[index, 'pred_Keq'])

    return df_res

def get_res_vec(in_array, files, solution, protein, resin):
    df_master = get_df_master(files)
    df_res = get_res_df(in_array, df_master, solution, protein, resin)
    residual = []
    for element in df_res['residual']:
        residual.append(element)
    return residual

def get_sum_sq(in_array, files, solution, protein, resin):
    residual = get_res_vec(in_array, files, solution, protein, resin)
    sum_sq = 0
    for element in residual:
        sum_sq += element**2
    return sum_sq

def get_res_vec_log(in_array, files, solution, protein, resin):
    df_master = get_df_master(files)
    df_res = get_res_df(in_array, df_master, solution, protein, resin)
    residual = []
    for element in df_res['log_residual']:
        residual.append(element)
    return residual

def get_sum_sq_log(in_array, files, solution, protein, resin):
    residual = get_res_vec_log(in_array, files, solution, protein, resin)
    sum_sq = 0
    for element in residual:
        sum_sq += element**2
    return sum_sq



def get_res_df_integrate(in_array, df_res, solution, protein, resin, x_up_bnd=1e-8):
    # in_array is [protein radius]
    protein.radius = in_array[0]
    protein.log_pK_surf_dens() #

    for index in range(len(df_res['Keq'])):
        solution.ion_str = df_res.at[index, 'IS(M)']
        solution.pH = df_res.at[index, 'pH']
        get_system_params(solution, protein, resin)

        df_res.at[index, 'pred_Keq'] = get_Keq_integrate(solution, protein, resin, x_up_bnd)
        df_res.at[index, 'residual'] = (df_res.at[index, 'Keq'] - df_res.at[index, 'pred_Keq'])
        df_res.at[index, 'log_residual'] = numpy.log10(df_res.at[index, 'Keq']/df_res.at[index, 'pred_Keq'])
    return df_res

def get_res_vec_log_integrate(in_array, files, solution, protein, resin, x_up_bnd=1e-8):
    df_master = get_df_master(files)
    df_res = get_res_df_integrate(in_array, df_master, solution, protein, resin, x_up_bnd)
    residual = []
    for element in df_res['log_residual']:
        residual.append(element)
    return residual

def get_sum_sq_log_integrate(in_array, files, solution, protein, resin, x_up_bnd=1e-8):
    residual = get_res_vec_log_integrate(in_array, files, solution, protein, resin, x_up_bnd)
    sum_sq = 0
    for element in residual:
        sum_sq += element**2
    return sum_sq
