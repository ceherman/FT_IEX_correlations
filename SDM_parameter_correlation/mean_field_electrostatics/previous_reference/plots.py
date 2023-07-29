def plot_titration(solution, charged_obj):
    pH_list = numpy.linspace(0, 14, 90)
    charge_dens_list = []

    for pH in pH_list:
        solution.pH = pH
        solution.get_ch()
        solution.get_kappa()
        solution.get_cap_dif()

        get_psi0(solution, charged_obj)
        get_cap_in(solution, charged_obj)
        get_p(solution, charged_obj)

        charge_dens_list.append(LHS_1([charged_obj.psi0], solution))

    # fig, axs = plt.subplots()
    fig.set_size_inches(8, 7, forward=True)
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 20}
    plt.rc('font', **font)
    axs.plot(pH_list, charge_dens_list)
    axs.set_xlabel('pH')
    axs.set_ylabel(r'Charge density [$C/m^2$]')
    axs.axhline(color='black', linewidth=1)
    return

def plot_titration_ion_str(solution, charged_obj):
    is_list = numpy.linspace(0.05, 0.3, 6)
    pH_list = numpy.linspace(0, 14, 90)

    fig, axs = plt.subplots()
    fig.set_size_inches(8, 7, forward=True)
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 20}
    plt.rc('font', **font)
    axs.set_xlabel('pH')
    axs.set_ylabel(r'Charge density [$C/m^2$]')
    axs.axhline(color='black', linewidth=1)

    for ion_str in is_list:
        sol.ion_str = ion_str
        charge_dens_list = []

        for pH in pH_list:
            solution.pH = pH
            solution.get_ch()
            solution.get_kappa()
            solution.get_cap_dif()

            get_psi0(solution, charged_obj)
            get_cap_in(solution, charged_obj)
            get_p(solution, charged_obj)

            charge_dens_list.append(LHS_1([charged_obj.psi0], solution))

        axs.plot(pH_list, charge_dens_list, label="IS = {:.2f} M".format(ion_str))
    axs.legend(loc='upper right', bbox_to_anchor=(1.0, 0.98), frameon=False)
    fig.savefig('Seph_titration_PB.png', bbox_inches='tight', dpi = 300)
    return

def plot_titration_var_obj(solution, charged_obj, leg_label):
    pH_list = numpy.linspace(0, 14, 90)
    charge_dens_list = []

    for pH in pH_list:
        solution.pH = pH
        solution.get_ch()
        solution.get_kappa()
        solution.get_cap_dif()

        get_psi0(solution, charged_obj)
        get_cap_in(solution, charged_obj)
        get_p(solution, charged_obj)

        charge_dens_list.append(LHS_1([charged_obj.psi0], solution))

    axs.plot(pH_list, charge_dens_list, label=leg_label)
    axs.set_xlabel('pH')
    axs.set_ylabel(r'Charge density [$C/m^2$]')
    axs.legend(loc='best', frameon=False)
    axs.axhline(color='black', linewidth=1)
    axs.axvline(x=7, color='black', linewidth=1)
    axs.set_xlim(6, 8)
    axs.set_ylim(-.02, 0.02)
    return

def plot_charge_dens(solution, protein, resin):
#     get_system_params(solution, protein, resin)
    psi = numpy.linspace(-0.5, 0.5, 2001)
    df = pandas.DataFrame(numpy.nan, index=range(len(psi)), columns=['psi', 'lhs', 'rhs_prot', 'rhs_res'])

    for index in range(len(df['psi'])):
        df.at[index, 'psi'] = psi[index]
        df.at[index, 'lhs'] = LHS_1([df.at[index, 'psi']], solution)
        df.at[index, 'rhs_prot'] = RHS_1([df.at[index, 'psi']], solution, protein)
        df.at[index, 'rhs_res'] = RHS_1([df.at[index, 'psi']], solution, resin)

    fig, axs = plt.subplots()
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)

    axs.plot(df['psi']*1e3, df['lhs'], label='Diffuse layer')
    axs.plot(df['psi']*1e3, df['rhs_prot'], label=protein.name)
    axs.plot(df['psi']*1e3, df['rhs_res'], label=resin.name)
    axs.axhline(color='black', linewidth=1)
    axs.legend(loc='center left', bbox_to_anchor=(-0.02, 0.65), frameon=False)

    axs.set_ylim(-0.4, 0.4)
    # axs.set_xlim(psi[0]*1e3, psi[len(psi)-1]*1e3)
    axs.set_xlim(-300, 300)


    axs.set_xlabel('Surface potential [mV]')
    axs.set_ylabel(r'Charge density [$C/m^2$]')
    fig.savefig('Seph_charge_dens_DH.png', bbox_inches='tight', dpi = 300)
    return

def plot_charge_dens_ion_str(solution, protein, resin):
#     get_system_params(solution, protein, resin)
    psi = numpy.linspace(-0.3, 0.3, 2001)
    ion_str_list = numpy.linspace(0.01, 0.1, 10)
    df = pandas.DataFrame(numpy.nan, index=range(len(psi)), columns=['psi', 'lhs', 'rhs_prot', 'rhs_res'])

    fig, axs = plt.subplots()
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)

    for ion_str in ion_str_list:
        sol.ion_str = ion_str
        get_system_params(solution, protein, resin)

        for index in range(len(df['psi'])):
            df.at[index, 'psi'] = psi[index]
            df.at[index, 'lhs'] = LHS_1([df.at[index, 'psi']], solution)
            df.at[index, 'rhs_prot'] = RHS_1([df.at[index, 'psi']], solution, protein)
            df.at[index, 'rhs_res'] = RHS_1([df.at[index, 'psi']], solution, resin)

        axs.plot(df['psi']*1e3, df['lhs'], label='{:.2f} M'.format(ion_str))
    axs.plot(df['psi']*1e3, df['rhs_prot'], label='Protein')
    axs.plot(df['psi']*1e3, df['rhs_res'], label='Resin')

    axs.axhline(color='black', linewidth=1)
    axs.legend(loc='center right', bbox_to_anchor=(1.5, 0.5), ncol=1, frameon=False)
    axs.set_ylim(-0.5, 0.3)
    # axs.set_xlim(psi[0]*1e3, psi[len(psi)-1]*1e3)
    axs.set_xlim(-300, 300)

    axs.set_xlabel('Surface potential [mV]')
    axs.set_ylabel(r'Charge density [$C/m^2$]')
    return

def plot_charge_dens_mult_obj(solution, protein1, protein2, resin, leg_label1, leg_label2):
#     get_system_params(solution, protein, resin)
    psi = numpy.linspace(-0.3, 0.3, 2001)
    ion_str_list = numpy.linspace(0.05, 0.2, 4)
    df = pandas.DataFrame(numpy.nan, index=range(len(psi)), columns=['psi', 'lhs', 'rhs_prot', 'rhs_res'])

    for ion_str in ion_str_list:
        sol.ion_str = ion_str
        get_system_params(solution, protein1, resin)

        for index in range(len(df['psi'])):
            df.at[index, 'psi'] = psi[index]
            df.at[index, 'lhs'] = LHS_1([df.at[index, 'psi']], solution)
            df.at[index, 'rhs_prot'] = RHS_1([df.at[index, 'psi']], solution, protein1)
            df.at[index, 'rhs_res'] = RHS_1([df.at[index, 'psi']], solution, resin)

        axs.plot(df['psi']*1e3, df['lhs'], label='{:.2f} M'.format(ion_str))

    # axs.plot(df['psi']*1e3, df['rhs_prot'], label=leg_label1)
    axs.plot(df['psi']*1e3, df['rhs_prot'])

    get_system_params(solution, protein2, resin)
    for index in range(len(df['psi'])):
        df.at[index, 'rhs_prot'] = RHS_1([df.at[index, 'psi']], solution, protein2)

    # axs.plot(df['psi']*1e3, df['rhs_prot'], label=leg_label2)
    axs.plot(df['psi']*1e3, df['rhs_prot'])


    axs.axhline(color='black', linewidth=1)
    # axs.legend(loc='center right', bbox_to_anchor=(1.5, 0.5), ncol=1, frameon=False)
    axs.legend(loc='center left', bbox_to_anchor=(-0.02, 0.65), frameon=False)

    axs.set_ylim(-0.4, 0.4)
    # axs.set_xlim(psi[0]*1e3, psi[len(psi)-1]*1e3)
    axs.set_xlim(-300, 300)
    # axs.axvline(x=0, color='black', linewidth=1)

    axs.set_xlabel('Surface potential [mV]')
    axs.set_ylabel(r'Charge density [$C/m^2$]')
    return

def plot_charge_dens_mult_resin(solution, protein, resin1, resin2):
#     get_system_params(solution, protein, resin)
    psi = numpy.linspace(-0.4, 0.2, 2001)
    ion_str_list = numpy.linspace(0.05, 0.5, 10)
    df = pandas.DataFrame(numpy.nan, index=range(len(psi)), columns=['psi', 'lhs', 'rhs_prot', 'rhs_res1', 'rhs_res2'])

    fig, axs = plt.subplots()
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)

    for ion_str in ion_str_list:
        sol.ion_str = ion_str
        get_system_params(solution, protein, resin1)
        get_system_params(solution, protein, resin2)

        for index in range(len(df['psi'])):
            df.at[index, 'psi'] = psi[index]
            df.at[index, 'lhs'] = LHS_1([df.at[index, 'psi']], solution)
            df.at[index, 'rhs_prot'] = RHS_1([df.at[index, 'psi']], solution, protein)
            df.at[index, 'rhs_res1'] = RHS_1([df.at[index, 'psi']], solution, resin1)
            df.at[index, 'rhs_res2'] = RHS_1([df.at[index, 'psi']], solution, resin2)

        axs.plot(df['psi']*1e3, df['lhs'], label='{:.2f} M'.format(ion_str))

    axs.plot(df['psi']*1e3, df['rhs_prot'], label=protein.name)
    axs.plot(df['psi']*1e3, df['rhs_res1'], label=resin1.name)
    axs.plot(df['psi']*1e3, df['rhs_res2'], label=resin2.name)

    axs.axhline(color='black', linewidth=1)
    axs.legend(loc='center right', bbox_to_anchor=(1.7, 0.5), ncol=1, frameon=False)
    axs.set_ylim(-0.5, 0.3)
    # axs.set_xlim(psi[0]*1e3, psi[len(psi)-1]*1e3)
    axs.set_xlim(-400, 200)

    axs.set_xlabel('Surface potential [mV]')
    axs.set_ylabel(r'Charge density [$C/m^2$]')
    # fig.savefig('Seph_charge_dens_IS_DH.png', bbox_inches='tight', dpi = 300)
    return

def plot_charge_dens_mult_resin_pH(solution, protein, resin1, resin2):
#     get_system_params(solution, protein, resin)
    psi = numpy.linspace(-0.3, 0.3, 2001)
    pH_list = numpy.linspace(6, 8, 3)
    df = pandas.DataFrame(numpy.nan, index=range(len(psi)), columns=['psi', 'lhs', 'rhs_prot', 'rhs_res1', 'rhs_res2'])

    fig, axs = plt.subplots()
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)

    for pH in pH_list:
        sol.pH = pH
        get_system_params(solution, protein, resin1)
        get_system_params(solution, protein, resin2)

        for index in range(len(df['psi'])):
            df.at[index, 'psi'] = psi[index]
            df.at[index, 'lhs'] = LHS_1([df.at[index, 'psi']], solution)
            df.at[index, 'rhs_prot'] = RHS_1([df.at[index, 'psi']], solution, protein)
            df.at[index, 'rhs_res1'] = RHS_1([df.at[index, 'psi']], solution, resin1)
            df.at[index, 'rhs_res2'] = RHS_1([df.at[index, 'psi']], solution, resin2)

        axs.plot(df['psi']*1e3, df['rhs_prot'], label='pH = {:.0f}'.format(pH))
        axs.plot(df['psi']*1e3, df['rhs_res1'], '--', color=plt.gca().lines[-1].get_color())

    axs.plot(df['psi']*1e3, df['lhs'])
    # axs.plot(df['psi']*1e3, df['rhs_res2'], label=resin2.name)

    axs.axhline(color='black', linewidth=1)
    # axs.legend(loc='center right', bbox_to_anchor=(1.8, 0.5), ncol=1, frameon=False)
    axs.legend(loc='center left', bbox_to_anchor=(-0.02, 0.63), frameon=False)

    axs.set_ylim(-0.4, 0.4)
    # axs.set_xlim(psi[0]*1e3, psi[len(psi)-1]*1e3)
    axs.set_xlim(-300, 300)

    axs.set_xlabel('Surface potential [mV]')
    axs.set_ylabel(r'Charge density [$C/m^2$]')
    fig.savefig('Seph_pH_effect_DH.png', bbox_inches='tight', dpi = 300)
    return

def plot_surf_IS(solution, protein, resin1, resin2):
#     get_system_params(solution, protein, resin)
    pH_list = numpy.linspace(5, 9, 5)
    ion_str_list = numpy.linspace(0.01, 0.5, 200)

    fig, axs = plt.subplots()
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)
    ax2 = axs.twinx()

    for pH in pH_list:
        df = pandas.DataFrame(numpy.nan, index=range(len(ion_str_list)), columns=['ion_str', 'psi1', 'lhs1', 'psi2', 'lhs2'])
        sol.pH = pH
        for index in range(len(ion_str_list)):
            sol.ion_str = ion_str_list[index]
            get_system_params(solution, protein, resin1)
            get_system_params(solution, protein, resin2)

            df.at[index, 'ion_str'] = ion_str_list[index]
            df.at[index, 'psi1'] = resin1.psi0
            df.at[index, 'lhs1'] = LHS_1([df.at[index, 'psi1']], solution)
            df.at[index, 'psi2'] = resin2.psi0
            df.at[index, 'lhs2'] = LHS_1([df.at[index, 'psi2']], solution)

        axs.plot(df['ion_str'], df['psi1']*1e3, label=r'$\psi_0$, pH = {:.0f}'.format(sol.pH))
        # axs.plot(df['ion_str'], df['psi2']*1e3, label=resin2.name + ' psi, pH = {:.0f}'.format(sol.pH))

        ax2.plot(df['ion_str'], df['lhs1'], '--', label=r'$\sigma_0$, pH = {:.0f}'.format(sol.pH))
        # ax2.plot(df['ion_str'], df['lhs2'], '--', label=resin2.name + ' sigma, pH = {:.0f}'.format(sol.pH))

    fig.legend(loc='center right', bbox_to_anchor=(1.5, 0.5), ncol=1, frameon=False)
    axs.set_ylim(-400, -100)
    # axs.set_xlim(psi[0]*1e3, psi[len(psi)-1]*1e3)
    ax2.set_ylim(-0.5, -0.1)

    axs.set_xlabel('Ionic strength [M]')
    axs.set_ylabel(r'$\psi_0$ [mV]')
    ax2.set_ylabel(r'$\sigma_0$ [$C/m^2$]')
    # fig.savefig('Seph_psi_IS_PB.png', bbox_inches='tight', dpi = 300)
    return

def plot_surf_pH(solution, protein, resin1, resin2):
    pH_list = numpy.linspace(0, 14, 400)
    ion_str_list = numpy.linspace(0.1, 0.5, 5)

    fig, axs = plt.subplots()
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)
    ax2 = axs.twinx()

    for ion_str in ion_str_list:
        df = pandas.DataFrame(numpy.nan, index=range(len(ion_str_list)), columns=['ion_str', 'psi1', 'lhs1', 'psi2', 'lhs2'])
        sol.ion_str = ion_str
        for index in range(len(pH_list)):
            sol.pH = pH_list[index]
            get_system_params(solution, protein, resin1)
            get_system_params(solution, protein, resin2)

            df.at[index, 'pH'] = pH_list[index]
            df.at[index, 'psi1'] = resin1.psi0
            df.at[index, 'lhs1'] = LHS_1([df.at[index, 'psi1']], solution)
            df.at[index, 'psi2'] = resin2.psi0
            df.at[index, 'lhs2'] = LHS_1([df.at[index, 'psi2']], solution)

        axs.plot(df['pH'], df['psi1']*1e3, label=r'$\psi_0$, {:.2f} M'.format(sol.ion_str))
        ax2.plot(df['pH'], df['lhs1'], '--', label=r'$\sigma_0$, {:.2f} M'.format(sol.ion_str))

    fig.legend(loc='center right', bbox_to_anchor=(1.5, 0.5), ncol=1, frameon=False)
    axs.set_ylim(-600, 0)
    ax2.set_ylim(-0.5, 0)
    axs.set_xlim(0, 15)

    axs.set_xlabel('pH')
    axs.set_ylabel(r'$\psi_0$ [mV]')
    ax2.set_ylabel(r'$\sigma_0$ [$C/m^2$]')
    # fig.savefig('Seph_psi_pH_DH.png', bbox_inches='tight', dpi = 300)
    return


def plot_ener(solution, protein, resin, leg_label=''):
    get_system_params(solution, protein, resin)
    z = numpy.linspace(0, 6e-9, 1001)
    df_ener = pandas.DataFrame(numpy.nan, index=range(len(z)), columns=['z', 'ener'])

    for index in range(len(df_ener['z'])):
        df_ener.at[index, 'z']    = z[index]
        df_ener.at[index, 'ener'] = get_ener(df_ener.at[index, 'z'], solution, protein, resin)/kT

    axs.plot(df_ener['z']*1e9, df_ener['ener'], label=leg_label)

    axs.set_xlim(z[0]*1e9, z[len(z)-1]*1e9)
    axs.set_ylabel('Energy [kT]')
    axs.set_xlabel('Surface separation [nm]')
    axs.axhline(color='black', linewidth=1)

    if leg_label=='None':
        pass
    else:
        axs.legend(loc='center right', bbox_to_anchor=(1.5, 0.5), frameon=False)

    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)
    return

def plot_effect(fit_array, solution, protein, resin, plt_label, leg_cols):
    protein.radius = fit_array[1]
    protein.log_pK_surf_dens()

    ion_str = numpy.logspace(-2, 0, 100)
    df_pred = pandas.DataFrame(numpy.nan, index=range(len(ion_str)), columns=['IS(M)', 'pred_Keq'])

    for index in range(len(df_pred['pred_Keq'])):
        df_pred.at[index, 'IS(M)'] = ion_str[index]
        sol.ion_str = df_pred.at[index, 'IS(M)']
        get_system_params(solution, protein, resin)

        u_min = find_ener_min(solution, protein, resin)
        df_pred.at[index, 'pred_Keq'] = get_Keq(fit_array[0], u_min)

    axs.loglog(df_pred['IS(M)'], df_pred['pred_Keq'], label=plt_label)
    if leg_cols==1:
        axs.legend(loc='center right', bbox_to_anchor=(1.5, 0.5), ncol=leg_cols, frameon=False)
    elif leg_cols==2:
        axs.legend(loc='center right', bbox_to_anchor=(2, 0.5), ncol=leg_cols, frameon=False)
#     axs.set_ylim(1e-8, 1e-5)

    axs.set_xlabel('Ionic strength [M]')
    axs.set_ylabel(r'$K_{eq}$ [m]')
    return

def plot_pred(fit_array, file, solution, protein, resin, plt_label):
    axs = plt.gca()

    df_exp = pandas.read_csv(file)
    protein.radius = fit_array[1]
    protein.log_pK_surf_dens()

    ion_str = numpy.linspace(0.05, 0.5, 200)
    df_pred = pandas.DataFrame(numpy.nan, index=range(len(ion_str)), columns=['IS(M)', 'pred_Keq'])

    for index in range(len(df_pred['pred_Keq'])):
        df_pred.at[index, 'IS(M)'] = ion_str[index]
        solution.ion_str = df_pred.at[index, 'IS(M)']
        solution.pH = df_exp.at[0, 'pH']
        get_system_params(solution, protein, resin)

        u_min = find_ener_min(solution, protein, resin)
        df_pred.at[index, 'pred_Keq'] = get_Keq(fit_array[0], u_min)

    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)

    if 'Keq_std' in df_exp.columns:
        axs.loglog(df_pred['IS(M)'], df_pred['pred_Keq']) # label=plt_label+' pred'
        axs.errorbar(df_exp['IS(M)'], df_exp['Keq'],  yerr=df_exp['Keq_std'], fmt='o', color=plt.gca().lines[-1].get_color(), label=plt_label, capsize=4)

    else:
        axs.loglog(df_pred['IS(M)'], df_pred['pred_Keq']) # label=plt_label+' pred'
        axs.scatter(df_exp['IS(M)'], df_exp['Keq'], color=plt.gca().lines[-1].get_color(), label=plt_label, marker='o')
        # axs.legend(loc='lower left', edgecolor='black', handletextpad=0.1)

    handles, labels = axs.get_legend_handles_labels()
    new_handles = []
    for h in handles:
        #only need to edit the errorbar legend entries
        if isinstance(h, container.ErrorbarContainer):
            new_handles.append(h[0])
        else:
            new_handles.append(h)
    axs.legend(new_handles, labels, loc='lower left', bbox_to_anchor=(-0.05, -0.03), frameon=False, handletextpad=0.1)  # loc='center right', bbox_to_anchor=(1.5, 0.5), frameon=False

    axs.set_ylim(8e-9, 1e-5)
    axs.set_xlim(0.1, 0.5)
    axs.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:.1f}'.format(x)))
    axs.xaxis.set_minor_formatter(ticker.FuncFormatter(lambda x, _: '{:.1f}'.format(x)))

    # # axs.set_ylim(7e-9, 1e-5)
    # axs.set_ylim(1e-8, 1e-4)
    # axs.set_xlim(0.05, 0.5)
    # axs.set_xscale('log')
    # axs.xaxis.set_major_formatter(NullFormatter())
    # axs.xaxis.set_minor_formatter(NullFormatter())
    # axs.set_xticks([0.05, 0.1, 0.5])
    # axs.get_xaxis().set_major_formatter(ticker.ScalarFormatter())

    axs.set_xlabel('Ionic strength [M]')
    axs.set_ylabel(r'$K_{eq}$ [m]')
    return


def plot_pred_integrate(fit_array, file, solution, protein, resin, plt_label, x_up_bnd=1e-8):
    df_exp = pandas.read_csv(file)
    protein.radius = fit_array[0]
    protein.log_pK_surf_dens()

    ion_str = numpy.linspace(0.05, 0.5, 500)
    df_pred = pandas.DataFrame(numpy.nan, index=range(len(ion_str)), columns=['IS(M)', 'pred_Keq'])

    for index in range(len(df_pred['pred_Keq'])):
        df_pred.at[index, 'IS(M)'] = ion_str[index]
        solution.ion_str = df_pred.at[index, 'IS(M)']
        solution.pH = df_exp.at[0, 'pH']
        get_system_params(solution, protein, resin)

        df_pred.at[index, 'pred_Keq'] = get_Keq_integrate(solution, protein, resin, x_up_bnd)

    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)

    if 'Keq_std' in df_exp.columns:
        axs.loglog(df_pred['IS(M)'], df_pred['pred_Keq']) # label=plt_label+' pred'
        axs.errorbar(df_exp['IS(M)'], df_exp['Keq'],  yerr=df_exp['Keq_std'], fmt='o', color=plt.gca().lines[-1].get_color(), label=plt_label, capsize=4)

    else:
        axs.loglog(df_pred['IS(M)'], df_pred['pred_Keq']) # label=plt_label+' pred'
        axs.scatter(df_exp['IS(M)'], df_exp['Keq'], color=plt.gca().lines[-1].get_color(), label=plt_label, marker='o')

    handles, labels = axs.get_legend_handles_labels()
    new_handles = []
    for h in handles:
        #only need to edit the errorbar legend entries
        if isinstance(h, container.ErrorbarContainer):
            new_handles.append(h[0])
        else:
            new_handles.append(h)
    axs.legend(new_handles, labels, loc='lower left', bbox_to_anchor=(-0.05, -0.03), frameon=False, handletextpad=0.1)  # loc='center right', bbox_to_anchor=(1.5, 0.5), frameon=False

    # axs.set_ylim(8e-9, 1e-5)
    # axs.set_xlim(0.1, 0.5)
    # axs.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:g}'.format(x)))
    # axs.xaxis.set_minor_formatter(ticker.ScalarFormatter())
    # axs.xaxis.set_major_formatter(ticker.ScalarFormatter())

    # axs.set_ylim(7e-9, 1e-5)
    axs.set_ylim(1e-8, 1e-4)
    axs.set_xlim(0.05, 0.5)
    axs.set_xscale('log')
    axs.xaxis.set_major_formatter(NullFormatter())
    axs.xaxis.set_minor_formatter(NullFormatter())
    axs.set_xticks([0.05, 0.1, 0.5])
    axs.get_xaxis().set_major_formatter(ticker.ScalarFormatter())

    axs.set_xlabel('Ionic strength [M]')
    axs.set_ylabel(r'$K_{eq}$ [m]')
    return

def m_plot_pred(in_array, file, apr, solution, protein, resin, plt_label):
    df_exp = pandas.read_csv(file)

    protein.radius = in_array[0]
    protein.log_pK_surf_dens()
    dpr = in_array[1]

    ion_str = numpy.linspace(0.05, 0.5, 500)
    df_pred = pandas.DataFrame(numpy.nan, index=range(len(ion_str)), columns=['IS(M)', 'pred_Keq'])

    for index in range(len(df_pred['pred_Keq'])):
        df_pred.at[index, 'IS(M)'] = ion_str[index]
        solution.ion_str = df_pred.at[index, 'IS(M)']
        solution.pH = df_exp.at[0, 'pH']
        solution.get_kappa()
        solution.get_ch()

        df_pred.at[index, 'pred_Keq'] = m_get_Keq(dpr, apr, solution, protein, resin)

    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    fig.set_size_inches(8, 7, forward=True)

    if 'Keq_std' in df_exp.columns:
        axs.loglog(df_pred['IS(M)'], df_pred['pred_Keq']) # label=plt_label+' pred'
        axs.errorbar(df_exp['IS(M)'], df_exp['Keq'],  yerr=df_exp['Keq_std'], fmt='o', color=plt.gca().lines[-1].get_color(), label=plt_label, capsize=4)

    else:
        axs.loglog(df_pred['IS(M)'], df_pred['pred_Keq']) # label=plt_label+' pred'
        axs.scatter(df_exp['IS(M)'], df_exp['Keq'], color=plt.gca().lines[-1].get_color(), label=plt_label, marker='o')

    handles, labels = axs.get_legend_handles_labels()
    new_handles = []
    for h in handles:
        #only need to edit the errorbar legend entries
        if isinstance(h, container.ErrorbarContainer):
            new_handles.append(h[0])
        else:
            new_handles.append(h)
    axs.legend(new_handles, labels, loc='lower left', bbox_to_anchor=(-0.05, -0.03), frameon=False, handletextpad=0.1)  # loc='center right', bbox_to_anchor=(1.5, 0.5), frameon=False

    # axs.set_ylim(8e-9, 1e-5)
    # axs.set_xlim(0.1, 0.5)
    # axs.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:g}'.format(x)))
    # axs.xaxis.set_minor_formatter(ticker.ScalarFormatter())
    # axs.xaxis.set_major_formatter(ticker.ScalarFormatter())

#     axs.set_ylim(7e-9, 1e-5)
    axs.set_ylim(1e-8, 1e-4)
    axs.set_xlim(0.05, 0.5)
    axs.set_xscale('log')
    axs.xaxis.set_major_formatter(NullFormatter())
    axs.xaxis.set_minor_formatter(NullFormatter())
    axs.set_xticks([0.05, 0.1, 0.5])
    axs.get_xaxis().set_major_formatter(ticker.ScalarFormatter())

    axs.set_xlabel('Ionic strength [M]')
    axs.set_ylabel(r'$K_{eq}$ [m]')
    return
