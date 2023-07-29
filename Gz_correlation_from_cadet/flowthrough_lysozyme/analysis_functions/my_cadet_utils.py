from analysis_functions.cadet_imports import *
from scipy.interpolate import interp1d

# Function for running CADET____________________________________________________

def run_simulation(cadet, file_name=None):
    if file_name is None:
        f = next(tempfile._get_candidate_names())
        cadet.filename = os.path.join(tempfile.tempdir, f + '.h5')
    else:
        cadet.filename = file_name
    # save the simulation
    cadet.save()

    # run the simulation
    data = cadet.run()

    if data.returncode == 0:
        # print("Simulation completed successfully")
        cadet.load()
    else:
        print(data)
        raise Exception("Simulation failed")

    if file_name is None:
        os.remove(os.path.join(tempfile.tempdir, f + '.h5'))

    return cadet

# General templates_____________________________________________________________

def get_cadet_template(n_units=3, split_components_data=False):
    cadet_template = Cadet()

    cadet_template.root.input.model.nunits = n_units

    # Store solution
    cadet_template.root.input['return'].split_components_data = split_components_data
    cadet_template.root.input['return'].split_ports_data = 0
    cadet_template.root.input['return'].unit_000.write_solution_inlet = 1
    cadet_template.root.input['return'].unit_000.write_solution_outlet = 1
    cadet_template.root.input['return'].unit_000.write_solution_bulk = 1
    cadet_template.root.input['return'].unit_000.write_solution_particle = 1
    cadet_template.root.input['return'].unit_000.write_solution_solid = 1
    cadet_template.root.input['return'].unit_000.write_solution_flux = 1
    cadet_template.root.input['return'].unit_000.write_solution_volume = 1
    cadet_template.root.input['return'].unit_000.write_coordinates = 1
    cadet_template.root.input['return'].unit_000.write_sens_outlet = 1

    for unit in range(n_units):
        cadet_template.root.input['return']['unit_{0:03d}'.format(unit)] =\
         cadet_template.root.input['return'].unit_000

    set_solver_settings(cadet_template)

    return cadet_template

def set_solver_settings(cadet_template):
    # Tolerances for the time integrator
    cadet_template.root.input.solver.time_integrator.abstol = 1e-6
    cadet_template.root.input.solver.time_integrator.algtol = 1e-10
    cadet_template.root.input.solver.time_integrator.reltol = 1e-6
    cadet_template.root.input.solver.time_integrator.init_step_size = 1e-6
    cadet_template.root.input.solver.time_integrator.max_steps = 1000000

    # Solver settings
    cadet_template.root.input.model.solver.gs_type = 1
    cadet_template.root.input.model.solver.max_krylov = 0
    cadet_template.root.input.model.solver.max_restarts = 10
    cadet_template.root.input.model.solver.schur_safety = 1e-8

    # Run the simulation on all threads
    cadet_template.root.input.solver.nthreads = -1
    return

def set_discretization(model, n_bound=None, n_col=150, n_par=30):
    columns = {'GENERAL_RATE_MODEL', 'LUMPED_RATE_MODEL_WITH_PORES',
    'LUMPED_RATE_MODEL_WITHOUT_PORES'}


    for unit_name, unit in model.root.input.model.items():
        if 'unit_' in unit_name and unit.unit_type in columns:
            unit.discretization.ncol = n_col
            unit.discretization.npar = n_par

            if n_bound is None:
                n_bound = unit.ncomp*[0]
            unit.discretization.nbound = n_bound

            unit.discretization.par_disc_type = 'EQUIDISTANT_PAR'
            unit.discretization.use_analytic_jacobian = 1
            unit.discretization.reconstruction = 'WENO'
            unit.discretization.gs_type = 1
            unit.discretization.max_krylov = 0
            unit.discretization.max_restarts = 10
            unit.discretization.schur_safety = 1.0e-8

            unit.discretization.weno.boundary_model = 0
            unit.discretization.weno.weno_eps = 1e-10
            unit.discretization.weno.weno_order = 3
    return

def set_discretization_settings(unit):
    unit.discretization.par_disc_type = 'EQUIDISTANT_PAR'
    unit.discretization.use_analytic_jacobian = 1
    unit.discretization.reconstruction = 'WENO'
    unit.discretization.gs_type = 1
    unit.discretization.max_krylov = 0
    unit.discretization.max_restarts = 10
    unit.discretization.schur_safety = 1.0e-8
    unit.discretization.weno.boundary_model = 0
    unit.discretization.weno.weno_eps = 1e-10
    unit.discretization.weno.weno_order = 3
    return

# Misc functions________________________________________________________________

def cond_2_tis(cond):
    tis = (0.0229207306*(cond**2) + 9.3999636424*cond + 21.2424910935)
    return tis

# Updated simulations - getting extra-column volume and fraction CSTR___________

def set_bypass_sim(cstr_frac, data_df, flow):
    fVolumetric = flow*1e-6/60 # m3/s
    t_load      = 0.1e-6/fVolumetric # s

    c_load = np.trapz(data_df.c_280nm_mM, data_df.UV2_280nm_ml)/0.1
    data_time = list(data_df['t_280nm_s'])

    tube_id  = 0.75 # mm
    area     = np.pi/4.0*(tube_id*1e-3)**2.0 # m2
    velocity = fVolumetric/area # m/s

    v_extra   = 0.44718051*1e-6 # m3

    nComp           = 1
    dpfr_dispersion = 1e-12;

    simulation = Cadet()

    simulation.root.input.model.nunits = 4
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,
     2, 3, -1, -1, fVolumetric]

    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    simulation.root.input.model.unit_000.sec_000.const_coeff = [c_load]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0,]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0,]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0,]

    simulation.root.input.model.unit_000.sec_001.const_coeff = [0.0]
    simulation.root.input.model.unit_000.sec_001.lin_coeff =   [0.0,]
    simulation.root.input.model.unit_000.sec_001.quad_coeff =  [0.0,]
    simulation.root.input.model.unit_000.sec_001.cube_coeff =  [0.0,]

    #create the mixer as a CSTR (for pre-column mixing)
    simulation.root.input.model.unit_001.unit_type = 'CSTR'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_volume = v_extra*cstr_frac
    simulation.root.input.model.unit_001.init_c = [0.0,]

    #for modeling tubing you can use the LUMPED_RATE_MODEL_WITHOUT_PORES (pre-column)
    simulation.root.input.model.unit_002.unit_type = 'LUMPED_RATE_MODEL_WITHOUT_PORES'
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.adsorption_model = 'NONE'
    simulation.root.input.model.unit_002.init_c = [0.0,]
    simulation.root.input.model.unit_002.init_q = [0.0,]
    simulation.root.input.model.unit_002.col_dispersion = dpfr_dispersion
    simulation.root.input.model.unit_002.col_length = v_extra*(1.0-cstr_frac)/area
    simulation.root.input.model.unit_002.total_porosity = 1.0
    simulation.root.input.model.unit_002.velocity = velocity
    simulation.root.input.model.unit_002.discretization.ncol = 25
    simulation.root.input.model.unit_002.discretization.nbound = [0.0,]
    set_discretization_settings(simulation.root.input.model.unit_002)

    #create an outlet
    simulation.root.input.model.unit_003.ncomp = nComp
    simulation.root.input.model.unit_003.unit_type = 'OUTLET'

    #set what values get saved
    simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_003.write_solution_outlet = 1

    #set the times that the simulator writes out data for
    simulation.root.input.solver.user_solution_times = data_time

    #solver settings
    set_solver_settings(simulation)
    simulation.root.input.solver.sections.nsec = 2
    simulation.root.input.solver.sections.section_times = [0.0, t_load, simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = [0.0]
    return simulation

def get_bypass_sim_residual(params_to_fit, data_df, flow, file_path):
    cstr_frac = params_to_fit[0]
    simulation = set_bypass_sim(cstr_frac, data_df, flow)
    run_simulation(simulation, file_path)
    simulation.load()
    c_out  = simulation.root.output.solution.unit_003.solution_outlet_comp_000
    residual_vec = np.array(data_df.c_280nm_mM) - c_out
    return list(residual_vec)

def get_cstr_frac(flow):
    fun = get_cstr_frac_interp()
    return fun(flow)

def get_cstr_frac_previous(flow):
    return -0.259*flow + 0.561

def get_cstr_frac_interp():
    x = [0.98, 0.49, 0.1, 0.3, 0.7]
    y = [0.32864428948848445, 0.4620417890992433, 0.5309155507994068, 0.5222044334365368, 0.3951717146519932]
    return interp1d(x, y, kind='linear', fill_value='extrapolate')




# Updated simulations - getting Dax, Dp, eps_c, and eps_p_______________________

def set_breakthrough_sim_pulse(Keq, Dax, Dp, Ds, eps_c, eps_p, data_df, flow):
    cstr_frac = get_cstr_frac(flow)

    fVolumetric = flow*1e-6/60 # m3/s
    t_load      = 0.1e-6/fVolumetric # s

    c_load = np.trapz(data_df.c_280nm_mM, data_df.UV2_280nm_ml)/0.1
    data_time = list(data_df['t_280nm_s'])

    tube_id  = 0.75 # mm
    area     = np.pi/4.0*(tube_id*1e-3)**2.0 # m2
    velocity = fVolumetric/area # m/s

    v_extra   = 0.3598*1e-6 # m3

    nComp           = 1
    dpfr_dispersion = 1e-12;

    simulation = Cadet()

    simulation.root.input.model.nunits = 5
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 3, -1, -1, fVolumetric,
     3, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,
     2, 4, -1, -1, fVolumetric]

    # [0, 1, -1, -1, fVolumetric,
    #  1, 2, -1, -1, fVolumetric,
    #  2, 3, -1, -1, fVolumetric,
    #  3, 4, -1, -1, fVolumetric]

    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    simulation.root.input.model.unit_000.sec_000.const_coeff = [c_load]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0,]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0,]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0,]

    simulation.root.input.model.unit_000.sec_001.const_coeff = [0.0]
    simulation.root.input.model.unit_000.sec_001.lin_coeff =   [0.0,]
    simulation.root.input.model.unit_000.sec_001.quad_coeff =  [0.0,]
    simulation.root.input.model.unit_000.sec_001.cube_coeff =  [0.0,]

    #create the mixer as a CSTR (for pre-column mixing)
    simulation.root.input.model.unit_001.unit_type = 'CSTR'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_volume = v_extra*cstr_frac
    simulation.root.input.model.unit_001.init_c = [0.0,]

    #for modeling tubing you can use the LUMPED_RATE_MODEL_WITHOUT_PORES (pre-column)
    simulation.root.input.model.unit_002.unit_type = 'LUMPED_RATE_MODEL_WITHOUT_PORES'
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.adsorption_model = 'NONE'
    simulation.root.input.model.unit_002.init_c = [0.0,]
    simulation.root.input.model.unit_002.init_q = [0.0,]
    simulation.root.input.model.unit_002.col_dispersion = dpfr_dispersion
    simulation.root.input.model.unit_002.col_length = v_extra*(1.0-cstr_frac)/area
    simulation.root.input.model.unit_002.total_porosity = 1.0
    simulation.root.input.model.unit_002.velocity = velocity
    simulation.root.input.model.unit_002.discretization.ncol = 25
    simulation.root.input.model.unit_002.discretization.nbound = [0.0,]
    set_discretization_settings(simulation.root.input.model.unit_002)

    # Column model with adsorption
    simulation.root.input.model.unit_003.unit_type = 'GENERAL_RATE_MODEL'
    simulation.root.input.model.unit_003.ncomp = nComp
    simulation.root.input.model.unit_003.init_c                  = [0.0,]
    simulation.root.input.model.unit_003.init_q                  = [0.0,]
    simulation.root.input.model.unit_003.col_dispersion          = Dax
    simulation.root.input.model.unit_003.col_length              = 0.10
    simulation.root.input.model.unit_003.col_porosity            = eps_c
    simulation.root.input.model.unit_003.film_diffusion          = [1e-4,]
    simulation.root.input.model.unit_003.par_porosity            = eps_p
    simulation.root.input.model.unit_003.par_radius              = 4.5e-5
    simulation.root.input.model.unit_003.par_diffusion           = [Dp]
    simulation.root.input.model.unit_003.par_surfdiffusion       = [Ds]
    simulation.root.input.model.unit_003.cross_section_area      = np.pi/4 * (0.5**2) * 1e-4
    simulation.root.input.model.unit_003.discretization.ncol     = 100
    simulation.root.input.model.unit_003.discretization.npar     = 30
    simulation.root.input.model.unit_003.discretization.nbound   = [1]
    set_discretization_settings(simulation.root.input.model.unit_003)
    # Adsorption model
    simulation.root.input.model.unit_003.adsorption_model = 'LINEAR'
    simulation.root.input.model.unit_003.adsorption.is_kinetic = 0
    simulation.root.input.model.unit_003.adsorption.lin_ka = [Keq]
    simulation.root.input.model.unit_003.adsorption.lin_kd = [1.0]

    #create an outlet
    simulation.root.input.model.unit_004.ncomp = nComp
    simulation.root.input.model.unit_004.unit_type = 'OUTLET'

    #set what values get saved
    simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_004.write_solution_outlet = 1

    #set the times that the simulator writes out data for
    simulation.root.input.solver.user_solution_times = data_time

    #solver settings
    set_solver_settings(simulation)
    simulation.root.input.solver.sections.nsec = 2
    simulation.root.input.solver.sections.section_times = [0.0, t_load, simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = [0.0]
    return simulation


def get_flowthrough_sim_residual(params_to_fit, Keq, Ds, eps_c, eps_p, data_df, flow, Dp):
    # a, b, c, Dp = params_to_fit[0], params_to_fit[1], params_to_fit[2], params_to_fit[3]
    # Dax = a*flow**2.0 + b*b + c
    # if Dax < 0:
    #     return 1e5
    # else:
    #     simulation = set_breakthrough_sim_pulse(Keq, Dax, Dp, Ds, eps_c, eps_p, data_df, flow)
    #     simulation = run_simulation(simulation)
    #     c_out  = simulation.root.output.solution.unit_004.solution_outlet_comp_000
    #     residual_vec = np.array(data_df.c_280nm_mM) - np.array(c_out)
    #     return np.linalg.norm(residual_vec)

    Dax = params_to_fit[0]
    simulation = set_breakthrough_sim_pulse(Keq, Dax, Dp, Ds, eps_c, eps_p, data_df, flow)
    simulation = run_simulation(simulation)
    c_out  = simulation.root.output.solution.unit_004.solution_outlet_comp_000
    residual_vec = np.array(data_df.c_280nm_mM) - np.array(c_out)
    return np.linalg.norm(residual_vec)



# Updated simulations - breakthrough simulations________________________________

def set_breakthrough_sim(Keq, c_load, Dax, Dp, Ds, eps_c, eps_p, data_df, flow):
    cstr_frac = get_cstr_frac(flow)

    fVolumetric = flow*1e-6/60 # m3/s

    data_time = list(data_df['t_215nm_s']) # Could improve by changing the df pass

    tube_id  = 0.75 # mm
    area     = np.pi/4.0*(tube_id*1e-3)**2.0 # m2
    velocity = fVolumetric/area # m/s

    v_extra   = 0.3098*1e-6 # m3


    nComp           = 1
    dpfr_dispersion = 1e-12;

    simulation = Cadet()

    simulation.root.input.model.nunits = 5
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,
     2, 3, -1, -1, fVolumetric,
     3, 4, -1, -1, fVolumetric]

    # [0, 3, -1, -1, fVolumetric,
    #  3, 1, -1, -1, fVolumetric,
    #  1, 2, -1, -1, fVolumetric,
    #  2, 4, -1, -1, fVolumetric]


    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    simulation.root.input.model.unit_000.sec_000.const_coeff = [c_load]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0,]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0,]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0,]

    #create the mixer as a CSTR (for pre-column mixing)
    simulation.root.input.model.unit_001.unit_type = 'CSTR'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_volume = v_extra*cstr_frac
    simulation.root.input.model.unit_001.init_c = [0.0,]

    #for modeling tubing you can use the LUMPED_RATE_MODEL_WITHOUT_PORES (pre-column)
    simulation.root.input.model.unit_002.unit_type = 'LUMPED_RATE_MODEL_WITHOUT_PORES'
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.adsorption_model = 'NONE'
    simulation.root.input.model.unit_002.init_c = [0.0,]
    simulation.root.input.model.unit_002.init_q = [0.0,]
    simulation.root.input.model.unit_002.col_dispersion = dpfr_dispersion
    simulation.root.input.model.unit_002.col_length = v_extra*(1.0-cstr_frac)/area
    simulation.root.input.model.unit_002.total_porosity = 1.0
    simulation.root.input.model.unit_002.velocity = velocity
    simulation.root.input.model.unit_002.discretization.ncol = 25
    simulation.root.input.model.unit_002.discretization.nbound = [0.0,]
    set_discretization_settings(simulation.root.input.model.unit_002)

    # Column model with adsorption
    simulation.root.input.model.unit_003.unit_type = 'GENERAL_RATE_MODEL'
    simulation.root.input.model.unit_003.ncomp = nComp
    simulation.root.input.model.unit_003.init_c                  = [0.0,]
    simulation.root.input.model.unit_003.init_q                  = [0.0,]
    simulation.root.input.model.unit_003.col_dispersion          = Dax
    simulation.root.input.model.unit_003.col_length              = 0.10
    simulation.root.input.model.unit_003.col_porosity            = eps_c
    simulation.root.input.model.unit_003.film_diffusion          = [1e-4,]
    simulation.root.input.model.unit_003.par_porosity            = eps_p
    simulation.root.input.model.unit_003.par_radius              = 4.5e-5
    simulation.root.input.model.unit_003.par_diffusion           = [Dp]
    simulation.root.input.model.unit_003.par_surfdiffusion       = [Ds]
    simulation.root.input.model.unit_003.cross_section_area      = np.pi/4 * (0.5**2) * 1e-4
    simulation.root.input.model.unit_003.discretization.ncol     = 100
    simulation.root.input.model.unit_003.discretization.npar     = 30
    simulation.root.input.model.unit_003.discretization.nbound   = [1]
    set_discretization_settings(simulation.root.input.model.unit_003)
    # Adsorption model
    simulation.root.input.model.unit_003.adsorption_model = 'LINEAR'
    simulation.root.input.model.unit_003.adsorption.is_kinetic = 0
    simulation.root.input.model.unit_003.adsorption.lin_ka = [Keq]
    simulation.root.input.model.unit_003.adsorption.lin_kd = [1.0]

    #create an outlet
    simulation.root.input.model.unit_004.ncomp = nComp
    simulation.root.input.model.unit_004.unit_type = 'OUTLET'

    #set what values get saved
    simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_004.write_solution_outlet = 1

    #set the times that the simulator writes out data for
    simulation.root.input.solver.user_solution_times = data_time

    #solver settings
    set_solver_settings(simulation)
    simulation.root.input.solver.sections.nsec = 1
    simulation.root.input.solver.sections.section_times = [0.0, simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = [0.0]
    return simulation


def get_breakthrough_sim_residual(params_to_fit, c_load, Dp, eps_c, eps_p, data_df, flow):
    Keq, Ds = params_to_fit[0], params_to_fit[1]
    Dax = get_Dax(flow)
    simulation = set_breakthrough_sim(Keq, c_load, Dax, Dp, Ds, eps_c, eps_p, data_df, flow)
    simulation = run_simulation(simulation)
    c_out  = simulation.root.output.solution.unit_004.solution_outlet_comp_000
    residual_vec = np.array(data_df.c_215nm_mM) - np.array(c_out)
    # print(np.linalg.norm(residual_vec), Keq, Ds)
    return np.linalg.norm(residual_vec)


def get_Dax(flow):
    # Dax_coeff = 1.75441212e-06 # flow in [ml/min]
    # return Dax_coeff * flow

    Dax_fun = np.poly1d(np.array([-3.20672717e-06,  5.68204258e-06, -4.11660273e-07,  3.32467066e-08]))
    return Dax_fun(flow)




































#
