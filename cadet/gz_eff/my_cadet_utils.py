from cadet_imports import *

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
        print("Simulation completed successfully")
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

    set_solver_settings_2(cadet_template)

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

def set_solver_settings_2(cadet_template):
    # Tolerances for the time integrator
    cadet_template.root.input.solver.time_integrator.abstol = 1e-7
    cadet_template.root.input.solver.time_integrator.algtol = 1e-11
    cadet_template.root.input.solver.time_integrator.reltol = 1e-7
    cadet_template.root.input.solver.time_integrator.init_step_size = 1e-7
    cadet_template.root.input.solver.time_integrator.max_steps = 1000000

    # Solver settings
    cadet_template.root.input.model.solver.gs_type = 1
    cadet_template.root.input.model.solver.max_krylov = 0
    cadet_template.root.input.model.solver.max_restarts = 10
    cadet_template.root.input.model.solver.schur_safety = 1e-9

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

# Template for Akta Explorer____________________________________________________

def get_explorer_system(salt_c=0.0, load_c=0.0, Dp=0.0, Ds=0.0, q_max=1.0, ka=1.0):
    simulation = Cadet()

    # My system parameters, 0.98 ml/min
    nComp = 2
    fVolumetric = 1.633e-8 # m**3/s, 0.98 ml/min
    velocity = 3.697e-2 # m/s
    area = fVolumetric/velocity # m**2
    cstr_frac_before = 0.5;
    cstr_frac_after  = 0.2;
    v_extra = 1.6e-7
    vol_before = v_extra;
    vol_after = v_extra;
    dpfr_dispersion = 1e-12;

    simulation.root.input.model.nunits = 7
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,
     2, 3, -1, -1, fVolumetric,
     3, 5, -1, -1, fVolumetric,
     5, 4, -1, -1, fVolumetric,
     4, 6, -1, -1, fVolumetric]

    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    #const_coeff + lin_coeff*t + quad_coeff*t**2 + cube_coeff*t**3
    simulation.root.input.model.unit_000.sec_000.const_coeff = [salt_c,load_c]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0,0.0,]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0,0.0,]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0,0.0,]

    #create the mixer as a CSTR (for pre-column mixing)
    simulation.root.input.model.unit_001.unit_type = 'CSTR'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_volume = vol_before*cstr_frac_before
    simulation.root.input.model.unit_001.init_c = [salt_c,0.0,]

    #for modeling tubing you can use the LUMPED_RATE_MODEL_WITHOUT_PORES (pre-column)
    simulation.root.input.model.unit_002.unit_type = 'LUMPED_RATE_MODEL_WITHOUT_PORES'
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.adsorption_model = 'NONE'
    simulation.root.input.model.unit_002.init_c = [salt_c,0.0,]
    simulation.root.input.model.unit_002.init_q = [0.0,0.0,]
    simulation.root.input.model.unit_002.col_dispersion = dpfr_dispersion
    simulation.root.input.model.unit_002.col_length = 4/np.pi*((1/(0.75*10**(-3)))**2)*vol_before*(1-cstr_frac_before)
    simulation.root.input.model.unit_002.total_porosity = 1.0
    simulation.root.input.model.unit_002.velocity = velocity
    simulation.root.input.model.unit_002.discretization.ncol = 25
    simulation.root.input.model.unit_002.discretization.nbound = [0,0,]
    set_discretization_settings(simulation.root.input.model.unit_002)

    # Column model with adsorption
    simulation.root.input.model.unit_003.unit_type = 'GENERAL_RATE_MODEL'
    simulation.root.input.model.unit_003.ncomp = nComp
    simulation.root.input.model.unit_003.init_c                  = [salt_c,0.0,]
    simulation.root.input.model.unit_003.init_q                  = [0.0,0.0,]
    simulation.root.input.model.unit_003.col_dispersion          = 5.0e-7
    simulation.root.input.model.unit_003.col_length              = 0.042
    simulation.root.input.model.unit_003.col_porosity            = 0.49
    simulation.root.input.model.unit_003.film_diffusion          = [6.9e-6, 1e-3,]
    simulation.root.input.model.unit_003.par_porosity            = 0.4
    simulation.root.input.model.unit_003.par_radius              = 2.5e-5
    simulation.root.input.model.unit_003.par_diffusion           = [7e-10, Dp]
    simulation.root.input.model.unit_003.par_surfdiffusion       = [0, Ds]
    simulation.root.input.model.unit_003.cross_section_area      = np.pi/4 * (0.5**2) * 1e-4
    simulation.root.input.model.unit_003.discretization.ncol     = 150
    simulation.root.input.model.unit_003.discretization.npar     = 75 # Used 50 for protein != BLG
    simulation.root.input.model.unit_003.discretization.nbound   = [0, 1]
    set_discretization_settings(simulation.root.input.model.unit_003)
    # Adsorption model
    simulation.root.input.model.unit_003.adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
    simulation.root.input.model.unit_003.adsorption.is_kinetic = 0
    simulation.root.input.model.unit_003.adsorption.mcl_qmax = [1, q_max]
    simulation.root.input.model.unit_003.adsorption.mcl_ka = [0, ka]
    simulation.root.input.model.unit_003.adsorption.mcl_kd = [1, 1]

    #create the mixer as a CSTR (for post-column mixing)
    simulation.root.input.model.unit_004 = copy.deepcopy(simulation.root.input.model.unit_001)
    simulation.root.input.model.unit_004.init_volume = vol_after*cstr_frac_after

    #for modeling tubing you can use the LUMPED_RATE_MODEL_WITHOUT_PORES (post-column)
    simulation.root.input.model.unit_005 = copy.deepcopy(simulation.root.input.model.unit_002)
    simulation.root.input.model.unit_005.col_length = 4/np.pi*((1/(0.75*10**(-3)))**2)*vol_after*(1-cstr_frac_after)

    #create an outlet
    simulation.root.input.model.unit_006.ncomp = nComp
    simulation.root.input.model.unit_006.unit_type = 'OUTLET'

    #set what values get saved
    simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_006.write_solution_outlet = 1

    simulation.root.input['return'].unit_003.write_solution_bulk =     1
    simulation.root.input['return'].unit_003.write_solution_particle = 1
    simulation.root.input['return'].unit_003.write_solution_solid =    1

    #set the times that the simulator writes out data for
    simulation.root.input.solver.user_solution_times = np.linspace(0, 5000, 5001)

    #solver settings
    set_solver_settings_2(simulation)
    simulation.root.input.solver.sections.nsec = 1
    simulation.root.input.solver.sections.section_times = [0.0, simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = []

    return simulation


# Template for column without extra-column effects______________________________

def get_column_system(salt_c=0.0, load_c=0.0, Dp=0.0, Ds=0.0, q_max=1.0, ka=1.0, linear=False):
    simulation = Cadet()

    # My system parameters, 0.98 ml/min
    nComp = 2
    fVolumetric = 1.633e-8 # m**3/s, 0.98 ml/min

    simulation.root.input.model.nunits = 3
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,]

    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    #const_coeff + lin_coeff*t + quad_coeff*t**2 + cube_coeff*t**3
    simulation.root.input.model.unit_000.sec_000.const_coeff = [salt_c,load_c]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0,0.0,]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0,0.0,]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0,0.0,]

    # Column model with adsorption
    simulation.root.input.model.unit_001.unit_type = 'GENERAL_RATE_MODEL'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_c                  = [salt_c,0.0,]
    simulation.root.input.model.unit_001.init_q                  = [0.0,0.0,]
    simulation.root.input.model.unit_001.col_dispersion          = 5.0e-7
    simulation.root.input.model.unit_001.col_length              = 0.042
    simulation.root.input.model.unit_001.col_porosity            = 0.49
    simulation.root.input.model.unit_001.film_diffusion          = [6.9e-6, 1e-3,]
    simulation.root.input.model.unit_001.par_porosity            = 0.4
    simulation.root.input.model.unit_001.par_radius              = 2.5e-5
    simulation.root.input.model.unit_001.par_diffusion           = [7e-10, Dp]
    simulation.root.input.model.unit_001.par_surfdiffusion       = [0, Ds]
    simulation.root.input.model.unit_001.cross_section_area      = np.pi/4 * (0.5**2) * 1e-4
    simulation.root.input.model.unit_001.discretization.ncol     = 150
    simulation.root.input.model.unit_001.discretization.npar     = 75 # Used 50 for protein != BLG
    simulation.root.input.model.unit_001.discretization.nbound   = [0, 1]
    set_discretization_settings(simulation.root.input.model.unit_001)
    # Adsorption model
    if linear == True:
        simulation.root.input.model.unit_001.adsorption_model = 'LINEAR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.lin_ka = [0, ka]
        simulation.root.input.model.unit_001.adsorption.lin_kd = [1, 1]
    else:
        simulation.root.input.model.unit_001.adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.mcl_qmax = [1, q_max]
        simulation.root.input.model.unit_001.adsorption.mcl_ka = [0, ka]
        simulation.root.input.model.unit_001.adsorption.mcl_kd = [1, 1]

    #create an outlet
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.unit_type = 'OUTLET'

    #set what values get saved
    simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_002.write_solution_outlet = 1
    # simulation.root.input['return'].unit_001.write_solution_bulk =     1
    # simulation.root.input['return'].unit_001.write_solution_particle = 1
    # simulation.root.input['return'].unit_001.write_solution_solid =    1

    #set the times that the simulator writes out data for
    simulation.root.input.solver.user_solution_times = np.linspace(0, 5000, 5001)

    #solver settings
    set_solver_settings_2(simulation)
    simulation.root.input.solver.sections.nsec = 1
    simulation.root.input.solver.sections.section_times = [0.0, simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = []

    return simulation


def get_column_system_1_comp(load_c=0.0, Dp=0.0, Ds=0.0, q_max=1.0, ka=1.0, t_max=5000,
                             use_linear=False, fVolumetric=None):
    simulation = Cadet()

    # My system parameters, 0.98 ml/min
    nComp = 1
    if fVolumetric is None:
        fVolumetric = 1.633e-8 #/4.0 # m**3/s, 0.98 ml/min

    simulation.root.input.model.nunits = 3
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,]

    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    #const_coeff + lin_coeff*t + quad_coeff*t**2 + cube_coeff*t**3
    simulation.root.input.model.unit_000.sec_000.const_coeff = [load_c]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0,]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0,]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0,]

    # Column model with adsorption
    simulation.root.input.model.unit_001.unit_type = 'GENERAL_RATE_MODEL'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_c                  = [0.0,]
    simulation.root.input.model.unit_001.init_q                  = [0.0,]
    simulation.root.input.model.unit_001.col_dispersion          = 5.0e-7/4.0
    simulation.root.input.model.unit_001.col_length              = 0.042
    simulation.root.input.model.unit_001.col_porosity            = 0.49
    simulation.root.input.model.unit_001.film_diffusion          = [1e-3,]
    simulation.root.input.model.unit_001.par_porosity            = 0.4
    simulation.root.input.model.unit_001.par_radius              = 2.5e-5
    simulation.root.input.model.unit_001.par_diffusion           = [Dp]
    simulation.root.input.model.unit_001.par_surfdiffusion       = [Ds]
    simulation.root.input.model.unit_001.cross_section_area      = np.pi/4 * (0.5**2) * 1e-4
    simulation.root.input.model.unit_001.discretization.ncol     = 150
    simulation.root.input.model.unit_001.discretization.npar     = 75 # Used 50 for protein != BLG
    simulation.root.input.model.unit_001.discretization.nbound   = [1]
    set_discretization_settings(simulation.root.input.model.unit_001)
    # Adsorption model
    if use_linear:
        simulation.root.input.model.unit_001.adsorption_model = 'LINEAR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.lin_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.lin_kd = [1]
    else:
        simulation.root.input.model.unit_001.adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.mcl_qmax = [q_max]
        simulation.root.input.model.unit_001.adsorption.mcl_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.mcl_kd = [1]

    #create an outlet
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.unit_type = 'OUTLET'

    #set what values get saved
    # simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_002.write_solution_outlet = 1
    # simulation.root.input['return'].unit_001.write_solution_bulk =     1
    # simulation.root.input['return'].unit_001.write_solution_particle = 1
    # simulation.root.input['return'].unit_001.write_solution_solid =    1

    #set the times that the simulator writes out data for
    simulation.root.input.solver.user_solution_times = np.linspace(0, int(t_max), int(t_max)+1)

    #solver settings
    set_solver_settings_2(simulation)
    simulation.root.input.solver.sections.nsec = 1
    simulation.root.input.solver.sections.section_times = [0.0, simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = []

    return simulation


def get_column_system_1_comp_with_Q(load_c=0.0, Dp=0.0, Ds=0.0, q_max=1.0, ka=1.0,
                                    t_max=5000, use_linear=False, fVolumetric=0.0):
    simulation = Cadet()

    # My system parameters, 0.98 ml/min
    nComp = 1

    simulation.root.input.model.nunits = 3
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,]

    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    #const_coeff + lin_coeff*t + quad_coeff*t**2 + cube_coeff*t**3
    simulation.root.input.model.unit_000.sec_000.const_coeff = [load_c]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0,]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0,]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0,]

    # Column model with adsorption
    simulation.root.input.model.unit_001.unit_type = 'GENERAL_RATE_MODEL'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_c                  = [0.0,]
    simulation.root.input.model.unit_001.init_q                  = [0.0,]
    simulation.root.input.model.unit_001.col_dispersion          = 5.0e-7*fVolumetric/1.633e-8
    simulation.root.input.model.unit_001.col_length              = 0.042
    simulation.root.input.model.unit_001.col_porosity            = 0.49
    simulation.root.input.model.unit_001.film_diffusion          = [1e-3,]
    simulation.root.input.model.unit_001.par_porosity            = 0.4
    simulation.root.input.model.unit_001.par_radius              = 2.5e-5
    simulation.root.input.model.unit_001.par_diffusion           = [Dp]
    simulation.root.input.model.unit_001.par_surfdiffusion       = [Ds]
    simulation.root.input.model.unit_001.cross_section_area      = np.pi/4 * (0.5**2) * 1e-4
    simulation.root.input.model.unit_001.discretization.ncol     = 150
    simulation.root.input.model.unit_001.discretization.npar     = 75 # Used 50 for protein != BLG
    simulation.root.input.model.unit_001.discretization.nbound   = [1]
    set_discretization_settings(simulation.root.input.model.unit_001)
    # Adsorption model
    if use_linear:
        simulation.root.input.model.unit_001.adsorption_model = 'LINEAR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.lin_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.lin_kd = [1]
    else:
        simulation.root.input.model.unit_001.adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.mcl_qmax = [q_max]
        simulation.root.input.model.unit_001.adsorption.mcl_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.mcl_kd = [1]

    #create an outlet
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.unit_type = 'OUTLET'

    #set what values get saved
    # simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_002.write_solution_outlet = 1
    # simulation.root.input['return'].unit_001.write_solution_bulk =     1
    # simulation.root.input['return'].unit_001.write_solution_particle = 1
    # simulation.root.input['return'].unit_001.write_solution_solid =    1

    #set the times that the simulator writes out data for
    simulation.root.input.solver.user_solution_times = np.linspace(0, int(t_max), int(t_max)+1)

    #solver settings
    set_solver_settings_2(simulation)
    simulation.root.input.solver.sections.nsec = 1
    simulation.root.input.solver.sections.section_times = [0.0, simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = []

    return simulation


def get_column_system_1_comp_pulse(load_c=0.0, Dp=0.0, Ds=0.0, q_max=0.0, ka=1.0, t_max=5000,
                             use_linear=True):
    simulation = Cadet()

    # My system parameters, 0.98 ml/min
    nComp = 1
    fVolumetric = 1.633e-8/4.0 # m**3/s, 0.98/4 ml/min

    simulation.root.input.model.nunits = 3
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,]

    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    #const_coeff + lin_coeff*t + quad_coeff*t**2 + cube_coeff*t**3
    simulation.root.input.model.unit_000.sec_000.const_coeff = [load_c]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0]

    simulation.root.input.model.unit_000.sec_001.const_coeff = [0.0]
    simulation.root.input.model.unit_000.sec_001.lin_coeff =   [0.0]
    simulation.root.input.model.unit_000.sec_001.quad_coeff =  [0.0]
    simulation.root.input.model.unit_000.sec_001.cube_coeff =  [0.0]

    # Column model with adsorption
    simulation.root.input.model.unit_001.unit_type = 'GENERAL_RATE_MODEL'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_c                  = [0.0,]
    simulation.root.input.model.unit_001.init_q                  = [0.0,]
    simulation.root.input.model.unit_001.col_dispersion          = 5.0e-7/4.0
    simulation.root.input.model.unit_001.col_length              = 0.042
    simulation.root.input.model.unit_001.col_porosity            = 0.49
    simulation.root.input.model.unit_001.film_diffusion          = [1e-3,]
    simulation.root.input.model.unit_001.par_porosity            = 0.4
    simulation.root.input.model.unit_001.par_radius              = 2.5e-5
    simulation.root.input.model.unit_001.par_diffusion           = [Dp]
    simulation.root.input.model.unit_001.par_surfdiffusion       = [Ds]
    simulation.root.input.model.unit_001.cross_section_area      = np.pi/4 * (0.5**2) * 1e-4
    simulation.root.input.model.unit_001.discretization.ncol     = 150
    simulation.root.input.model.unit_001.discretization.npar     = 75 # Used 50 for protein != BLG
    simulation.root.input.model.unit_001.discretization.nbound   = [1]
    set_discretization_settings(simulation.root.input.model.unit_001)
    # Adsorption model
    if use_linear:
        simulation.root.input.model.unit_001.adsorption_model = 'LINEAR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.lin_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.lin_kd = [1]
    else:
        simulation.root.input.model.unit_001.adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.mcl_qmax = [q_max]
        simulation.root.input.model.unit_001.adsorption.mcl_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.mcl_kd = [1]

    #create an outlet
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.unit_type = 'OUTLET'

    #set what values get saved
    # simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_002.write_solution_outlet = 1
    # simulation.root.input['return'].unit_001.write_solution_bulk =     1
    # simulation.root.input['return'].unit_001.write_solution_particle = 1
    # simulation.root.input['return'].unit_001.write_solution_solid =    1

    #set the times that the simulator writes out data for
    simulation.root.input.solver.user_solution_times = np.linspace(0, int(t_max), int(t_max)+1)

    #solver settings
    set_solver_settings_2(simulation)
    simulation.root.input.solver.sections.nsec = 2
    simulation.root.input.solver.sections.section_times = [0.0, 10.0,
                            simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = [0.0]

    return simulation


def get_more_general_system(load_c=0.0, Dp=0.0, Ds=0.0, q_max=1.0, ka=1.0,
                            t_max=5000, use_linear=False, fVolumetric=0.0,
                            d_part=5.0e-5, l_col=5.0e-2):
    simulation = Cadet()

    # My system parameters, 0.98 ml/min
    nComp = 1

    simulation.root.input.model.nunits = 3
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,]

    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    #const_coeff + lin_coeff*t + quad_coeff*t**2 + cube_coeff*t**3
    simulation.root.input.model.unit_000.sec_000.const_coeff = [load_c]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0,]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0,]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0,]

    # Column model with adsorption
    simulation.root.input.model.unit_001.unit_type = 'GENERAL_RATE_MODEL'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_c                  = [0.0,]
    simulation.root.input.model.unit_001.init_q                  = [0.0,]
    simulation.root.input.model.unit_001.col_dispersion          = 5.0e-7*fVolumetric/1.633e-8
    simulation.root.input.model.unit_001.col_length              = l_col
    simulation.root.input.model.unit_001.col_porosity            = 0.49
    simulation.root.input.model.unit_001.film_diffusion          = [1e-3,]
    simulation.root.input.model.unit_001.par_porosity            = 0.4
    simulation.root.input.model.unit_001.par_radius              = d_part/2.0
    simulation.root.input.model.unit_001.par_diffusion           = [Dp]
    simulation.root.input.model.unit_001.par_surfdiffusion       = [Ds]
    simulation.root.input.model.unit_001.cross_section_area      = np.pi/4 * (0.5**2) * 1e-4
    simulation.root.input.model.unit_001.discretization.ncol     = 150
    simulation.root.input.model.unit_001.discretization.npar     = 75 # Used 50 for protein != BLG
    simulation.root.input.model.unit_001.discretization.nbound   = [1]
    set_discretization_settings(simulation.root.input.model.unit_001)
    # Adsorption model
    if use_linear:
        simulation.root.input.model.unit_001.adsorption_model = 'LINEAR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.lin_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.lin_kd = [1]
    else:
        simulation.root.input.model.unit_001.adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.mcl_qmax = [q_max]
        simulation.root.input.model.unit_001.adsorption.mcl_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.mcl_kd = [1]

    #create an outlet
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.unit_type = 'OUTLET'

    #set what values get saved
    # simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_002.write_solution_outlet = 1
    # simulation.root.input['return'].unit_001.write_solution_bulk =     1
    # simulation.root.input['return'].unit_001.write_solution_particle = 1
    # simulation.root.input['return'].unit_001.write_solution_solid =    1

    #set the times that the simulator writes out data for
    simulation.root.input.solver.user_solution_times = np.linspace(0, int(t_max), int(t_max)+1)

    #solver settings
    set_solver_settings_2(simulation)
    simulation.root.input.solver.sections.nsec = 1
    simulation.root.input.solver.sections.section_times = [0.0, simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = []

    return simulation

def get_Dax(d_part, fVolum, area_col, eps_c, D=7.5e-11):
    logy_fun = np.poly1d([ 1.2267772, -0.1996732])
    x = d_part/D * fVolum/area_col * 1/(1-eps_c)
    logx = np.log10(x)
    logy = logy_fun(logx)
    Dax = D * 10.0**(logy)
    return Dax

def get_system_updated(load_c=0.0, Dp=0.0, Ds=0.0, q_max=1.0, ka=1.0,
                        t_max=5000, use_linear=False, fVolumetric=0.0,
                        d_part=5.0e-5, l_col=5.0e-2, eps_c=0.49, eps_p=0.4):
    simulation = Cadet()

    d_col = 0.5e-2
    area_col = np.pi/4 * d_col**2

    eps_t = eps_c + eps_p*(1.0 - eps_c)

    # My system parameters, 0.98 ml/min
    nComp = 1

    simulation.root.input.model.nunits = 3
    simulation.root.input.model.connections.nswitches = 1
    simulation.root.input.model.connections.switch_000.section = 0

    #this connects the system
    simulation.root.input.model.connections.switch_000.connections = \
    [0, 1, -1, -1, fVolumetric,
     1, 2, -1, -1, fVolumetric,]

    #create an inlet
    simulation.root.input.model.unit_000.unit_type = 'INLET'
    simulation.root.input.model.unit_000.ncomp = nComp
    simulation.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    #const_coeff + lin_coeff*t + quad_coeff*t**2 + cube_coeff*t**3
    simulation.root.input.model.unit_000.sec_000.const_coeff = [load_c]
    simulation.root.input.model.unit_000.sec_000.lin_coeff =   [0.0,]
    simulation.root.input.model.unit_000.sec_000.quad_coeff =  [0.0,]
    simulation.root.input.model.unit_000.sec_000.cube_coeff =  [0.0,]

    # Column model with adsorption
    simulation.root.input.model.unit_001.unit_type = 'GENERAL_RATE_MODEL'
    simulation.root.input.model.unit_001.ncomp = nComp
    simulation.root.input.model.unit_001.init_c                  = [0.0,]
    simulation.root.input.model.unit_001.init_q                  = [0.0,]
    simulation.root.input.model.unit_001.col_dispersion          = get_Dax(d_part, fVolumetric, area_col, eps_c, D=7.5e-11)
    simulation.root.input.model.unit_001.col_length              = l_col
    simulation.root.input.model.unit_001.col_porosity            = eps_c
    simulation.root.input.model.unit_001.film_diffusion          = [1e-3,]
    simulation.root.input.model.unit_001.par_porosity            = eps_p
    simulation.root.input.model.unit_001.par_radius              = d_part/2.0
    simulation.root.input.model.unit_001.par_diffusion           = [Dp]
    simulation.root.input.model.unit_001.par_surfdiffusion       = [Ds]
    simulation.root.input.model.unit_001.cross_section_area      = area_col
    simulation.root.input.model.unit_001.discretization.ncol     = int(2000*l_col)
    simulation.root.input.model.unit_001.discretization.npar     = max(25, int(1.5e6 * d_part)) # Used 50 for protein != BLG
    simulation.root.input.model.unit_001.discretization.nbound   = [1]
    set_discretization_settings(simulation.root.input.model.unit_001)
    # Adsorption model
    if use_linear:
        simulation.root.input.model.unit_001.adsorption_model = 'LINEAR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.lin_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.lin_kd = [1]
    else:
        simulation.root.input.model.unit_001.adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
        simulation.root.input.model.unit_001.adsorption.is_kinetic = 0
        simulation.root.input.model.unit_001.adsorption.mcl_qmax = [q_max]
        simulation.root.input.model.unit_001.adsorption.mcl_ka = [ka]
        simulation.root.input.model.unit_001.adsorption.mcl_kd = [1]

    #create an outlet
    simulation.root.input.model.unit_002.ncomp = nComp
    simulation.root.input.model.unit_002.unit_type = 'OUTLET'

    #set what values get saved
    # simulation.root.input['return'].unit_000.write_solution_inlet =  1
    simulation.root.input['return'].unit_002.write_solution_outlet = 1
    # simulation.root.input['return'].unit_001.write_solution_bulk =     1
    # simulation.root.input['return'].unit_001.write_solution_particle = 1
    # simulation.root.input['return'].unit_001.write_solution_solid =    1

    #set the times that the simulator writes out data for
    # if t_max > 1e5:
    #     n_t_pts = int(t_max/2) + 1
    # elif t_max > 1e6:
    #     n_t_pts = int(t_max/5) + 1
    # elif t_max > 1e7:
    #     n_t_pts = int(t_max/10) + 1
    # elif t_max > 1e8:
    #     n_t_pts = int(t_max/50) + 1
    # else:
    #     n_t_pts = int(t_max)+1
    # simulation.root.input.solver.user_solution_times = np.linspace(0, int(t_max), n_t_pts)
    simulation.root.input.solver.user_solution_times = np.linspace(0, int(t_max), int(t_max)+1)

    #solver settings
    set_solver_settings_2(simulation)
    simulation.root.input.solver.sections.nsec = 1
    simulation.root.input.solver.sections.section_times = [0.0, simulation.root.input.solver.user_solution_times[-1]]
    simulation.root.input.solver.sections.section_continuity = []

    return simulation





# Looking at isotherms__________________________________________________________

def get_mcl_isotherm(ka, q_max):
    linear_gradient_model = get_cadet_template(n_units=2)
    n_comp = 2
    Q = 1e-3

    # INLET
    linear_gradient_model.root.input.model.unit_000.unit_type = 'INLET'
    linear_gradient_model.root.input.model.unit_000.ncomp = n_comp
    linear_gradient_model.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    # CSTR
    linear_gradient_model.root.input.model.unit_001.unit_type = 'CSTR'
    linear_gradient_model.root.input.model.unit_001.ncomp = n_comp
    linear_gradient_model.root.input.model.unit_001.init_volume = 1e-3
    linear_gradient_model.root.input.model.unit_001.porosity = 0.694
    linear_gradient_model.root.input.model.unit_001.init_c = n_comp*[0]
    linear_gradient_model.root.input.model.unit_001.init_q = n_comp*[0]
    linear_gradient_model.root.input.model.unit_001.flow_rate_filter = Q

    # Sections and Switches
    linear_gradient_model.root.input.solver.sections.nsec = 1
    linear_gradient_model.root.input.solver.sections.section_times = [0.0, 50]

    linear_gradient_model.root.input.model.unit_000.sec_000.const_coeff = [0.0, 0.0]
    linear_gradient_model.root.input.model.unit_000.sec_000.lin_coeff = [0.0, 1]

    linear_gradient_model.root.input.model.connections.nswitches = 1
    linear_gradient_model.root.input.model.connections.switch_000.section = 0
    linear_gradient_model.root.input.model.connections.switch_000.connections = [0, 1, -1, -1, Q]

    adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
    adsorption_parameters = Dict()
    adsorption_parameters.is_kinetic = False
    adsorption_parameters.mcl_ka = [0, ka]
    adsorption_parameters.mcl_kd = [1, 1]
    adsorption_parameters.mcl_qmax = [1, q_max]
    linear_gradient_model.root.input.model.unit_001.nbound = [1, 1]
    linear_gradient_model.root.input.model.unit_001.adsorption_model = adsorption_model
    linear_gradient_model.root.input.model.unit_001.adsorption = adsorption_parameters

    run_simulation(linear_gradient_model)
    solution_bulk  = linear_gradient_model.root.output.solution.unit_001.solution_bulk
    solution_solid = linear_gradient_model.root.output.solution.unit_001.solution_solid
    return solution_bulk, solution_solid
