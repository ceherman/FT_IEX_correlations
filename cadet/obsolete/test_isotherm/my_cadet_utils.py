from cadet_imports import *

# Function for running CADET____________________________________________________

def run_simulation(cadet, file_name=None):
    try:
        os.mkdir('./temp_files')
    except:
        pass

    if file_name is None:
        f = next(tempfile._get_candidate_names())
        cadet.filename = os.path.join(Path().absolute(), 'temp_files', f + '.h5')
    else:
        cadet.filename = os.path.join(Path().absolute(), 'temp_files', file_name + '.h5')
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
        os.remove(os.path.join(Path().absolute(), 'temp_files', f + '.h5'))

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


# Misc functions________________________________________________________________

def cond_2_tis(cond):
    tis = (0.0229207306*(cond**2) + 9.3999636424*cond + 21.2424910935)
    return tis


# Looking at isotherms__________________________________________________________

def get_isotherm(ka, q_max):
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
