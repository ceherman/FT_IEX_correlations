from analysis_functions.cadet_imports import *

def get_ka(keq, q_max):
    return keq/q_max

def get_ds(keq, a, b):
    return a*keq**b


###########################################################################
# When will breakthrough occur?

# System parameters
fVolumetric = 1.633e-8/4.0 # m**3/s
V_col = np.pi/4 * (0.5**2) * 1e-4 * 0.042 # [m3]
eps_t = 0.694
eps_c = 0.49

def get_t0():
    return V_col*eps_t/fVolumetric

def get_min_t0():
    return V_col*eps_c/fVolumetric

def get_t_saturation(c_load, q_max):
    return V_col*(1-eps_t)*q_max/fVolumetric/c_load

def get_t_r(keq):
    return V_col/fVolumetric*(keq*(1-eps_t) + eps_t)

def get_t_r_min(keq):
    # Actually larger than t_r above for the given numbers
    return V_col/fVolumetric*(keq*(1-eps_c) + eps_c)

def get_keq_tr_tsat(q_max, c_load):
    return q_max/c_load - (eps_t/(1-eps_t))

def get_keq_lin(q_max, c_load, q_over_qlin):
    return q_max/c_load*(1/q_over_qlin - 1)

def get_min_t(c_load, keq, q_max):
    t_sat = get_t_saturation(c_load, q_max)
    t_r   = get_t_r(keq)
    t_min = min(t_sat, t_r)
    return t_min

def get_max_sim_time(c_load, keq, q_max):
    t_min = get_min_t(c_load, keq, q_max)
    if t_min == get_t_saturation(c_load, q_max):
        t_max = t_min + 500
    else:
        t_max = 4*t_min + 5000
    return t_max

##############################################################################
# Looking at results
def get_sim_results_file(results_folder, file):
    sim = Cadet()
    sim.filename = os.path.join(results_folder, file)
    sim.load()
    return sim

def get_sim_results(results_folder, file, flow_flag=False):
    sim    = get_sim_results_file(results_folder, file)
    c_load = sim.root.input.model.unit_000.sec_000.const_coeff[0]
    dp     = sim.root.input.model.unit_001.par_diffusion[0]
    t      = sim.root.output.solution.solution_times
    c_out  = sim.root.output.solution.unit_002.solution_outlet_comp_000
    fVolum = sim.root.input.model.connections.switch_000.connections[4]

    try:
        ka     = sim.root.input.model.unit_001.adsorption.mcl_ka[0]
        q_max  = sim.root.input.model.unit_001.adsorption.mcl_qmax[0]
        kd     = sim.root.input.model.unit_001.adsorption.mcl_kd[0]
        keq    = q_max*ka/kd
    except TypeError:
        ka    = sim.root.input.model.unit_001.adsorption.lin_ka[0]
        q_max = None
        kd    = sim.root.input.model.unit_001.adsorption.lin_kd[0]
        keq   = ka/kd
    if flow_flag:
        return c_load, dp, ka, q_max, kd, keq, t, c_out, fVolum
    else:
        return c_load, dp, ka, q_max, kd, keq, t, c_out,

def plot_outlet(t, c_out, c_load, keq, mass, ax=None):
    if ax == None:
        params = {'font.weight':'normal', 'font.size':20, 'figure.autolayout':True}
        plt.rcParams.update(params)
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 5, forward=True)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel(r'c$_{out}$/c$_{feed}$')
    ax.plot(t, c_out/c_load, label=f'K$_{{eq}}$ = {keq:.1e}, c = {mass*c_load:.1e} mg/ml')
    return

################################################################################
# Moments
def calc_first_moment(x, y):
    """Note:  no correction included for baseline drift."""
    x  = np.array(x)
    y  = np.array(y)
    y_num = x*y
    return np.trapz(y_num, x)/np.trapz(y, x)

def calc_k_prime(t_r, t_0):
    return (t_r - t_0)/t_0


###############################################################################

def get_linear_v(Q, D=5.0e-3):
    return Q/(np.pi/4.0*D**2.0)

def get_t_star(t, v_lin, l_col):
    return t*v_lin/l_col

def get_peclet(d_part, v_lin, Dp):
    return d_part*v_lin/Dp
