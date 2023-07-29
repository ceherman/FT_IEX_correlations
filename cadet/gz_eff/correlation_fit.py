from cadet_imports import *

import my_cadet_utils as cad_utils
import parameter_correlations as corr
import exploration_helper as explor

import multiprocessing as mp
import gc

from scipy import interpolate
from scipy import optimize
import time


def local_run_sim_fun(args):
    sim, res_file = args
    print(res_file[res_file.rfind('/'):])

#     run_time = time.time()
    cad_utils.run_simulation(sim, res_file)
#     rtime = (time.time()-run_time)/60.0
    return

def solve_time(time, frac, f):
    return f(time) - frac

def get_sim_data(args):
    file_name, results_folder = args
    c_load, c_out, t, fVolum, Dp, d_part, l_col, Ds =\
    explor.get_sim_results_more_general(results_folder, file_name)

    d_col = 0.5e-2
    area_col = np.pi/4 * d_col**2
    v_col = area_col * l_col
    norm_c = c_out/c_load
    f = interpolate.interp1d(t, norm_c, kind='cubic')

    these_res = [fVolum, Dp, d_part, l_col, Ds]
    fracs = np.arange(0.01, 0.81, 0.01)

    for frac in fracs:
        idx = np.searchsorted(norm_c, frac, side="left")
        t_guess = t[idx]
        ftime = optimize.fsolve(solve_time, t_guess, args=(frac, f))[0]
        v_break  = ftime * fVolum
        cv_break = v_break/v_col
        fun_pe   = (cv_break - eps_c)/keq
        these_res.append(fun_pe)
#     print(file_name)
    return these_res

def get_xval(fVolum, d_part, Dp, l_col, eps_c, area_col):
    return fVolum*d_part**2/(Dp*l_col)/(eps_c * area_col)

def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))

def get_sim_data_keq_slope(args):
    file_name, results_folder = args
    c_load, dp, ka, q_max, kd, keq, t, c_out, fVolum, Ds, d_part, l_col =\
    explor.get_sim_results(results_folder, file_name)

    d_col = 0.5e-2
    area_col = np.pi/4 * d_col**2
    v_col = area_col * l_col
    norm_c = c_out/c_load
    f = interpolate.interp1d(t, norm_c, kind='cubic')

    these_res = [fVolum, Dp, d_part, l_col, Ds, keq]
    fracs = np.arange(0.01, 0.81, 0.01)

    for frac in fracs:
        idx = np.searchsorted(norm_c, frac, side="left")
        t_guess = t[idx]
        ftime = optimize.fsolve(solve_time, t_guess, args=(frac, f))[0]
        v_break  = ftime * fVolum
        cv_break = v_break/v_col
        fun_pe   = (cv_break - eps_c)/keq
        these_res.append(fun_pe)
#     print(file_name)
    return these_res

def get_sim_data_keq(args):
    file_name, results_folder = args
    c_load, Dp, ka, q_max, kd, keq, t, c_out, fVolum, Ds, d_part, l_col =\
    explor.get_sim_results(results_folder, file_name)

    d_col = 0.5e-2
    area_col = np.pi/4 * d_col**2
    v_col = area_col * l_col
    norm_c = c_out/c_load
    f = interpolate.interp1d(t, norm_c, kind='cubic')

    these_res = [fVolum, Dp, d_part, l_col, Ds, keq]
    fracs = np.arange(0.01, 0.81, 0.01)

    for frac in fracs:
        idx = np.searchsorted(norm_c, frac, side="left")
        t_guess = t[idx]
        ftime = optimize.fsolve(solve_time, t_guess, args=(frac, f))[0]
        v_break  = ftime * fVolum
        cv_break = v_break/v_col
        these_res.append(cv_break)
#     print(file_name)
    return these_res

def get_D_combo_2(e_vals, Ds, Dp):
    w1 = 10.0**(e_vals[0])
    w2 = 10.0**(e_vals[1])
    return w1*Ds + w2*Dp

def mod_x_2(e_vals, df):
    df['D_combo'] = get_D_combo_2(e_vals, df['Ds'], df['Dp'])
    df['new_x'] = (df['v'] * df['d_part']**2) / (df['D_combo'] * df['l_col'])
    return df

def get_residuals_2(e_vals, df):
    df = mod_x_2(e_vals, df)
    df.sort_values(by='new_x', inplace=True)
    spl = interpolate.UnivariateSpline(df['new_x'], df['0.50']/(1.0-eps_t))
    spl.set_smoothing_factor(0.75)
    df['residual'] = df['0.50']/(1.0-eps_t) - spl(df['new_x'])
    return np.linalg.norm(df['residual'].values)

results_folder = Path().absolute() / 'sim_results_more_comprehensive'
if not os.path.exists( results_folder.as_posix() ):
    os.makedirs( results_folder.as_posix() )

image_folder = Path().absolute() / 'images_more_comprehensive_sims'
if not os.path.exists( image_folder.as_posix() ):
    os.makedirs( image_folder.as_posix() )

eps_c = 0.49
eps_p = 0.4
eps_t = eps_c + eps_p*(1.0 - eps_c)

mass       = 50
cap_mg     = 100 # [mg/ml column]
q_m_fac    = cap_mg/(1.0-eps_t)
q_max      = q_m_fac/mass # [mol m-3 resin]

d_col = 0.5e-2
area_col = np.pi/4 * d_col**2

a_pnas, b_pnas = (1.6598739258839034e-12, -0.23974534290342786)

fVolum_vals = np.linspace(1.633e-8/3, 1.633e-8*2/3, 3) # for d_c = 0.5 cm, corresponds to 100-200 cm/h
Dp_vals     = np.linspace(5.0e-12, 4.0e-11, 3)
d_part_vals = np.linspace(5.0e-6, 200.0e-6, 4)
l_col_vals  = np.linspace(5.0e-2, 50.0e-2, 4)
# keq_vals    = np.concatenate((np.linspace(0.01, 0.81, 5), np.logspace(0, 4, 50), np.linspace(2e4, 2e4, 1)))
# load_c_vals = np.array([2.0e-5, 2.0e-4, 2.0e-3, 2.0e-2, 2.0e-1])

keq_vals = np.array([6.86648845e+03])
load_c_vals = np.array([2.0e-5])

















def get_cv_fun(params, df):
    a, b = params
    df['keq_fun'] = df['keq'] + a * np.log(df['keq']) + b * np.log(df['keq'])**2
    return df

def get_D_combo_2(e_vals, Ds, Dp):
    w1 = 10.0**(e_vals[0])
    w2 = 10.0**(e_vals[1])
    return w1*Ds + w2*Dp

def mod_x_2(e_vals, v, d_part, Dp, Ds, l_col):
    df['D_combo'] = get_D_combo_2(e_vals, df['Ds'], df['Dp'])
    df['new_x'] = (df['v'] * df['d_part']**2) / (df['D_combo'] * df['l_col'])
    return df



df = pd.read_csv('./comprehensive_results_1.csv')

params_list = []
for fVolum in fVolum_vals:
    for Dp in Dp_vals:
        for d_part in d_part_vals:
            for l_col in l_col_vals:
                df_temp = df[(df.Dp - Dp < 1e-13) & (df.d_part == d_part) \
                                & (df.l_col == l_col) & (df.fVolum == fVolum)].copy()
                df_temp = get_cv_fun([-1.27139371, -1.4010829], df_temp)
                df_temp.sort_values(by='keq_fun', inplace=True)
                m = np.polyfit(df_temp['keq_fun'], df_temp['0.01'], 1)[0]

                v = list(set(df_temp['v']))
                assert len(v) == 1
                v = v[0]

                params_list.append([v, d_part, Dp, ])

# fit = optimize.minimize(get_cv_fun_res_group_2, x0=[1, 1], args=(dfs_list,))

st_time = time.time()
bounds = [(0.5, 2.0), (-10.0, 10.0)]
fit = optimize.differential_evolution(get_cv_fun_res_group_2,
                                      bounds, args=(dfs_list,), disp=True,
                                      workers=1, tol=1e-12)
print(fit)
run_time = (time.time() - st_time)/60.0
print(f'{run_time:.2f} min')














#
