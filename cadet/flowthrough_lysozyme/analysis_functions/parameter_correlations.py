import numpy as np
from scipy import optimize
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import transforms

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


# Functions for parameter correlation___________________________________________

def my_fit_function_keq(tis, a, b, c):
    return c*np.exp(a*(tis)**b)

def my_fit_function_keq_2(tis, a, b):
    return a*tis**b

def my_fit_function_ds(keq, a, b):
    return np.log(a*keq**b)

def get_fit_param_inputs():
    data_file    = pd.ExcelFile('aex_pH_7_all_data.xlsx')
    protein_list = ['adh', 'blg', 'cat', 'ova']
    return data_file, protein_list

def get_correlation(resin):
    data_file, protein_list = get_fit_param_inputs()

    keq_fit_params = []
    Ds_fit_params = []
    if resin != 'phq':
        print('Need to inspect data to see if ovalbumin series must be treated \
        separately')

    for protein in protein_list:
        tab_name = protein + '_' + resin
        df = pd.read_excel(data_file, tab_name)
        df.dropna(inplace=True)
        df_ds_fit = df.copy()

        if protein != 'ova':
            fit_keq = optimize.curve_fit(my_fit_function_keq, df['IS (M)']*1e3,
                      df['Keq'], p0=(0.0, -1.0, 1.0e-2), maxfev=10000)
        else:
            fit_keq = optimize.curve_fit(my_fit_function_keq_2, df['IS (M)']*1e3,
                      df['Keq'], p0=(1.0e-2, 0.0), maxfev=10000)
            df_ds_fit.drop(2, inplace=True)

        keq_fit_params.append(fit_keq[0])
        fit_ds = optimize.curve_fit(my_fit_function_ds, df_ds_fit['Keq'],
                 np.log(df_ds_fit['Ds']), p0=(np.exp(-11), -1.5), maxfev=10000)
        Ds_fit_params.append(fit_ds[0])

    return keq_fit_params, Ds_fit_params

def get_column_params(resin, salt_c):
    keq_fit_params, Ds_fit_params = get_correlation(resin)

    keq_vals = [my_fit_function_keq(salt_c, keq_fit_params[i][0],
    keq_fit_params[i][1], keq_fit_params[i][2]) for i in [0, 1, 2]]

    keq_vals.append(my_fit_function_keq_2(salt_c, keq_fit_params[3][0],
    keq_fit_params[3][1])) # for ovalbumin

    ds_vals  = [np.exp(my_fit_function_ds(keq_vals[i], Ds_fit_params[i][0],
    Ds_fit_params[i][1])) for i in range(len(keq_fit_params))]

    return keq_vals, ds_vals


# # Plots_______________________________________________________________________
def plot_keq(resin):
    keq_fit_params, Ds_fit_params = get_correlation(resin)
    data_file, protein_list = get_fit_param_inputs()
    tis_new = np.linspace(50, 500, 1000)

    params = {'font.weight':'normal', 'font.size':20, 'figure.autolayout':True}
    plt.rcParams.update(params)
    fig, ax = plt.subplots()
    fig.set_size_inches(7, 6, forward=True)
    ax.set_xlabel('Ionic Strength [mM]')
    ax.set_ylabel(r'K$_{eq}$ [-]')
    ax.set_ylim(0.1, 100)

    for i, protein in enumerate(protein_list):
        tab_name = protein + '_' + resin
        df = pd.read_excel(data_file, tab_name)
        df.dropna(inplace=True)
        fit_params = keq_fit_params[i]

        if protein != 'ova':
            keq_new = [my_fit_function_keq(tis, fit_params[0], fit_params[1],
            fit_params[2]) for tis in tis_new]
        else:
            keq_new = [my_fit_function_keq_2(tis, fit_params[0],
            fit_params[1]) for tis in tis_new]

        ax.loglog(df['IS (M)']*1e3, df['Keq'], 'o', label=protein)
        ax.loglog(tis_new, keq_new, color=plt.gca().lines[-1].get_color())

    ax.legend(loc='best', frameon=False, handlelength=1.0)
    plt.show()
    return

def plot_ds(resin, fig=None, ax=None):
    keq_fit_params, Ds_fit_params = get_correlation(resin)
    data_file, protein_list = get_fit_param_inputs()
    keq_new = np.logspace(-1, 2)

    params = {'font.weight':'normal', 'font.size':20, 'figure.autolayout':True}
    plt.rcParams.update(params)

    if fig == None and ax == None:
        fig, ax = plt.subplots()

    fig.set_size_inches(7, 6, forward=True)
    ax.set_ylabel(r'D$_s$ [m$^2$/s]')
    ax.set_xlabel(r'K$_{eq}$ [-]')
    ax.set_ylim(1e-18, 1e-8)

    for i, protein in enumerate(protein_list):
        tab_name = protein + '_' + resin
        df = pd.read_excel(data_file, tab_name)
        df.dropna(inplace=True)
        fit_params = Ds_fit_params[i]
        ds_new = [np.exp(my_fit_function_ds(keq, fit_params[0],
        fit_params[1])) for keq in keq_new]

        ax.loglog(df['Keq'],  df['Ds'], 'o', label=protein)
        ax.loglog(keq_new, ds_new, color=plt.gca().lines[-1].get_color())

    ax.legend(loc='best', frameon=False, handlelength=1.0)
    # plt.show()
    return
