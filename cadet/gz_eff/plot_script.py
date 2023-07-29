import platform
from pathlib import Path
import os
import subprocess
import time
import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import exploration_helper as explor
import gc


results_folder = Path().absolute() / 'sim_results_updated_PREP_2021_transport_effect'
image_folder = Path().absolute() / 'images_updated_PREP_2021'

sims = os.listdir(results_folder)
sims.sort()

mass       = 50
mg_per_ml_load = np.logspace(-3, 1, 5)
c_load_vals    = mg_per_ml_load/mass # [mol m-3]

v_col = 0.042 * np.pi/4 * (0.5**2) * 1e-4

dp_vals = [6.0e-12, 2.0e-11]
fVolum_vals = [1.633e-8, 1.633e-8/4]

for dp in dp_vals:
    for fVolum in fVolum_vals:
        params = {'font.weight':'normal', 'font.size':22, 'figure.autolayout':True}
        plt.rcParams.update(params)
        fig, ax = plt.subplots()
        ax.set_xlabel('CV [-]')
        ax.set_ylabel(r'c$_{out}$/c$_{load}$ [-]')
        for file in sims:
            try:
                c_load, this_dp, ka, q_max, kd, keq, t, c_out, this_fVolum =\
                explor.get_sim_results(results_folder, file, True)

                if this_dp == dp and this_fVolum == fVolum:
                    ax.plot(t*fVolum/v_col, c_out/c_load)
            except:
                print(file)

        ax.set_xscale('log')
        fig.patch.set_facecolor('white')
        fig.set_size_inches(10, 6, forward=True)
        ax.set_xlim(0.1, 1e5)

        fig.savefig(os.path.join(image_folder, f'cv_transport_effect_dp_{dp:.1e}_Q_{fVolum:.2e}.png'), dpi=300)

        plt.show()
        fig.clf()
        plt.close()
        gc.collect()
