import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class data():
    def __init__(self, directory, file):
        self.file = file
        self.directory = directory
        self.name = self.directory + self.file
        self.load_and_clean()
        return

    def load_and_clean(self):
        self.df = pd.read_excel(self.name)
        units = list(self.df.loc[1])

        if list(self.df.iloc[1]).count('min') == 1:
            properties = self.get_names_normalized(['min'])
        elif list(self.df.iloc[1]).count('ml') == 1:
            properties = self.get_names_normalized(['ml'])
        else:
            properties = list(self.df.loc[0])
            for ind, prop in enumerate(properties):
                if ind%2 == 0:
                    name        = properties[ind][(properties[ind].rfind(':')+1):]
                    abbrev_name = name[(name.find('_')+1):]
                    units       = self.df.iat[1, ind]
                else:
                    units       = self.df.iat[1, ind][(self.df.iat[1, ind].find(' ')+1):]
                properties[ind] = abbrev_name + '_' + units

        self.df.columns = properties
        self.df.drop(labels=[0, 1], axis='index', inplace=True)
        self.df.reset_index(drop=True, inplace=True)
        return

    def get_names_normalized(self, names):
        for i in range(len(self.df.iloc[0])):
            if i != 0:
                names.append(self.df.iloc[0][i][self.df.iloc[0][i].find(':')+4:]\
                            + '_' + self.df.iloc[1][i][1:])
        return names

    def add_velocity_and_cv(self, v_col, d_col, flow_ind=500, flow=None):
        # Could maybe improve by distinguishing between cases and separating
        self.v_col = v_col # [ml]
        self.d_col = d_col # [cm]

        if flow is None:
            try:
                if 'P960_Flow_ml/min' in self.df.columns:
                    self.flow = self.df.at[flow_ind, 'P960_Flow_ml/min']
                elif 'Flow_ml/min' in self.df.columns:
                    self.flow = self.df.at[flow_ind, 'Flow_ml/min']
            except:
                print(f'Make sure the flow rate is provided in {self.name}, or enter the flow rate manually.')
        else:
            self.flow = flow # [ml/min]

        self.velocity = self.flow/(np.pi/4.0*self.d_col**2.0)*60.0
        for col_name in self.df.columns:
            if 'min' in col_name:
                self.df[col_name[:-4] + '_cv'] = self.df[col_name]*self.flow/self.v_col
            elif 'ml' in col_name:
                self.df[col_name[:-3] + '_cv'] = self.df[col_name]/self.v_col
        return




# Plotting functions ########################3

def align_yaxis_multiple(ax_list):
    y_lims = np.array([ax.get_ylim() for ax in ax_list])

    # force 0 to appear on both axes, comment if don't need
    y_lims[:, 0] = y_lims[:, 0].clip(None, 0)
    y_lims[:, 1] = y_lims[:, 1].clip(0, None)

    # normalize both axes
    y_mags = (y_lims[:,1] - y_lims[:,0]).reshape(len(y_lims),1)
    y_lims_normalized = y_lims / y_mags

    # find combined range
    y_new_lims_normalized = np.array([np.min(y_lims_normalized),
                                      np.max(y_lims_normalized)])

    # denormalize combined range to get new axes
    new_lim = y_new_lims_normalized * y_mags
    for ax, lim in zip(ax_list, new_lim):
        ax.set_ylim(lim)
    ax.axhline(color='k')
    return

def add_fractions(ax, df, x_conversion=1, text=False, xlim_times=None):
    if 'Fractions_min' in df.columns:
        frac_times = [x*x_conversion for x in df['Fractions_min'] if str(x) != 'nan']
    elif 'Fractions_ml' in df.columns:
        frac_times = [x*x_conversion for x in df['Fractions_ml'] if str(x) != 'nan']
    else:
        print('No x-coordinates found')
        return

    frac_names = [x for x in df['Fractions_(Fractions)'] if str(x) != 'nan']

    if xlim_times is not None:
        xlim_min = xlim_times[0]
        xlim_max = xlim_times[1]

    for (time, name) in zip(frac_times, frac_names):
        if name=='Waste': # hack to remove overlapping fraction cuts
            pass
        else:
            ax.axvline(x=time, ymin=0, ymax=0.05, color='red')
            if text==True and xlim_min < time < xlim_max:
                plt.text(time, 0.02, name, rotation=90, color='red', size=12, transform=ax.get_xaxis_transform())
    return
















##
