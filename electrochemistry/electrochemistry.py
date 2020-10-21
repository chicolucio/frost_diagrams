import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ast import literal_eval
from collections import namedtuple
from helpers.plot_annotations import Annotations
from helpers.core import columns_df, coef_ac, coef_bas

DATA = pd.read_csv("data/frost.csv")
DATA = DATA.set_index('charge')

GAS_CONSTANT = 8.3144598
FARADAY_CONSTANT = 96485.33289


class Electrochemistry:

    def __init__(self, elem_symbol):
        self.elem_symbol = elem_symbol
        self._column_names = columns_df(elem_symbol)
        self._df = self.element_data

    @property
    def element_data(self):
        df = DATA.filter(like=self.elem_symbol).dropna(how='all')

        acid = self._column_names['acid']
        basic = self._column_names['basic']
        acid_bal = self._column_names['acid_bal']
        basic_bal = self._column_names['basic_bal']

        df[acid + '_E0'] = np.where(df[acid] != 0, df[acid] / df[acid].index, 0)
        df[basic + '_E0'] = np.where(df[basic] != 0, df[basic] / df[basic].index, 0)

        acid_mask = pd.notnull(df[acid_bal])
        df_acid = df.filter(like=acid)[acid_mask]
        df_acid[acid_bal] = df_acid[acid_bal].apply(literal_eval)

        basic_mask = pd.notnull(df[basic_bal])
        df_basic = df.filter(like=basic)[basic_mask]
        df_basic[basic_bal] = df_basic[basic_bal].apply(literal_eval)

        Data = namedtuple("element_data", ["full", "acid", "basic"])

        return Data(df, df_acid, df_basic)

    def _nernst_acid(self, pH, conc_ion, temperature):
        hydronium = 10 ** (-pH)
        acid = self._column_names['acid']
        acid_bal = self._column_names['acid_bal']
        acid_state = self._column_names['acid_state']
        E_ac = []
        df = self._df.acid

        for i in df[acid + '_E0'].index:
            if i == 0:
                E_ac.append(0)
            else:
                x, y, n, w, z = coef_ac(df[acid_bal][i])
                if df[acid_state][0] == '(s)' or df[acid_state][0] == '(l)':
                    w = 0
                if df[acid_state][i] == '(s)' or df[acid_state][i] == '(l)':
                    x = 0
                nernst = ((GAS_CONSTANT * temperature) / (n * FARADAY_CONSTANT)) * np.log(conc_ion ** w / (
                        conc_ion ** x * hydronium ** y))
                E0_ac = df[acid + '_E0'][i]
                E_ac.append(E0_ac + nernst)

        df[acid + '_E'] = E_ac
        df[acid + '_new'] = df[acid + '_E'] * df.index
        return df

    def _nernst_basic(self, pH, conc_ion, temperature):
        hydroxide = 10 ** (-(14 - pH))
        basic = self._column_names['basic']
        basic_bal = self._column_names['basic_bal']
        basic_state = self._column_names['basic_state']
        df = self._df.basic
        E_bas = []
        for i in df[basic + '_E0'].index:
            if i == 0:
                E_bas.append(0)
            else:
                x, r, n, w, z, s = coef_bas(df[basic_bal][i])
                if df[basic_state][0] == '(s)' or df[basic_state][0] == '(l)':
                    w = 0
                if df[basic_state][i] == '(s)' or df[basic_state][i] == '(l)':
                    x = 0
                nernst = ((GAS_CONSTANT * temperature) / (n * FARADAY_CONSTANT)) * np.log(
                    ((conc_ion ** w) * (hydroxide ** s)) / (conc_ion ** x))
                E0_bas = df[basic + '_E0'][i]
                E_bas.append(E0_bas + nernst)

        df[basic + '_E'] = E_bas
        df[basic + '_new'] = df[basic + '_E'] * df.index
        return df

    def nernst(self, pH=0, conc_ion=1, temperature=298.15):

        if pH <= 7:
            df = self._nernst_acid(pH, conc_ion, temperature)
        else:
            df = self._nernst_basic(pH, conc_ion, temperature)
        return df

    def _plot_param(self, ax=None):
        ax.axhline(color='gray', zorder=-1)
        ax.axvline(color='gray', zorder=-1)
        ax.grid(b=True, which='major', linestyle=':', linewidth=2)
        ax.minorticks_on()
        ax.grid(b=True, which='minor', axis='y', linestyle=':', linewidth=1.0)
        ax.tick_params(which='both', labelsize=16)
        ax.tick_params(which='minor', axis='x', bottom=False)
        ax.set_aspect('equal')  # so that no scale adjust is needed
        ax.set_xlabel('Oxidation number', size=18)
        ax.set_ylabel(r'$\Delta G/F$', size=18)
        ax.tick_params(axis='both', length=6, which='major', width=1.5)
        ax.tick_params(axis='both', length=3, which='minor', width=1.0)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.5)

    def _plot_standard_curves(self, solution=None, chem_labels=True):
        ax = plt.gca()
        df_acid = self.nernst()
        df_basic = self.nernst(14)
        acid_plot = {'linewidth': 3,
                     'label': fr'$\mathrm{{pH = 0, [ions] = 1.0 \, mol \cdot L^{{-1}}}}$',
                     'color': "blue",
                     'marker': 'o',
                     'markersize': 7, }
        basic_plot = {'linewidth': 3,
                      'label': fr'$\mathrm{{pH = 14, [ions] = 1.0 \, mol \cdot L^{{-1}}}}$',
                      'color': "red",
                      'marker': 's',
                      'markersize': 7,
                      'linestyle': '-.', }

        acid = self._column_names['acid']
        acid_label = self._column_names['acid_label']
        basic = self._column_names['basic']
        basic_label = self._column_names['basic_label']

        if solution == 'both':
            ax.plot(df_acid.index, df_acid[acid], **acid_plot)
            ax.plot(df_basic.index, df_basic[basic], **basic_plot)
        elif solution == 'acid':
            ax.plot(df_acid.index, df_acid[acid], **acid_plot)
        elif solution == 'basic':
            ax.plot(df_basic.index, df_basic[basic], **basic_plot)
        else:
            raise ValueError('Solution must be acid, basic or both')

        if chem_labels and (solution == 'acid' or solution == 'both'):
            Annotations(df_acid.index, df_acid[acid], df_acid[acid_label],
                        color='blue', facecolor='dodgerblue').text_plotter()
        if chem_labels and (solution == 'basic' or solution == 'both'):
            Annotations(df_basic.index, df_basic[basic], df_basic[basic_label],
                        color='red', facecolor='darkorange').text_plotter()

    def _plot_not_standard_curves(self, pH=0, conc_ion=1, temperature=298.15, chem_labels=True):
        ax = plt.gca()
        df = self.nernst(temperature=temperature, conc_ion=conc_ion, pH=pH)

        acid_label = self._column_names['acid_label']
        basic_label = self._column_names['basic_label']

        if pH <= 7:
            color = 'cyan'
            column = self._column_names['acid_new']
        else:
            color = 'orange'
            column = self._column_names['basic_new']

        ax.plot(df.index, df[column], linewidth=3,
                label=fr'$\mathrm{{pH = {pH}, [ions] = {{{conc_ion:1.1E}}} \, mol \cdot L^{{-1}}}}$',
                color=color, marker='s', markersize=7)

        if chem_labels and pH <= 7:
            Annotations(df.index, df[column], df[acid_label], color='blue', facecolor='dodgerblue').text_plotter()
        if chem_labels and pH > 7:
            Annotations(df.index, df[column], df[basic_label], color='red', facecolor='darkorange').text_plotter()

    def plot_frost(self, ax=None, chemical_labels=True, plot_standards=True, pH=0, conc_ion=1, temperature=298.15,
                   legend=True):

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8), facecolor=(1.0, 1.0, 1.0))
        self._plot_param(ax)

        if plot_standards:
            self._plot_standard_curves(solution='both', chem_labels=chemical_labels)

        elif pH <= 7:
            self._plot_standard_curves(solution='acid', chem_labels=chemical_labels)
            if pH != 0 or conc_ion != 1 or temperature != 298.15:
                self._plot_not_standard_curves(pH, conc_ion, temperature, chem_labels=False)

        elif pH > 7:
            self._plot_standard_curves(solution='basic', chem_labels=chemical_labels)
            if pH != 14 or conc_ion != 1 or temperature != 298.15:
                self._plot_not_standard_curves(pH, conc_ion, temperature, chem_labels=False)

        if legend:
            xbbox = 0.5
            ybbox = -0.15
            ax.legend(fontsize=14, loc='upper center', shadow=True, fancybox=True,
                      bbox_to_anchor=(xbbox, ybbox), ncol=1)
        ax.set_title('Frost diagram - {0}'.format(self.elem_symbol), size=20)
        ax.set_xticks(np.arange(round(ax.get_xlim()[0]), round(ax.get_xlim()[1]) + 1, 1))
