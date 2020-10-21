import pandas as pd
import numpy as np
from ast import literal_eval


def labels_df(elem_symbol):
    acid = elem_symbol + '_acid'
    basic = elem_symbol + '_basic'
    acid_state = acid + '_state'
    basic_state = basic + '_state'
    acid_label = elem_symbol + '_acid_label'
    basic_label = elem_symbol + '_basic_label'
    acid_bal = elem_symbol + '_acid_bal'
    basic_bal = elem_symbol + '_basic_bal'
    acid_E0 = acid + '_E0'
    basic_E0 = basic + '_E0'
    acid_E = acid + '_E'
    basic_E = basic + '_E'
    acid_new = acid + '_new'
    basic_new = basic + '_new'
    return (acid, basic, acid_state, basic_state, acid_label, basic_label,
            acid_bal, basic_bal, acid_E0, basic_E0, acid_E, basic_E, acid_new,
            basic_new)


def coef_ac(lst):
    a, b, c, d, m = lst[:]
    x = 1
    z = c
    w = b / d
    y = 2 * z - a
    n = m + y
    return [x, y, n, w, z]


def coef_bas(lst):
    a, b, c, d, m = lst[:]
    x = 1
    w = b / d
    z = c
    r = 2 * z - a
    s = r
    n = m + s
    return [x, r, n, w, z, s]


def potentials(elem_symbol, data, temperature=298.15, conc_ion=1, pH=0):
    gas_constant = 8.3144598
    faraday_constant = 96485.33289

    df = data.filter(like=elem_symbol).dropna(how='all')

    acid, basic, acid_state, basic_state, acid_label, basic_label, acid_bal, \
        basic_bal = labels_df(elem_symbol)[:8]

    df[acid + '_E0'] = np.where(df[acid] != 0, df[acid] / df[acid].index, 0)
    df[basic + '_E0'] = np.where(df[basic] != 0, df[basic] / df[basic].index, 0)

    acid_mask = pd.notnull(df[acid_bal])
    df_acid = df.filter(like=acid)[acid_mask]
    df_acid[acid_bal] = df_acid[acid_bal].apply(literal_eval)

    basic_mask = pd.notnull(df[basic_bal])
    df_basic = df.filter(like=basic)[basic_mask]
    df_basic[basic_bal] = df_basic[basic_bal].apply(literal_eval)

    hydronium = 10 ** (-pH)
    hydroxide = 10 ** (-(14 - pH))

    E_ac = []
    for i in df_acid[acid + '_E0'].index:
        if i == 0:
            E_ac.append(0)
        else:
            x, y, n, w, z = coef_ac(df_acid[acid_bal][i])
            if df_acid[acid_state][0] == '(s)' or df_acid[acid_state][0] == '(l)':
                w = 0
            if df_acid[acid_state][i] == '(s)' or df_acid[acid_state][i] == '(l)':
                x = 0
            nernst = ((gas_constant * temperature) / (n * faraday_constant)) * np.log((conc_ion) ** w / (conc_ion ** x * hydronium ** y))
            E0_ac = df_acid[acid + '_E0'][i]
            E_ac.append(E0_ac + nernst)

    E_bas = []
    for i in df_basic[basic + '_E0'].index:
        if i == 0:
            E_bas.append(0)
        else:
            x, r, n, w, z, s = coef_bas(df_basic[basic_bal][i])
            if df_basic[basic_state][0] == '(s)' or df_basic[basic_state][0] == '(l)':
                w = 0
            if df_basic[basic_state][i] == '(s)' or df_basic[basic_state][i] == '(l)':
                x = 0
            nernst = ((gas_constant * temperature) / (n * faraday_constant)) * np.log(((conc_ion ** w) * (hydroxide ** s)) / (conc_ion ** x))
            E0_bas = df_basic[basic + '_E0'][i]
            E_bas.append(E0_bas + nernst)

    df_acid[acid + '_E'] = E_ac
    df_basic[basic + '_E'] = E_bas

    df_acid[acid + '_new'] = df_acid[acid + '_E'] * df_acid.index
    df_basic[basic + '_new'] = df_basic[basic + '_E'] * df_basic.index

    if pH <= 7:
        return df_acid
    return df_basic


def get_text_positions(x_data, y_data, txt_width, txt_height):
    # code from https://stackoverflow.com/questions/8850142/matplotlib-overlapping-annotations/10739207
    a = zip(y_data, x_data)
    text_positions = y_data.copy()
    for index, (y, x) in enumerate(a):
        local_text_positions = [i for i in a if i[0] > (y - txt_height)
                                and (abs(i[1] - x) < txt_width * 2) and i != (y, x)]
        if local_text_positions:
            sorted_ltp = sorted(local_text_positions)
            if abs(sorted_ltp[0][0] - y) < txt_height:  # True == collision
                differ = np.diff(sorted_ltp, axis=0)
                a[index] = (sorted_ltp[-1][0] + txt_height, a[index][1])
                text_positions[index] = sorted_ltp[-1][0] + txt_height
                for k, (j, m) in enumerate(differ):
                    # j is the vertical distance between words
                    if j > txt_height * 2:  # if True then room to fit a word in
                        a[index] = (sorted_ltp[k][0] + txt_height, a[index][1])
                        text_positions[index] = sorted_ltp[k][0] + txt_height
                        break
    return text_positions


def text_plotter(x_data, y_data, label, text_positions, axis, txt_width,
                 txt_height, solution, xmove=0, ymove=-0.05):
    # adapted from https://stackoverflow.com/questions/8850142/matplotlib-overlapping-annotations/10739207
    if solution == 'acid':
        color = 'blue'
        facecolor = 'dodgerblue'
    else:
        color = 'red'
        facecolor = 'darkorange'

    for x, y, l, t in zip(x_data, y_data, label, text_positions):
        axis.text(x + txt_width + xmove, t + ymove, '$\mathrm{{{0}}}$'.format(label[x]),
                  rotation=0, color=color, fontsize=20,
                  bbox=dict(facecolor=facecolor, alpha=0.1))
        if y != t:
            axis.arrow(x, t, 0, y - t, color='red', alpha=0.3, width=txt_width * 0.1,
                       head_width=txt_width, head_length=txt_height * 0.5,
                       zorder=0, length_includes_head=True)


def plot_param(ax=None):
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


def plot_frost(elem_symbol, data, ax=None, ac_xmove=0, ac_ymove=-0.05,
               bas_xmove=0, bas_ymove=-0.05, xbbox=1.2, ybbox=0.5, label=True,
               plotboth=True, pH=0, conc_ion=1, temperature=298.15, legend=True):
    plot_param(ax)

    acid, basic, acid_state, basic_state, acid_label, basic_label, acid_bal, \
        basic_bal, acid_E0, basic_E0, acid_E, basic_E, acid_new, basic_new = labels_df(elem_symbol)

    if pH <= 7 or plotboth:
        if plotboth:
            pH = 0
            conc_ion = 1
            df = potentials(elem_symbol, data, temperature=298.15, conc_ion=conc_ion, pH=pH)
            ax.plot(df.index, df[acid], linewidth=3,
                    label='$\mathrm{{pH = {{{0}}}, [ions] = {{{1:1.1f}}} \, mol \cdot L^{{-1}}}}$'.format(pH, conc_ion),
                    color="blue", marker='o', markersize=7)
        else:
            df = potentials(elem_symbol, data, temperature=298.15, conc_ion=1, pH=0)
            ax.plot(df.index, df[acid], linewidth=3,
                    label='$\mathrm{{pH = 0, [ions] = 1 \, mol \cdot L^{{-1}}}}$'.format(pH, conc_ion),
                    color="blue", marker='o', markersize=7, linestyle='-.')
            if pH != 0 or conc_ion != 1 or temperature != 298.15:
                df = potentials(elem_symbol, data, temperature=temperature, conc_ion=conc_ion, pH=pH)
                ax.plot(df.index, df[acid_new], linewidth=3,
                        label='$\mathrm{{pH = {{{0}}}, [ions] = {{{1:1.1e}}} \, mol \cdot L^{{-1}}}}$'.format(pH,
                                                                                                              conc_ion),
                        color="cyan", marker='s', markersize=7)

        if label:
            txt_height = 0.04 * (ax.get_ylim()[1] - ax.get_ylim()[0])
            txt_width = 0.03 * (ax.get_ylim()[1] - ax.get_ylim()[0])

            text_positions = get_text_positions(df.index, df[acid], txt_width, txt_height)
            text_plotter(df.index, df[acid], df[acid_label], text_positions,
                         ax, txt_width, txt_height, 'acid', xmove=ac_xmove, ymove=ac_ymove)

    if pH >= 7 or plotboth:
        if plotboth:
            pH = 14
            conc_ion = 1
            df = potentials(elem_symbol, data, temperature=298.15, conc_ion=conc_ion, pH=pH)
            ax.plot(df.index, df[basic], linewidth=3,
                    label='$\mathrm{{pH = {{{0}}}, [ions] = {{{1:1.1f}}} \, mol \cdot L^{{-1}}}}$'.format(pH, conc_ion),
                    color="red", marker='s', markersize=7, linestyle='-.')
        else:
            df = potentials(elem_symbol, data, temperature=298.15, conc_ion=1, pH=14)
            ax.plot(df.index, df[basic], linewidth=3,
                    label='$\mathrm{{pH = 14, [ions] = 1 \, mol \cdot L^{{-1}}}}$'.format(pH, conc_ion),
                    color="red", marker='o', markersize=7, linestyle='-.')
            if pH != 14 or conc_ion != 1 or temperature != 298.15:
                df = potentials(elem_symbol, data, temperature=temperature, conc_ion=conc_ion, pH=pH)
                ax.plot(df.index, df[basic_new], linewidth=3,
                        label='$\mathrm{{pH = {{{0}}}, [ions] = {{{1:1.1e}}} \, mol \cdot L^{{-1}}}}$'.format(pH,
                                                                                                              conc_ion),
                        color="orange", marker='s', markersize=7)

        if label:
            txt_height = 0.04 * (ax.get_ylim()[1] - ax.get_ylim()[0])
            txt_width = 0.03 * (ax.get_ylim()[1] - ax.get_ylim()[0])

            text_positions = get_text_positions(df.index, df[basic], txt_width, txt_height)
            text_plotter(df.index, df[basic], df[basic_label], text_positions,
                         ax, txt_width, txt_height, 'basic', xmove=bas_xmove, ymove=bas_ymove)

    if legend:
        ax.legend(fontsize=14, loc='upper center', shadow=True, fancybox=True,
                  bbox_to_anchor=(xbbox, ybbox), ncol=1)
    ax.set_title('Frost diagram - {0}'.format(elem_symbol), size=20)
    ax.set_xticks(np.arange(round(ax.get_xlim()[0]), round(ax.get_xlim()[1]) + 1, 1))
