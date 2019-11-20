import csv
import os
from collections import defaultdict
import itertools

from matplotlib import pyplot as plt
from matplotlib.transforms import BboxBase
from sklearn.linear_model import LinearRegression
from data import fe, column, element, cohesive_energies, N2_engs, N2H_engs, d_band, electronegativity,d_cohesive, s_cohesive, plus_4_fe
import matplotlib
from data import fe_dict, d_band
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
import numpy as np

from ase.units import kJ, mol


plt.rcParams["figure.figsize"] = (4.5,3.5)
plt.rcParams["font.size"] = 10
plt.rcParams['legend.fontsize'] =  9
plt.rcParams['legend.handlelength'] = 1


metal_dict = defaultdict(dict)
path_counter = 0
compound = []
for pathway in os.listdir('../data/corrected_data'):
    if '_' in pathway:
        if 'sa' not in pathway:
            continue
    with open('../data/corrected_data/' + pathway, 'r') as f:
        compound.append(pathway.split('.')[0])
        data = csv.reader(f)
        for row in data:
            metal_dict[pathway.split('.')[0]][row[0]] = float(row[1])
        # figure out how long one of them is

print(metal_dict.keys())
# make a dummy comparison entry
#all_metals = set(list(metal_dict['NH2'].keys()) + list(metal_dict['NH3'].keys()))
#shared_metals = [a for a in all_metals if a in metal_dict[s1].keys() and a in metal_dict[s2].keys()]

#reaction_data = [float(a) for a in reaction_data]
# build a matrix 

def build_lists(species, return_electronegativity=False):
    d_band_center = []
    cohesive_energy = []
    bindings = []
    metals = []
    electro_n = []
    for metal, binding in metal_dict[species].items():
        #if metal in ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'zn']:
        #    continue
        if binding is None:
            continue
        d_band_center.append(d_band[metal])
        cohesive_energy.append(d_cohesive[metal])
        bindings.append(float(binding))
        metals.append(metal)
        electro_n.append(electronegativity[metal])
    if return_electronegativity:
        return bindings, cohesive_energy, d_band_center, metals, electro_n
    return bindings, cohesive_energy, d_band_center, metals


def compare_species_bindings(s1, s2):
    all_metals = set(list(metal_dict[s1].keys()) + list(metal_dict[s2].keys()))
    shared_metals = [a for a in all_metals if a in metal_dict[s1].keys() and a in metal_dict[s2].keys()]
    s1_binding = []
    s2_binding = []
    metals = []
    for metal in shared_metals:
        if metal_dict[s1][metal] is None or metal_dict[s2][metal] is None:
            continue
        if abs(metal_dict[s1][metal]) > 13 or abs(metal_dict[s2][metal]) > 13:
            continue
        s1_binding.append(metal_dict[s1][metal])
        s2_binding.append(metal_dict[s2][metal])
        metals.append(metal)
    return s1_binding, s2_binding, metals

def make_complete_set(data_dict):
    metals = []
    for key, value in data_dict.items():
        for metal in value.keys():
            metals.append(metal)
    metals = list(set(metals))
    remove_metals = []
    for key, value in data_dict.items():
        metals_for_compound =  list(value.keys())
        for metal in metals:
            if metal not in metals_for_compound:
                remove_metals.append(metal)
    for metal in list(set(remove_metals)):
        metals.remove(metal)

    data_list = []
    compound_order = []
    for compound, value in data_dict.items():
        compound_order.append(compound)
        for metal in metals:
            data_list.append(value[metal])
    data_matrix = np.array(data_list).reshape((len(metals), len(data_dict.keys())))
    return data_matrix, compound_order, metals

#data_matrix, compound_order, metals = make_complete_set(metal_dict)

#combinations = itertools.combinations(range(len(compound_order)), 2)

"""
for i, j in combinations:
    tmp_matrix = data_matrix.copy()
    feature_1 = tmp_matrix[:, i]
    feature_2 = tmp_matrix[:, j]
    feature_1 = feature_1.reshape(len(feature_1),1)
    feature_2 = feature_2.reshape(len(feature_2),1)
    feature = np.hstack((feature_1, feature_2))
    tmp_matrix = np.delete(tmp_matrix, [i, j], 1)
    lin = LinearRegression().fit(feature, tmp_matrix)
    #print(lin.score(feature, tmp_matrix))
"""

############# d-band vs fe
common_elements = list(set(list(d_band.keys())+ list(fe_dict.keys())))
fe_s = []
d_band_s = []
syms = []
for symbol in common_elements:
    if symbol not in d_band.keys() or symbol not in fe_dict.keys():
        continue
    if fe_dict[symbol] == None or d_band[symbol] == None:
        continue
    fe_s.append(fe_dict[symbol])
    d_band_s.append(d_band[symbol])
    syms.append(symbol)

slope, intercept, r_value, p_value, std_err = linregress(d_band_s, fe_s)
fig = plt.figure()
ax = fig.add_axes([0.14,0.14,0.76,0.76])
ax.scatter(d_band_s, fe_s)
for i, j, metal in zip(d_band_s, fe_s, syms):
    ax.text(i + 0.05, j + 0.05, metal)
x_buffered_loc = (max(d_band_s) - min(d_band_s)) * 0.82 + min(d_band_s)
ax.text(x_buffered_loc, max(fe_s) + 0.5, 'R$^2$ = {}'.format(round(r_value ** 2, 2)))
plt_data = np.array([min(d_band_s), max(d_band_s)])
ax.plot(plt_data, plt_data * slope + intercept)
ax.set_title('d-band Center vs Site 2+ Formation Energy')
ax.set_ylabel('Site Formation Energy (eV)', labelpad = -0.1)
ax.set_xlabel('d-band Center (eV)')
plt.savefig('../Images/2+_d_band_vs_formation.pdf')
plt.show()

############ d-band 4+ fe
o_fe_s = fe_s.copy()
o_d_band_s = d_band_s.copy()
o_syms = syms.copy()

common_elements = list(set(list(d_band.keys())+ list(plus_4_fe.keys())))
fe_s = []
d_band_s = []
syms = []
for symbol in common_elements:
    if symbol not in d_band.keys() or symbol not in plus_4_fe.keys():
        continue
    if plus_4_fe[symbol] == None or d_band[symbol] == None:
        continue
    fe_s.append(plus_4_fe[symbol] +  1.54)
    d_band_s.append(d_band[symbol])
    syms.append(symbol)

#d_band_s += o_d_band_s
#syms += o_syms
#fe_s += o_fe_s

slope, intercept, r_value, p_value, std_err = linregress(d_band_s, fe_s)
fig = plt.figure()
ax = fig.add_axes([0.14,0.14,0.76,0.76])
ax.scatter(d_band_s, fe_s)
for i, j, metal in zip(d_band_s, fe_s, syms):
    ax.text(i + 0.05, j + 0.05, metal)
x_buffered_loc = (max(d_band_s) - min(d_band_s)) * 0.82 + min(d_band_s)
ax.text(x_buffered_loc, max(fe_s) + 0.5, 'R$^2$ = {}'.format(round(r_value ** 2, 2)))
plt_data = np.array([min(d_band_s), max(d_band_s)])
ax.plot(plt_data, plt_data * slope + intercept)
ax.set_title('d-band Center vs Site 4+ Formation Energy')
ax.set_ylabel('Site Formation Energy (eV)', labelpad = -0.1)
ax.set_xlabel('d-band Center (eV)')
plt.savefig('../Images/4+_d_band_vs_formation.pdf')
plt.show()

################ combined d-band vs fe


#d_band_s += o_d_band_s
#syms += o_syms
#fe_s += o_fe_s

slope_4, intercept_4, r_value_4, p_value, std_err = linregress(d_band_s, fe_s)
slope_2, intercept_2, r_value_2, p_value, std_err = linregress(o_d_band_s, o_fe_s)
fig = plt.figure()
ax = fig.add_axes([0.14,0.14,0.76,0.76])
ax.scatter(d_band_s, fe_s, label='4+ slabs')
ax.scatter(o_d_band_s, o_fe_s, marker='s', label='2+ slabs')
for i, j, metal in zip(d_band_s+o_d_band_s, fe_s+o_fe_s, syms+o_syms):
    ax.text(i + 0.07, j + 0.07, metal)
x_buffered_loc = (max(d_band_s) - min(d_band_s)) * 0.82 + min(d_band_s)
#ax.text(x_buffered_loc, max(fe_s) + 0.5, 'R$^2$ = {}'.format(round(r_value_4 ** 2, 2)))
plt_data = np.array([min(d_band_s), max(d_band_s)])
ax.plot(plt_data, plt_data * slope_4 + intercept_4, label='4+ fit, R$^2$={}'.format(round(r_value_4 ** 2, 2)))
ax.plot(plt_data, plt_data * slope_2 + intercept_2, '--', label='2+ fit, R$^2$={}'.format(round(r_value_2 ** 2, 2)))
plt.legend(prop={'size': 4}, bbox_to_anchor=BboxBase(), fontsize=6)
ax.set_title('d-band Center vs Site Formation Energy')
ax.set_ylabel('Site Formation Energy (eV)', labelpad = -0.1)
ax.set_xlabel('d-band Center (eV)')
#ax.set_ylim([-2, 17])
plt.legend()
plt.savefig('../Images/combined_d_band_vs_formation.pdf')
plt.show()


############# N2H vs NH2
NH2_bindings, cohesive_energy, d_band_center, NH2_metals = build_lists('NH2')
NH2, N2H, metals = compare_species_bindings('NH2', 'N2H')
slope, intercept, r_value, p_value, std_err = linregress(NH2, N2H)
fig = plt.figure()
ax = fig.add_axes([0.14,0.14,0.76,0.76])
ax.scatter(NH2, N2H)
for i, j, metal in zip(NH2, N2H, metals):
    ax.text(i + 0.01, j + 0.01, metal)
x_buffered_loc = (max(NH2_bindings) - min(NH2_bindings)) * 0 + min(NH2_bindings)
ax.text(x_buffered_loc, max(N2H) - 0.1, 'R$^2$ = {}'.format(round(r_value ** 2, 2)))
plt_data = np.array([min(NH2), max(NH2)])
ax.plot(plt_data, plt_data * slope + intercept)
ax.set_title('$\Delta E_{NH_2}$ vs $\Delta E_{N_2H}$')
ax.set_ylabel('$\Delta E_{N_2H}$ (eV)', labelpad = -0.1)
ax.set_xlabel('$\Delta E_{NH_2}$ (eV)')
plt.savefig('NH2_N2H.pdf')
plt.show()


############# N2H vs NH2
NH2, N2H, metals = compare_species_bindings('N2H', 'sa_N2H')
slope, intercept, r_value, p_value, std_err = linregress(NH2, N2H)
fig = plt.figure()
ax = fig.add_axes([0.14,0.14,0.76,0.76])
ax.scatter(NH2, N2H)
for i, j, metal in zip(NH2, N2H, metals):
    ax.text(i + 0.01, j + 0.01, metal)
x_buffered_loc = (max(NH2_bindings) - min(NH2_bindings)) * 0 + min(NH2_bindings)
ax.text(x_buffered_loc, max(N2H) - 0.1, 'R$^2$ = {}'.format(round(r_value ** 2, 2)))
plt_data = np.array([min(NH2), max(NH2)])
ax.plot(plt_data, plt_data * slope + intercept)
ax.set_title('$\Delta E_{NH_2}$ vs $\Delta E_{N_2H}$')
ax.set_ylabel('$\Delta E_{N_2H}$ (eV)', labelpad = -0.1)
ax.set_xlabel('$\Delta E_{NH_2}$ (eV)')
plt.savefig('NH2_N2H.pdf')
plt.show()

############## cohesive energy vs NH2
#sq_d_band_center = [a ** 2 for a in d_band_center]
cohesive_energy = [a * kJ / mol for a in cohesive_energy]
slope, intercept, r_value, p_value, std_err = linregress(cohesive_energy, NH2_bindings)
fig = plt.figure()
ax = fig.add_axes([0.14,0.14,0.76,0.76])
ax.scatter(cohesive_energy, NH2_bindings)
for i, j, metal in zip(cohesive_energy, NH2_bindings, NH2_metals):
    ax.text(i + 0.01, j + 0.01, metal)
plt_data = np.array([min(cohesive_energy), max(cohesive_energy)])
x_buffered_loc = (max(cohesive_energy) - min(cohesive_energy)) * 0.8 + min(cohesive_energy)
ax.text(x_buffered_loc, max(NH2_bindings) - 0.1, 'R$^2$ = {}'.format(round(r_value ** 2, 2)))
ax.plot(plt_data, plt_data * slope + intercept)
ax.set_title('$\Delta E_{NH_2}$ vs d Band Contribution of Cohesive Energy')
ax.set_ylabel('$\Delta E_{N_2H}$ (eV)', labelpad=-0.1)
ax.set_xlabel('d-Band Contribution of Cohesive Energy (eV)')
ax.yaxis.labelpad = 0
plt.savefig('../Images/cohesive_eng_vs_N2H.pdf')
plt.show()

############## cohesive energy vs N2H
N2H_bindings, cohesive_energy, d_band_center, N2H_metals = build_lists('N2H')
cohesive_energy = [a * kJ / mol for a in cohesive_energy]
slope, intercept, r_value, p_value, std_err = linregress(cohesive_energy, N2H_bindings)
fig = plt.figure()
ax = fig.add_axes([0.14,0.14,0.76,0.76])
ax.scatter(cohesive_energy, N2H_bindings)
for i, j, metal in zip(cohesive_energy, N2H_bindings, N2H_metals):
    ax.text(i + 0.01, j + 0.01, metal)
x_buffered_loc = (max(cohesive_energy) - min(cohesive_energy)) * 0.85 + min(cohesive_energy)
ax.text(x_buffered_loc, max(N2H_bindings) - 0.1, 'R$^2$ = {}'.format(round(r_value ** 2, 2)))
plt_data = np.array([min(cohesive_energy), max(cohesive_energy)])
ax.plot(plt_data, plt_data * slope + intercept)
ax.set_title('$\Delta E_{N_2H}$ vs d Band Contribution of Cohesive Energy (eV)')
ax.set_xlabel('d-band Contribution of Cohesive Energy (eV)')
ax.set_ylabel('$\Delta E_{N_2H} (eV)$', labelpad = -0.1)
plt.savefig('../Images/cohesive_eng_vs_NH2.pdf')
plt.show()


############## cohesive energy vs N2
N2H_bindings, cohesive_energy, d_band_center, N2H_metals = build_lists('N2')
cohesive_energy = [a * kJ / mol for a in cohesive_energy]
slope, intercept, r_value, p_value, std_err = linregress(cohesive_energy, N2H_bindings)
fig = plt.figure()
ax = fig.add_axes([0.14,0.14,0.76,0.76])
ax.scatter(cohesive_energy, N2H_bindings)
for i, j, metal in zip(cohesive_energy, N2H_bindings, N2H_metals):
    ax.text(i + 0.01, j + 0.01, metal)
x_buffered_loc = (max(cohesive_energy) - min(cohesive_energy)) * 0.85 + min(cohesive_energy)
ax.text(x_buffered_loc, max(N2H_bindings) - 0.1, 'R$^2$ = {}'.format(round(r_value ** 2, 2)))
plt_data = np.array([min(cohesive_energy), max(cohesive_energy)])
ax.plot(plt_data, plt_data * slope + intercept)
ax.set_title('$\Delta E_{N_2}$ vs d Band Contribution of Cohesive Energy (eV)')
ax.set_xlabel('d-band Contribution of Cohesive Energy (eV)')
ax.set_ylabel('$\Delta E_{N_2} (eV)$', labelpad = -0.1)
plt.savefig('../Images/cohesive_eng_vs_N2.pdf')
plt.show()


############## electronegativity vs formation
common_elements = list(set(list(electronegativity.keys())+ list(fe_dict.keys())))
fe_s = []
electro_n = []
syms = []
for symbol in common_elements:
    if symbol not in electronegativity.keys() or symbol not in fe_dict.keys():
        continue
    if fe_dict[symbol] == None or electronegativity[symbol] == None:
        continue
    fe_s.append(fe_dict[symbol])
    electro_n.append(electronegativity[symbol])
    syms.append(symbol)

o_fe_s = fe_s.copy()
o_electro_n = electro_n.copy()
o_syms = syms.copy()

common_elements = list(set(list(d_band.keys())+ list(plus_4_fe.keys())))
fe_s = []
electro_n = []
syms = []
for symbol in common_elements:
    if symbol not in d_band.keys() or symbol not in plus_4_fe.keys():
        continue
    if plus_4_fe[symbol] == None or d_band[symbol] == None:
        continue
    fe_s.append(plus_4_fe[symbol] + 1.54)
    electro_n.append(electronegativity[symbol])
    syms.append(symbol)

#electro_n += o_electro_n
#syms += o_syms
#fe_s += o_fe_s

print(electro_n, fe_s)
slope_4, intercept_4, r_value_4, p_value, std_err = linregress(electro_n, fe_s)
slope_2, intercept_2, r_value_2, p_value, std_err = linregress(o_electro_n, o_fe_s)
fig = plt.figure()
ax = fig.add_axes([0.14,0.14,0.76,0.76])
ax.scatter(electro_n, fe_s, label='4+ slabs')
ax.scatter(o_electro_n, o_fe_s, marker='s', label='2+ slabs')
for i, j, metal in zip(electro_n+o_electro_n, fe_s+o_fe_s, syms+o_syms):
    ax.text(i + 0.01, j - 0.07, metal)
#x_buffered_loc = (max(d_band_s) - min(d_band_s)) * 0.82 + min(d_band_s)
#ax.text(x_buffered_loc, max(fe_s) + 0.5, 'R$^2$ = {}'.format(round(r_value_4 ** 2, 2)))
plt_data = np.array([min(electro_n), max(electro_n)])
ax.plot(plt_data, plt_data * slope_4 + intercept_4, label='4+ fit, R$^2$={}'.format(round(r_value_4 ** 2, 2)))
ax.plot(plt_data, plt_data * slope_2 + intercept_2, '--', label='2+ fit, R$^2$={}'.format(round(r_value_2 ** 2, 2)))
ax.set_title('Electonegativity vs Site Formation Energy')
ax.set_ylabel('Site Formation Energy (eV)', labelpad = -0.1)
ax.set_xlabel('Metal Electronegativity')
plt.legend()
plt.savefig('../Images/electronegativity_vs_formation.pdf')
plt.show()

