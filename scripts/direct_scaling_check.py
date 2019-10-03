import csv
import os
from collections import defaultdict

from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
from data import fe, column, element, cohesive_energies, N2_engs, N2H_engs, d_band, electronegativity,d_cohesive, s_cohesive
from scipy.stats import linregress
import numpy as np



metal_dict = defaultdict(dict)
path_counter = 0
compound = []
for pathway in os.listdir('../data/data'):
    if '_' in pathway:
        continue
    with open('../data/data/' + pathway, 'r') as f:
        compound.append(pathway.split('.')[0])
        data = csv.reader(f)
        for row in data:
            metal_dict[pathway.split('.')[0]][row[0]] = float(row[1])
        # figure out how long one of them is

# make a dummy comparison entry
#all_metals = set(list(metal_dict['NH2'].keys()) + list(metal_dict['NH3'].keys()))
#shared_metals = [a for a in all_metals if a in metal_dict[s1].keys() and a in metal_dict[s2].keys()]

#reaction_data = [float(a) for a in reaction_data]
# build a matrix 

def build_lists(species):
    d_band_center = []
    cohesive_energy = []
    bindings = []
    for metal, binding in metal_dict[species].items():
        if binding is None:
            continue
        d_band_center.append(d_band[metal])
        cohesive_energy.append(d_cohesive[metal])
        bindings.append(float(binding))
    return bindings, cohesive_energy, d_band_center

def compare_species_bindings(s1, s2):
    all_metals = set(list(metal_dict[s1].keys()) + list(metal_dict[s2].keys()))
    shared_metals = [a for a in all_metals if a in metal_dict[s1].keys() and a in metal_dict[s2].keys()]
    s1_binding = []
    s2_binding = []
    for metal in shared_metals:
        if metal_dict[s1][metal] is None or metal_dict[s2][metal] is None:
            continue
        s1_binding.append(metal_dict[s1][metal])
        s2_binding.append(metal_dict[s2][metal])
    return s1_binding, s2_binding

NH2_bindings, cohesive_energy, d_band_center = build_lists('NH2')
N2H_bindings, cohesive_energy, d_band_center = build_lists('NH3')
NH2, N2H = compare_species_bindings('NH3', 'N2H')
slope, intercept, r_value, p_value, std_err = linregress(NH2, N2H)
fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.scatter(NH2, N2H)
plt_data = np.array([min(NH2), max(NH2)])
ax.plot(plt_data, plt_data * slope + intercept)
print(r_value ** 2)
plt.show()
