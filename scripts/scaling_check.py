import csv
import os
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
from data import fe, column, element, cohesive_energies, N2_engs, N2H_engs, d_band, electronegativity,d_cohesive, s_cohesive
import numpy as np

reaction_data =[]
path_counter = 0
element_symbol = []
for pathway in os.listdir('../data/zp_pathway_data'):
    if 'associative.' not in pathway or 'rate' in pathway:
        continue
    with open('../data/pathway_data/' + pathway, 'r') as f:
        element_symbol.append(pathway.split('_')[0])
        data = csv.reader(f)
        #reaction_data += data[-3]
        path_counter += 1 # can't think of a better way quickly
        for row in data:
            reaction_data += [float(a) for a in row]
        # figure out how long one of them is
        if path_counter == 1:
            path_length = len(reaction_data)

reaction_data = [float(a) for a in reaction_data]
# build a matrix 
data_matrix = np.array(reaction_data).reshape(path_length, path_counter)

nh2 = data_matrix[-3] - data_matrix[-1]
nh3 = data_matrix[-2] - data_matrix[-1]
d_band_center = [d_band[a] for a in element_symbol]
cohesive_energy = [d_cohesive[a] for a in element_symbol]

plt.scatter(cohesive_energy, data_matrix[2] + 0.15)
plt.show()
