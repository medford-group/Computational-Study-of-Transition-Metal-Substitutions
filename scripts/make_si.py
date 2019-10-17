import os
from collections import defaultdict
import csv

metal_dict = defaultdict(dict)
all_metals = []
all_species = []

for pathway in os.listdir('../data/corrected_data'):
    if '_' in pathway:
        if pathway != 'formation_energy.csv':
            continue
    with open('../data/corrected_data/' + pathway, 'r') as f:
        data = csv.reader(f)
        for row in data:
            metal_dict[pathway.split('.')[0]][row[0]] = float(row[1])


for key, value in metal_dict.items():
    for metal in value.keys():
        all_metals.append(metal)
    all_metals = list(set(all_metals))

all_species = list(metal_dict.keys())

g = open('../SI.tex', 'w')

g.write('Supplementary Information\n')

g.write('\onecolumn\n')

g.write('\\begin{center}\n\\begin{tabular}{| c | c | c | c | c | c | c | c | c | c | c | c | c | c |}\n')

g.write('\hline\n')
g.write('element & ')
all_species = sorted(all_species)
for i, species in enumerate(all_species):
    subscripted = ''
    for char in species:
        if char.isdecimal():
            subscripted += '$_' + char + '$'
        elif char == '_':
            subscripted += ' '
        else:
            subscripted += char
    g.write(subscripted)
    if i != len(all_species) - 1:
        g.write(' & ')
g.write('\\\\\n\hline\n')
g.write('\n')

for metal in all_metals:
    g.write(metal + ' & ')
    for i, species in enumerate(all_species):
    #for i, species, values in zip(range(len(metal_dict.keys())), metal_dict.keys(),\
    #                              metal_dict.values()):
        if metal in metal_dict[species].keys():
            g.write(str(round(metal_dict[species][metal], 2)))
        if i != len(metal_dict.keys()) - 1:
            g.write(' & ')
    g.write(' \\\\\n')

g.write('\hline\n')
g.write('\\end{tabular}\n')
g.write('\\caption{The calculated relative energies of all surface species on all metal substituents at standard state. All energies are referenced with respect to N$_2$ gas and H$_2$ gas at 300K and 1 bar of pressure. Blank spaces represent calculations that could not be converged}\n')
g.write('\\label{table:energies}\n')
g.write('\\end{center}\n\n\n\n')

g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/electronegativity_vs_formation.pdf}\n\\caption{Electronegativity vs formation energy of 2+ dopant site}\n\\end{figure}\n\n')


g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/Valence_vs_formation_energy.pdf}\n\\caption{Valence number vs formation energy of 2+ dopant site}\n\\end{figure}\n\n')

for plot in os.listdir('../data/plots/'):
    g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{data/plots/')
    g.write(plot)
    g.write('}\n\\end{figure}\n\n')


g.close()
