import os
from collections import defaultdict
from data import plus_4_fe, plus_4_N2, plus_4_N2H 
import csv

metal_dict = defaultdict(dict)
all_metals = []
all_species = []


N2_corr = 0.0356539561469
N2H_corr = 0.482259647646

N2_gas_corr = -0.351226774491
H2_gas_corr = -0.0394892748343


def subscipt(species):
    subscripted = ''
    for i, char in enumerate(species):
        if i == 0:
            subscripted += char
        elif char.isdecimal():
            subscripted += '$_' + char + '$'
        elif char == '_':
            subscripted += ' '
        else:
            subscripted += char
    return subscripted


for pathway in os.listdir('../data/corrected_data'):
    if '_' in pathway:
        if pathway != 'formation_energy.csv':
            continue
    with open('../data/corrected_data/' + pathway, 'r') as f:
        if pathway == 'scaling.py':
            continue
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

g.write('\\begin{table}\n')
g.write('\\begin{center}\n\\begin{tabular}{| c | c | c | c | c | c | c | c | c | c | c | c | c | c |}\n')

g.write('\hline\n')
g.write('Element & ')
all_species = sorted(all_species)
for i, species in enumerate(all_species):
    subscripted = ''
    for j, char in enumerate(species):
        if j == 0:
            subscripted += char
        elif char.isdecimal():
            subscripted += '$_' + char + '$'
        elif char == '_':
            subscripted += ' '
        else:
            subscripted += char
        if subscripted == 'formation energy':
            subscripted = 'Formation Energy'
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
g.write('\\end{center}\n')
g.write('\\caption{The calculated relative energies of all 2+ surface species on all metal substituents at standard state. All energies are referenced with respect to N$_2$ gas and H$_2$ gas at 300K and 1 bar of pressure. Blank spaces represent calculations that could not be converged}\n')
g.write('\\label{table:energies}\n')
g.write('\\end{table}\n\n')

g.write('\\begin{table}\n')
g.write('\\begin{center}\n\\begin{tabular}{| c | c | c | c |}\n')
g.write('\hline\n')
g.write('Element & N$_2$ & N$_2$H & Formation Energy \\\\\n')
g.write('\\hline\n')
available_elements = list(set(list(plus_4_fe.keys())+ list(plus_4_N2.keys())+ list(plus_4_N2H.keys())))
for element in available_elements:
    g.write(element + ' & ')
    if element in plus_4_N2.keys():
        g.write(str(round(plus_4_N2[element] + N2_corr - N2_gas_corr, 2)) + ' & ')
    else:
        g.write(' & ')
    if element in plus_4_N2H.keys():
        g.write(str(round(plus_4_N2H[element] + N2H_corr - N2_gas_corr - H2_gas_corr/2 , 2)) + ' & ')
    else:
        g.write(' & ')
    if element in plus_4_fe.keys():
        g.write(str(round(plus_4_fe[element], 2)))
    else:
        pass
    g.write(' \\\\\n')



g.write('\hline\n')
g.write('\\end{tabular}\n')
g.write('\\end{center}\n')
g.write('\\caption{The calculated relative energies of all 4+ surface species on all metal substituents at standard state. All energies are referenced with respect to N$_2$ gas and H$_2$ gas at 300K and 1 bar of pressure. Blank spaces represent calculations that could not be converged}\n')

g.write('\hline\n')
g.write('\\end{table}\n\n')


g.write('\\begin{table}\n\\begin{center}\n\\begin{tabular}{| c | c |c |}\n')

g.write('\hline\n')
g.write('Element & Limiting Potential & Limiting Step \\\\\n')
g.write('\\hline\n')
txt = csv.reader(open('../data/pathway_data/limiting_potential_associative_2.csv', 'r'))
txt = list(txt)
for metal in txt:
    g.write(metal[0] + ' & ' + str(round(float(metal[1]), 2)) + ' & ' + subscipt(metal[2]) + ' $\\rightarrow$ ' + subscipt(metal[3]))

    g.write('\\\\\n')
g.write('\\hline\n')
g.write('\\end{tabular}\n\\end{center}\n')
g.write('\\caption{The limiting potentials and limiting steps for each dopant metal}')
g.write('\\label{table:limiting_steps}')
g.write('\\end{table}')


g.write('\hline\n')
g.write('\\end{table}\n\n')


g.write('\\begin{table}\n\\begin{center}\n\\begin{tabular}{| c | c |c |}\n')
g.write('\hline\n')
g.write('Element & Limiting Potential & Limiting Step \\\\\n')
g.write('\\hline\n')
txt = csv.reader(open('../data/pathway_data/thermo_limiting_potential_associative_2.csv', 'r'))
txt = list(txt)
for metal in txt:
    g.write(metal[0] + ' & ' + str(round(float(metal[1]), 2)) + ' & ' + subscipt(metal[2]) + ' $\\rightarrow$ ' + subscipt(metal[3]))

    g.write('\\\\\n')
g.write('\\hline\n')
g.write('\\end{tabular}\n\\end{center}\n')
g.write('\\caption{The limiting potentials and limiting steps for each dopant metal}')
g.write('\\label{table:limiting_steps}')
g.write('\\end{table}')

g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/scaling_species.pdf}\n\\caption{The calculated scaling relations between the binding energies of various species and the binding energies of N$_2$H and NH$_2$ on 2+ dopant sites}\n\\end{figure}\n\n')

g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/scaling_reactions.pdf}\n\\caption{The calculated scaling relations between the binding energies of various species and the binding energies of N$_2$H and NH$_2$ on 2+ dopant sites}\n\\end{figure}\n\n')

g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/electronegativity_vs_formation.pdf}\n\\caption{Electronegativity vs formation energy of 2+ dopant site}\n\\end{figure}\n\n')

g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/N2_vs_N2H.pdf}\n\\caption{The binding energy of N$_2$ vs N$_2$H}\n\\end{figure}\n\n')

g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/Valence_vs_formation_energy.pdf}\n\\caption{Valence number vs formation energy of 2+ dopant site}\n\\end{figure}\n\n')

for plot in os.listdir('../data/plots/'):
    g.write('\\begin{figure}\n\\includegraphics[width=0.8\\linewidth]{data/plots/')
    g.write(plot)
    g.write('}\n\\end{figure}\n\n')


g.close()
