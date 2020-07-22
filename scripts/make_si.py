import os
from collections import defaultdict
from data import plus_4_fe, plus_4_N2, plus_4_N2H 
import json
import csv

metal_dict = defaultdict(dict)
mag_dict = defaultdict(dict)
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
    if pathway == 'formation_energy.csv':
        continue
    if '_' in pathway:
        if pathway != 'formation_energy.csv':
            continue
    with open('../data/corrected_data/' + pathway, 'r') as f:
        if pathway == 'scaling.py':
            continue
        data = csv.reader(f)
        for row in data:
            metal_dict[pathway.split('.')[0]][row[0]] = float(row[1])

for pathway in os.listdir('../data/magnetizations'):
    with open('../data/magnetizations/' + pathway, 'r') as f:
        data = csv.reader(f)
        for row in data:
            species = pathway.split('_magnetizations.')[0].split('_')[1]
            mag_dict[species][row[0]] = float(row[1])


fes = json.load(open('../data/corrected_data/formation_energy.csv', 'r'))
metal_dict['formation_energy'] = fes

for key, value in metal_dict.items():
    for metal in value.keys():
        all_metals.append(metal)
    all_metals = list(set(all_metals))

all_species = list(metal_dict.keys())


g = open('../SI.tex', 'w')

header = r"""\documentclass{article}

\usepackage{graphicx}
\usepackage{longtable}
\usepackage{lastpage}
\usepackage{array}
\setlength\extrarowheight{2pt}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\thefigure}{S\arabic{figure}}

\usepackage{geometry}
\geometry{
letterpaper,
total={170mm,257mm},
left=20mm,
top=20mm,
}

\newcommand*{\TitleFont}{%
      \usefont{\encodingdefault}{\rmdefault}{b}{n}%
            \fontsize{25}{30}%
                  \selectfont}
\title{\TitleFont Supplementary Information for Computational Study of Transition-Metal Substitutions in Rutile TiO$_2$ (110) for Photoelectrocatalytic Ammonia Synthesis}

\author{Benjamin M. Comer, Max H. Lenk, Aradhya P. Rajanala, Emma L. Flynn, Andrew J. Medford}
\begin{document}

\maketitle"""

g.write(header)

g.write('\\begin{table}\n')
g.write('\\setlength\\tabcolsep{2pt}\n')
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
g.write('\\label{table:4+_energies}\n')
g.write('\\caption{The calculated relative energies of all 4+ surface species on all metal substituents at standard state. All energies are referenced with respect to N$_2$ gas and H$_2$ gas at 300K and 1 bar of pressure. Blank spaces represent calculations that could not be converged}\n')

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
g.write('\\caption{The limiting potentials and limiting steps for each dopant metal on 2+ surfaces}')
g.write('\\label{table:pot_limiting_steps}')
g.write('\\end{table}')



g.write('\\begin{table}\n\\begin{center}\n\\begin{tabular}{| c | c |c |}\n')
g.write('\hline\n')
g.write('Element & Largest Thermodynamic Step & Limiting Step \\\\\n')
g.write('\\hline\n')
txt = csv.reader(open('../data/pathway_data/thermo_limiting_potential_associative_2.csv', 'r'))
txt = list(txt)
for metal in txt:
    g.write(metal[0] + ' & ' + str(round(float(metal[1]), 2)) + ' & ' + subscipt(metal[2]) + ' $\\rightarrow$ ' + subscipt(metal[3]))

    g.write('\\\\\n')
g.write('\\hline\n')
g.write('\\end{tabular}\n\\end{center}\n')
g.write('\\caption{The largest barrier for thermochemical steps and corresponding steps for each dopant metal on 2+ surfaces}')
g.write('\\label{table:thermo_limiting_steps}')
g.write('\\end{table}')


g.write('\\begin{table}\n\\begin{center}\n\\begin{tabular}{| c | c |c |}\n')
g.write('\hline\n')
g.write('Element & Rate Limiting Step & Limiting Step \\\\\n')
g.write('\\hline\n')
txt = csv.reader(open('../data/pathway_data/rate_limiting_steps_associative_2.csv', 'r'))
txt = list(txt)
for metal in txt:
    g.write(metal[0] + ' & ' + str(round(float(metal[1]), 2)) + ' & ' + subscipt(metal[2]) + ' $\\rightarrow$ ' + subscipt(metal[3]))

    g.write('\\\\\n')
g.write('\\hline\n')
g.write('\\end{tabular}\n\\end{center}\n')
g.write('\\caption{The largest thermodynamic barrier and corresponding steps for each dopant metal on 2+ surfaces when set at the band edge of rutile, -0.142V}')
g.write('\\label{table:rate_limiting_steps}')
g.write('\\end{table}')


g.write('\\begin{table}\n')
g.write('\\setlength\\tabcolsep{2pt}\n')
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
        if metal in mag_dict[species].keys():
            g.write(str(round(mag_dict[species][metal], 2)))
        if i != len(metal_dict.keys()) - 1:
            g.write(' & ')
    g.write(' \\\\\n')

g.write('\hline\n')

g.write('\\end{tabular}\n')
g.write('\\end{center}\n')
g.write('\\caption{The calculated total magnetic moments for the unit cells of each 2+ species. Blank spaces represent calculations that could not be converged}\n')
g.write('\\label{table:mags}\n')
g.write('\\end{table}\n\n')

g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/scaling_species.pdf}\n\\caption{The calculated scaling relations between the binding energies of various species and the binding energies of N$_2$H and NH$_2$ on 2+ dopant sites}\n\\label{fig:scaling_species}\n\\end{figure}\n\n')

g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/scaling_reactions.pdf}\n\\caption{The calculated scaling relations between the reaction energies energies of all electrochemical reations and the binding energies of N$_2$H and NH$_2$ on 2+ dopant sites}\n\\label{fig:scaling_reactions}\n\\end{figure}\n\n')

g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/electronegativity_vs_formation.pdf}\n\\caption{Electronegativity vs formation energy of 2+ dopant site}\n\\label{fig:electronegativity}\n\\end{figure}\n\n')

#g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/N2_vs_N2H.pdf}\n\\caption{The binding energy of N$_2$ vs N$_2$H}\n\\end{figure}\n\n')

#g.write('\\begin{figure}\n\\centering\n\\includegraphics[width=0.8\\linewidth]{Images/Valence_vs_formation_energy.pdf}\n\\caption{Valence number vs formation energy of 2+ dopant site}\n\\end{figure}\n\n')

g.write('\\twocolumn\n')

for i, plot in enumerate(os.listdir('../data/plots/')):
    if not i % 6:
        g.write('\\clearpage\n')
    g.write('\\begin{figure}\n\\includegraphics[width=1\\linewidth]{data/plots/')
    g.write(plot)
    g.write('}\n\\label{fig:' + plot.split('.')[0] + '}\n\\caption{Free energy diagram for '+plot.split('_')[0]+'}\n\\end{figure}\n\n')

g.write('\\end{document}')

g.close()
