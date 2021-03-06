from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress
from data import fe_dict, column_dict, element, cohesive_energies, d_band, electronegativity,d_cohesive, s_cohesive, plus_4_N2H, plus_4_N2, NH2_engs_dict, N2H_engs_dict, N2_engs_dict 

plt.rcParams["figure.figsize"] = (4.5,3.5)
plt.rcParams["font.size"] = 10


#cohesive_energies = s_cohesive

n_col = []
n_fe = []
n_ele = []
n_N2_engs = []
n_N2H_engs = []
nn_N2H_engs = []
N2H_engs = []
N2_engs = []


for key, e in fe_dict.items():
    #print(element[i], column[i], e)
    if e is None:
        continue
    if key in ['Fe', 'Y']:
        continue
    n_col.append(column_dict[key])
    n_fe.append(fe_dict[key])
    n_ele.append(key)
    n_N2H_engs.append(N2H_engs_dict[key])
    

n_ele_N2 = []
for e in element:
    if e in ['Fe', 'Mn', 'Y']:
        continue
    n_N2_engs.append(N2_engs_dict[e])
    nn_N2H_engs.append(N2H_engs_dict[e])
    n_ele_N2.append(e)
    N2_engs.append(N2_engs_dict[e])
    N2_engs.append(N2H_engs_dict[e])

########################### periodic table position

slope, intercept, r_value, p_value, std_err = linregress(n_col, n_fe)

plt.scatter(n_col, n_fe)
#plt.scatter([cohesive_energies[a] for a in n_ele], n_fe)
plt.plot([min(n_col),max(n_col)],np.array([min(n_col),max(n_col)]) * slope + intercept)
plt.title('Valence Number vs Formation Energy')
plt.ylabel('2+ Substitution Formation Energy')
plt.xlabel('Valence Number')
plt.text(min(n_col,),max(n_fe)-0.25, 'R$^2$ = {}'.format(round(r_value ** 2, 2)))
for c, e, t in zip(n_col, n_fe, n_ele):
    plt.text(c + 0.1, e + 0.1, t)

plt.savefig('../Images/Valence_vs_formation_energy.pdf')
plt.show()

############################ cohesive energy N2H

slope, intercept, r_value, p_value, std_err = linregress([cohesive_energies[a] for a in n_ele_N2],
                                                         nn_N2H_engs)

plt.title('Metal Cohesive Energy vs N$_2$H Binding Energy')
plt.ylabel('N$_2$H Binding Energy')
plt.xlabel('Metal Cohesive Energy (eV/atom)')
co_e = [cohesive_energies[a] for a in n_ele_N2]
plt.text(max(co_e)-0.75,max(nn_N2H_engs), 'R$^2$ = {}'.format(round(r_value ** 2,2)))
plt.scatter(co_e, nn_N2H_engs)
plt.plot([min(co_e),max(co_e)],np.array([min(co_e),max(co_e)]) * slope + intercept)
for c, e, t in zip(co_e, nn_N2H_engs, n_ele_N2):
    plt.text(c + 0.05, e + 0.05, t)

plt.savefig('cohesive_eng_vs_N2H.pdf')
plt.show()


################################ cohesive energy N2
slope, intercept, r_value, p_value, std_err = linregress([cohesive_energies[a] for a in n_ele_N2],
                                                         n_N2_engs)

plt.title('Metal Cohesive Energy vs N$_2$ Binding Energy')
plt.ylabel('N$_2$ Binding Energy')
plt.xlabel('Metal Cohesive Energy (eV/atom)')
co_e = [cohesive_energies[a] for a in n_ele_N2]
plt.text(max(co_e)-0.75,max(n_N2_engs), 'R$^2$ = {}'.format(round(r_value ** 2,2)))
plt.scatter(co_e, n_N2_engs)
plt.plot([min(co_e),max(co_e)],np.array([min(co_e),max(co_e)]) * slope + intercept)
for c, e, t in zip(co_e, n_N2_engs, n_ele_N2):
    plt.text(c + 0.05, e + 0.05, t)

plt.savefig('cohesive_eng_vs_N2.pdf')
plt.show()

################################ cohvesive energy N2H

common_elements = list(set(list(d_cohesive.keys())+ list(plus_4_N2H.keys())))
N2H_s = []
d_coh_s = []
syms = []
for symbol in common_elements:
    if symbol not in d_cohesive.keys() or symbol not in plus_4_N2H.keys():
        continue
    if plus_4_N2H[symbol] == None or d_cohesive[symbol] == None:
        continue
    if symbol in ['Cu','Ag']:
        continue
    N2H_s.append(plus_4_N2H[symbol])
    d_coh_s.append(d_cohesive[symbol])
    syms.append(symbol)

slope, intercept, r_value, p_value, std_err = linregress(d_coh_s, N2H_s)
plt.title('Metal Cohesive Energy vs 4+ N$_2$H Binding Energy')
plt.ylabel('N$_2$H Binding Energy')
plt.xlabel('Metal Cohesive Energy (eV/atom)')
co_e = d_coh_s
#plt.text(max(co_e)-0.75,max(n_N2_engs), 'R$^2$ = {}'.format(round(r_value ** 2,2)))
plt.scatter(d_coh_s, N2H_s)
plt.plot([min(co_e),max(co_e)],np.array([min(co_e),max(co_e)]) * slope + intercept)
for e, c, t in zip(co_e, d_coh_s, N2H_s):
    plt.text(c + 0.05, e + 0.05, t)

plt.savefig('cohesive_eng_vs_4+_N2H.pdf')
plt.show()


########################### N2 vs N2H

slope, intercept, r_value, p_value, std_err = linregress(n_N2_engs, nn_N2H_engs)

plt.title('N$_2$ Binding Energy vs N$_2$H Binding Energy')
plt.ylabel('N$_2$H Binding Energy (eV)')
plt.xlabel('N$_2$ Binding Energy (eV)')
plt.text(max(n_N2_engs)-0.75,max(nn_N2H_engs), 'R$^2$ = {}'.format(round(r_value ** 2,2)))
plt.scatter(n_N2_engs, nn_N2H_engs)
plt.plot([min(n_N2_engs),max(n_N2_engs)],np.array([min(n_N2_engs),max(n_N2_engs,)]) * slope + intercept)
for c, e, t in zip(n_N2_engs, nn_N2H_engs, n_ele_N2):
    plt.text(c + 0.05, e + 0.05, t)

plt.savefig('../Images/N2_vs_N2H.pdf')
plt.show()

###########################  d band center vs formation

d_band_fe_list = [a for a, b in zip(n_fe, n_ele) if d_band[b] is not None]
d_band_list = [d_band[b] for b in n_ele if d_band[b] is not None]
d_band_element_list = [b for a, b in zip(n_fe, n_ele) if d_band[b] is not None]
d_band_col_list = [a for a, b in zip(n_col, n_ele) if d_band[b] is not None] 

slope, intercept, r_value, p_value, std_err = linregress(d_band_list, d_band_fe_list)

plt.title('d Band Center vs 2+ Formation Energy')
plt.xlabel('d Band Center (eV)')
plt.ylabel('2+ Formation Energy (eV)')
plt.text(max(d_band_list)-1,max(d_band_fe_list)+.5, 'R$^2$ = {}'.format(round(r_value ** 2,2)))
plt.scatter(d_band_list, d_band_fe_list)
plt.plot([max(d_band_list),min(d_band_list)],np.array([max(d_band_list),min(d_band_list)]) * slope + intercept)
for c, e, t in zip(d_band_list, d_band_fe_list, d_band_element_list):
    plt.text(c + 0.05, e + 0.05, t)

plt.savefig('d_band_vs_formation.pdf')
plt.show()

########################### row bar charts


plt.rcParams["figure.figsize"] = (5.5,8)

fig, _axs = plt.subplots(nrows=2, ncols=1)
axs = _axs.flatten()

row_1_elements = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
row_2_elements = ['Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag']
row_3_elements = ['', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au']

positive_row_energies = [[],[],[]]
negative_row_energies = [[],[],[]]

rows = [row_1_elements, row_2_elements]#, row_3_elements]

for i, row in enumerate(rows):
    for element_sym in row:
        if element_sym == '':
            eng = 0
        elif element_sym not in element:
            eng = 0
        else:
            index = element.index(element_sym)
            eng = N2_engs[index] 
            if eng is None:
                eng = 0
        if eng >= 0:
            positive_row_energies[i].append(eng)
            negative_row_energies[i].append(0)
        elif eng < 0:
            positive_row_energies[i].append(0)
            negative_row_energies[i].append(eng)
#del rows[0]
#del positive_row_energies[0]
#del negative_row_energies[0]

for i, row in enumerate(rows):
    axs[i].bar(row, positive_row_energies[i], color='#0390fc')
    axs[i].bar(row, negative_row_energies[i], color='#0390fc')
    #axs[i].scatter(row, positive_row_energies[i], color='#0390fc')
    #axs[i].scatter(row, negative_row_energies[i], color='#0390fc')

    #axs[i].set_title('Row {}'.format(i + 4))
    axs[i].set_ylim([-1.6,0])
    if i == 1:
        axs[i]. set_ylabel('N$_2$ Adsorption Energy (eV)')
    if i == 1:
        axs[i].set_xlabel('Element')
    if i == 0:
        axs[i].set_title('N$_2$ Adsorption Energy by Periodic Row')
plt.savefig('N2_adsorption_rows.pdf')
plt.show()
#####################

plt.rcParams["figure.figsize"] = (5.5,8)

fig, _axs = plt.subplots(nrows=2, ncols=1)
axs = _axs.flatten()

row_1_elements = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
row_2_elements = ['Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag']
row_3_elements = ['', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au']

positive_row_energies = [[],[],[]]
negative_row_energies = [[],[],[]]

rows = [row_1_elements, row_2_elements, row_3_elements]

for i, row in enumerate(rows):
    for element_sym in row:
        if element_sym in ['', 'Mn', 'Fe', 'Y']:
            eng = 0
        elif element_sym not in element:
            eng = 0
        else:
            index = element.index(element_sym)
            eng = N2H_engs_dict[element_sym]
            if eng is None:
                eng = 0
        if eng >= 0:
            positive_row_energies[i].append(eng)
            negative_row_energies[i].append(0)
        elif eng < 0:
            positive_row_energies[i].append(0)
            negative_row_energies[i].append(eng)
del rows[0]
del positive_row_energies[0]
del negative_row_energies[0]

for i, row in enumerate(rows):
    axs[i].bar(row, positive_row_energies[i], color='#0390fc')
    axs[i].bar(row, negative_row_energies[i], color='#0390fc')
    #axs[i].set_title('Row {}'.format(i + 4))
    axs[i].set_ylim([-1.5,2])
    if i == 1:
        axs[i]. set_ylabel('N$_2$H Adsorption Energy (eV)')
    if i == 1:
        axs[i].set_xlabel('Element')
    if i == 0:
        axs[i].set_title('N$_2$H Adsorption Energy by Periodic Row')
plt.savefig('../Images/N2H_adsorption_rows.pdf')
plt.show()


#################################### 4+

plt.rcParams["figure.figsize"] = (5.5,8)

fig, _axs = plt.subplots(nrows=3, ncols=1)
axs = _axs.flatten()

row_1_elements = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
row_2_elements = ['Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag']
row_3_elements = ['', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au']

positive_row_energies = [[],[],[]]
negative_row_energies = [[],[],[]]

rows = [row_1_elements, row_2_elements, row_3_elements]

for i, row in enumerate(rows):
    for element_sym in row:
        if element_sym == '':
            eng = 0
        elif element_sym not in plus_4_N2.keys():
            eng = 0
        else:
            eng = plus_4_N2[element_sym] 
            if eng is None:
                eng = 0
        if eng >= 0:
            positive_row_energies[i].append(eng)
            negative_row_energies[i].append(0)
        elif eng < 0:
            positive_row_energies[i].append(0)
            negative_row_energies[i].append(eng)

for i, row in enumerate(rows):
    axs[i].bar(row, positive_row_energies[i], color='#0390fc')
    axs[i].bar(row, negative_row_energies[i], color='#0390fc')
    #axs[i].set_title('Row {}'.format(i + 4))
    axs[i].set_ylim([-1.6,0])
    if i == 1:
        axs[i]. set_ylabel('N$_2$ Adsorption Energy (eV)')
    if i == 2:
        axs[i].set_xlabel('Element')
    if i == 0:
        axs[i].set_title('N$_2$ Adsorption Energy by Periodic Row')
plt.savefig('N2_adsorption_rows.pdf')
plt.show()


########################################### N2H

fig, _axs = plt.subplots(nrows=3, ncols=1)
axs = _axs.flatten()
N2H_corr = 0.482259647646
NH2_corr = 0.681757082272
N2_corr = 0.0356539561469
N2_gas_corr = -0.351226774491
H2_gas_corr = -0.0394892748343


row_1_elements = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
row_2_elements = ['Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag']
row_3_elements = ['', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au']

positive_row_energies = [[],[],[]]
negative_row_energies = [[],[],[]]

rows = [row_1_elements, row_2_elements, row_3_elements]

for i, row in enumerate(rows):
    for element_sym in row:
        if element_sym == '':
            eng = 0
        elif element_sym not in plus_4_N2H.keys():
            eng = 0
        else:
            eng = plus_4_N2H[element_sym] + N2H_corr - N2_gas_corr - H2_gas_corr/2 
            if eng is None:
                eng = 0
        if eng >= 0:
            positive_row_energies[i].append(eng)
            negative_row_energies[i].append(0)
        elif eng < 0:
            positive_row_energies[i].append(0)
            negative_row_energies[i].append(eng)

for i, row in enumerate(rows):
    axs[i].bar(row, positive_row_energies[i], color='#0390fc')
    axs[i].bar(row, negative_row_energies[i], color='#0390fc')
    #axs[i].set_title('Row {}'.format(i + 4))
    axs[i].set_ylim([-1.5,2.5])
    if i == 1:
        axs[i].set_ylabel('N$_2$H Adsorption Energy (eV)')
    if i == 2:
        axs[i].set_xlabel('Element')
    if i == 0:
        axs[i].set_title('N$_2$H Adsorption Energy by Periodic Row')

plt.savefig('N2H_adsorption_rows')
plt.show()


########################################### Combined scatters

plt.rcParams["figure.figsize"] = (4.8, 12.5)
fig, _axs = plt.subplots(nrows=3, ncols=1)
axs = _axs
axs = _axs.flatten()
N2H_corr = 0.482259647646
N2_gas_corr = -0.351226774491
H2_gas_corr = -0.0394892748343


row_1_elements = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
row_2_elements = ['Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag']
row_3_elements = ['', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au']

row_nums = [3,4,5,6,7,8,9,10,11]
markers = ['o', 's', '^']


positive_row_energies = [[],[],[]]
negative_row_energies = [[],[],[]]
numbered_row_energies = [[],[],[]]

rows = [row_1_elements, row_2_elements, row_3_elements]

for j, species_dict in enumerate([N2_engs_dict, N2H_engs_dict, NH2_engs_dict]):
    positive_row_energies = [[],[],[]]
    negative_row_energies = [[],[],[]]
    numbered_row_energies = [[],[],[]]

    for i, row in enumerate(rows):
        for element_sym in row:
            if element_sym == '':
                eng = None
            elif element_sym not in species_dict.keys():
                eng = None
            # iron is messed up for N2H
            elif  j == 1 and element_sym == 'Fe':
                eng == None
            else:
                if species_dict[element_sym] is None:
                    eng = None
                elif j == 0:
                    eng = species_dict[element_sym] + N2_corr - N2_gas_corr
                elif j == 1:
                    eng = species_dict[element_sym] + N2H_corr - N2_gas_corr - H2_gas_corr/2
                elif j ==2:
                    eng = species_dict[element_sym] + NH2_corr - N2_gas_corr/2 - H2_gas_corr
            if eng == None:
                numbered_row_energies[i].append(None)
            else:
                numbered_row_energies[i].append(eng)
            if eng == None:
                eng = 0
            if eng >= 0:
                positive_row_energies[i].append(eng)
                negative_row_energies[i].append(0)
            elif eng < 0:
                positive_row_energies[i].append(0)
                negative_row_energies[i].append(eng)

    for i, row in enumerate(rows):
        #axs[i].bar(row, positive_row_energies[i], color='#0390fc')
        #axs[i].bar(row, negative_row_energies[i], color='#0390fc')
        #axs[i].set_title('Row {}'.format(i + 4))
        axs[j].scatter(row_nums, numbered_row_energies[i], label='row {}'.format(i + 4),
                       marker=markers[i])
        for c, e, t in zip(row_nums, numbered_row_energies[i], row):
            if c is None or e is None:
                continue
            axs[j].text(c + 0.05, e + 0.05, t)
        if j == 0:
            axs[j].set_ylabel('N$_2$ Adsorption Energy (eV)')
            axs[j].set_ylim([-1.5,1.2])
            axs[j].set_xlabel('(a)')
        if j == 1:
             axs[j].set_ylabel('N$_2$H Adsorption Energy (eV)')
             axs[j].set_ylim([-0.5,2.9])
             axs[j].set_xlabel('(b)')
        if j == 2:
             axs[j].set_ylabel('NH$_2$ Adsorption Energy (eV)')
             axs[j].set_ylim([-1.5,2.5])
        if j == 2:
            axs[j].set_xlabel('(c)\nColumn')
        if j == 0:
            axs[j].set_title('Binding Energy vs Periodic Column')
        axs[0].legend()

plt.tight_layout()
plt.savefig('../Images/adsorption_rows.pdf')
plt.show()

