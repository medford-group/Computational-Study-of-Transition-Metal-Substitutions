from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress
from data import fe, column, element, cohesive_energies, N2_engs, N2H_engs, d_band, electronegativity,d_cohesive, s_cohesive
plt.rcParams["figure.figsize"] = (5.2,4.2)

#cohesive_energies = s_cohesive

n_col = []
n_fe = []
n_ele = []
n_N2_engs = []
n_N2H_engs = []
nn_N2H_engs = []

for i, e in enumerate(fe):
    if e is None:
        continue
    n_col.append(column[i])
    n_fe.append(fe[i])
    n_ele.append(element[i])
    n_N2H_engs.append(N2H_engs[i])

n_ele_N2 = []
for e, n2,n2h in zip(element,N2_engs,N2H_engs):
    if n2 is None:
        continue
    if e in ['Fe', 'Ni', 'Co', 'Mn', 'Cr']:
        continue
    n_N2_engs.append(n2)
    nn_N2H_engs.append(n2h)
    n_ele_N2.append(e)

########################### periodic table position

slope, intercept, r_value, p_value, std_err = linregress(n_col, n_fe)

plt.scatter(n_col, n_fe)
#plt.scatter([cohesive_energies[a] for a in n_ele], n_fe)
plt.plot([min(n_col),max(n_col)],np.array([min(n_col),max(n_col)]) * slope + intercept)
plt.title('Valence Number vs Formation Energy')
plt.ylabel('2+ Substitution Formation Energy')
plt.xlabel('Valence Number')
plt.text(min(n_col,),max(n_fe), 'R$^2$ = {}'.format(r_value ** 2))
for c, e, t in zip(n_col, n_fe, n_ele):
    plt.text(c + 0.1, e + 0.1, t)

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

plt.savefig('N2_vs_N2H.pdf')
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

rows = [row_1_elements, row_2_elements, row_3_elements]

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
del rows[0]
del positive_row_energies[0]
del negative_row_energies[0]

for i, row in enumerate(rows):
    axs[i].bar(row, positive_row_energies[i], color='#0390fc')
    axs[i].bar(row, negative_row_energies[i], color='#0390fc')
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
        if element_sym == '':
            eng = 0
        elif element_sym not in element:
            eng = 0
        else:
            index = element.index(element_sym)
            eng = N2H_engs[index]
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
plt.savefig('N2H_adsorption_rows.pdf')
plt.show()


####################################

"""
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
        elif element_sym not in element:
            eng = 0
        else:
            index = element.index(element_sym)
            eng = N2H_engs[index] 
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
    axs[i].set_ylim([-1.5,2])
    if i == 1:
        axs[i].set_ylabel('N$_2$H Adsorption Energy (eV)')
    if i == 2:
        axs[i].set_xlabel('Element')
    if i == 0:
        axs[i].set_title('N$_2$H Adsorption Energy by Periodic Row')

plt.savefig('N2H_adsorption_rows')
plt.show()
"""
