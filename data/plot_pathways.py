import io, re
import csv
import numpy as np

from catmap.analyze.analysis_base import MechanismPlot
import matplotlib.pyplot as plt

data_folder = 'data/'
output_folder = 'plots/'
output_energies_folder = 'pathway_data/'

slab_list = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Hf','Ta','W','Re','Os','Ir','Pt','Au']


#adsorption energies calculated for Ti spots. Used when the N-N bond dissociates and one molecule hops to another adsorption spot.
#Ti_N = 4.602056151
#Ti_NH = 4.244323642
#Ti_NH2 = 1.532267866
#Ti_NH3 = -0.8112892265
Ti_N = 4.783985488
Ti_NH = 4.029732435
Ti_NH2 = 0.9095197761
Ti_NH3 = -1.665674942

NH3 = -0.7394391845

#electric potential correction for adding hydrogen. Approximating as valence band energy of pure titania
#hydrogen_potential = -0.15
#not including this effect currently
hydrogen_potential = 0

#free energy corrections for different adsorbates
H2NNH2_corr = 1.41847290209
HNNH_corr = 0.724334755451
N_corr = 0.0627455534737
N2_corr = 0.0356539561469
N2H_corr = 0.482259647646
N2H2_corr = 0.755080777255
N2H3_corr = 1.06805895958
NH_corr = 1.0859324152
NH2_corr = 0.681757082272
NH3_corr = 0.926362617529

#free energy corrections of gas phase molecules
N2_gas_corr = -0.351226774491
H2_gas_corr = -0.0394892748343
NH3_gas_corr = 0.426183909776

#correct NH3 gas
#NH3 = NH3 + NH3_corr + 3*hydrogen_potential
NH3 = NH3 + NH3_gas_corr + 3*hydrogen_potential


#associative_adsorbates = ['N2','N2*','N2H*','N2H2*','HNNH2*','H2NNH2*','NH2*+NH3*','2NH3*','2NH3']
#associative_columns = [0,1,2,3,4,7,8]
associative_adsorbates = ['N2','N2*','N2H*','HNNH*','HNNH2*','H2NNH2*','2NH2*','NH2*+NH3','NH3*+NH3','2NH3']
associative_columns = [0,1,2,21,4,13,19,20]

dissociative_adsorbates = ['N2','N2*','2N*','2NH*','2NH2*','2NH3*','2NH3']
dissociative_columns = [0,11,12,13,8]
distal_1_adsorbates = ['N2','N2*','N2H*','N2H2*','N*+NH2*','N*+NH3*','NH*+NH3*','NH2*+NH3*','2NH3*','2NH3']
distal_1_columns = [0,1,2,14,15,16,7,8]
distal_2_adsorbates = ['N2','N2*','N2H*','N2H2*','HNNH2*','NH*+NH2*','NH*+NH3*','NH2*+NH3*','2NH3*','2NH3']
distal_2_columns = [0,1,2,3,17,16,7,8]

#associative_2_adsorbates = ['N2','N2*','N2H*','N2H2*','HNNH2*','H2NNH2*','NH*+NH3','NH2*+NH3','NH3*+NH3','2NH3']
associative_2_adsorbates = ['N2','N2*','N2H*','N2H2*','HNNH2*','H2NNH2*','2NH2*','NH2*+NH3','NH3*+NH3','2NH3']
#associative_2_columns = [0,1,2,3,4,18,19,20]
associative_2_columns = [0,1,2,3,4,13,19,20]

with open(data_folder+'H2NNH2.csv') as csvfile:
    H2NNH2_energies = list(csv.reader(csvfile))
with open(data_folder+'HNNH.csv') as csvfile:
    HNNH_energies = list(csv.reader(csvfile))
with open(data_folder+'N.csv') as csvfile:
    N_energies = list(csv.reader(csvfile))
with open(data_folder+'N2.csv') as csvfile:
    N2_energies = list(csv.reader(csvfile))
with open(data_folder+'N2H.csv') as csvfile:
    N2H_energies = list(csv.reader(csvfile))
with open(data_folder+'N2H2.csv') as csvfile:
    N2H2_energies = list(csv.reader(csvfile))
with open(data_folder+'N2H3.csv') as csvfile:
    HNNH2_energies = list(csv.reader(csvfile))
with open(data_folder+'NH.csv') as csvfile:
    NH_energies = list(csv.reader(csvfile))
with open(data_folder+'NH2.csv') as csvfile:
    NH2_energies = list(csv.reader(csvfile))
with open(data_folder+'NH3.csv') as csvfile:
    NH3_energies = list(csv.reader(csvfile))

def make_float(num):
    num = num.replace(' ','').replace(',','.').replace("-", "-")
    return float(num)

def plot_pathway(energies,adsorbates,title,file_name):
    ax = plt.subplot()

    barriers = []
    for i in energies:
        barriers.append(i-1)
    barriers.pop(0)
    plot = MechanismPlot(energies, labels = adsorbates, barriers = barriers)

    plot.energy_mode = 'absolute'

    plot.draw(ax = ax)

    plt.title(title)

    plt.ylabel('Free Energy (eV)')

    plt.xlabel('Reaction Coordinate')

    plt.gca().axes.get_xaxis().set_ticks([])
    plt.ylim(top = max(energies)+1)

    plt.savefig(output_folder+file_name)
    plt.close()
    return

pathway_array = [[0 for col in range(0,22)] for row in range(0,26)]
#array, each row is a substitute, each column is the energy of an adsorbate combination.
#adding corrections for free energy and number of hydrogens attached.
#Indexed rather arbitrarily:
#column 0: N2*
for energy in N2_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][0] = make_float(energy[1]) + N2_corr - N2_gas_corr
#1: N2H*
for energy in N2H_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][1] = make_float(energy[1]) + N2H_corr - N2_gas_corr - H2_gas_corr/2 + hydrogen_potential
#2: N2H2*
for energy in N2H2_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][2] = make_float(energy[1]) + N2H2_corr - N2_gas_corr - H2_gas_corr + 2*hydrogen_potential
#3: HNNH2*
for energy in HNNH2_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][3] = make_float(energy[1]) + N2H3_corr - N2_gas_corr - H2_gas_corr*3/2 + 3*hydrogen_potential
#4: H2NNH2*
for energy in H2NNH2_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][4] = make_float(energy[1]) + H2NNH2_corr - N2_gas_corr - 2*H2_gas_corr + 4*hydrogen_potential
#5: NH2*
for energy in NH2_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][5] = make_float(energy[1]) + NH2_corr - N2_gas_corr/2 - H2_gas_corr + 2*hydrogen_potential
#6: NH3*
for energy in NH3_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][6] = make_float(energy[1]) + NH3_corr - N2_gas_corr/2 - H2_gas_corr*3/2 + 3*hydrogen_potential
#7: NH2* + NH3*
for energy in NH2_energies:
    for energy_2 in NH3_energies:
        for i in range(0,26):
            if energy[0] == slab_list[i] and energy_2[0] == slab_list[i]:
                pathway_array[i][7] = min(make_float(energy[1])+Ti_NH3,make_float(energy_2[1])+Ti_NH2) + NH2_corr + NH3_corr - N2_gas_corr - H2_gas_corr*5/2 + 5*hydrogen_potential
#8: 2NH3*
for energy in NH3_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][8] = make_float(energy[1])+Ti_NH3 + NH3_corr*2 - N2_gas_corr - H2_gas_corr*3 + 6*hydrogen_potential
#9: N*
for energy in N_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][9] = make_float(energy[1]) + N_corr - N2_gas_corr/2
#10: NH*
for energy in NH_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][10] = make_float(energy[1]) + NH_corr - N2_gas_corr/2 - H2_gas_corr/2 + hydrogen_potential
#11: 2N*
for energy in N_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][11] = make_float(energy[1])+Ti_N + N_corr*2 - N2_gas_corr
#12: 2NH*
for energy in NH_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][12] = make_float(energy[1])+Ti_NH + NH_corr*2 - N2_gas_corr - H2_gas_corr + 2*hydrogen_potential
#13: 2NH2*
for energy in NH2_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][13] = make_float(energy[1])+Ti_NH2 + NH2_corr*2 - N2_gas_corr - H2_gas_corr*2 + 4*hydrogen_potential
            #pathway_array[i][13] = 2*make_float(energy[1])+Ti_NH2 + NH2_corr*2 - N2_gas_corr - H2_gas_corr*2 + 4*hydrogen_potential
#14: N* + NH2*
for energy in N_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][14] = make_float(energy[1])+Ti_NH2 + N_corr + NH2_corr - N2_gas_corr - H2_gas_corr + 2*hydrogen_potential
#15: N* + NH3*
for energy in N_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][15] = make_float(energy[1])+Ti_NH3 + N_corr + NH3_corr - N2_gas_corr - H2_gas_corr*3/2 + 3*hydrogen_potential
#16: NH* + NH3*
for energy in NH_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][16] = make_float(energy[1])+Ti_NH3 + NH_corr + NH3_corr - N2_gas_corr - H2_gas_corr*2 + 4*hydrogen_potential
#17: NH* + NH2*
for energy in NH_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][17] = make_float(energy[1])+Ti_NH2 + NH_corr + NH2_corr - N2_gas_corr - H2_gas_corr*3/2 + 3*hydrogen_potential
#18: NH* + NH3
for energy in NH_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][18] = make_float(energy[1])+NH_corr - N2_gas_corr/2 - H2_gas_corr*1/2 + hydrogen_potential + NH3

#19: NH2* + NH3
for energy in NH2_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][19] = make_float(energy[1])+NH2_corr - N2_gas_corr/2 - H2_gas_corr*2/2 + 2*hydrogen_potential + NH3

#20: NH3* + NH3
for energy in NH3_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][20] = make_float(energy[1])+NH3_corr - N2_gas_corr/2 - H2_gas_corr*3/2 + 3*hydrogen_potential + NH3

#21: HNNH*
for energy in HNNH_energies:
    for i in range(0,26):
        if energy[0] == slab_list[i]:
            pathway_array[i][21] = make_float(energy[1]) + HNNH_corr - N2_gas_corr - H2_gas_corr*2 + 2*hydrogen_potential




#make plots and write associative data to file
for i in range(0,26):
    slab = slab_list[i]
    #associative
    current_energies = [0]
    missing_step = 0
    for column in associative_columns:
        current_energies.append(pathway_array[i][column])
        if pathway_array[i][column] == 0 or abs(pathway_array[i][column]) > 100:
            missing_step = 1
    if missing_step == 0:
        current_energies.append(NH3*2)
        plot_pathway(current_energies,associative_adsorbates,slab+' Associative Pathway',slab+'_associative.pdf')
        np.savetxt(output_energies_folder+slab+"_associative.csv",current_energies,delimiter=",",fmt='%s')

    #associative 2
    current_energies = [0]
    missing_step = 0
    for column in associative_2_columns:
        current_energies.append(pathway_array[i][column])
        if pathway_array[i][column] == 0 or abs(pathway_array[i][column]) > 100:
            missing_step = 1
    if missing_step == 0:
        current_energies.append(NH3*2)
        plot_pathway(current_energies,associative_2_adsorbates,slab+' Second Associative Pathway',slab+'_associative_2.pdf')
        np.savetxt(output_energies_folder+slab+"_associative_2.csv",current_energies,delimiter=",",fmt='%s')

    #dissociative
    current_energies = [0]
    missing_step = 0
    for column in dissociative_columns:
        current_energies.append(pathway_array[i][column])
        if pathway_array[i][column] == 0 or abs(pathway_array[i][column]) > 100:
            missing_step = 1
    if missing_step == 0:
        current_energies.append(NH3*2)
        plot_pathway(current_energies,dissociative_adsorbates,slab+' Dissociative Pathway',slab+'_dissociative.pdf')
    #distal 1
    current_energies = [0]
    missing_step = 0
    for column in distal_1_columns:
        current_energies.append(pathway_array[i][column])
        if pathway_array[i][column] == 0 or abs(pathway_array[i][column]) > 100:
            missing_step = 1
    if missing_step == 0:
        current_energies.append(NH3*2)
        plot_pathway(current_energies,distal_1_adsorbates,slab+' First Distal Pathway',slab+'_distal_1.pdf')
    #distal 2
    current_energies = [0]
    missing_step = 0
    for column in distal_2_columns:
        current_energies.append(pathway_array[i][column])
        if pathway_array[i][column] == 0 or abs(pathway_array[i][column]) > 100:
            missing_step = 1
    if missing_step == 0:
        current_energies.append(NH3*2)
        plot_pathway(current_energies,distal_2_adsorbates,slab+' Second Distal Pathway',slab+'_distal_2.pdf')



