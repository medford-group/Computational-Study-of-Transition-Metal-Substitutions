import io, re
import csv
import numpy
import os

slab_list = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Hf','Ta','W','Re','Os','Ir','Pt','Au']

pathway_data_folder = 'pathway_data/'

associative_adsorbates = ['N2','N2*','N2H*','N2H2*','HNNH2*','H2NNH2*','NH2*+NH3*','2NH3*','2NH3']
associative_2_adsorbates = ['N2','N2*','N2H*','N2H2*','HNNH2*','H2NNH2*','2NH2*','NH2*+NH3','NH3*+NH3','2NH3']

associative_subs = [0,0,1,2,3,4,5,6,6]
associative_2_subs = [0,0,1,2,3,4,4,5,6,6]

e_transfers = [associative_subs, associative_2_subs]


path_names = ['_associative.csv','_associative_2.csv']
adsorbates = [associative_adsorbates,associative_2_adsorbates]
output_names = ["limiting_potential_associative.csv","limiting_potential_associative_2.csv"]

def make_float(num):
    num = num.replace(' ','').replace(',','.').replace("-", "-")
    return float(num)


for counter in range(len(path_names)):
    cur_path_name = path_names[counter]
    cur_adsorbates = adsorbates[counter]
    cur_output_name = output_names[counter]
    output_array = []
    for slab in slab_list:
        cur_path = pathway_data_folder+slab+cur_path_name
        cur_e_transfers = e_transfers[counter]
        if os.path.exists(cur_path):
            with open(cur_path) as csvfile:
                current_pathway_data = list(csv.reader(csvfile))
            
            largest_energy = 0
            start_step = "nothing"
            end_step = "nothing"
            previous_positive = 0
            cur_step_sum = 0
            largest_pot = 0
            temp_end_i = 0
            temp_start_i = 0
            prev_diff = 1
            for i in range(len(cur_adsorbates)-1):
                cur_step = make_float(current_pathway_data[i+1][0]) - make_float(current_pathway_data[i][0])
                if previous_positive == 1:
                    if cur_step > 0:
                        cur_step += cur_step_sum
                        cur_step_sum = cur_step
                        temp_end_step = cur_adsorbates[i+1]
                        temp_end_i = i + 1
                        previous_positive = 1
                    else:
                        previous_positive = 0
                        cur_step_sum = 0
                else:
                    temp_start_step = cur_adsorbates[i]
                    temp_start_i = i
                    temp_end_step = cur_adsorbates[i+1]
                    temp_end_i = i + 1
                    if cur_step > 0:
                        previous_positive = 1
                        cur_step_sum = cur_step

                if  cur_e_transfers[temp_end_i] - cur_e_transfers[temp_start_i] != 0:
                    cur_pot_step = cur_step / (cur_e_transfers[temp_end_i] - cur_e_transfers[temp_start_i])
                    prev_diff = cur_e_transfers[temp_end_i] - cur_e_transfers[temp_start_i]
                else:
                    cur_pot_step = cur_step / prev_diff

                if cur_pot_step > largest_pot:
                    largest_pot = cur_pot_step
                    start_step = temp_start_step
                    end_step = temp_end_step

            output_array.append([slab,str(-1 * largest_pot),start_step,end_step])
                
            

    numpy.savetxt(pathway_data_folder+cur_output_name,output_array,delimiter=",",fmt='%s')
