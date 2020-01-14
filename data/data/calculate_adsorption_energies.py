import io, re
import csv
import numpy

slab_list = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Hf','Ta','W','Re','Os','Ir','Pt','Au']

N2_energy = -554.5670507
H2_energy = -32.88420862

def make_float(num):
    num = num.replace(' ','').replace(',','.').replace("-", "-")
    return float(num)

def get_energies(data_file, gas_energy):
    with open('2+_Clean_energies.csv') as csvfile:
        pristine_energies = list(csv.reader(csvfile))
    with open(data_file) as csvfile:
        adsorbed_energies = list(csv.reader(csvfile))

    energies = []

    found = 0
    for name in slab_list:
            found = 0
            for energy in pristine_energies:
                    if energy[0] == name:
                            for energy_2 in adsorbed_energies:
                                    if energy_2[0] == name:
                                            current_energy_1 = make_float(energy_2[1])-make_float(energy[1])- gas_energy
                                            found = found + 1
                            if found == 1:
                                    energies.append([name,current_energy_1])
    return energies


def sa_get_energies(data_file, gas_energy):
    with open('single_atom.csv') as csvfile:
        pristine_energies = list(csv.reader(csvfile))
    with open(data_file) as csvfile:
        adsorbed_energies = list(csv.reader(csvfile))

    energies = []

    found = 0
    for name in slab_list:
            found = 0
            for energy in pristine_energies:
                    if energy[0] == name:
                            for energy_2 in adsorbed_energies:
                                    if energy_2[0] == name:
                                            current_energy_1 = make_float(energy_2[1])-make_float(energy[1])- gas_energy
                                            found = found + 1
                            if found == 1:
                                    energies.append([name,current_energy_1])
    return energies

energies = get_energies('2+_N2_energies.csv',N2_energy)
numpy.savetxt("N2.csv",energies,delimiter=",",fmt='%s')

energies = get_energies('2+_N2H_energies.csv',N2_energy+H2_energy/2)
numpy.savetxt("N2H.csv",energies,delimiter=",",fmt='%s')

energies = get_energies('2+_H2NNH2_energies.csv',N2_energy+2*H2_energy)
numpy.savetxt("H2NNH2.csv",energies,delimiter=",",fmt='%s')

energies = get_energies('2+_HNNH_energies.csv',N2_energy+H2_energy)
numpy.savetxt("HNNH.csv",energies,delimiter=",",fmt='%s')

energies = get_energies('2+_N_energies.csv',N2_energy/2)
numpy.savetxt("N.csv",energies,delimiter=",",fmt='%s')

energies = get_energies('2+_N2H2_energies.csv',N2_energy+H2_energy)
numpy.savetxt("N2H2.csv",energies,delimiter=",",fmt='%s')

energies = get_energies('2+_N2H3_energies.csv',N2_energy+H2_energy*3/2)
numpy.savetxt("N2H3.csv",energies,delimiter=",",fmt='%s')

energies = get_energies('2+_NH_energies.csv',N2_energy/2 + H2_energy/2)
numpy.savetxt("NH.csv",energies,delimiter=",",fmt='%s')

energies = get_energies('2+_NH2_energies.csv',N2_energy/2 + H2_energy)
numpy.savetxt("NH2.csv",energies,delimiter=",",fmt='%s')

energies = get_energies('2+_NH3_energies.csv',N2_energy/2 + H2_energy*3/2)
numpy.savetxt("NH3.csv",energies,delimiter=",",fmt='%s')
