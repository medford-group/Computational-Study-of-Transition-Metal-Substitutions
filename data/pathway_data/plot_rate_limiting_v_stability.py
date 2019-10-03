import io, re
import matplotlib.pyplot as plt
import csv
import numpy as np

slab_list = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Hf','Ta','W','Re','Os','Ir','Pt','Au']


def make_float(num):
    num = num.replace(' ','').replace(',','.').replace("-", "-")
    return float(num)

with open('stability_energy.csv') as csvfile:
    stability_energies = list(csv.reader(csvfile))
with open('rate_limiting_steps_associative_2.csv') as csvfile:
    rate_limiting_energies = list(csv.reader(csvfile))

plot_stability_energies = []
plot_rate_limiting_energies = []
name_list = []

found = 0
for name in slab_list:
        found = 0
        for energy in stability_energies:
                if energy[0] == name:
                        current_energy_1 = make_float(energy[1])
                        for energy_2 in rate_limiting_energies:
                                if energy_2[0] == name:
                                        current_energy_2 = make_float(energy_2[1])
                                        current_name = name+" "+energy_2[2]+"->"+energy_2[3]
                                        found = found + 1
                        if found == 1:
                                plot_stability_energies.append(current_energy_1)
                                plot_rate_limiting_energies.append(current_energy_2)
                                name_list.append(current_name)


plt.scatter(plot_stability_energies, plot_rate_limiting_energies, alpha=0.5)
#plt.show()

#plt.xticks(y_pos,names)
plt.ylabel('Rate Limiting Step Energies')
plt.xlabel('2+ Formation Energies')
plt.title('Rate Limiting Step Energies v. 2+ Formation Energies')
for i, txt in enumerate(name_list):
    plt.annotate(txt, (plot_stability_energies[i], plot_rate_limiting_energies[i]))

plt.show()
