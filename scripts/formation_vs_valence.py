from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress
from data import fe, column, element, cohesive_energies, N2_engs, N2H_engs

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
    n_N2_engs.append(n2)
    nn_N2H_engs.append(n2h)
    n_ele_N2.append(e)


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

slope, intercept, r_value, p_value, std_err = linregress([cohesive_energies[a] for a in n_ele], n_N2H_engs)

plt.title('Metal Cohesive Energy vs N$_2$H Binding Energy')
plt.ylabel('N$_2$H Binding Energy')
plt.xlabel('Metal Cohesive Energy (eV/atom)')
co_e = [cohesive_energies[a] for a in n_ele]
plt.text(max(co_e)-0.75,max( n_N2H_engs), 'R$^2$ = {}'.format(round(r_value ** 2,2)))
plt.scatter(co_e, n_N2H_engs)
plt.plot([min(co_e),max(co_e)],np.array([min(co_e),max(co_e)]) * slope + intercept)
for c, e, t in zip(co_e, n_N2H_engs, n_ele):
    plt.text(c + 0.05, e + 0.05, t)

plt.savefig('cohesive_eng_vs_N2H.pdf')
plt.show()



slope, intercept, r_value, p_value, std_err = linregress(n_N2_engs, nn_N2H_engs)

plt.title('N$_2$ Binding Energy vs N$_2$H Binding Energy')
plt.ylabel('N$_2$H Binding Energy (eV)')
plt.xlabel('N$_2$ Binding Energy (eV)')
plt.text(max(n_N2_engs)-0.75,max(nn_N2H_engs), 'R$^2$ = {}'.format(round(r_value ** 2,2)))
plt.scatter(n_N2_engs, nn_N2H_engs)
plt.plot([min(n_N2_engs),max(n_N2_engs)],np.array([min(n_N2_engs),max(n_N2_engs,)]) * slope + intercept)
for c, e, t in zip(n_N2_engs, nn_N2H_engs, n_ele):
    plt.text(c + 0.05, e + 0.05, t)

plt.savefig('N2_vs_N2H.pdf')
plt.show()
