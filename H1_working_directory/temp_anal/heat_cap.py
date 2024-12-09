import numpy as np
import matplotlib.pyplot as plt

#data parse
data1 = open("anal_files/energy_evol_t4_temp.txt").read().split('\n')
data2 = open("anal_files/energy_evol_t3_temp.txt").read().split('\n')

data1.pop()
data2.pop()

for i in range(len(data1)):
    data1[i] = data1[i].split('\t')
    data2[i] = data2[i].split('\t')

    data1[i] = [float(item) for item in data1[i]]
    data2[i] = [float(item) for item in data2[i]]

data1 = np.array(data1)
data2 = np.array(data2)

#simulation constants
k_B = f #boltzmann constant in eV/K
n_atoms = 256

pot_e_liquid = data1[:,0]
kin_e_liquid = data1[:,1]

pot_e_solid = data2[:,0]
kin_e_solid = data2[:,1]

#average temperature from kin_e
T_liquid = 2/(3 * n_atoms * k_B) * kin_e_liquid
T_avg_liquid = np.mean(T_liquid)

T_solid = 2/(3 * n_atoms * k_B) * kin_e_solid
T_avg_solid = np.mean(T_solid)

#variance of energies
pot_e_var_liquid = np.mean(np.pow(pot_e_liquid,2)) - np.pow(np.mean(pot_e_liquid),2)
kin_e_var_liquid = np.mean(np.pow(kin_e_liquid,2)) - np.pow(np.mean(kin_e_liquid),2)

pot_e_var_solid = np.mean(np.pow(pot_e_solid,2)) - np.pow(np.mean(pot_e_solid),2)
kin_e_var_solid = np.mean(np.pow(kin_e_solid,2)) - np.pow(np.mean(kin_e_solid),2)

#heat capacity
cv_pot_liquid = 1.5 * n_atoms * k_B / (1 - 2/(3 * n_atoms * k_B * k_B * T_avg_liquid * T_avg_liquid) * pot_e_var_liquid)
cv_kin_liquid = 1.5 * n_atoms * k_B / (1 - 2/(3 * n_atoms * k_B * k_B * T_avg_liquid * T_avg_liquid) * kin_e_var_liquid)

cv_pot_solid = 1.5 * n_atoms * k_B / (1 - 2/(3 * n_atoms * k_B * k_B * T_avg_solid * T_avg_solid) * pot_e_var_solid)
cv_kin_solid = 1.5 * n_atoms * k_B / (1 - 2/(3 * n_atoms * k_B * k_B * T_avg_solid * T_avg_solid) * kin_e_var_solid)

print(f"Heat Capacity for liquid (pot e):\t{cv_pot_liquid} eV/K")
print(f"Heat Capacity for liquid (kin e):\t{cv_kin_liquid} eV/K")
print()
print(f"Heat Capacity for solid (pot e):\t{cv_pot_solid} eV/K")
print(f"Heat Capacity for solid (kin e):\t{cv_kin_solid} eV/K")