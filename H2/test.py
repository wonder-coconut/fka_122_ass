import numpy as np

energy = open("op_text/energy_evol_alpha_t4.txt","r").read().split('\n')
energy.pop()
energy = [float(item) for item in energy]
energy = np.array(energy)

variance = np.mean(energy*energy) - (np.mean(energy) * np.mean(energy))
error = np.sqrt(18 * variance/len(energy))
print(variance)
print(error)