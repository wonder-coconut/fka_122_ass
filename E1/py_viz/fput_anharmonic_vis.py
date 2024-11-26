import numpy as np
import matplotlib.pyplot as plt

data = open('../op_text/fput_energies.txt').read().split('\n')
params = data.pop(0).split('\t')

timestep = float(params[0])
alpha = float(params[1])

data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

x = np.arange(0,len(data)) * timestep * 1000

energies = [[] for item in range(data.shape[1])]

for i in range(len(data)):
    energies[0].append((np.sum(data[:i,0]))/(i+1))
    energies[1].append((np.sum(data[:i,1]))/(i+1))
    energies[2].append((np.sum(data[:i,2]))/(i+1))
    energies[3].append((np.sum(data[:i,3]))/(i+1))
    energies[4].append((np.sum(data[:i,4]))/(i+1))

plt.yscale('log')
plt.plot(x,energies[0])
plt.plot(x,energies[1])
plt.plot(x,energies[2])
plt.plot(x,energies[3])
plt.plot(x,energies[4])

plt.xlabel('Time (ps)')
plt.ylabel('Energy (a.u.)')

plt.title(f'Energy evolution\nalpha = {alpha}')

plt.savefig(f'../op_figs/anharmonic_equipartition_{alpha}_res.png')