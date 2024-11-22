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

x = np.arange(0,len(data)) * timestep

plt.plot(x,data[:,0],label = 'Normal mode #1')
plt.plot(x,data[:,1],label = 'Normal mode #2')
plt.plot(x,data[:,2],label = 'Normal mode #4')
plt.plot(x,data[:,3],label = 'Normal mode #4')
plt.plot(x,data[:,4],label = 'Normal mode #5')

plt.xlabel('Time (ps)')
plt.ylabel('Energy (a.u.)')
plt.legend(loc = 'center right')
plt.title(f'FPUT Normal Mode Energy Distribution\n alpha = {alpha}')

plt.savefig(f'../op_figs/fput_lin_energy_{alpha}.png')