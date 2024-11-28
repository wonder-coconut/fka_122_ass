import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/energy_evol_0.000100.txt","r").read().split('\n')
data.pop()
timestep = float(data.pop(0))

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)
x = np.arange(len(data)) * timestep

kb = 1.380649e-23
joule_in_ev = 6.242e+18
kb_in_ev = kb * joule_in_ev

temp = data[:,1] * 2/(3 * 256 * kb_in_ev)
temp_avg = temp.mean()

fig,axes = plt.subplots(2,1,figsize = (10,10))

axes[0].plot(x,data[:,0], label='Potential energy')
axes[0].plot(x,data[:,1], label='Kinetic energy')
axes[0].plot(x,(data[:,0] + data[:,1]), label = 'Total energy')
axes[0].set_ylabel('Energy (eV)')
axes[0].legend(loc = 'center right')

axes[1].plot(x,temp)
axes[1].axhline(temp_avg,label=f'Average Temp = {temp_avg}',color = 'red')
axes[1].set_ylim(600,1000)
axes[1].legend(loc='upper right')
axes[1].set_ylabel('Temperature (K)')

axes[1].set_xlabel('Time (ps)')

axes[0].set_title(f'Energy and temperature evolution\n(timestep = {timestep})')

plt.savefig(f"../op_fig/e_evol_{timestep}.png")