import numpy as np
import matplotlib.pyplot as plt

data = open('../op_text/output.txt','r').read().split('\n')

params = data.pop(0).split()

n_samples= int(params[0])
timestep = float(params[1])

data.pop(0)
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)
x = np.arange(0,len(data)) * timestep

fig,axes = plt.subplots(2,1,figsize=(20,7))

axes[0].plot(x,data[:,0],label = f'Particle 0')
axes[0].plot(x,data[:,1],label = f'Particle 1')
axes[0].plot(x,data[:,2],label = f'Particle 2')

axes[0].set_xlabel('Time (ps)')
axes[0].set_ylabel('Displacement (Å)')
axes[0].legend(loc = 'upper right')

axes[0].set_ylim(-0.008,0.016)
axes[0].set_xlim(0,x[-1])

axes[0].set_title(f"N = {n_samples}, timestep = {timestep}")

axes[1].plot(x,data[:,3] * 1000,label = f'Potential energy')
axes[1].plot(x,data[:,4] * 1000,label = f'Kinetic Energy')
axes[1].plot(x,(data[:,3] + data[:,4]) * 1000,label = f'Total Energy')

axes[1].set_xlabel('Time (ps)')
axes[1].set_ylabel('Energy (meV)')
axes[1].legend(loc = 'upper right')

axes[1].set_ylim(-2,8)
axes[1].set_xlim(0,x[-1])

plt.savefig('../op_figs/coupled_harm_osc_char.png')