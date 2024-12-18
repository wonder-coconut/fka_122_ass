import numpy as np
import matplotlib.pyplot as plt

dt = 0.001
param_set = 'low'

filename = f'../op_text/traj_out_{param_set}_{dt}.txt'
data = open(filename,"r").read().split('\n')
pressure = float(data.pop(0))
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

t = data[:,0]
x = data[:,1]
v = data[:,2]

fig,axes = plt.subplots(2,1,figsize = (10,10))

axes[0].scatter(t,x,s=10)
axes[0].plot(t,x)

axes[0].set_xlabel('Time (ms)')
axes[0].set_ylabel(r'Position ($\mu m$)')

axes[1].scatter(t,v,s=10)
axes[1].plot(t,v)

axes[1].set_xlabel('Time (ms)')
axes[1].set_ylabel(r'Velocity ($\mu m/ms$)')

axes[0].set_title(f'Trajectory for {pressure} KPa at ' + r'$\Delta t$ = ' + f'{dt}')

plt.savefig(f'../op_fig/traj_viz_{param_set}_{dt}.png')