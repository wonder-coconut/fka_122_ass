import numpy as np
import matplotlib.pyplot as plt

data = open('../op_text/power_spectrum.txt').read().split('\n')

params = data.pop(0).split()
N = float(params[0])
dt = float(params[1])

data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)
fig,axes = plt.subplots(3,1,figsize=(10,10))

axes[0].plot(data[:,0],data[:,1],label = 'Particle 0')
axes[1].plot(data[:,0],data[:,2],label = 'Particle 1')
axes[2].plot(data[:,0],data[:,3],label = 'Particle 2')

axes[0].legend(loc = 'upper right')
axes[1].legend(loc = 'upper right')
axes[2].legend(loc = 'upper right')

plt.xlabel("Frequency (rad/s)")

axes[0].set_ylabel("Power (a.u.)")
axes[1].set_ylabel("Power (a.u.)")
axes[2].set_ylabel("Power (a.u.)")

axes[0].set_xlim(-600,600)
axes[1].set_xlim(-600,600)
axes[2].set_xlim(-600,600)

axes[0].set_title(f'Fourier Transform Power Spectrum\nN = {N}, dt = {dt}')

plt.savefig('../op_figs/coupled_harm_osc_fft.png')