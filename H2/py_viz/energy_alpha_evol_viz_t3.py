import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/energy_evol_alpha_t3.txt").read().split('\n')
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

plt.scatter(data[:,0],data[:,1])
plt.axhline(np.min(data[:,1]), color = 'red', label = f'Min E = {np.min(data[:,1])} a.u.' + r'at $\alpha$ =' + f'{data[np.argmin(data[:,1])][0]}')

#error plots
plt.vlines(data[:,0],data[:,1] + data[:,2],data[:,1] - data[:,2],color='black')
plt.scatter(data[:,0],data[:,1] + data[:,2],facecolors='none', edgecolors='black')
plt.scatter(data[:,0],data[:,1] - data[:,2],facecolors='none', edgecolors='black')

plt.legend(loc = 'upper right')
plt.xlabel(r"$\alpha$")
plt.ylabel("Local energy (a.u.)")

plt.title(r"Energy evolution w.r.t. variational parameter $\alpha$")

plt.savefig('../op_fig/energy_evol_t3.png')
#plt.show()