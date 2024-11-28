import numpy as np
import matplotlib.pyplot as plt

data = open("../op_text/equilibration_lp_t3.txt","r").read().split('\n')

timestep = float(data.pop(0))

data.pop()

data = [float(item) for item in data]
data = np.array(data)

x = np.arange(len(data)) * timestep * 50

plt.plot(x,data)
plt.xlabel('Time (ps)')
plt.ylabel('Lattice Constant (Å)')
plt.title('Time evolution of lattice constant during equilibration')

plt.savefig('../op_fig/lp_evol_t3.png')