import matplotlib.pyplot as plt
import numpy as np

#data parse
data = open("../op_text/lp_evol_t3.txt","r").read().split('\n')

timestep = float(data.pop(0))
data.pop()

data = np.array([float(item) for item in data])

#runtime
x = np.arange(len(data)) * timestep

plt.plot(x,data)
plt.show()