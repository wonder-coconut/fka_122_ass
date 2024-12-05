import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/lp_evol_t4.txt","r").read().split('\n')

timestep = float(data.pop(0))
data.pop()

data = np.array([float(item) for item in data])

x = np.arange(len(data))

plt.plot(x,data)
plt.show()