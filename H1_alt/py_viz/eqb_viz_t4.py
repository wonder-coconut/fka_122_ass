import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/eqb_evol_t4.txt").read().split('\n')

timestep = float(data.pop(0))
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)
x = np.arange(len(data)) * 0.001

fig,axes = plt.subplots(2,1,figsize = (10,10))
axes[0].plot(x,data[:,0])
axes[1].plot(x,data[:,1])

plt.show()