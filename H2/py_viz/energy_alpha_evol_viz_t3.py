import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/energy_evol_alpha_t3.txt").read().split('\n')
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)


plt.scatter(data[:,0],data[:,1])
#plt.vlines(data[:1,0],data[:1,1] + data[:1,2],data[:1,1] - data[:1,2],color='black')
plt.show()