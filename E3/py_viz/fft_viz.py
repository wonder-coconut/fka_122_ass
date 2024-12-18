import numpy as np
import matplotlib.pyplot as plt

dt = 0.025
param_set = 'high'

filename = f'../op_text/pspec_{param_set}_{dt}.txt'
data = open(filename,"r").read().split('\n')
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

plt.xlim(0,10)
plt.plot(data[:,0],data[:,1])
plt.xlabel("frequency")
plt.ylabel("amplitude")
plt.savefig(f'../op_fig/pspec_{param_set}_{dt}.png')