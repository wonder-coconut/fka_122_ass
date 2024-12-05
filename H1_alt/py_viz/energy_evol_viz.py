import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/energy_evol.txt").read().split('\n')

data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)
x = np.arange(len(data))

kb = 1.380649e-23
joule_in_ev = 6.242e+18
kb_in_ev = kb * joule_in_ev

temp = (data[:,0] * 2/(3 * 256 * kb_in_ev))

print(np.mean(temp))

fig,axes = plt.subplots(2,1,figsize=(10,10))

axes[0].plot(x,data)
axes[0].plot(x,data[:,0] + data[:,1])

axes[1].plot(x,temp)
axes[1].axhline(np.mean(temp),color='red')

plt.show()
