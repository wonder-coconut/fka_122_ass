import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/pos_evol_t3.txt","r").read().split('\n')

timestep = float(data.pop(0))
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i].pop()
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

n = 5
pos = []
i = 0
while(i < n*3):
    pos.append([data[:,i],data[:,i+1],data[:,i+2]])
    i += 3

pos = np.array(pos)

x = np.arange(len(data))

fig,axes = plt.subplots(3,1,figsize=(10,10))

axes[0].plot(x,pos[0][0])
axes[0].plot(x,pos[1][0])
axes[0].plot(x,pos[2][0])
axes[0].plot(x,pos[3][0])
axes[0].plot(x,pos[4][0])

axes[1].plot(x,pos[0][1])
axes[1].plot(x,pos[1][1])
axes[1].plot(x,pos[2][1])
axes[1].plot(x,pos[3][1])
axes[1].plot(x,pos[4][1])

axes[2].plot(x,pos[0][2])
axes[2].plot(x,pos[1][2])
axes[2].plot(x,pos[2][2])
axes[2].plot(x,pos[3][2])
axes[2].plot(x,pos[4][2])

plt.show()