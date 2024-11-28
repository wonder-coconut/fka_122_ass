import numpy as np
import matplotlib.pyplot as plt

data = open("../op_text/equilibration_pos_t4.txt").read().split('\n')
timestep = float(data.pop(0))

data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i].pop()
    data[i] = [float(item) for item in data[i]]
    res = []
    n = len(data[i])
    r = 3
    for j in range(int(n/r)):
        res.append([data[i][j*3], data[i][j*3+1], data[i][j*3+2]])
    data[i] = res

data = np.array(data)
x = np.arange(len(data)) * timestep

p0 = data[:,0]
p1 = data[:,1]
p2 = data[:,2]
p3 = data[:,3]
p4 = data[:,4]

fig,axes = plt.subplots(3,1,figsize=(10,15))

axes[0].plot(x, p0[:,0],label='p0')
axes[1].plot(x, p0[:,1],label='p0')
axes[2].plot(x, p0[:,2],label='p0')

axes[0].plot(x, p1[:,0],label='p1')
axes[1].plot(x, p1[:,1],label='p1')
axes[2].plot(x, p1[:,2],label='p1')

axes[0].plot(x, p2[:,0],label='p2')
axes[1].plot(x, p2[:,1],label='p2')
axes[2].plot(x, p2[:,2],label='p2')

axes[0].plot(x, p3[:,0],label='p3')
axes[1].plot(x, p3[:,1],label='p3')
axes[2].plot(x, p3[:,2],label='p3')

axes[0].plot(x, p4[:,0],label='p4')
axes[1].plot(x, p4[:,1],label='p4')
axes[2].plot(x, p4[:,2],label='p4')

axes[0].legend(loc = 'upper right')
axes[1].legend(loc = 'upper right')
axes[2].legend(loc = 'upper right')

axes[0].set_title('x axis')
axes[1].set_title('y axis')
axes[2].set_title('z axis')

plt.xlabel('Time (ps)')
axes[1].set_ylabel('Position (Å)')

axes[0].set_title('Position evolution with time')

plt.savefig('../op_fig/pos_viz_liquid.png')