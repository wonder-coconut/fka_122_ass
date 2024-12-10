import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/samples.txt").read().split('\n')
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(data[:,0],data[:,1],data[:,2])
ax.plot3D(data[:,3],data[:,4],data[:,5])
plt.show()