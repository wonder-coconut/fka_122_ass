import numpy as np
import matplotlib.pyplot as plt

data = open("../op_text/integral_stats.txt","r").read().split('\n')
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)
x = np.pow(10,np.arange(1,len(data) + 1))

fig,axes = plt.subplots(2,1, figsize=(6.5,9.5))

y_true = 1/2 - 1/3

axes[0].set_title("Monte Carlo Integration Comparison\nUniform sampling vs Importance Sampling")

axes[0].set_xscale('log')
axes[0].axhline(y_true, color='red', label = 'True value')

axes[0].scatter(x,data[:,0] + -data[:,1],marker = '_', color = 'black')
axes[0].scatter(x,data[:,0] + data[:,1],marker = '_', color = 'black')
axes[0].vlines(x, data[:,0] + -data[:,1], data[:,0] + data[:,1], color = 'black', label = 'Error Margin')
axes[0].scatter(x,data[:,0], label = 'MC Integral value')

axes[0].legend(loc = 'upper right')
axes[0].set_xlabel('Number of samples')

axes[0].set_ylim(-0.05 + y_true,0.05 + y_true)
axes[0].set_xlim(1,1e6)


axes[1].set_xscale('log')
axes[1].axhline(y_true, color='red',label = 'True value')

axes[1].scatter(x,data[:,2] + -data[:,3],marker = '_', color = 'black')
axes[1].scatter(x,data[:,2] + data[:,3],marker = '_', color = 'black')
axes[1].vlines(x, data[:,2] + -data[:,3], data[:,2] + data[:,3], color = 'black', label = 'Error Margin')
axes[1].scatter(x,data[:,2], label = 'MC Integral value')

axes[1].legend(loc = 'upper right')
axes[1].set_xlabel('Number of samples')

axes[1].set_ylim(-0.05 + y_true,0.05 + y_true)
axes[1].set_xlim(1,1e6)

plt.savefig("../op_fig/mc_integral_stats.png")
