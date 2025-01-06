import matplotlib.pyplot as plt
import numpy as np

data1 = open("../op_text/vmc_energy_0.50_t4.txt","r").read().split('\n')
data1.pop()
data2 = open("../op_text/vmc_energy_0.75_t4.txt","r").read().split('\n')
data2.pop()
data3 = open("../op_text/vmc_energy_0.90_t4.txt","r").read().split('\n')
data3.pop()
data4 = open("../op_text/vmc_energy_1.00_t4.txt","r").read().split('\n')
data4.pop()

for i in range(len(data1)):
    data1[i] = data1[i].split('\t')
    data1[i] = [float(item) for item in data1[i]]
    data2[i] = data2[i].split('\t')
    data2[i] = [float(item) for item in data2[i]]
    data3[i] = data3[i].split('\t')
    data3[i] = [float(item) for item in data3[i]]
    data4[i] = data4[i].split('\t')
    data4[i] = [float(item) for item in data4[i]]

data1 = np.array(data1)
data2 = np.array(data2)
data3 = np.array(data3)
data4 = np.array(data4)

x = np.arange(len(data1))

fig,axes = plt.subplots(2,1,figsize=(10,10))

axes[0].plot(x,data1[:,0], label = r'$\beta = $' + f'0.50')
axes[0].plot(x,data2[:,0], label = r'$\beta = $' + f'0.75')
axes[0].plot(x,data3[:,0], label = r'$\beta = $' + f'0.90')
axes[0].plot(x,data4[:,0], label = r'$\beta = $' + f'1.00')
axes[0].axhline(data4[-1][0], color = 'black', label = r'$\alpha_{opt} = $' + f'{data4[-1][0]:.4f}')
axes[0].legend(loc='lower right',fontsize='11')
axes[0].set_ylabel(r'$\alpha$',fontsize='15')

axes[1].plot(x,data1[:,1], label = r'$\beta = $' + f'0.50')
axes[1].plot(x,data2[:,1], label = r'$\beta = $' + f'0.75')
axes[1].plot(x,data3[:,1], label = r'$\beta = $' + f'0.90')
axes[1].plot(x,data4[:,1], label = r'$\beta = $' + f'1.00')
axes[1].axhline(data4[-1][1], color='black', label = r'$E_{opt} = $' + f'{data4[-1][1]:.4f}')
axes[1].set_ylim(np.mean(data4[:,1])-0.01,np.mean(data1[:,1])+0.01)
axes[1].legend(loc='lower right',fontsize='11')
axes[1].set_ylabel('Mean Energy ' + r'$E_h$',fontsize='15')
axes[1].set_xlabel('VMC step number',fontsize='15')

axes[0].set_title(r"VMC optimization of $\alpha$",fontsize='22')

plt.savefig('../op_fig/vmc_evol_t4.png')