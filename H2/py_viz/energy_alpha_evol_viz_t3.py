import matplotlib.pyplot as plt
import scipy.optimize as sp
import numpy as np

def cubic_model(x,a,b,c,d):
    return d*x*x*x + a*x*x + b*x + c

data = open("../op_text/energy_evol_alpha_t3.txt").read().split('\n')
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

#curve fitting
x = np.linspace(min(data[:,0]),max(data[:,0]),len(data))
res = sp.curve_fit(cubic_model, x, data[:,1])

params = [0,0,0,0]
params[0] = res[0][0]
params[1] = res[0][1]
params[2] = res[0][2]
params[3] = res[0][3]

y = cubic_model(x,params[0],params[1],params[2],params[3])
y_res = min(y)

plt.figure(figsize=(8,5))
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)

plt.scatter(data[:,0],data[:,1], label = 'Average Energy')
plt.plot(x,y,color='red', label='Cubic curve fit')
plt.axhline(min(y), color='purple',label=f'Minimum E ={min(y):.3f}' + r' $E_h$' + '\n'+r'at $\alpha$ = ' + f'{x[np.argmin(y)]:.3f}')

plt.vlines(data[:,0],data[:,1] + data[:,2],data[:,1] - data[:,2],color='black', label='Error range',alpha = 0.35)
plt.scatter(data[:,0],data[:,1] + data[:,2],marker='s', color ='black',alpha=0.5)
plt.scatter(data[:,0],data[:,1] - data[:,2],marker='s', color ='black',alpha=0.5)


plt.legend(loc = 'upper right', fontsize = '12')
plt.xlabel(r"$\alpha$", fontsize = '17')
plt.ylabel(r"Local energy ($E_h$)", fontsize = '17')

plt.title(r"Energy evolution w.r.t. variational parameter $\alpha$",fontsize='18')

plt.savefig('../op_fig/energy_evol_t3.png')
#plt.show()