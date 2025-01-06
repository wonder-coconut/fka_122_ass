import matplotlib.pyplot as plt
import scipy.optimize as sp
import numpy as np

def exp_decay(x,a,b,c):
    return -a * np.exp(-b*x) + c

data = open("../op_text/block_avg_ineff_t2.txt","r").read().split('\n')
data.pop()
data = [float(item) for item in data]
data = np.array(data)

x = np.arange(len(data))

res = sp.curve_fit(exp_decay,x,data)

params = [0,0,0]
params[0] = res[0][0]
params[1] = res[0][1]
params[2] = res[0][2]

y = exp_decay(x,params[0],params[1],params[2])
y_res = max(y)

plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

plt.scatter(x,data,color ='black',s=1)
plt.plot(x,y,color='red',label='Exponential decay fit')
plt.axhline(y_res,label=f'Convergence at {y_res:.4f}')

plt.legend(loc='lower right', fontsize='15')
plt.xlabel('Block size', fontsize='18')
plt.ylabel('Statistical Inefficiency', fontsize='18')
plt.title('Statistical Inefficiency:\nBlock Averaging Method', fontsize='19')

plt.savefig('../op_fig/block_avg_t2.png')
#plt.show()