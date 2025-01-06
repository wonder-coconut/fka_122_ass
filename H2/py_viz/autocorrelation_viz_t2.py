import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp

#exponential model for curve fitting
def exp_decay(x,a,b,c):
    return a * np.exp(-b*x) + c

#data read
data = open("../op_text/autocorrelation_ineff_t2.txt","r").read().split('\n')
data.pop()
data = [float(item) for item in data]

#curve fitting
x = np.arange(len(data))
res = sp.curve_fit(exp_decay, x, data)

params = [0,0,0]

params[0] = res[0][0]
params[1] = res[0][1]
params[2] = res[0][2]

y = exp_decay(x,params[0],params[1],params[2])

#finding stat inefficiency using correlation = 0.135
y_res = min(y, key=lambda x:abs(x - 0.135))
x_res = np.where(y == y_res)[0][0]

plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

plt.plot(x,data)
plt.plot(x,y)
plt.axhline(y_res,color='red',label = f'Autocorrelation = 0.135\nfor k = {x_res}')
plt.axvline(x_res,color='red')

plt.xlabel('K',fontsize='18')
plt.ylabel(r'Autocorrelation($r_k$)',fontsize='18')
plt.title('Statistical Inefficiency:\nAutocorrelation method',fontsize='19')
plt.legend(loc='upper right',fontsize='15')

plt.savefig('../op_fig/stat_ineff_t2.png')