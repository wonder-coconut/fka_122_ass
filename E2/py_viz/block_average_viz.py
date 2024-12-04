import matplotlib.pyplot as plt
import scipy.optimize as sp
import numpy as np

def exp_decay(x,a,b,c):
    return -a * np.exp(-b*x) + c

data = open("../OUT","r").read().split('\n')
data.pop()
data = [float(item) for item in data]
data = np.array(data)

data = data[:10000]

x = np.arange(len(data))

res = sp.curve_fit(exp_decay,x,data)

params = [0,0,0]
params[0] = res[0][0]
params[1] = res[0][1]
params[2] = res[0][2]

y = exp_decay(x,params[0],params[1],params[2])
y_res = max(y)

plt.plot(x,data,color ='black')
plt.plot(x,y,color='red')
plt.axhline(y_res,label=f'Convergence at {y_res}')

plt.legend(loc='lower right')
plt.xlabel('Block size')
plt.ylabel('Statistical Inefficiency')
plt.title('Block Averaged Statistical Inefficiency evolution')

plt.savefig('../op_fig/block_avg.png')
plt.show()