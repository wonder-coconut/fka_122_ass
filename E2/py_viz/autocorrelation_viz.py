import matplotlib.pyplot as plt
import scipy.optimize as sp
import numpy as np

def exp_decay(x,a,b,c):
    return a * np.exp(-b*x) + c

data = open("../OUT","r").read().split('\n')
data.pop()
data = [float(item) for item in data]
data = np.array(data)


x = np.arange(len(data)) * 1000

#res = sp.curve_fit(exp_decay,x, data )

#params = [0,0,0]
#params[0] = res[0][0]
#params[1] = res[0][1]
#params[2] = res[0][2]

#y = exp_decay(x,params[0],params[1],params[2])

#y_res = min(y, key=lambda x:abs(x - 0.135))
#x_res = np.where(y == y_res)[0][0]

plt.plot(x,data)
#plt.plot(x,y)
#plt.axhline(y_res, color = 'red')
#plt.axvline(x_res, color = 'red', label = f'Autocorrelation = 0.135 for k = {x_res}')

plt.title('Autocorrelation evolution for increasing time lag')
plt.xlabel('K')
plt.ylabel('Autocorrelation')
#plt.legend(loc = 'upper right')


plt.savefig('../op_fig/autocorrelation_maxiterations.png')
plt.show()

