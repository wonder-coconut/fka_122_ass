import matplotlib.pyplot as plt
import scipy.optimize as sp
import numpy as np

def quad(x,a,b,c):
    return a*x*x + b*x + c

data = open("../op_text/pot_e_lp.txt").read().split('\n')
out = open("../op_text/lp_0_k.txt","w")

data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

data[:,0] = data[:,0]*data[:,0]*data[:,0]  
data[:,1] /= 4*4*4

params = sp.curve_fit(quad, data[:,0], data[:,1])[0]

x = np.arange(data[:,0][0],data[:,0][-1],0.1)

y = quad(x, params[0],params[1],params[2])
lp_opt = x[np.argmin(y)]

out.write(str(pow(lp_opt,1/3)))

plt.scatter(data[:,0],data[:,1])
plt.plot(x,y)
plt.axvline(lp_opt, color='red')
plt.show()