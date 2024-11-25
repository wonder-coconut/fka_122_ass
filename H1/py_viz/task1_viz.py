import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp

def quad_model(x,a,b,c):
    return a*x*x + b*x + c

data = open("../op_text/pot_e_dist.txt","r").read().split('\n')
op_write = open("../op_text/lp_0k.txt","w")

data.pop()
for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)
data[:,0] = np.power(data[:,0],3)
data[:,1] /= 4*4*4

res = sp.curve_fit(quad_model,data[:,0], data[:,1])

a_opt = res[0][0]
b_opt = res[0][1]
c_opt = res[0][2]

x = np.arange(min(data[:,0]),max(data[:,0]),0.001)
y = quad_model(x,a_opt,b_opt,c_opt)

y_min = min(y)
x_min = x[np.argmin(y)]

op_write.write(str(np.power(x_min,1/3)))

plt.figure(figsize = (10,7))

plt.axvline(x_min,color='red')
plt.axhline(y_min,color='red')

plt.scatter(data[:,0],data[:,1])
plt.plot(x,y)

plt.scatter(x_min,y_min, label=f'Å = {np.power(x_min,1/3)}')

plt.xlabel(r"Volume ($Å^3$)")
plt.ylabel(r"Energy per unit cell volume (eV/$Å^3$)")
plt.legend(loc = 'upper left')
plt.savefig('../op_fig/pot_e_lp_dist.png')


op_write.close()
