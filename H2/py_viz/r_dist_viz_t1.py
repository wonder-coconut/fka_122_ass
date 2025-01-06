import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

def prob(Z,r):
    return Z*Z*Z * 4*r*r * np.exp(-2*Z*r)

data = open("../op_text/samples.txt","r").read().split('\n')
data.pop()
n = int(data.pop(0))

for i in range(n):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

r1_samples = data[:,0:3]
r2_samples = data[:,3:6]

#distance from origin
r1 = np.sqrt(np.sum(r1_samples**2, axis=1))
r2 = np.sqrt(np.sum(r2_samples**2, axis=1))

#radial distance attribute for defined pdf
x = np.arange(0,6,0.1)

p1 = prob(2,x) #z = 2
p2 = prob(27/16,x) #z = 27/16

#plotting the two z cases

plt.rc('xtick',labelsize=18)
plt.rc('ytick',labelsize=18)

plt.figure(figsize=(8,6))

plt.plot(x,p1, label = 'Z = 2')
plt.plot(x,p2,color='red', label = 'Z = 27/16')

#sample distribution as a histogram
r_res = np.concatenate((r1,r2))
plt.hist(r_res,bins=100,density=True, color='orange', label = 'Sampling distribution')
plt.legend(loc = 'upper right',fontsize='18')
plt.xlabel(r'Radial distance r ($a_0$)',fontsize='18') #bohr radius
plt.ylabel('Radial distribution probability P(r)',fontsize='18')
plt.title(f'Comparison of Radial Distribution Probabilites\n(MCMC Samples={n})',fontsize='20')

plt.savefig('../op_fig/radial_dist_prob_t1.png')