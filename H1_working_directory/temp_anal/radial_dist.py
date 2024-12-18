import numpy as np
import matplotlib.pyplot as plt

lp = 4.236854996336203
n_atoms = 256
N = 4
L = N*lp
V = L*L*L
#file read
data = open("anal_files/pos_evol_md copy.txt").read().split('\n')
data.pop()
for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i].pop()
    data[i] = [[float(data[i][j*3]),float(data[i][j*3+1]),float(data[i][j*3+2])] for j in range(n_atoms)]

data = np.array(data)
#bin width
del_r = 0.032
#k => max number of bins
k_max = 250
k = np.arange(1,k_max+1)
#radial bins
r_k = (k - 0.5) * del_r
r_max = np.max(r_k)
#histogram
bin_counts = np.zeros(k_max)

#ideal density
n_ideal = (n_atoms-1)/V * 4*np.pi/3 * (3*k*k - 3*k + 1)*del_r*del_r*del_r

time = len(data)
#histogram calc

#time iteration
for i in range(time):
    print(i)
    #particle 1 iteration
    for j in range(n_atoms):
        #particle 2 iteration
        for k in range(j + 1, n_atoms):
            distance = np.linalg.norm(data[i][j] - data[i][k])
            #print(distance)
            if(distance < r_max):
                #sort in histogram
                bin_index = int(distance/r_max * k_max)
                bin_counts[bin_index] += 1

np.savetxt("anal_files/bins.txt",bin_counts)
np.savetxt("anal_files/r_k.txt",r_k)
np.savetxt("anal_files/n_ideal.txt",n_ideal)
np.savetxt("anal_files/time.txt",[time])

#plt.plot(r_k,bin_counts/(time*n_ideal))
#plt.savefig('radial_dist.png')