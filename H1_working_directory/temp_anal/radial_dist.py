import numpy as np
import matplotlib.pyplot as plt

#sim constants
n_atoms = 256
lp = 4.236854996336203

#read and parse data
data = open("anal_files/pos_evol_t4_temp.txt").read().split('\n')
f_out = open("anal_files/bins.txt")
data.pop()
for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i].pop()
    data[i] = [float(item) for item in data[i]]

data = np.array(data)
positions = []

for i in range(len(data)):
    temp = []
    j = 0
    while(j < 256*3):
        temp.append([data[i][j],data[i][j+1],data[i][j+2]])
        j += 3
    positions.append(temp)

positions = np.array(positions)

#simulation box characteristics setup
L = 4 * lp
V = L * L * L

del_r = 0.1 #bin width

k_max = (L)/del_r #maximum k => max distance
k = np.arange(1,k_max)

r_k = (k-0.5)*del_r #radius values => pairwise distance

num_bins = len(k)
bin_counts = np.zeros(num_bins) #bins for pairwise distance collection


max_distance = L/2

#ideal calc
n_ideal = (255)/V * 4*np.pi/3 *(3*k*2 - 3*k + 1) * del_r**3 

np.savetxt('anal_files/n_ideal.txt',n_ideal)
np.savetxt('anal_files/r_k.txt',r_k)

#distance_calc
for t in range(len(data)):
    print(t)
    for i in range(256):
        for j in range(i + 1, 256):
            distance = np.linalg.norm(positions[t][i] - positions[t][j])
            if(distance < max_distance):
                bin_index = int(distance/max_distance*num_bins)
                bin_counts[bin_index] += 2

f_out.write("\n")
f_out.write(str(bin_counts))

#bin_counts = np.loadtxt('anal_files/bins.txt')

bin_counts = np.array(bin_counts)
bin_counts = bin_counts/999

plt.title("Radial Distribution Function for T=973 K, P=1 bar")
plt.xlabel("r (Å)")
plt.ylabel("g(r)")
plt.plot(r_k,bin_counts/n_ideal)
plt.savefig("radial_dist.png")