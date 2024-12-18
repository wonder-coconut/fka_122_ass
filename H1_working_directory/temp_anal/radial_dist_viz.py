import numpy as np
import matplotlib.pyplot as plt

bin_counts = np.loadtxt("anal_files/bins.txt")
r_k = np.loadtxt("anal_files/r_k.txt")
n_ideal = np.loadtxt("anal_files/n_ideal.txt")
time = np.loadtxt("anal_files/time.txt")

max_idx = np.argmax(bin_counts/n_ideal)
min_r_k = r_k[max_idx + np.argmin(bin_counts[max_idx:]/n_ideal[max_idx:])+4]


#plt.plot(r_k,bin_counts/n_ideal/(time))
#plt.axvline(min_r_k,color='r', label=f'First minimum at r = {min_r_k:.3f}')
#plt.axhline(1,linestyle='--', alpha = 0.5, color='black', label='Function convergence')
#plt.title("Radial Distribution Function for T=973 K, P=1 bar")
#plt.xlabel("r (Å)")
#plt.ylabel("g(r)")
#plt.legend(loc = 'upper right')
#plt.savefig('radial_dist.png')
#plt.show()