import matplotlib.pyplot as plt
import numpy as np

n_atoms = 256
lp = 4.236854996336203
L = 4*lp
V = L*L*L
rho =  n_atoms/V

n_ideal = open('anal_files/n_ideal.txt').read().split('\n')
n_ideal.pop()
n_ideal = [float(item) for item in n_ideal]


r_k = open('anal_files/r_k.txt').read().split('\n')
r_k.pop()
r_k = [float(item) for item in r_k]

hist = open('anal_files/bins.txt').read().split(' ')
hist = [float(item) for item in hist]

n_ideal = np.array(n_ideal)
hist = np.array(hist)
r_k = np.array(r_k)

#print(data)

g_r = hist/(1000*n_ideal)
min_idx = np.argmin(g_r[60:]) + 60


#calculate integral

valid_indices = r_k <= r_k[min_idx]

# Compute integral using the trapezoidal rule
dr = np.diff(r_k)  # Spacing between r values

g_r_integral = np.sum(g_r[:-1][valid_indices[:-1]] * r_k[:-1][valid_indices[:-1]]**2 * dr[valid_indices[:-1]])

coordination_number = 4 * np.pi * rho * g_r_integral

print(r_k[min_idx])
print(coordination_number)