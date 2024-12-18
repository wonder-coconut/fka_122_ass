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

bin_counts = open('anal_files/bins.txt').read().split('\n')
bin_counts.pop()
bin_counts = [float(item) for item in bin_counts]

time = np.loadtxt("anal_files/time.txt")

n_ideal = np.array(n_ideal)
bin_counts = np.array(bin_counts)
r_k = np.array(r_k)

#print(data)

max_idx = np.argmax(bin_counts/n_ideal)
min_idx = max_idx + np.argmin(bin_counts[max_idx:]/n_ideal[max_idx:])+4

#calculate integral
g_r = bin_counts/n_ideal/(time)

# Compute integral using the trapezoidal rule
dr = np.diff(r_k)  # Spacing between r values

g_r_integral = np.sum(g_r[:min_idx]*3/3.5 * r_k[:min_idx] * r_k[:min_idx] * dr[:min_idx]) * 4 * np.pi * rho

print(f'R_min = {r_k[min_idx]}')
print(f'Coordination number = {g_r_integral}')