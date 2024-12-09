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

#wave vector q
q_vector = np.zeros(3)
n_max = 3
L = 4*lp

i = j = k = -n_max
q_grid = []

while(i <= n_max):
    j = -n_max
    while(j <= n_max):
        k = -n_max
        while(k <= n_max):
            q_grid.append([i,j,k])
            k += 1
        j += 1
    i += 1

#3d grid containing q vectors
q_grid = np.array(q_grid)
q_grid = q_grid * 2*np.pi/L


s_q = 0
#static structure computation (naive)
for q in q_grid:
    print(f"q = {q}")
    s_q = 0
    for position_t in positions:
        temp = 0
        for px in position_t:
            temp += np.exp(1j * np.dot(px,q))
        s_q += np.abs(temp)**2
    print(s_q/999)