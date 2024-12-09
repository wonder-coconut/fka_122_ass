import numpy as np
import matplotlib as pyplot

data = open("anal_files/t4_lattice.xyz").read().split('\n')
file_out = open("anal_files/positions_liquid.txt","w")

data.pop(0)
data.pop(0)
data.pop()

for i in range(len(data)):
    data[i] = data[i].split(' ')
    data[i].pop(0)
    data[i].pop()    
    data[i].pop()
    data[i].pop()
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

for pos in data:
    for i in pos:
        file_out.write(str(i)+'\t')
    file_out.write('\n')