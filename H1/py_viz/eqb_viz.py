import numpy as np
import matplotlib.pyplot as plt

data = open("../op_text/equilibration_t_p.txt","r").read().split('\n')

data.pop()

t_eq = 0
p_eq = 0

for i in range(len(data)):
    data[i] = data[i].split('\t')
    if(i == 0):
        t_eq = float(data[i][1])
        p_eq = float(data[i][3])
    
    j = 0
    res = []
    while(j < len(data[i])):
        res.append(float(data[i][j]))
        j += 2
    data[i] = res

data = np.array(data)
p_eq/=1000000
data[:,1]/=1000000

x = np.arange(len(data))

fig,axes = plt.subplots(2,1,figsize=(10,10))
axes[0].plot(x,data[:,1])
axes[0].axhline(p_eq, color = 'red')

axes[1].plot(x,data[:,0])
axes[1].axhline(t_eq, color = 'red')
plt.savefig('../op_fig/eqb.png')