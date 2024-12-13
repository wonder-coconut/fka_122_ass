import numpy as np

data = open('../op_text/mcmc_energy_eqb_t2.txt',"r").read().split('\n')

n = float(data.pop(0))
data.pop
data = [float(item) for item in data]
data = np.array(data)
print(np.std(data))