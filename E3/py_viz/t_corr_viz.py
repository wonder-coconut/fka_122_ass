import numpy as np
import matplotlib.pyplot as plt

data = open('../op_text/t_corr_low.txt',"r").read().split('\n')
data.pop()
data = [float(item) for item in data]

data = np.array(data)
plt.plot(data)
plt.show()