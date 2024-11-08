import matplotlib.pyplot as plt
import numpy as np

data = np.array(list(map(float,open("rng_dist.txt","r").read().split('\n'))))

x = np.arange(len(data))

nbins = 50
bin_width = (max(data) - min(data))/nbins
bins = [bin_width * i for i in range(nbins)]

plt.hist(data,bins,edgecolor='black')
plt.xlabel('Value')
plt.ylabel('Tally')
plt.title('Uniform distribution')
plt.show()