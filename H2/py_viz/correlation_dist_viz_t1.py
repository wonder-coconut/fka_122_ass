import numpy as np
import matplotlib.pyplot as plt

def uncorrelated_prob(theta):
    return 0.5 * np.sin(theta)

#distribution for both theta and x
data_t = open("../op_text/theta_t1.txt","r").read().split('\n')
data_x = open("../op_text/x_uniform_t1.txt","r").read().split('\n')
data_t.pop()
data_x.pop()

data_t = [float(item) for item in data_t]
data_x = [float(item) for item in data_x]
data_t = np.array(data_t)
data_x = np.array(data_x)

#uncorrelated distributions
theta = np.arange(0,np.pi,0.1)
x = np.arange(-1,1,0.1)

p_theta = 0.5 * np.sin(theta)
p_x = 1/2 * np.ones(len(x))

fig,axes = plt.subplots(2,1,figsize = (8,10))

axes[0].hist(data_t,bins = 100, density=True,label='Sample distribution')
axes[0].plot(theta,p_theta, label='Uncorrelated distribution')
axes[0].set_xlabel(r"$\theta$")
axes[0].set_ylabel(r"P($\theta$)")
axes[0].legend(loc='upper right')

axes[1].hist(data_x,bins = 100, density=True,label='Sample distribution')
axes[1].plot(x,p_x, label='Uncorrelated distribution')
axes[1].set_xlabel("x")
axes[1].set_ylabel("P(x)")
axes[1].legend(loc='upper right')

axes[0].set_title("Electron correlation")

plt.savefig('../op_fig/correlation_t1.png')