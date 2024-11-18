import numpy as np
import matplotlib.pyplot as plt

def h(a,f,t,psi):
    return a * np.cos(2*np.pi * f * t + psi)

a = 1
psi = 0
f = 2
dt = 0.1
N = 250

t_k = np.arange(0,N) * dt
h_t_1 = h(a,f,t_k,psi)
h_t_2 = h(a,1,t_k,psi)
h_t_3 = h(a,f,t_k,np.pi/2)

fig,axes = plt.subplots(3,1,figsize=(10,10))
axes[0].plot(t_k,h_t_1)
axes[0].set_title('a = 1, f = 2, psi = 0')
axes[0].set_ylabel('signal (a.u.)')
axes[1].plot(t_k,h_t_2)
axes[1].set_title('a = 1, f = 1, psi = 0')
axes[1].set_ylabel('signal (a.u.)')
axes[2].plot(t_k,h_t_3)
axes[2].set_title('a = 1, f = 2, psi = pi/2')
axes[2].set_ylabel('signal (a.u.)')
plt.xlabel('time (a.u.)')
plt.savefig('signal_viz.png')