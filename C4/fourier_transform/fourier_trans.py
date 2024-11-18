import numpy as np
import matplotlib.pyplot as plt

def h(a,f,t,psi):
    return a * np.cos(2*np.pi * f * t + psi)

def h2(a,f,t,psi):
    return a[0] * np.cos(2*np.pi * f[0] * t + psi[0]) + a[1] * np.cos(2*np.pi * f[1] * t + psi[1])

a = [1,1]
psi = [0,0]
f = [2,6]
dt = 0.08
N = 250

t_k = np.arange(0,N) * dt #time domain

h_t = h2(a,f,t_k,psi) #signal

H_n = np.fft.fft(h_t) #fourier transformed signal
H_n = np.fft.fftshift(H_n)

#f = np.arange(0,N)/(dt * N) #asymmetric freq domain
f2 = np.fft.fftfreq(N, dt)  #symmetric Frequency domain 
f2 = np.fft.fftshift(f2)

P_n = (np.abs(H_n) * np.abs(H_n))/N

fig,axes = plt.subplots(2,1,figsize = (10,10))
axes[0].set_title('Symmetric Power Spectrum (sampling = 0.075)')

axes[0].plot(t_k,h_t)
axes[0].set_ylabel('signal (a.u.)')
axes[0].set_xlabel('time (a.u.)')

axes[1].plot(f2,P_n)
axes[1].set_ylabel('power (a.u.)')
axes[1].set_xlabel('frequency (a.u.)')

#plt.show()
plt.savefig('sym_pow_2signals_altered_sampling.png')