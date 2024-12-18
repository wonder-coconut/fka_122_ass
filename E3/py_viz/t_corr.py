import numpy as np
import matplotlib.pyplot as plt

def psi(t):
    tp = 147.3 * 1E-3
    w0 = 19.477874452256717
    w1 = np.sqrt(w0*w0 - 1/((2*tp)*(2*tp)))

    psi = np.exp(-t/(2*tp)) * (np.cos(w1*t) - np.sin(w1*t)/(2*w1*tp))
    return psi

t = np.arange(0,1,0.001)
res = psi(t)

plt.plot(t,res)
plt.show()