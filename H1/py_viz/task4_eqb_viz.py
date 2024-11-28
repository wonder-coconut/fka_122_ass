import numpy as np
import matplotlib.pyplot as plt

data = open("../op_text/equilibration_t_p_t4.txt","r").read().split('\n')

params = data.pop(0).split('\t')

timestep = float(params[0])
t_eq = float(params[1])
p_eq = float(params[2])

data.pop()


for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)
gpa_in_asu = 6.2415076486555486E-03
data[:,1] /= gpa_in_asu
p_eq /= gpa_in_asu

x = np.arange(len(data)) * timestep

fig,axes = plt.subplots(2,1,figsize=(10,10))

axes[0].plot(x,data[:,0], label = 'Temperature fluctuation')
axes[0].axhline(t_eq, color = 'red', label = f'Equillibrium Temperature T = {t_eq} K')
axes[0].legend(loc = 'upper right')
axes[0].set_xlabel('Time (ps)')
axes[0].set_ylabel('Temperature (K)')

axes[1].plot(x,data[:,1],label = 'Pressure fluctuation')
axes[1].axhline(p_eq, color = 'red', label = f'Equilibrium Pressure P = {p_eq:1.4f} GPa')
axes[1].legend(loc = 'upper right')
axes[1].set_xlabel('Time (ps)')
axes[1].set_ylabel('Pressure (GPa)')

axes[0].set_title(f"Equilibriation for {t_eq} K and {p_eq:1.4f} GPa")

plt.savefig('../op_fig/eqb_t4.png')