import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/mcmc_energy_t2.txt","r").read().split('\n')
data.pop()
data = [float(item) for item in data]
data = np.array(data)

#visually inferred eqb
i = 3250

#average energy post eq
e_avg = np.mean(data[i:])

x = np.arange(0,len(data))

plt.figure(figsize=(20,10))
#plt.plot(x,data)

plt.rc('xtick',labelsize=23)
plt.rc('ytick',labelsize=23)

plt.scatter(x,data, s=5,alpha=1)
plt.axvline(i,color='black',label = r'N_eq ~ ' + f'{i}')
plt.axhline(e_avg,color='red',label=r'E_avg = ' f"{e_avg:.3f}" + r' E_h')

plt.xlabel("Sample number",fontsize='27')
plt.ylabel(r"Energy ($E_h$)",fontsize='27') #hartree energy
plt.title("MCMC Equilibration (Energy)",fontsize='30')

plt.legend(loc='upper right',fontsize='25')

#plt.show()
plt.savefig("../op_fig/mcmc_eqb_t2.png")
