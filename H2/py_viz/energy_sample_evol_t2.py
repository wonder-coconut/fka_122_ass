import matplotlib.pyplot as plt
import numpy as np

data = open("../op_text/mcmc_energy_t2.txt","r").read().split('\n')
data.pop()
data = [float(item) for item in data]
data = np.array(data)

#finding n_eq
for i in range(len(data) - 1):
    if (data[i+1] - data[i] > 1E-2):
        break

#average energy post eq
e_avg = np.mean(data[i:])

x = np.arange(0,len(data))

plt.figure(figsize=(20,10))
plt.scatter(x,data, s=1,alpha=0.05)
plt.axvline(i,color='red',label = r'N_eq = ' + f'{i}')
plt.axhline(e_avg,color='red',label=r'E_avg = ' f"{e_avg}")

plt.xlabel("Sample number")
plt.ylabel("Energy (a.u.)")
plt.title("MCMC Equilibration (Energy)")

plt.legend(loc='upper right')

#plt.show()
plt.savefig("../op_fig/mcmc_eqb_t2.png")

#print eq samples to file
file = open("../op_text/mcmc_energy_eqb_t2.txt","w")

#length
file.write(str(len(data[i:])) + '\n')

for j in range(len(data[i:])):
    if(j != len(data[i:]) - 1):
        file.write(str(data[j+i])+"\n")
    else:
        file.write(str(data[j+i]))

