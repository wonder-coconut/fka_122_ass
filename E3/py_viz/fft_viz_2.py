import numpy as np
import matplotlib.pyplot as plt

dt = 0.001
param_set = 'high'

filename = f'../op_text/traj_out_{param_set}_{dt}.txt'
data = open(filename,"r").read().split('\n')
pressure = float(data.pop(0))
n = int(data.pop(0))
data.pop()

for i in range(len(data)):
    data[i] = data[i].split('\t')
    data[i] = [float(item) for item in data[i]]

data = np.array(data)

t = data[:,0]
x = data[:,1]
v = data[:,2]

#split fft
dt_fft = 50*dt
M = 30
signal_length = len(v)
segment_length = int(signal_length/M)

p_sum = np.zeros(segment_length)
f = np.fft.fftfreq(segment_length,dt_fft) * 2 * np.pi

i = 0
while(i < signal_length):
    p = np.abs(np.fft.fft(v[i:i+segment_length]))**2
    if(p_sum.shape == p.shape):
        p_sum += p
    i += int(segment_length)

plt.xlim(0,10)
plt.plot(f,p_sum/M * (dt_fft/segment_length))
plt.xlabel('frequency')
plt.ylabel('power')
plt.title(f'averaged fft for {pressure}KPa at {dt_fft}ms sampling\n{M} segments')
plt.savefig(f'../op_fig/fft_avg_{param_set}_{dt_fft}.png')