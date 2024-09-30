# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

ffttypes = ['npfft', 'scfft', 'pyfftw']

# collect data
data = []
nn = []
energies = []
ttimes = []
for fft in ffttypes:
    data.append(np.loadtxt('ch4_%s.dat' % fft, usecols=[0,2,5]))
    nn.append(data[-1][:,0])
    energies.append(data[-1][:,1])
    ttimes.append(data[-1][:,2])

fig, ax = plt.subplots(2, 1, dpi=144, figsize=(8,6))

# plot energy scaling
ax[0].semilogx(nn[0]**3, energies[0], '--o')
ax[0].set_xlabel('Number of grid points per cartesian direction [-]')
ax[0].set_ylabel('Energy [Ht]')
ax[0].grid(linestyle='--')
ax[0].set_title('Total electronic energy')

labels = ['NumPy FFT',
          'SciPy FFT',
          'pyFFTW']
colors = ['tab:blue',
          'tab:orange',
          'tab:green']

# plot timings
for i,l in enumerate(labels):
    # plot result
    ax[1].bar(np.arange(3,8)+(i-1)*0.3, ttimes[i][3:8], zorder=2, 
              width=0.3, label=l, alpha=0.5, color=colors[i])
    
ax[1].set_xticks(np.arange(3,8), [r'%i' % i for i in nn[0][3:8]])
ax[1].set_xlabel('Number of grid points per cartesian direction [-]')
ax[1].set_ylabel('Computation time [s]')
ax[1].legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.5))
ax[1].set_title('Comparison of FFT methods')
ax[1].grid(linestyle='--', alpha=0.5, color='black')

plt.tight_layout()
plt.show()