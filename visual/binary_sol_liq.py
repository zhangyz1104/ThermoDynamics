import numpy as np
import matplotlib.pyplot as plt

filename = '../data/binary_sol_liq'
[c0, Tsol, Tliq] = np.loadtxt(filename).transpose()

''' plot figure '''
plt.plot(c0, Tsol, linewidth=3, c='k', linestyle='-')
plt.plot(c0, Tliq, linewidth=3, c='k', linestyle='-')
plt.plot([0, 1], [2050, 2050], linewidth=3, c='gray', linestyle='--')
plt.plot([0, 1], [1350, 1350], linewidth=3, c='gray', linestyle='--')
# plt.plot(Tm_guess, radius, linewidth=2, c='k', linestyle='--', label='Tm(guess)')
plt.xlim(0, 1)
plt.xlabel('Olivine concentration', fontsize=14)
plt.ylabel('Temperature(k)', fontsize=14)
plt.show()