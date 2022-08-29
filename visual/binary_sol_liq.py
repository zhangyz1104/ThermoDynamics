import numpy as np
import matplotlib.pyplot as plt

filename = '../data/binary_sol_liq'
[c0, Tsol, Tliq] = np.loadtxt(filename).transpose()

''' plot figure '''
plt.figure(figsize=(5, 4))
plt.plot(c0, Tsol, linewidth=3, c='k', linestyle='-')
plt.plot(c0, Tliq, linewidth=3, c='k', linestyle='-')
plt.plot([0, 1], [2050, 2050], linewidth=3, c='r', linestyle='--', label='Tm(olv)')
plt.plot([0, 1], [1350, 1350], linewidth=3, c='gray', linestyle='--', label='Tm(fsp)')
# plt.plot(Tm_guess, radius, linewidth=2, c='k', linestyle='--', label='Tm(guess)')
plt.xlim(0, 1)
plt.xlabel('Olivine concentration', fontsize=14)
plt.ylabel('Temperature(k)', fontsize=14)
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('../figure/binary_sol_liq.png')
plt.show()