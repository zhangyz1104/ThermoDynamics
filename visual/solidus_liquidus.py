import numpy as np
import matplotlib.pyplot as plt

filename = '../data/solidus_liquidus'
[radius, Tsol, Tliq] = np.loadtxt(filename).transpose()
P = 2 * 3.14 * 6.67e-11 / 3 * 3400 ** 2 * 1e-3 * (1740 ** 2 - radius ** 2)
Tm_ol = 2050 + 60 * P
Tm_px = 1500 + 100 * P
Tm_fs = 1350 + 120 * P
Tm_guess = 0.5 * Tm_ol + 0.35 * Tm_px + 0.15 * Tm_fs

''' plot figure '''
plt.plot(Tsol, radius, linewidth=3, c='k', linestyle='-')
plt.plot(Tliq, radius, linewidth=3, c='k', linestyle='-')
plt.plot(Tm_ol, radius, linewidth=2, c='g', linestyle='--', label='Tm(olv)')
plt.plot(Tm_px, radius, linewidth=2, c='b', linestyle='--', label='Tm(pxn)')
plt.plot(Tm_fs, radius, linewidth=2, c='gray', linestyle='--', label='Tm(fsp)')
# plt.plot(Tm_guess, radius, linewidth=2, c='k', linestyle='--', label='Tm(guess)')
plt.legend()
plt.ylim(radius[0], radius[-1])
plt.show()