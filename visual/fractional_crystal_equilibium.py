import matplotlib.pyplot as plt
import numpy as np

''' read data '''
filename = '../data/fractional_crystal'
[T, Tliq, Tsol, cs1, cl1, cs2, cl2, cs3, cl3, rc] = np.loadtxt(filename)[::-1].transpose()


''' plot figure '''
figure, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(5, 6), sharex=True)
ax1.plot(T, 1 - rc, linewidth=3, c='k', linestyle='-')
ax1.set_ylim(0, 1)
ax1.set_ylabel('Melt fraction')

ax2.plot(T, cl1, linewidth=3, c='g', linestyle='-', label='olv')
ax2.plot(T, cl2, linewidth=3, c='b', linestyle='-', label='pxn')
ax2.plot(T, cl3, linewidth=3, c='gray', linestyle='-', label='fsp')
ax2.set_ylim(0, 1)
ax2.set_ylabel('Melt composition')
ax2.legend()

ax3.plot(T, cs1, linewidth=3, c='g', linestyle='-', label='olv')
ax3.plot(T, cs2, linewidth=3, c='b', linestyle='-', label='pxn')
ax3.plot(T, cs3, linewidth=3, c='gray', linestyle='-', label='fsp')
ax3.set_ylim(0, 1)
ax3.set_ylabel('Rock composition')
ax3.legend()

plt.xlim(T[0], T[-1])
# plt.xticks([1500, 1600, 1700, 1800], ['1500', '1600', '1700', '1800'])
plt.xlabel('Temperature(K)')
plt.tight_layout()
plt.savefig('../figure/fractional_crystal_equilibrium.png')
plt.show()
