import numpy as np
import matplotlib.pyplot as plt

filename = '../data/fractional_crystal'
[T, Tliq, Tsol, cs1, cl1, cs2, cl2, cs3, cl3, rc] = np.loadtxt(filename).transpose()

zero = np.zeros(len(T))
one = np.ones(len(T))

# plt.plot(cs1, rc, linewidth=3, c='g', linestyle='-', label='olv')
# plt.plot(cs2, rc, linewidth=3, c='b', linestyle='-', label='pxn')
# plt.plot(cs3, rc, linewidth=3, c='gray', linestyle='-', label='fsp')
# plt.ylim(0, 1)
# plt.ylabel('Radius')
# plt.legend()
# plt.show()

plt.figure(figsize=(1.5, 5))
plt.fill_betweenx(rc, zero, cs1, color='g', alpha=1)
plt.fill_betweenx(rc, cs1, cs1 + cs2, color='b', alpha=1)
plt.fill_betweenx(rc, cs1 + cs2, one, color='gray', alpha=1)
# plt.axis('off')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xticks([])
plt.tight_layout()
plt.savefig('../figure/fractional_crystal_colume.png')
plt.show()