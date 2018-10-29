import numpy as np
import matplotlib.pyplot as plt

T = 26.872
x = np.linspace(-T, T, 64)
phi1 = np.load('p1.npy')
phi2 = np.load('p2.npy')

phi1 = phi1 / np.sum(phi1)
phi2 = phi2 / np.sum(phi2)

fig = plt.figure(0)
ax = fig.add_subplot(111)
ax.set_xlabel('x (cm)')
ax.set_ylabel('Flux Profile')
ax.plot(x, phi1, 'k', label='Monte Carlo')
ax.plot(x, phi2, 'k--', label='Power Iteration')

ax.legend()
fig.savefig('comparison', dpi=250)
