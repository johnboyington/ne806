import numpy as np
import matplotlib.pyplot as plt


def f(w, psi):
    return (1 / (2*np.pi)) * (1 + w)



fig = plt.figure(0)
ax = fig.add_subplot(111)
ax.set_xlabel('$\psi$')
ax.set_ylabel('$f(\psi)$')

w = np.zeros(100)
psi = np.linspace(0, 2*np.pi, 100)
ax.plot(psi, f(w, psi), color='k')


fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.set_xlabel('$\omega_c$')
ax.set_ylabel('$f(\omega_c)$')

w = np.linspace(-1, 1, 100)
ax.plot(w, f(w, psi), color='k')