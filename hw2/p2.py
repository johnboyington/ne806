import numpy as np
from numpy.random import rand


def phi(x, mu):
    return 4e4 * np.exp(-(x/3)) * (2 + np.sin((x * mu)/3))


N = int(1e6)

x_0 = 5
s = 0
s2 = 0

for i in range(N):
    rho = rand() * 2 - 1
    v = phi(x_0, rho)
    s += v
    s2 += v**2


I = (s / N)
err = np.sqrt(((s2 / N) - I**2) * (N / (N - 1)))
err = err / np.sqrt(N)
total_flux_density = I * 2

print('(a) Total Flux Density:  {:4.2e}'.format(total_flux_density))
print('(b) Error:  {:4.2e}'.format(err))
