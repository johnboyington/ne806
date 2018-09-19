import numpy as np
import matplotlib.pyplot as plt


def w_c(alpha, ratio):
    return (1/(1-alpha)) * (2 * ratio - (1 + alpha))


def P(alpha):
    x = np.linspace(alpha, 0.999, 100)
    y = 0.5 * (1 + w_c(alpha, x))
    return x, y


fig = plt.figure(0)
ax = fig.add_subplot(111)

alphas = [0, 0.2, 0.5, 0.9]
for a in alphas:
    ax.plot(*P(a))
