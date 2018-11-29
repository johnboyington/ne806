import numpy as np
import matplotlib.pyplot as plt


def diffusion(n, T, D, Sig_a, S_0, source, geometry='slab'):
    """A finite difference diffusion solver."""

    # calculate L
    L = np.sqrt(D / Sig_a)

    # calculate x points and the step between them
    x = np.linspace(0, T, n)
    delta = x[1] - x[0]

    # set up A matrix
    A = np.zeros((len(x), len(x)))

    # set up source vector
    print(source)
    s = source(x, S_0) / D

    # loop through x's and update A matrix
    for i in range(len(x)):

        # first boundary condition
        if i == 0:
            A[i, i] = (1/4) - (D / (2 * delta))
            A[i, i + 1] = (D / (2 * delta))
            s[i] = 0

        # second boundary condition
        elif i == len(x) - 1:
            A[i, i] = (1/4) - (D / (2 * delta))
            A[i, i - 1] = (D / (2 * delta))
            s[i] = 0

        # for every other point
        else:
            A[i, i - 1] = (-1) / delta**2
            A[i, i] = (2 / delta**2) + (1 / L**2)
            A[i, i + 1] = (-1) / delta**2

    # solve the system of equations
    phi = np.linalg.solve(A, s)
    return phi


def source_a(x, S_0):
    """The distributed source for part a."""
    return S_0 * x**2


def source_b(x, S_0):
    """The distributed source for parts b and c."""
    return S_0


if __name__ == '__main__':
    phi = diffusion(100, 45, 0.600, 0.005, 10**8, source_a)
    plt.figure(0)
    plt.plot(phi)
