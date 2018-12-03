import numpy as np
import matplotlib.pyplot as plt


def diffusion(n, T, D, Sig_a, S_0, source, geometry='slab'):
    """A finite difference diffusion solver."""

    # guarantee geometry in correct format
    assert geometry in ['slab', 'cylinder', 'sphere'], "Geometry must either be 'slab', 'cylinder' or 'sphere'."

    # calculate L
    L = np.sqrt(D / Sig_a)

    # calculate x points and the step between them
    x = np.linspace(0, T, n)
    delta = x[1] - x[0]

    # set up A matrix
    A = np.zeros((len(x), len(x)))

    # set up source vector
    s = source(x, S_0) / D

    # loop through x's and update A matrix
    for i in range(len(x)):

        # first boundary condition
        if i == 0:
            if geometry == 'slab':
                A[i, i] = (1/4) - (D / (2 * delta))
                A[i, i + 1] = (D / (2 * delta))
                s[i] = 0

            elif geometry in ['sphere', 'cylinder']:
                A[i, i] = (-1) / delta
                A[i, i + 1] = 1 / delta
                s[i] = 0

        # second boundary condition
        elif i == len(x) - 1:
            if geometry == 'slab':
                A[i, i] = (1/4) - (D / (2 * delta))
                A[i, i - 1] = (D / (2 * delta))
                s[i] = 0

            elif geometry in ['sphere', 'cylinder']:
                A[i, i - 1] = (1/4) - (D / (2 * delta))
                A[i, i] = (D / (2 * delta))
                s[i] = 0

        # for every other point
        else:
            if geometry == 'slab':
                A[i, i - 1] = (-1) / delta**2
                A[i, i] = (2 / delta**2) + (1 / L**2)
                A[i, i + 1] = (-1) / delta**2

            elif geometry == 'sphere':
                A[i, i - 1] = (-1) / delta**2
                A[i, i] = (x[i+1]**2/(x[i]**2 * delta**2) + 1 / delta**2 + 1 / L**2)
                A[i, i + 1] = (-x[i+1]**2) / (x[i]**2 * delta**2)

            elif geometry == 'cylinder':
                A[i, i - 1] = (-1) / delta**2
                A[i, i] = (x[i+1]/(x[i] * delta**2) + 1 / delta**2 + 1 / L**2)
                A[i, i + 1] = (-x[i+1]) / (x[i] * delta**2)

    # solve the system of equations
    return np.linalg.solve(A, s)


def source_a(x, S_0):
    """The distributed source for part a."""
    return S_0 * x**2


def source_b(x, S_0):
    """The distributed source for parts b and c."""
    return S_0 * np.ones(len(x))


def plot_it(phi, label, savename):
    """Utility for plotting."""
    x = np.linspace(0, 45, 100)
    plt.figure(0)
    plt.xlabel('$\hat{r}$')
    plt.ylabel('$\Phi$')
    plt.plot(x, phi, 'k', label=label)
    plt.legend()
    plt.savefig(savename + '.png', dpi=300)
    plt.clf()

if __name__ == '__main__':
    phi_slab = diffusion(100, 45, 0.600, 0.005, 10**8, source_a, geometry='slab')
    phi_cyl = diffusion(100, 45, 0.600, 0.005, 10**8, source_b, geometry='cylinder')
    phi_sph = diffusion(100, 45, 0.600, 0.005, 10**8, source_b, geometry='sphere')

    # plot the slab
    plot_it(phi_slab, 'Slab', 'slab')
    plot_it(phi_cyl, 'Cylinder', 'cyl')
    plot_it(phi_sph, 'Sphere', 'sph')
