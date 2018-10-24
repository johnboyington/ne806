import numpy as np
import matplotlib.pyplot as plt
import mpmath
from gauss_factors import gaussian_integration_factors


class Problem_Parameters(object):
    """A container for all of the NTE problem parameters."""
    def __init__(self, N, T, tol):
        self.N = N
        self.T = T
        self.tol = tol


def initialize_flux(N):
    """A function that produces an initial flux guess."""
    return np.ones(N)


def E(n, z):
    return mpmath.expint(n, z)


def y(x, T):
    return (T/2) * x + (T/2)


def create_B(N, T):
    """Generates the B matrix for power iteration."""
    # initialize B and get guassian integration factors
    B = np.zeros((N, N))
    X, W = gaussian_integration_factors(N)

    # loop through B and calculate elements
    for i in range(N):
        for j in range(N):
            B[i, j] = W[j] * E(1, abs(y(X[j], T) - y(X[i], T)))

    # set diagonal back to zero
    np.fill_diagonal(B, 0)

    # calculate diagonal
    for i in range(N):
        B[i, i] = 2 - E(2, y(X[i], T)) - E(2, T - y(X[i], T)) - np.sum(B[i])
    return B


def power_iteration(phi, B, tol):
    """Function that returns c_crit using the power iteration method."""
    diff = 1E5
    lam_1 = lam_2 = 0

    # begin iterative convergence
    while diff > tol:
        lam_1 = lam_2
        x = phi
        phi = np.matmul(B, phi)
        lam_2 = np.linalg.norm(phi)/np.linalg.norm(x)
        diff = lam_2-lam_1
        c = 2 / lam_2

    # return c_crit
    return phi, c


def do_problem_two():
    """Does problem 2."""
    params = Problem_Parameters(64, 2, 0.00000001)
    phi_i = initialize_flux(params.N)
    B = create_B(params.N, params.T)
    phi, c_crit = power_iteration(phi_i, B, params.tol)
    return phi, c_crit, B


if __name__ == '__main__':
    phi, c, B = do_problem_two()
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    ax.plot(phi)
    print(c)
