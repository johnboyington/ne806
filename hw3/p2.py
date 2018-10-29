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


def do_problem_two(params):
    """Does problem 2."""
    phi_i = initialize_flux(params.N)
    B = create_B(params.N, params.T)
    phi, c_crit = power_iteration(phi_i, B, params.tol)
    return phi, c_crit, B


def loop_problem_two():
    """Loops problem two."""
    c_crits = []
    t_vals = np.linspace(0.01, 5, 5)
    for t in t_vals:
        params = Problem_Parameters(64, t, 0.00000001)
        phi, c, B = do_problem_two(params)
        c_crits.append(c)

    # plotting c_crit
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    ax.set_xlabel('Thickness (cm)')
    ax.set_ylabel('$\ln (c_{crit} - 1)$')
    c_crits = np.array(c_crits)
    ax.plot(t_vals, np.log(c_crits - 1), 'k')
    fig.savefig('c_crit.png', dpi=250)

    # plotting phi
    T = 2
    x = np.linspace(-T, T, 64)
    params = Problem_Parameters(64, 2, 0.00000001)
    phi, c, B = do_problem_two(params)
    np.save('p2.npy', phi)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.set_xlabel('$x$ (cm)')
    ax.set_ylabel('$\Phi$')
    ax.plot(x, phi, 'k')
    fig.savefig('phi.png', dpi=250)

    # plotting phi
    T = 0.5
    x = np.linspace(-T, T, 64)
    params = Problem_Parameters(64, T, 0.00000001)
    phi, c, B = do_problem_two(params)
    np.save('p2_05.npy', phi)

    # plotting phi
    T = 10
    x = np.linspace(-T, T, 64)
    params = Problem_Parameters(64, T, 0.00000001)
    phi, c, B = do_problem_two(params)
    np.save('p2_10.npy', phi)
    return

if __name__ == '__main__':
    loop_problem_two()
