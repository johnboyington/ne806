import numpy as np


def power_iteration(A, num_simulations):
    # Ideally choose a random vector
    # To decrease the chance that our vector
    # Is orthogonal to the eigenvector
    b_k = np.random.rand(A.shape[1])
    num_simulations = int(num_simulations)

    for _ in range(num_simulations):
        # calculate the matrix-by-vector product Ab
        b_k1 = np.dot(A, b_k)

        # calculate the norm
        b_k1_norm = np.max(b_k1)
        c = 2 / b_k1_norm

        # re normalize the vector
        b_k = b_k1 / b_k1_norm

    return b_k, c


B = np.loadtxt('B.txt', delimiter='&')
phi = np.ones(len(B))
c = 1
F = 2

phi_k, c_crit = power_iteration(B, 1e3)
print('c_crit:  {}'.format(c_crit))
