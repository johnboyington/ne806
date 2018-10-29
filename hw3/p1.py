import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt


class Slab_Parameters(object):
    """Container for all problem parameters."""
    def __init__(self, T, Sig_c, Sig_s, Sig_f, nu_bar, n_t):
        self.T = T
        self.Sig_c = Sig_c
        self.Sig_s = Sig_s
        self.Sig_f = Sig_f
        self.Sig_t = Sig_f + Sig_c + Sig_s
        self.nu_bar = nu_bar
        self.n_t = n_t
        self.edges = np.linspace(-self.T, self.T, self.n_t + 1)
        return


class Counter(object):
    """Container for problem tally information."""
    def __init__(self, n_t):
        self.N = 0
        self.N_l = 0
        self.N_c = 0
        self.fission_spec = np.zeros(n_t)
        self.N_f = 0
        self.N_s = 0
        self.D = 0
        return


class Score(object):
    """Container for problem score information."""
    def __init__(self):
        self.N_cumulative = 0
        self.D_cumulative = 0
        self.S_l = 0
        self.S2_l = 0
        self.S_c = 0
        self.S2_c = 0
        self.S_f = 0
        self.S2_f = 0
        self.S_s = 0
        self.S2_s = 0
        self.S_t = 0
        self.S2_t = 0
        return


class Particle(object):
    """Container for an individual particle's information."""
    def __init__(self):
        self.position = None
        self.direction = None


class Source_PDF(object):
    """Container for source distribution."""
    def __init__(self, n_t):
        self.pdf = np.ones(n_t) / n_t
        self.cdf = np.cumsum(self.pdf)

    def update_dist(self, dist):
        self.pdf = dist / np.sum(dist)
        self.cdf = np.cumsum(self.pdf)
        return


def leaked(position, params):
    """Returns whether or not a particle has leaked."""
    return False if position is None else (abs(position) > params.T)


def sample_source_position(q, params):
    """Chooses a source position from the source pdf."""
    i = np.searchsorted(q.cdf, rand())
    l, r = params.edges[i:i+2]
    return ((r - l) * rand()) + l


def choose_direction():
    """Chooses a direction for the particle."""
    return 1 - 2 * rand()


def choose_path_length(params):
    """Determines the path length for the particle to travel."""
    return -(1 / params.Sig_t) * np.log(rand())


def determine_reaction_type_and_update_counter(particle, params, counter):
    """Chooses a reaction type for the particle."""
    rho = rand()
    if rho <= params.Sig_f / params.Sig_t:
        i = np.searchsorted(params.edges, particle.position)
        counter.fission_spec[i - 1] += 1
        counter.N_f = np.sum(counter.fission_spec)
        return 'fission'
    elif rho < (params.Sig_f + params.Sig_c) / params.Sig_t:
        counter.N_c += 1
        return 'capture'
    else:
        counter.N_s += 1
        return 'scatter'


def update_score(score, counter, N_b):
    """Tallies the particle information whenever a history ends."""
    score.N_cumulative += N_b
    score.S_f += counter.N_f
    score.S_s += counter.N_s
    score.S_c += counter.N_c
    score.S_l += counter.N_l
    score.S_t += counter.N_c + counter.N_s + counter.N_f
    score.D_cumulative += counter.D


def estimate_keff(score, params):
    """Estimates the k effective for the slab."""
    # fission estimator
    k_eff_f = (params.nu_bar * score.S_f) / score.N_cumulative

    # collision estimator
    k_eff_c = ((score.S_t) * (params.nu_bar * (params.Sig_f / params.Sig_t))) / score.N_cumulative

    # absorption estimator
    k_eff_a = ((score.S_c + score.S_f) / score.N_cumulative) * (params.nu_bar * (params.Sig_f / (params.Sig_c + params.Sig_f)))

    # track length estimator
    k_eff_t = (score.D_cumulative * params.Sig_f * params.nu_bar) / score.N_cumulative

    # format printing
    s = 'k_eff_f: {:10.8f}\tk_eff_c: {:10.8f}\tk_eff_a: {:10.8f}\tk_eff_t: {:10.8f}'.format(k_eff_f, k_eff_c, k_eff_a, k_eff_t)
    print(s)
    return k_eff_t


def run_batch(N_b, q, params):
    """Runs a single batch of particles."""
    # initialize counter
    counter = Counter(params.n_t)

    # loop through N_b histories
    for i in range(N_b):
        par = Particle()
        while not leaked(par.position, params):
            par.position = sample_source_position(q, params)
            par.direction = choose_direction()
            d = choose_path_length(params)
            new_position = par.position + d * par.direction

            # if it leaks, add it to leaked particles, else pick a reaction
            if leaked(new_position, params):
                # add to leaked particles
                if par.direction < 0:
                    counter.D += (params.T + par.position) / abs(par.direction)
                else:
                    counter.D += (params.T - par.position) / abs(par.direction)
                counter.N_l += 1
                par.position = new_position
            else:
                # update track length and pick a reaction
                counter.D += d
                par.position = new_position
                reaction = determine_reaction_type_and_update_counter(par, params, counter)
                if reaction in ('fission', 'capture'):
                    break
    return counter


def run(N_b, n_b, q, params, plot_fiss=False):
    """Runs a series of batches to estimate the slab k effective."""
    # initialize scores and problem parameters and plotting
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    score = Score()

    # loop through n_b batches
    print('\nRunning Thickness:\t{}'.format(params.T))
    for i in range(n_b):
        counter = run_batch(N_b, q, params)
        q.update_dist(counter.fission_spec)
        update_score(score, counter, N_b)
        k_eff = estimate_keff(score, params)

        # plot updated fission probabilities
        if plot_fiss:
            x = (params.edges[1:] + params.edges[:-1]) / 2
            ax.set_xlabel('$x (cm)$')
            ax.set_ylabel('$\Phi (cm^{-2}s^{-1})$')
            phi = (counter.fission_spec / params.Sig_f) / N_b
            np.save('p1.npy', phi)
            ax.plot(x, phi, 'k', alpha=i/n_b)
            fig.savefig('t_crit.png', dpi=250)

    return k_eff


if __name__ == '__main__':
    ts = np.linspace(15, 40, 6)
    k_eff = []
    for T in ts:
        params = Slab_Parameters(T, 0.011437, 0.05, 0.013, 2.5, 64)
        q = Source_PDF(params.n_t)
        k_eff.append(run(10000, 10, q, params))

    # t_crit case
    params = Slab_Parameters(26.872, 0.011437, 0.05, 0.013, 2.5, 64)
    q = Source_PDF(params.n_t)
    run(1000000, 30, q, params, plot_fiss=True)

    # plotting
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.set_xlabel('T (cm)')
    ax.set_ylabel('$k_{eff}$')

    ax.plot(ts, k_eff, 'k')
    fig.savefig('k_vs_T.png', dpi=250)
