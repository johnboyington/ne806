import numpy as np
from numpy.random import rand


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
    def __init__(self):
        self.N = 0
        self.N_l = 0
        self.N_c = 0
        self.N_f = 0
        self.N_s = 0
        self.D = 0
        return


class Score(object):
    """Container for problem score information."""
    def __init__(self):
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
    return False if position is None else (abs(position) < params.T)


def sample_source_position(q, params):
    """Chooses a source position from the source pdf."""
    i = np.searchsorted(q.cdf, rand())
    l, r = params.edges[i:i+2]
    return ((r - l) * rand()) + l


def choose_direction():
    """Chooses a direction for the particle."""
    return 1 if rand() > 0.5 else -1


def choose_path_length(params):
    """Determines the path length for the particle to travel."""
    return -(1 / params.Sig_t) * np.log(rand())


def determine_reaction_type_and_update_counter(params, counter):
    """Chooses a reaction type for the particle."""
    rho = rand()
    if rho <= params.Sig_f / params.Sig_t:
        counter.N_f += 1
        return 'fission'
    elif rho < (params.Sig_f + params.Sig_c) / params.Sig_t:
        counter.N_c += 1
        return 'capture'
    else:
        counter.N_s += 1
        return 'scatter'


def update_scores():
    """Tallies the particle information whenever a history ends."""
    pass


def estimate_keff():
    """Estimates the k effective for the slab."""
    pass


def run_batch(N_b, q, params):
    """Runs a single batch of particles."""
    # initialize counter
    counter = Counter()

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
                counter.D += (params.T - par.position) * par.direction
                counter.N_l += 1
            else:
                # update track length and pick a reaction
                counter.D += d
                par.position = new_position
                reaction = determine_reaction_type_and_update_counter(params, counter)
                if reaction in ('fission', 'capture'):
                    break


def run(N_b, n_b, q, params):
    """Runs a series of batches to estimate the slab k effective."""
    # initialize scores and problem parameters
    
    score = Score()

    # loop through n_b batches
    for i in range(n_b):
        run_batch(N_b, q, params)


if __name__ == '__main__':
    params = Slab_Parameters(40, 0.011437, 0.05, 0.013, 2.5, 10)
    q = Source_PDF(params.n_t)
    run(1000, 10, q, params)
