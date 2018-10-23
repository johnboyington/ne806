from numpy.random import rand


class Slab_Parameters(object):
    """Container for all problem parameters."""
    def __init__(self, T, Sig_f, Sig_c, Sig_s):
        self.T = T
        self.Sig_f = Sig_f
        self.Sig_c = Sig_c
        self.Sig_s = Sig_s
        self.Sig_t = Sig_f + Sig_c + Sig_s
        self.nu_bar = 2.54
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


class Source_PDF(object):
    """Container for source distribution."""
    def __init__(self):
        pass


def initialize():
    """Initializes problem tally and source pdf."""
    pass


def sample_source_position():
    """Chooses a source position from the source pdf."""
    pass


def choose_direction():
    """Chooses a direction for the particle."""
    pass


def choose_path_length():
    """Determines the path length for the particle to travel."""
    pass


def determine_reaction_type():
    """Chooses a reaction type for the particle."""
    pass


def update_scores():
    """Tallies the particle information whenever a history ends."""
    pass


def estimate_keff():
    """Estimates the k effective for the slab."""
    pass


def run_batch(N_b, q):
    """Runs a single batch of particles."""
    # initialize counter
    counter = Counter()

    # loop through N_b histories
    for i in range(N_b):
        par = Particle()
        par.position = sample_source_position()


def run(N_b, n_b, q, params):
    """Runs a series of batches to estimate the slab k effective."""
    # initialize scores and problem parameters
    
    score = Score()

    # loop through n_b batches
    for i in range(n_b):
        run_batch(N_b, q)


if __name__ == '__main__':
    params = Slab_Parameters()
    q = Source_PDF()
    run(1000, 10, params)
