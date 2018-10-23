
class Slab_Parameters(object):
    """Container for all problem parameters."""
    def __init__(self):
        pass


class Tally(object):
    """Container for problem tally information."""
    def __init__(self):
        pass


class Particle(object):
    """Container for an individual particle's information."""
    def __init__(self):
        pass


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


def run_batch():
    """Runs a single batch of particles."""
    pass


def run():
    """Runs a series of batches to estimate the slab k effective."""
    pass
