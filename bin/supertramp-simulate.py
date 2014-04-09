#! /usr/bin/env python

##############################################################################
##
##  Copyright 2010-2014 Jeet Sukumaran.
##  All rights reserved.
##
##  Redistribution and use in source and binary forms, with or without
##  modification, are permitted provided that the following conditions are met:
##
##      * Redistributions of source code must retain the above copyright
##        notice, this list of conditions and the following disclaimer.
##      * Redistributions in binary form must reproduce the above copyright
##        notice, this list of conditions and the following disclaimer in the
##        documentation and/or other materials provided with the distribution.
##      * The names of its contributors may not be used to endorse or promote
##        products derived from this software without specific prior written
##        permission.
##
##  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
##  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
##  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
##  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN OR MARK T. HOLDER
##  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
##  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
##  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
##  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
##  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
##  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
##  POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import sys
import argparse
import random
import logging
import collections
import inspect
from supertramp.BitVector import BitVector
import dendropy

_LOGGING_LEVEL_ENVAR = "SUPERTRAMP_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "SUPERTRAMP_LOGGING_FORMAT"
_DEBUG_MODE = False

def dump_stack():
    for frame, filename, line_num, func, source_code, source_index in inspect.stack()[2:]:
        if source_code is None:
            print("{}: {}".format(filename, line_num))
        else:
            print("{}: {}: {}".format(filename, line_num, source_code[source_index].strip()))

def weighted_choice(seq, weights, rng=None):
    """
    Selects an element out of seq, with probabilities of each element
    given by the list `weights` (which must be at least as long as the
    length of `seq` - 1).
    """
    if weights is None:
        weights = [1.0/len(seq) for count in range(len(seq))]
    else:
        weights = list(weights)
    if len(weights) < len(seq) - 1:
        raise Exception("Insufficient number of weights specified")
    if len(weights) == len(seq) - 1:
        weights.append(1 - sum(weights))
    return seq[weighted_index_choice(weights, rng)]

def weighted_index_choice(weights, rng=None):
    """
    (From: http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/)
    The following is a simple function to implement weighted random choice in
    Python. Given a list of weights, it returns an index randomly, according
    to these weights [1].
    For example, given [2, 3, 5] it returns 0 (the index of the first element)
    with probability 0.2, 1 with probability 0.3 and 2 with probability 0.5.
    The weights need not sum up to anything in particular, and can actually be
    arbitrary Python floating point numbers.
    If we manage to sort the weights in descending order before passing them
    to weighted_choice_sub, it will run even faster, since the random call
    returns a uniformly distributed value and larger chunks of the total
    weight will be skipped in the beginning.
    """
    rnd = rng.uniform(0, 1) * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

class RunLogger(object):

    def __init__(self, **kwargs):
        self.name = kwargs.get("name", "RunLog")
        self._log = logging.getLogger(self.name)
        self._log.setLevel(logging.DEBUG)
        if kwargs.get("log_to_stderr", True):
            ch1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", logging.INFO))
            ch1.setLevel(stderr_logging_level)
            ch1.setFormatter(self.get_default_formatter())
            self._log.addHandler(ch1)
        if kwargs.get("log_to_file", True):
            log_stream = kwargs.get("log_stream", \
                open(kwargs.get("log_path", self.name + ".log"), "w"))
            ch2 = logging.StreamHandler(log_stream)
            file_logging_level = self.get_logging_level(kwargs.get("file_logging_level", logging.DEBUG))
            ch2.setLevel(file_logging_level)
            ch2.setFormatter(self.get_default_formatter())
            self._log.addHandler(ch2)

    def get_logging_level(self, level=None):
        if level in [logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING,
            logging.ERROR, logging.CRITICAL]:
            return level
        elif level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        if level_name == "NOTSET":
            level = logging.NOTSET
        elif level_name == "DEBUG":
            level = logging.DEBUG
        elif level_name == "INFO":
            level = logging.INFO
        elif level_name == "WARNING":
            level = logging.WARNING
        elif level_name == "ERROR":
            level = logging.ERROR
        elif level_name == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
        return level

    def get_default_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_rich_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(filename)s (%(lineno)d): %(levelname) 8s: %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_simple_formatter(self):
        return logging.Formatter("%(levelname) 8s: %(message)s")

    def get_raw_formatter(self):
        return logging.Formatter("%(message)s")

    def get_logging_formatter(self, format=None):
        if format is not None:
            format = format.upper()
        elif _LOGGING_FORMAT_ENVAR in os.environ:
            format = os.environ[_LOGGING_FORMAT_ENVAR].upper()
        if format == "RICH":
            logging_formatter = self.get_rich_formatter()
        elif format == "SIMPLE":
            logging_formatter = self.get_simple_formatter()
        elif format == "NONE":
            logging_formatter = self.get_raw_formatter()
        else:
            logging_formatter = self.get_default_formatter()
        if logging_formatter is not None:
            logging_formatter.datefmt='%H:%M:%S'

    def debug(self, msg, *args, **kwargs):
        self._log.debug(msg, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self._log.info(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self._log.warning(msg, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self._log.error(msg, *args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self._log.critical(msg, *args, **kwargs)

class HabitatType(object):

    counter = 0

    def reset_counter(cls):
        cls.counter = 0
    reset_counter = classmethod(reset_counter)

    def __init__(self, label):
        self.index = self.__class__.counter
        self.__class__.counter += 1
        self.label = label

    def __str__(self):
        return self.label

class Habitat(object):

    counter = 0

    def reset_counter(cls):
        cls.counter = 0
    reset_counter = classmethod(reset_counter)

    def __init__(self,
            habitat_type,
            carrying_capacity,
            island):
        self.index = self.__class__.counter
        self.__class__.counter += 1
        self.habitat_type = habitat_type
        self.carrying_capacity = carrying_capacity
        self.island = island
        self.lineages = set()
        self.migrants = set()

    def process_migrants(self):
        for lineage in self.migrants:
            self.add_lineage(lineage)
        self.migrants.clear()

    def receive_migrant(self, lineage):
        self.migrants.add(lineage)

    def add_lineage(self, lineage):
        self.lineages.add(lineage)
        lineage.register_habitat(self)

    def remove_lineage(self, lineage):
        self.lineages.remove(lineage)
        lineage.deregister_habitat(self)

    def __str__(self):
        return "{}-{}".format(self.island.label, self.habitat_type.label)

class Island(object):

    counter = 0

    def reset_counter(cls):
        cls.counter = 0
    reset_counter = classmethod(reset_counter)

    def __init__(self,
            rng,
            label,
            habitat_types,
            habitat_carrying_capacities=None,
            dispersal_rates=None,
            dispersal_source_habitat_types=None):
        self.index = self.__class__.counter
        self.__class__.counter += 1
        self.rng = rng
        self.label = label
        self.habitat_types = habitat_types
        if habitat_carrying_capacities is None:
            self.habitat_carrying_capacities = [None for h in self.habitat_types]
        else:
            if len(habitat_carrying_capacities) != len(self.habitat_types):
                raise ValueError("Require exactly {} carrying capacities, but given {}".format(len(self.habitat_types),
                    len(habitat_carrying_capacities)))
            self.habitat_carrying_capacities = habitat_carrying_capacities
        self.habitat_list = []
        self.habitats_by_type = {}
        if dispersal_source_habitat_types is None:
            self.dispersal_source_habitat_types = []
        else:
            self.dispersal_source_habitat_types = dispersal_source_habitat_types
        self.dispersal_source_habitat_list = []

        # Rate of dispersal to other habitats:
        #   - keys = destination (Habitat object)
        #   - values = rate (float)
        self._dispersal_rates = collections.OrderedDict()
        self._dispersal_dest_list = []
        self._dispersal_rate_list = []
        self._aggregate_rate_of_dispersal = 0.0
        if dispersal_rates is not None:
            for dest in dispersal_rates:
                self._dispersal_rates[dest] = dispersal_rates[dest]
            self.compile_dispersal_rates()

        # construct habitats
        for ht_idx, ht in enumerate(self.habitat_types):
            h = Habitat(
                    habitat_type=ht,
                    carrying_capacity=self.habitat_carrying_capacities[ht_idx],
                    island=self)
            self.habitat_list.append(h)
            self.habitats_by_type[ht] = h
            if ht in self.dispersal_source_habitat_types:
                self.dispersal_source_habitat_list.append(h)

    def __str__(self):
        return self.label

    def set_dispersal_rate(self, dest, rate):
        """
        Set a specific dispersal.
        """
        self._dispersal_rates[dest] = rate
        self._dispersal_dest_list = []
        self._dispersal_rate_list = []
        self._aggregate_rate_of_dispersal = 0.0

    def compile_dispersal_rates(self):
        """
        Split dictionary into list of destinations and rates.
        """
        self._dispersal_dest_list = []
        self._dispersal_rate_list = []
        sum_of_rates = 0.0
        for dest in self._dispersal_rates:
            rate = self._dispersal_rates[dest]
            self._dispersal_dest_list.append(dest)
            self._dispersal_rate_list.append(rate)
            sum_of_rates += rate
        self._aggregate_rate_of_dispersal = sum_of_rates

    def run_dispersals(self):
        if self.rng.uniform(0, 1) <= self._aggregate_rate_of_dispersal:
            dest = weighted_choice(self._dispersal_dest_list, self._dispersal_rate_list, rng=self.rng)
            self.disperse_to(dest)

    def process_migrants(self):
        for habitat in self.habitat_list:
            habitat.process_migrants()

    def disperse_to(self, destination_island):
        hx = [h for h in self.dispersal_source_habitat_list if h.lineages]
        if not hx:
            return
        habitat = self.rng.choice(hx)
        lineage = self.rng.choice(list(habitat.lineages))
        destination_island.receive_migrant(lineage=lineage, habitat_type=lineage.habitat_type)

    def receive_migrant(self, lineage, habitat_type):
        assert lineage.habitat_type is habitat_type
        self.habitats_by_type[habitat_type].receive_migrant(lineage)

    def add_lineage(self, lineage, habitat_type):
        assert lineage.habitat_type is habitat_type
        self.habitats_by_type[habitat_type].add_lineage(lineage)

class Lineage(dendropy.Node):

    counter = 0

    def reset_counter(cls):
        cls.counter = 0
    reset_counter = classmethod(reset_counter)

    def __init__(self,
            habitat_type=None,
            system=None):
        self.index = self.__class__.counter
        self.__class__.counter += 1
        super(Lineage, self).__init__()
        self.habitat_type = habitat_type
        self.system = system
        self.habitat_types = None
        self.island_localities = None
        self.habitats = None
        self.final_distribution_label = None
        self.edge.length = 0
        if self.system is not None:
            self.bootstrap()

    def bootstrap(self):
        # Note that all islands and habitat types need to be defined for this
        # to work (or at least, the maximum number of habitat types and islands
        # must be known.
        assert self.habitat_types is None
        assert self.island_localities is None
        assert self.habitats is None
        self.habitat_types = BitVector(size=self.system.num_habitat_types)
        self.island_localities = BitVector(size=self.system.num_islands)
        self.habitats = BitVector(size=self.system.num_islands * self.system.num_habitat_types)

    def register_habitat(self, habitat):
        self.habitats[habitat.index] = 1
        self.island_localities[habitat.island.index] = 1
        self.habitat_types[habitat.habitat_type.index] = 1

    def deregister_habitat(self, habitat):
        self.habitats[habitat.index] = 0
        self.island_localities[habitat.island.index] = 0
        # TODO: simply because it is removed from one particular habitat on one
        # particular island, does not mean that it is no longer associated with
        # this habitat type!!!
        # self.habitat_types[habitat.habitat_type.index] = 0

    def _get_label(self):
        return "S{:d}.{}".format(self.index, self.distribution_label)
    def _set_label(self, v):
        self._label = v
    label = property(_get_label, _set_label)

    @property
    def distribution_label(self):
        # this label gets locked to `final_distribution_label` when the species
        # diversifies
        if self.final_distribution_label is not None:
            return self.final_distribution_label
        return "{}.{}".format(self.island_localities, self.habitat_types)

    def add_age_to_tips(self, ngens=1):
        """
        Grows tree by adding ``ngens`` time unit(s) to all tips.
        """
        if self._child_nodes:
            for nd in self.leaf_iter():
                nd.edge.length += ngens
        else:
            self.edge.length += ngens

    def diversify(self, finalize_distribution_label=True):
        """
        Spawns two child lineages with self as parent.
        Returns tuple consisting of these two lineages.
        """
        if self._child_nodes:
            raise Exception("Trying to diversify internal node: {}: {}".format(self.label, ", ".join(c.label for c in self._child_nodes)))
        if finalize_distribution_label:
            self.final_distribution_label = self.distribution_label
        c1 = Lineage(habitat_type=self.habitat_type, system=self.system)
        c2 = Lineage(habitat_type=self.habitat_type, system=self.system)
        self.add_child(c1)
        self.add_child(c2)
        assert c1.parent_node is self
        assert c2.parent_node is self
        return (c1, c2)

    def _debug_check_dump_biogeography(self, out):
        out.write("[{}:{}:{}:  ".format(id(self), self.index, self.label))
        out.write("islands='{}'  ".format(self.island_localities))
        out.write("habitat_types='{}'  ".format(self.habitat_types))
        out.write("habitats='{}'".format(self.habitats))
        out.write("]\n")

    def num_child_nodes(self):
        try:
            return super(Lineage, self).num_child_nodes()
        except AttributeError:
            return len(self._child_nodes)


class Phylogeny(dendropy.Tree):

    def node_factory(cls, **kwargs):
        return Lineage(**kwargs)
    node_factory = classmethod(node_factory)

    def add_age_to_tips(self, ngens=1):
        """
        Grows tree by adding ``ngens`` time unit(s) to all tips.
        """
        self.seed_node.add_age_to_tips(ngens)

class TotalExtinctionException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class System(object):

    def __init__(self,
            dispersal_model="unconstrained",
            global_per_lineage_birth_prob=0.01,
            density_dependent_lineage_birth_model=False,
            global_per_lineage_death_prob=0.01,
            density_dependent_lineage_death_model=True,
            global_lineage_niche_evolution_probability=0.01,
            global_dispersal_rate=0.01,
            rng=None,
            output_prefix=None,
            tree_log=None,
            random_seed=None,
            logger=None,
            log_frequency=100,
            ):
        if output_prefix is None:
            self.output_prefix = "supertramp"
        else:
            self.output_prefix = output_prefix
        if logger is None:
            self.logger = RunLogger(name="supertramp",
                    log_path=self.output_prefix + ".log")
        else:
            self.logger = logger
        if tree_log is None:
            self.tree_log = open(self.output_prefix + ".trees", "w")
        else:
            self.tree_log = tree_log
        if rng is None:
            if random_seed is None:
                self.random_seed = random.randint(0, sys.maxsize)
            else:
                self.random_seed = random_seed
            self.logger.info("Initializing with random seed {}".format(self.random_seed))
            # self.rng = numpy.random.RandomState(seed=[self.random_seed])
            self.rng = random.Random(self.random_seed)
        else:
            self.rng = rng

        self.island_labels = ["A", "B", "C", "D"]
        self.habitat_type_labels = ["coastal", "interior", "deep"]
        self.island_habitat_carrying_capacity = 10
        self.dispersal_model = dispersal_model
        self.habitat_types = []
        self.islands = []
        self.phylogeny = None

        self.global_per_lineage_birth_prob = global_per_lineage_birth_prob
        self.density_dependent_lineage_birth_model = density_dependent_lineage_birth_model
        self.single_birth_per_generation = True
        self.global_per_lineage_death_prob = global_per_lineage_death_prob
        self.single_death_per_generation = True
        self.density_dependent_lineage_death_model = density_dependent_lineage_death_model
        self.global_lineage_niche_evolution_probability = global_lineage_niche_evolution_probability
        self.global_dispersal_rate = global_dispersal_rate

        self.current_gen = 0
        self.log_frequency = log_frequency

    def reset_system_globals(self):
        HabitatType.reset_counter()
        Habitat.reset_counter()
        Island.reset_counter()
        Lineage.reset_counter()

    def bootstrap(self):

        # reset
        self.reset_system_globals()

        # create habitat types
        for ht_label in self.habitat_type_labels:
            h = HabitatType(label=ht_label)
            self.habitat_types.append(h)
        self.all_habitat_types_bitmask = (1 << len(self.habitat_types)) - 1

        # set up dispersal regime
        if self.dispersal_model == "unconstrained":
            self.dispersal_source_habitat_types = list(self.habitat_types)
        else:
            self.dispersal_source_habitat_types = [self.habitat_types[0]]

        # create islands
        cc = [self.island_habitat_carrying_capacity] * len(self.habitat_types)
        for island_label in self.island_labels:
            island = Island(
                    rng=self.rng,
                    label=island_label,
                    habitat_types=self.habitat_types,
                    habitat_carrying_capacities=cc,
                    dispersal_source_habitat_types=self.dispersal_source_habitat_types)
            self.islands.append(island)
        self.all_islands_bitmask = (1 << len(self.islands)) - 1

        # sum of rates of dispersing out of any island == global dispersal rate
        island_dispersal_rate = self.global_dispersal_rate / ((len(self.islands) ** 2) - 1)
        for isl1 in self.islands:
            for isl2 in self.islands:
                if isl1 is not isl2:
                    isl1.set_dispersal_rate(isl2, island_dispersal_rate)
        for isl1 in self.islands:
            isl1.compile_dispersal_rates()

        # initialize lineages
        self.seed_habitat = self.dispersal_source_habitat_types[0]
        seed_node = Lineage(habitat_type=self.seed_habitat, system=self)
        self.phylogeny = Phylogeny(seed_node=seed_node)

        # seed lineage
        self.islands[0].habitat_list[0].add_lineage(self.phylogeny.seed_node)

    @property
    def num_islands(self):
        return len(self.islands)

    @property
    def num_habitat_types(self):
        return len(self.habitat_types)

    def run(self, ngens, repeat_on_total_extinction=True):
        for x in range(ngens):
            self.execute_life_cycle()
        return True

    def execute_life_cycle(self):
        self.current_gen += 1
        self.phylogeny.add_age_to_tips(1)
        if self.current_gen % self.log_frequency == 0:
            self.logger.info("Executing life-cycle {}".format(self.current_gen))
        for island in self.islands:
            island.run_dispersals()
        for island in self.islands:
            island.process_migrants()
        self.run_diversification()

    def run_diversification(self):
        if self.density_dependent_lineage_birth_model:
            self.run_density_dependent_birth()
        else:
            self.run_simple_birth()
        if self.density_dependent_lineage_death_model:
            self.run_density_dependent_death()
        else:
            self.run_simple_death()

    def _debug_check_habitat_000(self):
        for nd in self.phylogeny:
            if nd.is_leaf():
                if str(nd.habitat_types) == "000":
                    print("[{}]\n[{}]\n[{}]\n[{}]\n[{}]".format(
                            nd.label,
                            nd.distribution_label,
                            nd.island_localities,
                            nd.habitat_types,
                            nd.habitats))
                    print(self.phylogeny._as_newick_string())
                assert str(nd.habitat_types) != "000"

    def run_simple_birth(self):
        tips = self.phylogeny.leaf_nodes()
        if not tips:
            raise TotalExtinctionException()
        if self.rng.uniform(0, 1) > self.global_per_lineage_birth_prob * len(tips):
            return
        lineage_localities = collections.defaultdict(list)
        for island in self.islands:
            for habitat in island.habitat_list:
                for lineage in habitat.lineages:
                    lineage_localities[lineage].append(habitat)
        splitting_lineage = self.rng.choice(list(lineage_localities.keys()))
        c0, c1 = splitting_lineage.diversify(finalize_distribution_label=True)
        if _DEBUG_MODE:
            try:
                self.phylogeny._debug_check_tree()
            except AttributeError:
                self.phylogeny.debug_check_tree()
        if self.rng.uniform(0, 1) <= self.global_lineage_niche_evolution_probability:
            c1.habitat_type = self.rng.choice([ h for h in self.habitat_types if h is not c1.habitat_type ])
        splitting_lineage_localities = lineage_localities[splitting_lineage]
        assert splitting_lineage_localities
        if len(splitting_lineage_localities) == 1:
            target = splitting_lineage_localities[0]
        else:
            target = self.rng.choice(splitting_lineage_localities)
        for habitat in splitting_lineage_localities:
            habitat.remove_lineage(splitting_lineage)
            if habitat is target:
                # sympatric speciation: "old" species retained in original habitat on island
                # new species added to new habitat on island
                habitat.island.add_lineage(lineage=c0, habitat_type=c0.habitat_type)
                habitat.island.add_lineage(lineage=c1, habitat_type=c1.habitat_type)
            else:
                habitat.island.add_lineage(lineage=c0, habitat_type=c0.habitat_type)

    def run_density_dependent_birth(self):
        self.run_simple_birth() # TODO!!!

    def run_simple_death(self):
        tips = self.phylogeny.leaf_nodes()
        if not tips:
            raise TotalExtinctionException()
        if self.rng.uniform(0, 1) > self.global_per_lineage_death_prob * len(tips):
            return
        extinguishing_lineage = self.rng.choice(tips)
        lineage_localities = []
        for island in self.islands:
            for habitat in island.habitat_list:
                if extinguishing_lineage in habitat.lineages:
                    lineage_localities.append(habitat)
        for loc in lineage_localities:
            loc.remove_lineage(extinguishing_lineage)
        if extinguishing_lineage is self.phylogeny.seed_node:
            raise TotalExtinctionException()
        else:
            self.phylogeny.prune_subtree(node=extinguishing_lineage,
                    update_splits=False, delete_outdegree_one=True)
            if self.phylogeny.seed_node.num_child_nodes() == 0:
                raise TotalExtinctionException()

    def run_density_dependent_death(self):
        lineage_counts = collections.Counter()
        for island in self.islands:
            for habitat in island.habitat_list:
                lineage_counts.update(habitat.lineages)
                n = len(habitat.lineages)
                if n <= 0:
                    continue
                assert habitat.carrying_capacity is not None
                K = 4 # habitat.carrying_capacity
                if n > K:
                    while n > K:
                        lineage = self.rng.choice(list(habitat.lineages))
                        habitat.remove_lineage(lineage)
                        lineage_counts.subtract([lineage])
                        n -= 1
                else:
                    weight = 1.0 - ((K-n)//K)
                    prob = self.global_per_lineage_death_prob * weight
                    # print("n={:2d}, K={:2d}, weight={:8.6f}, prob={:8.6f}".format(n, K, weight, prob))
                    if self.rng.uniform(0, 1) <= prob:
                        lineage = self.rng.choice(list(habitat.lineages))
                        habitat.remove_lineage(lineage)
                        lineage_counts.subtract([lineage])
        for lineage in lineage_counts:
            count = lineage_counts[lineage]
            if count == 0:
                if lineage is self.phylogeny.seed_node:
                    raise TotalExtinctionException()
                else:
                    self.phylogeny.prune_subtree(node=lineage,
                            update_splits=False, delete_outdegree_one=True)
                    if self.phylogeny.seed_node.num_child_nodes() == 0:
                        raise TotalExtinctionException()
            elif _DEBUG_MODE:
                ## sanity checking ...
                found = True
                for island in self.islands:
                    for habitat in island.habitat_list:
                        if lineage in habitat.lineages:
                            found = True
                            break
                    if found:
                        break
                assert found, lineage

    def report(self):
        # report_prefix = self.output_prefix + ".T{:08}d".format(self.current_gen)
        self.tree_log.write("[&R] ")
        try:
            self.tree_log.write(self.phylogeny._as_newick_string())
        except AttributeError:
            self.tree_log.write(self.phylogeny.as_newick_string())
        self.tree_log.write(";\n")
        self.tree_log.flush()

def main():
    parser = argparse.ArgumentParser(description="Biogeographical simulator")

    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument("dispersal_model",
            choices=["constrained", "unconstrained"],
            help="Dispersal model: constrained or unconstrained by habitat")
    run_options.add_argument("-z", "--random-seed",
            default=None,
            help="Seed for random number generator engine.")
    run_options.add_argument("-n", "--nreps",
            type=int,
            default=10,
            help="number of replicates (default = %(default)s).")
    run_options.add_argument("-g", "--log-frequency",
            default=100,
            type=int,
            help="Frequency that background progress messages get written to the log (default = %(default)s).")
    run_options.add_argument("--debug-mode",
            action="store_true",
            default=False,
            help="Run in debugging mode.")
    simulation_param_options = parser.add_argument_group("Simulation Parameters")
    simulation_param_options.add_argument("-b", "--birth-probability",
            dest="global_per_lineage_birth_prob",
            default=0.01,
            type=float,
            help="Global per-lineage birth probability (default = %(default)s).")
    simulation_param_options.add_argument("-d", "--death-probability",
            dest="global_per_lineage_death_prob",
            default=0.00,
            type=float,
            help="Global per-lineage death probability (default = %(default)s).")
    simulation_param_options.add_argument("--fixed-death-rate-model", "--disable-density-dependent-lineage-death-model",
            dest="density_dependent_lineage_death_model",
            action="store_false",
            default=True,
            help="By default, the lineage extinction follows a density-dependent "
            "model. Under this default model, once the carrying capacity of "
            "a habitat is exceeded, lineages will be removed at random until the "
            "number of lineages equals the carrying capacity, while if the "
            "carrying capacity of a habitat is not exceeded, lineages will die "
            "with probability given by the global per-lineage death probability "
            "multiplied by 1-(n-K)/K, where where 'n' is the number of lineages "
            "in the habitat and 'K' is the maximum number of lineages (carrying "
            "capacity) of the habitat. Select this option to replace this default "
            "of density-dependent lineage death model with a fixed-death rate model "
            "where the lineage death rate will be fixed to the global per-lineage "
            "death probability through the course of the simulation.")
    simulation_param_options.add_argument("-y", "--niche-evolution-probability",
            default=0.00,
            type=float,
            help="Lineage (post-splitting) niche evolution probability (default = %(default)s).")
    simulation_param_options.add_argument("-D", "--dispersal-rate",
            default=0.01,
            type=float,
            help="Dispersal rate (default = %(default)s).")
    termination_condition_options = parser.add_argument_group("Termination Condition Options")
    termination_condition_options.add_argument("--ngens",
            type=int,
            default=10000,
            help="Number of generations to run (default = %(default)s).")
    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument('-o', '--output-prefix',
        action='store',
        dest='output_prefix',
        type=str,
        default='supertramp_run',
        metavar='OUTPUT-FILE-PREFIX',
        help="Prefix for output files (default='%(default)s').")
    args = parser.parse_args()

    _DEBUG_MODE = args.debug_mode
    logger = RunLogger(name="supertramp",
            log_path=args.output_prefix + ".log")
    if args.random_seed is None:
        args.random_seed = random.randint(0, sys.maxsize)
    logger.info("Initializing with random seed: {}".format(args.random_seed))
    # rng = numpy.random.RandomState(seed=[args.random_seed])
    rng = random.Random(args.random_seed)
    rep = 0
    tree_log = open(args.output_prefix + ".trees", "w")
    args.density_dependent_lineage_birth_model = False
    while rep < args.nreps:
        run_output_prefix = "{}.R{:04d}".format(args.output_prefix, rep+1)
        logger.info("Run {} of {}: starting".format(rep+1, args.nreps))
        supertramp_system = System(
                dispersal_model="unconstrained",
                random_seed=args.random_seed,
                global_per_lineage_birth_prob=args.global_per_lineage_birth_prob,
                density_dependent_lineage_birth_model=args.density_dependent_lineage_birth_model,
                global_per_lineage_death_prob=args.global_per_lineage_death_prob,
                density_dependent_lineage_death_model=args.density_dependent_lineage_death_model,
                global_lineage_niche_evolution_probability=args.niche_evolution_probability,
                global_dispersal_rate=args.dispersal_rate,
                rng=rng,
                output_prefix=run_output_prefix,
                tree_log=tree_log,
                logger=logger,
                log_frequency=args.log_frequency)
        supertramp_system.bootstrap()
        success = False
        while not success:
            try:
                success = supertramp_system.run(args.ngens)
            except TotalExtinctionException:
                logger.info("Run {} of {}: [t={}] total extinction of all lineages before termination condition".format(rep+1, args.nreps, supertramp_system.current_gen))
                logger.info("Run {} of {}: restarting".format(rep+1, args.nreps))
            else:
                logger.info("Run {} of {}: completed to termination condition".format(rep+1, args.nreps))
                supertramp_system.report()
                break
        rep += 1

if __name__ == "__main__":
    main()