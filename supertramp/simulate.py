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
import random
import collections
import argparse
import supertramp
from supertramp import utility
from supertramp.BitVector import BitVector
import dendropy

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
            island):
        self.index = self.__class__.counter
        self.__class__.counter += 1
        self.habitat_type = habitat_type
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
            run_logger=None):
        self.index = self.__class__.counter
        self.__class__.counter += 1
        self.rng = rng
        self.label = label
        self.habitat_types = habitat_types
        self.habitat_list = []
        self.habitats_by_type = {}
        self.run_logger = run_logger

        # construct habitats
        for ht_idx, ht in enumerate(self.habitat_types):
            h = Habitat(
                    habitat_type=ht,
                    island=self)
            self.habitat_list.append(h)
            self.habitats_by_type[ht] = h

        # initialize dispersal regime
        self._dispersal_rates = {}
        for ht in self.habitat_types:
            self._dispersal_rates[ht] = {}

    def __str__(self):
        return self.label

    def set_dispersal_rate(self, habitat_type, dest_island, rate):
        """
        Set a specific dispersal.
        """
        self._dispersal_rates[habitat_type][dest_island] = rate

    def run_dispersals(self):
        for habitat_type in self._dispersal_rates:
            for dest_island in self._dispersal_rates[habitat_type]:
                habitat = self.habitats_by_type[habitat_type]
                rate = self._dispersal_rates[habitat_type][dest_island]
                if not habitat.lineages or rate <= 0.0:
                    continue
                if self.rng.uniform(0, 1) <= rate:
                    lineage = self.rng.choice(list(habitat.lineages))
                    self.run_logger.debug("{lineage}, with habitat type '{habitat_type}', dispersing from island {island1} to {island2}, ".format(
                        island1=self.label,
                        island2=dest_island.label,
                        habitat_type=lineage.habitat_type,
                        lineage=lineage.logging_label))
                    dest_island.receive_migrant(lineage=lineage, habitat_type=lineage.habitat_type)

    def process_migrants(self):
        for habitat in self.habitat_list:
            habitat.process_migrants()

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
        self.island_habitat_localities = None
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
        assert self.island_habitat_localities is None
        assert self.habitats is None
        self.habitat_types = BitVector(size=self.system.num_habitat_types)
        self.island_habitat_localities = BitVector(size=self.system.num_islands)
        self.habitats = BitVector(size=self.system.num_islands * self.system.num_habitat_types)

    def register_habitat(self, habitat):
        self.habitats[habitat.index] = 1
        self.island_habitat_localities[habitat.island.index] = 1
        self.habitat_types[habitat.habitat_type.index] = 1

    def deregister_habitat(self, habitat):
        self.habitats[habitat.index] = 0
        self.island_habitat_localities[habitat.island.index] = 0
        # TODO: simply because it is removed from one particular habitat on one
        # particular island, does not mean that it is no longer associated with
        # this habitat type!!!
        # self.habitat_types[habitat.habitat_type.index] = 0

    def _get_label(self):
        return "S{:d}.{}".format(self.index, self.distribution_label)
    def _set_label(self, v):
        self._label = v
    label = property(_get_label, _set_label)

    def _get_logging_label(self):
        return "<Lineage S{:d}: {}.{}>".format(
                self.index,
                self.island_habitat_localities,
                self.habitat_types)
    logging_label = property(_get_logging_label)

    @property
    def distribution_label(self):
        # this label gets locked to `final_distribution_label` when the species
        # diversifies
        if self.final_distribution_label is not None:
            return self.final_distribution_label
        return "{}.{}".format(self.island_habitat_localities, self.habitat_types)

    def add_age_to_tips(self, ngens=1):
        """
        Grows tree by adding ``ngens`` time unit(s) to all tips.
        """
        if self._child_nodes:
            for nd in self.leaf_iter():
                nd.edge.length += ngens
        else:
            self.edge.length += ngens

    def diversify(self, finalize_distribution_label=True, nsplits=1):
        """
        Spawns two child lineages with self as parent.
        Returns tuple consisting of these two lineages.
        """
        if self._child_nodes:
            raise Exception("Trying to diversify internal node: {}: {}".format(self.label, ", ".join(c.label for c in self._child_nodes)))
        if finalize_distribution_label:
            self.final_distribution_label = self.distribution_label
        children = []
        for i in range(nsplits+1):
            c1 = Lineage(habitat_type=self.habitat_type, system=self.system)
            children.append(c1)
            self.add_child(c1)
            assert c1.parent_node is self
        return children

    def _debug_check_dump_biogeography(self, out):
        out.write("[{}:{}:{}:  ".format(id(self), self.index, self.label))
        out.write("islands='{}'  ".format(self.island_habitat_localities))
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

class SupertrampSimulator(object):

    @staticmethod
    def simulation_model_arg_parser():
        parser = argparse.ArgumentParser(add_help=False)
        model_landscape_options = parser.add_argument_group("MODEL: Landscape Configuration")
        model_landscape_options.add_argument("--num-islands",
                type=int,
                default=4,
                help="number of islands (default = %(default)s).")
        model_landscape_options.add_argument("--num-habitat-types",
                type=int,
                default=3,
                help="number of habitat types per island (default = %(default)s).")
        model_diversification_submodel_params = parser.add_argument_group("MODEL: Diversification Submodel Parameters")
        model_diversification_submodel_params.add_argument("-a", "--diversification-model-a",
                type=float,
                default=-0.5,
                help="'a' parameter of the diversfication model (default: %(default)s).")
        model_diversification_submodel_params.add_argument("-b", "--diversification-model-b",
                type=float,
                default=0.5,
                help="'b' parameter of the diversfication model (default: %(default)s).")
        model_diversification_submodel_params.add_argument("-s", "--diversification-model-s0", "--s0",
                type=float,
                default=0.001,
                help="'s' parameter of the diversfication model (default: %(default)s).")
        model_diversification_submodel_params.add_argument("-e", "--diversification-model-e0", "--e0",
                type=float,
                default=0.0001,
                help="'e' parameter of the diversfication model (default: %(default)s).")
        model_dispersal_submodel_params = parser.add_argument_group("MODEL: Dispersal Submodel Parameters")
        model_dispersal_submodel_params.add_argument("--dispersal-model",
                type=str,
                default="unconstrained",
                choices=["constrained", "unconstrained"],
                help="Dispersal model: constrained or unconstrained by habitat")
        model_dispersal_submodel_params.add_argument("-d", "--dispersal-rate",
                default=0.01,
                type=float,
                help="Dispersal rate (default = %(default)s).")
        lineage_evolution_submodel_params = parser.add_argument_group("MODEL: Lineage Evolution Submodel Parameters")
        lineage_evolution_submodel_params.add_argument("-y", "--niche-evolution-probability",
                default=0.01,
                type=float,
                help="Lineage (post-splitting) niche evolution probability (default = %(default)s).")
        return parser

    def __init__(self, **kwargs):
        self.reset_system_globals()
        self.configure_simulator(kwargs)
        self.set_model(kwargs)
        self.sympatric_speciation = False
        if kwargs:
            raise TypeError("Unsupported configuration keywords: {}".format(kwargs))
        self.bootstrap()

    def reset_system_globals(self):
        HabitatType.reset_counter()
        Habitat.reset_counter()
        Island.reset_counter()
        Lineage.reset_counter()
        self.current_gen = 0
        self.habitat_types = []
        self.islands = []
        self.phylogeny = None

    def set_model(self, model_params_d):

        # Landscape: islands
        self.island_labels = model_params_d.pop("island_labels", None)
        if self.island_labels is None:
            self.island_labels = []
            num_islands = model_params_d.pop("num_islands", 4)
            for i in range(num_islands):
                label = "I{}".format(i+1)
                self.island_labels.append(label)
        else:
            if num_islands in model_params_d and num_islands != len(self.island_labels):
                raise ValueError("Number of islands requested ({}) does not match number of islands specified ({})".format(
                    model_params_d["num_islands"], len(self.island_labels)))
        self.run_logger.info("Configuring {} islands: {}".format(
            len(self.island_labels), self.island_labels))

        # Landscape: habitats
        self.habitat_type_labels = model_params_d.pop("habitat_type_labels", None)
        if self.habitat_type_labels is None:
            num_habitat_types = model_params_d.pop("num_habitat_types", 3)
            self.habitat_type_labels = []
            for i in range(num_habitat_types):
                label = "H{}".format(i+1)
                self.habitat_type_labels.append(label)
        else:
            if num_habitat_types in model_params_d and num_habitat_types != len(self.habitat_type_labels):
                raise ValueError("Number of habitat_types requested ({}) does not match number of habitat_types specified ({})".format(
                    model_params_d["num_habitat_types"], len(self.habitat_type_labels)))
        self.run_logger.info("Configuring {} habitat types per island: {}".format(
            len(self.habitat_type_labels), self.habitat_type_labels))

        # Dispersal submodel
        self.dispersal_model = model_params_d.pop("dispersal_model", "unconstrained")
        self.run_logger.info("Dispersal model category: '{}'".format(self.dispersal_model))
        self.global_dispersal_rate = model_params_d.pop("dispersal_rate", 0.01)
        self.run_logger.info("Dispersal rate, d: {}".format(self.global_dispersal_rate))

        # Diversification submodel
        self.diversification_model_a = model_params_d.pop("diversification_model_a", -0.5)
        self.run_logger.info("Diversification model, a: {}".format(self.diversification_model_a))
        self.diversification_model_b = model_params_d.pop("diversification_model_b", 0.5)
        self.run_logger.info("Diversification model, b: {}".format(self.diversification_model_b))
        self.diversification_model_s0 = model_params_d.pop("diversification_model_s0", 0.1)
        self.run_logger.info("Diversification model, s0: {}".format(self.diversification_model_s0))
        self.diversification_model_e0 = model_params_d.pop("diversification_model_e0", 0.001)
        self.run_logger.info("Diversification model, e0: {}".format(self.diversification_model_e0))
        self.run_logger.info("Projected habitat species richness (s0/e0): {}".format(self.diversification_model_s0/self.diversification_model_e0))

        # Nice Shift/Evolution submodel
        self.global_lineage_niche_evolution_probability = model_params_d.pop("niche_evolution_probability", 0.01)
        self.run_logger.info("Niche evolution probability: {}".format(self.global_lineage_niche_evolution_probability))

    def configure_simulator(self, configd):

        self.output_prefix = configd.pop("output_prefix", "supertramp")

        self.run_logger = configd.pop("run_logger", None)
        if self.run_logger is None:
            self.run_logger = utility.RunLogger(name="supertramp",
                    log_path=self.output_prefix + ".log")
        self.run_logger.system = self

        self.debug_mode = configd.pop("debug_mode", False)
        if self.debug_mode:
            self.run_logger.info("Running in DEBUG mode")

        self.name = configd.pop("name", None)
        if self.name is None:
            self.name = str(id(self))
        self.run_logger.info("Configuring simulation '{}'".format(self.name))

        self.tree_log = configd.pop("tree_log", None)
        if self.tree_log is None:
            self.tree_log = open(self.output_prefix + ".trees", "w")
        self.run_logger.info("Tree log filepath: {}".format(self.tree_log.name))

        self.general_stats_log = configd.pop("general_stats_log", None)
        if self.general_stats_log is None:
            self.general_stats_log = open(self.output_prefix + ".trees", "w")
        self.run_logger.info("Statistics log filepath: {}".format(self.general_stats_log.name))

        self.rng = configd.pop("rng", None)
        if self.rng is None:
            self.random_seed = configd.pop("random_seed", None)
            if self.random_seed is None:
                self.random_seed = random.randint(0, sys.maxsize)
            self.run_logger.info("Initializing with random seed {}".format(self.random_seed))
            # self.rng = numpy.random.RandomState(seed=[self.random_seed])
            self.rng = random.Random(self.random_seed)
        else:
            if "random_seed" in configd:
                raise TypeError("Cannot specify both 'rng' and 'random_seed'")
            self.run_logger.info("Using existing random number generator")

        self.log_frequency = configd.pop("log_frequency", 1000)
        self.report_frequency = configd.pop("report_frequency", None)

    def poisson_rv(rate):
        """
        Returns a random number from a Poisson distribution with rate of
        `rate` (mean of 1/rate).
        """
        MAX_EXPECTATION = 64.0 # larger than this and we have underflow issues
        if rate > MAX_EXPECTATION:
            r = rate/2.0
            return poisson_rv(r) + poisson_rv(r)
        L = math.exp(-1.0 * rate)
        p = 1.0
        k = 0.0
        while p >= L:
            k = k + 1.0
            u = self.rng.random()
            p = p * u
        return int(k - 1.0)

    def bootstrap(self):

        # create habitat types
        for ht_label in self.habitat_type_labels:
            h = HabitatType(label=ht_label)
            self.habitat_types.append(h)
        self.all_habitat_types_bitmask = (1 << len(self.habitat_types)) - 1

        # create islands
        for island_label in self.island_labels:
            island = Island(
                    rng=self.rng,
                    label=island_label,
                    habitat_types=self.habitat_types,
                    run_logger=self.run_logger)
            self.islands.append(island)
        self.all_islands_bitmask = (1 << len(self.islands)) - 1

        # set up dispersal regime
        if self.dispersal_model == "unconstrained":
            self.dispersal_source_habitat_types = list(self.habitat_types)
        else:
            self.dispersal_source_habitat_types = [self.habitat_types[0]]

        # sum of rates of dispersing out of any island == global dispersal rate
        dispersal_rates = []
        if len(self.islands) <= 1:
            if self.global_dispersal_rate > 0:
                self.run_logger.info("Only {} island: forcing dispersal rate to 0.0".format(len(self.islands)))
            island_dispersal_rate = 0
        else:
            island_dispersal_rate = float(self.global_dispersal_rate) / ((len(self.islands) * (len(self.islands) - 1)))
        habitat_dispersal_rates = {}
        for idx, habitat_type in enumerate(self.habitat_types):
            if habitat_type in self.dispersal_source_habitat_types:
                disp_rate = island_dispersal_rate / len(self.dispersal_source_habitat_types)
            else:
                disp_rate = 0.0
            habitat_dispersal_rates[habitat_type] = disp_rate
        for isl1 in self.islands:
            for isl2 in self.islands:
                if isl1 is not isl2:
                    for habitat_type in self.habitat_types:
                        disp_rate = habitat_dispersal_rates[habitat_type]
                        isl1.set_dispersal_rate(habitat_type, isl2, disp_rate)
                        self.run_logger.info("Island {} to {} dispersal rate for habitat {} = {}".format(
                            isl1, isl2, habitat_type, disp_rate))
                        dispersal_rates.append(disp_rate)
        total_dispersal_rate = sum(dispersal_rates)
        if len(self.islands) > 1 and abs(total_dispersal_rate - self.global_dispersal_rate) > 1e-8:
            self.run_logger.critical("Error in dispersal rate distribution: {} != {}: {}".format(
                total_dispersal_rate, self.global_dispersal_rate, dispersal_rates))
            sys.exit(1)

        # initialize lineages
        self.seed_habitat = self.dispersal_source_habitat_types[0]
        seed_node = Lineage(habitat_type=self.seed_habitat, system=self)
        self.phylogeny = Phylogeny(seed_node=seed_node)

        # seed lineage
        self.islands[0].habitat_list[0].add_lineage(self.phylogeny.seed_node)

        # begin logging generations
        self.run_logger.system = self

    @property
    def num_islands(self):
        return len(self.islands)

    @property
    def num_habitat_types(self):
        return len(self.habitat_types)

    def run(self, ngens, repeat_on_total_extinction=True):
        for x in range(ngens):
            if ( (self.report_frequency is not None)
                and ( (self.current_gen % self.report_frequency) == 0 )
                ):
                self.report()
            self.execute_life_cycle()
        return True

    def execute_life_cycle(self):
        self.current_gen += 1
        self.phylogeny.add_age_to_tips(1)
        if self.current_gen % self.log_frequency == 0:
            self.run_logger.info("Executing life-cycle {}".format(self.current_gen))
        for island in self.islands:
            island.run_dispersals()
        for island in self.islands:
            island.process_migrants()
        self.run_diversification()

    def run_diversification(self):
        if self.rng.uniform(0, 1) <= 0.5:
            self.run_lineage_birth()
            self.run_lineage_death()
        else:
            self.run_lineage_death()
            self.run_lineage_birth()

    def _debug_check_habitat_000(self):
        for nd in self.phylogeny:
            if nd.is_leaf():
                if str(nd.habitat_types) == "000":
                    print("[{}]\n[{}]\n[{}]\n[{}]\n[{}]".format(
                            nd.label,
                            nd.distribution_label,
                            nd.island_habitat_localities,
                            nd.habitat_types,
                            nd.habitats))
                    print(self.phylogeny._as_newick_string())
                assert str(nd.habitat_types) != "000"

    def run_lineage_birth(self):
        # - Each lineage in each habitat has an (independent) probability of
        #   splitting given by the habitat-specific birth rate.
        # - A particular lineage's probability of splitting in a particular
        #   habitat is independent of the any other lineage splitting in
        #   the same or any other habitat.
        # - A particular lineage's probability of splitting in a particular
        #   habitat is also independent of the *same* lineage going extinct
        #   in any other habitat on any other island.
        # - A particular lineage may speciate in multiple habitats
        #   simultaneously in the same generation.
        lineage_splitting_habitat_localities = collections.defaultdict(set)
        lineage_habitats = collections.defaultdict(set)
        for island in self.islands:
            for habitat in island.habitat_list:
                if not habitat.lineages:
                    continue
                splitting_rate = self.diversification_model_s0 * (len(habitat.lineages) ** self.diversification_model_a)
                for lineage in habitat.lineages:
                    lineage_habitats[lineage].add(habitat)
                    if self.rng.uniform(0, 1) <= splitting_rate:
                        lineage_splitting_habitat_localities[lineage].add(habitat)
        for lineage in lineage_splitting_habitat_localities:
            splitting_habitats = lineage_splitting_habitat_localities[lineage]
            if not self.sympatric_speciation:
                if len(lineage_habitats[lineage]) == 1:
                    # lineage occurs only in a habitat of only one island;
                    # since sympatric speciation is disallowed: skipped
                    continue
                if len(lineage_habitats[lineage]) == len(splitting_habitats):
                    # lineage in splitting in every island in which it occurs;
                    # since sympatric speciation is disallowed: drop one
                    splitting_habitats.remove(self.rng.choice(splitting_habitats))
            children = lineage.diversify(finalize_distribution_label=True,
                    nsplits=len(splitting_habitats))
            self.run_logger.debug("{splitting_lineage} speciating in {num_islands} islands: {islands}".format(
                splitting_lineage=lineage.logging_label,
                num_islands=len(splitting_habitats),
                islands=",".join([habitat.island.label for habitat in splitting_habitats]),
                ))
            assert len(children) == len(splitting_habitats) + 1
            if self.debug_mode:
                try:
                    self.phylogeny._debug_check_tree()
                except AttributeError:
                    self.phylogeny.debug_check_tree()
                self.run_logger.debug("DEBUG MODE: phylogeny structure is valid")
            c0 = children[0]
            c_remaining = set(children[1:])
            if len(self.habitat_types) > 1:
                # assumes all children have the same habitat type
                habitats_to_evolve_into = [ h for h in self.habitat_types if h is not c0.habitat_type ]
                for c1 in c_remaining:
                    if self.rng.uniform(0, 1) <= self.global_lineage_niche_evolution_probability:
                        c1.habitat_type = self.rng.choice(habitats_to_evolve_into)
            c0_placed = False
            for habitat in lineage_habitats[lineage]:
                habitat.remove_lineage(lineage)
                if habitat in splitting_habitats:
                    c1 = c_remaining.pop()
                    habitat.island.add_lineage(lineage=c1, habitat_type=c1.habitat_type)
                    self.run_logger.debug("{splitting_lineage} (with habitat type '{splitting_lineage_habitat_type}') speciating to {daughter_lineage1} (with habitat type '{daughter_lineage1_habitat_type}') in island {island}".format(
                        splitting_lineage=lineage.logging_label,
                        splitting_lineage_habitat_type=lineage.habitat_type.label,
                        daughter_lineage0=c0.logging_label,
                        daughter_lineage1=c1.logging_label,
                        daughter_lineage1_habitat_type=c1.habitat_type.label,
                        island=habitat.island.label,
                        ))
                    if self.sympatric_speciation:
                        habitat.island.add_lineage(lineage=c0, habitat_type=c0.habitat_type)
                        self.run_logger.debug("{splitting_lineage} (with habitat type '{splitting_lineage_habitat_type}') continuing as {daughter_lineage0} in island {island}".format(
                            splitting_lineage=lineage.logging_label,
                            splitting_lineage_habitat_type=lineage.habitat_type.label,
                            daughter_lineage0=c0.logging_label,
                            daughter_lineage1=c1.logging_label,
                            daughter_lineage1_habitat_type=c1.habitat_type.label,
                            island=habitat.island.label,
                            ))
                else:
                    habitat.island.add_lineage(lineage=c0, habitat_type=c0.habitat_type)
                    c0_placed = True
                    self.run_logger.debug("{splitting_lineage} (with habitat type '{splitting_lineage_habitat_type}') continuing as {daughter_lineage0} in island {island}".format(
                        splitting_lineage=lineage.logging_label,
                        splitting_lineage_habitat_type=lineage.habitat_type.label,
                        daughter_lineage0=c0.logging_label,
                        daughter_lineage1=c1.logging_label,
                        daughter_lineage1_habitat_type=c1.habitat_type.label,
                        island=habitat.island.label,
                        ))
            assert c0_placed
            assert len(c_remaining) == 0

    def run_lineage_death(self):
        # - Each lineage in each habitat has an (independent) probability of
        #   going locally extinct (i.e., extirpated from a particular
        #   habitat) given by the local, habitat-specific death rate.
        # - A particular lineage's probability of going extinct in a
        #   particular habitat is independent of the any other lineage going
        #   extinct in the same or any other habitat.
        # - A particular lineage's probability of going extinct in a
        #   particular habitat is also independent of the *same* lineage
        #   going extinct in any other habitat on any other island.
        # - A particular lineage may be extirpated from multiple habitats
        #   simultaneously in the same generation.
        # - Once a lineage is no longer present in any habitat across all
        #   islands, it is considered to have gone globally-extinct and
        #   removed from the system.
        lineage_counts = collections.Counter()
        for island in self.islands:
            for habitat in island.habitat_list:
                lineage_counts.update(habitat.lineages)
                n = len(habitat.lineages)
                if n <= 0:
                    continue
                death_rate = self.diversification_model_e0 * (len(habitat.lineages) ** self.diversification_model_b)
                # print(self.diversification_model_e0 , (len(habitat.lineages) , self.diversification_model_b))
                # print(death_rate)
                to_remove = []
                for lineage in habitat.lineages:
                    if self.rng.uniform(0, 1) <= death_rate:
                        to_remove.append(lineage)
                for lineage in to_remove:
                    self.run_logger.debug("{lineage} extirpated from island {island}".format(
                        lineage=lineage.logging_label,
                        island=habitat.island.label,
                        ))
                    habitat.remove_lineage(lineage)
                    lineage_counts.subtract([lineage])
        for lineage in lineage_counts:
            count = lineage_counts[lineage]
            if count == 0:
                self.run_logger.debug("{lineage} extirpated from all islands and is now globally extinct".format(
                    lineage=lineage.logging_label,
                    ))
                if lineage is self.phylogeny.seed_node:
                    raise TotalExtinctionException()
                else:
                    self.phylogeny.prune_subtree(node=lineage,
                            update_splits=False, delete_outdegree_one=True)
                    if self.phylogeny.seed_node.num_child_nodes() == 0:
                        raise TotalExtinctionException()
            elif self.debug_mode:
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
                self.run_logger.debug("DEBUG MODE: Lineage extinction book-keeping is OK")

    def report_trees(self):
        self.tree_log.write("[&R][simulation={},generation={}]".format(self.name, self.current_gen))
        try:
            self.tree_log.write(self.phylogeny._as_newick_string())
        except AttributeError:
            self.tree_log.write(self.phylogeny.as_newick_string())
        self.tree_log.write(";\n")
        self.tree_log.flush()

    def write_general_stats_header(self, out):
        header = []
        header.append("name")
        header.append("generation")
        for island in self.islands:
            # number of lineages per island, across all habitats in island i
            header.append("island.{}.richness".format(island.label))
        for habitat_type in self.habitat_types:
            # number of lineages per habitat type, across all islands for each island i, I_{i}
            header.append("habitat.type.{}.richness".format(habitat_type.label))
        for island_idx, island in enumerate(self.islands):
            for habitat_type_idx, habitat in enumerate(island.habitat_list):
                # number of lineage in each habitat j of each island i, H_{i,j}
                header.append("{}.{}.richness".format(island.label, habitat.habitat_type.label))
        # mean number of lineages per island, across all habitat types
        header.append("mean.island.richness")
        # mean number of lineages per habitat type, across all islands
        header.append("mean.habitat.type.richness")
        # mean number of lineages per habitat, across all habitat types and islands
        header.append("mean.habitat.richness")
        header = "\t".join(header)
        out.write("{}\n".format(header))
        out.flush()

    def write_general_stats(self, out):
        island_habitat_richness = collections.OrderedDict()
        island_lineages = collections.defaultdict(set)
        habitat_type_lineages = collections.defaultdict(set)
        num_habitats_counted = 0
        sum_of_richness_in_each_habitat = 0
        for island in self.islands:
            for habitat in island.habitat_list:
                n = len(habitat.lineages)
                island_habitat_richness[ (island, habitat) ] = n
                island_lineages[island].update(habitat.lineages)
                habitat_type_lineages[habitat.habitat_type].update(habitat.lineages)
                sum_of_richness_in_each_habitat += n
                num_habitats_counted += 1
        island_richness = collections.OrderedDict()
        for island in self.islands:
            island_richness[island] = len(island_lineages[island])
        mean_island_richness = sum(island_richness.values())/len(island_richness)
        habitat_type_richness = collections.OrderedDict()
        for habitat_type in self.habitat_types:
            habitat_type_richness[habitat_type] = len(habitat_type_lineages[habitat_type])
        mean_habitat_type_richness = sum(habitat_type_richness.values())/len(habitat_type_richness)
        mean_habitat_richness = sum_of_richness_in_each_habitat/num_habitats_counted

        int_template = "{}"
        float_template = "{}"
        stat_values = collections.OrderedDict()
        stat_values["name"] = self.name
        stat_values["generation"] = int_template.format(self.current_gen)
        for island in self.islands:
            stat_values["island.{}.richness".format(island.label)] = int_template.format(island_richness[island])
        for habitat_type in self.habitat_types:
            stat_values["habitat.type.{}.richness".format(habitat_type.label)] = int_template.format(habitat_type_richness[habitat_type])
        for island_idx, island in enumerate(self.islands):
            for habitat_type_idx, habitat in enumerate(island.habitat_list):
                stat_values["{}.{}.richness".format(island.label, habitat.habitat_type.label)] = int_template.format(island_habitat_richness[ (island, habitat) ])
        stat_values["mean.island.richness"] = float_template.format(mean_island_richness)
        stat_values["mean.habitat_type.richness"] = float_template.format(mean_habitat_type_richness)
        stat_values["mean.habitat.richness"] = float_template.format(mean_habitat_richness)

        values = stat_values.values()
        out.write("\t".join(values))
        out.write("\n")

    def report_general_stats(self):
        if not hasattr(self.general_stats_log, "header_written") or not self.general_stats_log.header_written:
            self.write_general_stats_header(self.general_stats_log)
            self.general_stats_log.header_written = True
        self.write_general_stats(self.general_stats_log)

    def report(self):
        self.report_trees()
        self.report_general_stats()

def repeat_run_supertramp(
        model_params_d,
        ngens,
        nreps,
        output_prefix,
        random_seed="None",
        stderr_logging_level="info",
        file_logging_level="debug"):
    """
    Executes multiple runs of the Supertramp simulator under identical
    parameters to produce the specified number of replicates, discarding failed
    runs.

    Parameters
    ----------
    model_params_d : dict
        Simulator model parameters as keyword-value pairs. To be re-used for
        each replicate.
    ngens : integer
        Number of generations for which to run each individual replicate.
    nreps : integer
        Number of replicates to produce. f
    output_prefix : string
        Path prefix for output files.
    random_seed : integer
        Random seed to be used (for single random number generator across all
        replicates).
    stderr_logging_level : string or None
        Message level threshold for screen logs; if 'none' or `None`, screen
        logs will be supprsed.
    file_logging_level : string or None
        Message level threshold for file logs; if 'none' or `None`, file
        logs will be supprsed.
    """
    configd = dict(model_params_d)
    if stderr_logging_level is None or stderr_logging_level.lower() == "none":
        log_to_stderr = False
    else:
        log_to_stderr = True
    if file_logging_level is None or file_logging_level.lower() == "none":
        log_to_file = False
    else:
        log_to_file = True
    configd["run_logger"] = utility.RunLogger(
            name="supertramp",
            log_path=output_prefix + ".log",
            log_to_stderr=log_to_stderr,
            stderr_logging_level=stderr_logging_level,
            log_to_file=log_to_file,
            file_logging_level=file_logging_level,
            )
    run_logger = configd["run_logger"]
    run_logger.info("Starting: {}".format(supertramp.description()))
    if random_seed is None:
        random_seed = random.randint(0, sys.maxsize)
    run_logger.info("Initializing with random seed: {}".format(random_seed))
    configd["rng"] = random.Random(random_seed)
    configd["tree_log"] = open(output_prefix + ".trees",
            "w")
    configd["general_stats_log"] = open(output_prefix + ".general_stats.txt",
            "w")
    configd["general_stats_log"].header_written = False
    header_written = False
    rep = 0
    while rep < nreps:
        simulation_name="Run{}".format((rep+1))
        run_output_prefix = "{}.R{:04d}".format(output_prefix, rep+1)
        run_logger.info("Run {} of {}: starting".format(rep+1, nreps))
        supertramp_simulator = SupertrampSimulator(
                name=simulation_name,
                **configd)
        success = False
        while not success:
            try:
                success = supertramp_simulator.run(ngens)
            except TotalExtinctionException:
                run_logger.info("Run {} of {}: [t={}] total extinction of all lineages before termination condition".format(rep+1, nreps, supertramp_simulator.current_gen))
                run_logger.info("Run {} of {}: restarting".format(rep+1, nreps))
                supertramp_simulator = SupertrampSimulator(
                        name=simulation_name,
                        **configd)
            else:
                run_logger.info("Run {} of {}: completed to termination condition of {} generations".format(rep+1, nreps, ngens))
                supertramp_simulator.report()
                break
        rep += 1
