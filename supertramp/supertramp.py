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

import sys
import argparse
import random
import logging
import collections
import numpy

_LOGGING_LEVEL_ENVAR = "SUPERTRAMP_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "SUPERTRAMP_LOGGING_FORMAT"

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


class Habitat(object):
    """
    A single node in the landscape graph.
        - Supports populations of individuals from different lineages up to
          ``carrying_capacity``.
        - Fitness of each individual is given by the fitness function evaluated
          on the lineage phenotype.
        - Population constituted through multinomial selection with total
          population size given by ``carrying_capacity``, selection classes
          given by all lineages present in habitat in previous generation, and
          selection weight of each class (lineage) given by fitness of that
          lineage in current habitat.
        - ``dispersal_rates`` gives the probability of dispersal to other
          habitats (either on the same island or other islands). A successful
          dispersal event involves:
                - Selection of a lineage with probability proportional to its
                  representation in the population of this habitat.
                - If the destination habitat does not have any members of the
                  selected lineage present, then a single individual of this
                  lineage is added to the population of the destination
                  habitat.
                - If the destination habitat already has a member of the
                  lineage present, then the ``isolation_factor[s]`` value for
                  the destination island, where ``s`` is the lineage, is
                  reduced to 0, indicating that gene flow has occurred.
    """

    def __init__(self,
            label,
            system,
            island,
            fitness_function=None,
            dispersal_rates=None,
            carrying_capacity=10000):
        """
            ``system``
                Handle to System object.

            ``island``
                Host island.

            ``fitness_function``
                Function object that takes a Lineage object as an argument and
                returns a fitness score. The fitness score can be any floating
                point value. It will serve as a weight for multinomial
                selection.

            ``carrying_capacity``
                Total number of individuals (of any lineage) that can be
                supported by this habitat.
        """

        self.label = label
        self.system = system
        self.rng = system.rng
        self.island = island
        self.carrying_capacity = carrying_capacity

        self._fitness_function = None
        self._cached_lineage_fitnesses = {}
        self.fitness_function = fitness_function

        # How many individuals of each lineage are present in this habitat.
        #   - keys = Lineage
        #   - values = number of individuals (int)
        self.populations = collections.OrderedDict()

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

        self.migrants = []

    def execute_lifecycle(self):
        """
        Execute a single step of the lifecycle.
        """
        self.process_migrants()
        self.reproduce_populations()
        self.run_dispersals()
        self.clean_up()

    def process_migrants(self):
        """
        Integrate migrants from last round of dispersals into population.
        """
        for lineage in self.migrants:
            try:
                self.populations[lineage] += 1
            except KeyError:
                self.populations[lineage] = 1
        self.migrants = []

    def reproduce_populations(self):
        """
        Constitutes next generation.
        """
        lineages, probs = self._get_lineage_selection_probabilities()
        if not lineages:
            return
        self.populations.clear()
        try:
            pop_sizes = numpy.random.multinomial(self.carrying_capacity, probs)
        except ValueError:
            print("::: {}".format(lineages))
            print("::: {}".format(probs))
            raise
        assert len(lineages) == len(pop_sizes)
        for lineage, pop_size in zip(lineages, pop_sizes):
            if pop_size > 0:
                self.populations[lineage] = pop_size

    def run_dispersals(self):
        """
        Disperse.
        """
        if not self._dispersal_dest_list:
            self.compile_dispersal_rates()
        if not self.populations.keys():
            return
        if self.rng.uniform(0, 1) <= self._aggregate_rate_of_dispersal:
            dest = weighted_choice(self._dispersal_dest_list, self._dispersal_rate_list, rng=self.rng)
            lineage = weighted_choice(list(self.populations.keys()), list(self.populations.values()), rng=self.rng)
            dest.migrants.append(lineage)

    def clean_up(self):
        lineages_to_delete = [lineage for lineage in self.populations if self.populations[lineage] <= 0]
        for lineage in lineages_to_delete:
            del(self.populations[lineage])

    def get_fitness(self, lineage):
        """
        Returns fitness score of ``lineage`` in this habitat.
        Fitness is memoized for lazy evaluation.
        """
        try:
            return self._cached_lineage_fitnesses[lineage]
        except KeyError:
            if self._fitness_function is None:
                f = 1.0
            else:
                f = self._fitness_function(lineage)
            self._cached_lineage_fitnesses[lineage] = f
            return f

    def get_fitness(self, lineage):
        """
        Returns fitness score of ``lineage`` in this habitat.
        Fitness is memoized for lazy evaluation.
        """
        try:
            return self._cached_lineage_fitnesses[lineage]
        except KeyError:
            if self._fitness_function is None:
                f = 1.0
            else:
                f = self._fitness_function(lineage)
            self._cached_lineage_fitnesses[lineage] = f
            return f

    def clear_dispersal_rates(self):
        """
        Reset dispersals.
        """
        self._dispersal_rates.clear()

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

    def _get_fitness_function(self):
        """
        Returns the current fitness function.
        """
        return self._fitness_function

    def _set_fitness_function(self, fitness_func):
        """
        Sets the fitness function: a function that takes a Lineage object as an
        argument and returns the fitness score of the lineage in this habitat.
        Can be set to ``None`` to return an arbitrary fitness score of 1.0 for
        all and any lineages.
        Setting the fitness function clears the cached fitness values, forcing
        recalculation of fitnesses on next requests.
        """
        self._fitness_function = fitness_func
        self._lineage_fitnesses = {}

    fitness_function = property(_get_fitness_function, _set_fitness_function)

    def _get_lineage_selection_probabilities(self):
        """
        Returns a tuple consisting of:
            - List of Lineages present in this habitat (population size > 1)
            - List of normalized fitness scores: fitness of each Lineage
              present in this habitat, normalized such that they sum to 1.
        Lineages are indexed in the same order in both lists.
        """
        lineages = []
        fitnesses = []
        total_fitness = 0.0
        for lineage in self.populations:
            pop_size = self.populations[lineage]
            if pop_size > 0:
                lineages.append(lineage)
                fitness = self.get_fitness(lineage)
                fitnesses.append(fitness)
                total_fitness += fitness
        for idx, f in enumerate(fitnesses):
            fitnesses[idx] = fitnesses[idx]/total_fitness
        return lineages, fitnesses

class Island(object):
    """
    An Island is a suite of habitats.
    """

    def __init__(self):

        # How long individuals of each lineage on this island have been
        # isolated from other individuals of the same lineage in other
        # island.
        #   - keys = Lineage
        #   - values = number of generations of isolation (int)
        self.isolation_factor = collections.defaultdict(int)

class Lineage(object):

    counter = 0

    def __init__(self, parent=None):
        Lineage.counter += 1
        self.index = Lineage.counter
        self.age = 0
        self.parent = parent
        self.child_nodes = []

    def add_age_to_tips(self, ngens=1):
        """
        Grows tree by adding ``ngens`` time unit(s) to all tips.
        """
        if self.parent is None:
            self.age += 1
        else:
            for nd in self.leaf_iter():
                nd.age += 1

    def leaf_iter(self, filter_fn=None):
        """
        Returns an iterator over the leaf_nodes that are descendants of self
        (with leaves returned in same order as a post-order traversal of the
        tree).
        """
        if filter_fn:
            ff = lambda x: x.parent is None and filter_fn(x) or None
        else:
            ff = lambda x: x.parent is None and x or None
        for node in self.postorder_iter(ff):
            yield node

    def postorder_iter(self, filter_fn=None):
        """
        Postorder traversal of the self and its child_nodes.  Returns self
        and all descendants such that a node's child_nodes (and their
        child_nodes) are visited before node.  Filtered by filter_fn:
        node is only returned if no filter_fn is given or if filter_fn
        returns True.
        """
        stack = [(self, False)]
        while stack:
            node, state = stack.pop(0)
            if state:
                if filter_fn is None or filter_fn(node):
                    yield node
            else:
                stack.insert(0, (node, True))
                child_nodes = [(n, False) for n in node.child_nodes]
                child_nodes.extend(stack)
                stack = child_nodes

    def diversify(self):
        """
        Spawns two child lineages with self as parent.
        Returns tuple consisting of these two lineages.
        """
        c1 = Lineage(parent=self)
        c2 = Lineage(parent=self)
        self.child_nodes.append(c1)
        self.child_nodes.append(c2)
        return (c1, c2)

    @property
    def label(self):
        return "s{}".format(self.index)

    def __str__(self):
        return self.label

    def __repr__(self):
        return "<Lineage {}>".format(self.label)

class System(object):

    def __init__(self, random_seed=None):
        self.logger = RunLogger(name="supertramp")

        if random_seed is None:
            self.random_seed = random.randint(0, sys.maxsize)
        else:
            self.random_seed = random_seed
        self.log_frequency = 1
        self.current_gen = 0

        self.global_lineage_birth_rate = 0.01
        self.global_dispersal_rate = 1.0

        self.logger.info("Initializing with random seed {}".format(self.random_seed))
        self.rng = numpy.random.RandomState(seed=[self.random_seed])
        self.habitats = []
        self.seed_lineage = Lineage()

    def bootstrap(self):
        self.logger.info("Bootstrapping ...")
        for i in range(4):
            h = Habitat(
                    label="{}".format(i+1),
                    system=self,
                    island=None)
            self.habitats.append(h)
        for h1 in self.habitats:
            for h2 in self.habitats:
                if h1 is not h2:
                    h1.set_dispersal_rate(h2, self.global_dispersal_rate)
        self.habitats[0].populations[self.seed_lineage] = 1

    def run(self, ngens):
        for i in range(ngens):
            self.current_gen += 1
            self.seed_lineage.add_age_to_tips(1)
            if self.current_gen % self.log_frequency == 0:
                self.logger.info("Executing life-cycle {}".format(self.current_gen))

            # if self.rng.uniform(0, 1) <= self.global_lineage_birth_rate:
            for h in self.habitats:
                h.execute_lifecycle()

            #########################################################################
            ## Pure-birth (habitat-based) diversification process
            ## (hacked in here for now)
            lineage_habitats = {}
            for h in self.habitats:
                for lineage in h.populations:
                    if h.populations[lineage] <= 0:
                        continue
                    try:
                        lineage_habitats[lineage].append(h)
                    except KeyError:
                        lineage_habitats[lineage] = [h]
            birth_rate = len(lineage_habitats) * self.global_lineage_birth_rate
            if self.rng.uniform(0, 1) <= birth_rate:
                target_lineage = self.rng.choice(list(lineage_habitats.keys()))
                c1, c2 = target_lineage.diversify()
                target_habitat = self.rng.choice(lineage_habitats[target_lineage])
                for h in self.habitats:
                    if h is target_habitat:
                        # note that if populations need sync/management, this
                        # will not do
                        h.populations[c2] = h.populations[target_lineage]
                        # print(">>>> {}".format(h.populations[c2]))
                        del h.populations[target_lineage]
                    else:
                        if target_lineage in h.populations:
                            # note that if populations need sync/management, this
                            # will not do
                            h.populations[c1] = h.populations[target_lineage]
                            # print("---- {}".format(h.populations[c1]))
                            del h.populations[target_lineage]
            #########################################################################

def main():
    parser = argparse.ArgumentParser(description="Biogeographical simulator")

    parser.add_argument("-z", "--random-seed",
            default=None,
            help="Seed for random number generator engine.")
    args = parser.parse_args()
    sys = System(random_seed=args.random_seed)
    sys.bootstrap()
    sys.run(100)

if __name__ == "__main__":
    main()



