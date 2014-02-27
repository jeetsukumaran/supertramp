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

import random
import collections
import numpy

def weighted_choice(seq, weights, rng=None):
    """
    Selects an element out of seq, with probabilities of each element
    given by the list `weights` (which must be at least as long as the
    length of `seq` - 1).
    """
    if rng is None:
        rng = GLOBAL_RNG
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
    if rng is None:
        rng = GLOBAL_RNG
    rnd = rng.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

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
            self.populations[lineage] += 1
        self.migrants = []

    def reproduce_populations(self):
        """
        Constitutes next generation.
        """
        self.populations.clear()
        lineages, probs = self._get_lineage_selection_probabilities()
        pop_sizes = numpy.random.multinomial(self.carrying_capacity, probs)
        assert len(lineages) == len(pop_sizes)
        for lineage, pop_size in zip(lineages, pop_sizes):
            self.populations[lineage] = pop_size

    def run_dispersals(self):
        """
        Disperse.
        """
        if not self._dispersal_dest_list:
            self.compile_dispersal_rates()
        if self.rng.random() <= self._aggregate_rate_of_dispersal:
            dest = weighted_choice(self._dispersal_dest_list, self._dispersal_rate_list, rng=self.rng)
            lineage = weighted_choice(self.populations.keys(), self.populations.values(), rng=self.rng)
            dist.migrants.append(lineage)

    def clean_up(self):
        lineages_to_delete = [lineage for lineage in self.populations if self.populations[lineages] <= 0]
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
        self._dispersal_dest_list.clear()
        self._dispersal_rate_list.clear()
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
            fitnesses[idx] = fitnesses[idx]/f
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
    pass

class System(object):

    def __init__(self):
        self.rng = random.Random()
        self.habitats = []

    def bootstrap(self):
        for i in range(4):
            h = Habitat(
                    system=self,
                    island=None)
            self.habitats.append(h)
        for h1 in self.habitats:
            for h2 in self.habitats:
                if h1 is not h2:
                    h1.set_dispersal_rate(h2, 0.001)

    def run(self, ngens):
        for i in range(ngens):
            for h in self.habitats:
                h.execute_lifecycle()

def main():
    sys = System()
    sys.bootstrap()
    sys.run(100)

if __name__ == "__main__":
    main()



