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

import collections
import numpy


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
            island,
            carrying_capacity=10000,
            fitness_function=None):
        """
            ``island``
                Host island.

            ``carrying_capacity``
                Total number of individuals (of any lineage) that can be
                supported by this habitat.

            ``fitness_function``
                Function object that takes a Lineage object as an argument and
                returns a fitness score. The fitness score can be any floating
                point value. It will serve as a weight for multinomial
                selection.
        """

        self.island = island
        self.carrying_capacity = carrying_capacity

        self._fitness_function = None
        self._cached_lineage_fitnesses = {}
        self.fitness_function = fitness_function

        # How many individuals of each lineage are present in this habitat.
        #   - keys = Lineage
        #   - values = number of individuals (int)
        self.populations = collections.defaultdict(int)

        # Rate of dispersal to other habitats:
        #   - keys = destination (Habitat object)
        #   - values = rate (float)
        self.dispersal_rates = collections.defaultdict(float)

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
    pass

def main():
    pass

if __name__ == "__main__":
    main()



