#! /usr/bin/env python


    # def run_lineage_birth_0(self):
    #     # - Each lineage has a probability of splitting given by
    #     #   the sum of the local birth rates for the lineage across all
    #     #   habitats in which it occurs.
    #     # - The probability that a particular lineage in a particular habitat
    #     #   speciates, given that the lineage has split in that generation,
    #     #   is given by the local (habitat-specific) birth rate normalized by
    #     #   the sum of birth rates across all habitats.
    #     # - In any particular generation, a particular lineage splits at most
    #     #   once.
    #     lineage_splitting_rates = collections.defaultdict(list)
    #     lineage_splitting_habitat_localities = collections.defaultdict(list)
    #     for island in self.islands:
    #         for habitat in island.habitat_list:
    #             if not habitat.lineages:
    #                 continue
    #             splitting_rate = self.diversification_model_s0 * (len(habitat.lineages) ** self.diversification_model_a)
    #             for lineage in habitat.lineages:
    #                 lineage_splitting_rates[lineage].append(splitting_rate)
    #                 lineage_splitting_habitat_localities[lineage].append(habitat)
    #     for lineage in lineage_splitting_rates:
    #         total_rate = sum(lineage_splitting_rates[lineage])
    #         if self.rng.uniform(0, 1) <= total_rate:
    #             # print(">>> {}: splitting: {}:{}".format(self.current_gen, lineage, lineage.label))
    #             c0, c1 = lineage.diversify(finalize_distribution_label=True)
    #             if _DEBUG_MODE:
    #                 try:
    #                     self.phylogeny._debug_check_tree()
    #                 except AttributeError:
    #                     self.phylogeny.debug_check_tree()
    #             if len(self.habitat_types) > 1 and self.rng.uniform(0, 1) <= self.global_lineage_niche_evolution_probability:
    #                 c1.habitat_type = self.rng.choice([ h for h in self.habitat_types if h is not c1.habitat_type ])
    #             selected_habitat_idx = weighted_index_choice(lineage_splitting_rates[lineage],
    #                     self.rng)
    #             # selected_habitat = lineage_split_habitat_localities[selected_habitat_idx]
    #             for habitat_idx, habitat in enumerate(lineage_splitting_habitat_localities[lineage]):
    #                 habitat.remove_lineage(lineage)
    #                 if habitat_idx == selected_habitat_idx:
    #                     # sympatric speciation: "old" species retained in original habitat on island
    #                     # new species added to new habitat on island
    #                     habitat.island.add_lineage(lineage=c0, habitat_type=c0.habitat_type)
    #                     habitat.island.add_lineage(lineage=c1, habitat_type=c1.habitat_type)
    #                 else:
    #                     habitat.island.add_lineage(lineage=c0, habitat_type=c0.habitat_type)

    # def run_lineage_death_2(self):
    #     # - Each lineage in each habitat has a probability of going
    #     #   locally-extinct (i.e., extirpated from a particular habitat)
    #     #   based on its local (habitat-specific) death rate.
    #     # - A particular lineage's probability of going extinct in a particular
    #     #   habitat is independent of the any other lineage going extinct in
    #     #   the same or any other habitat.
    #     # - A lineage can only be extirpated from a single habitat
    #     #   in any given generation.
    #     # - Once a lineage is no longer present in any habitat across all
    #     #   islands, it is considered to have gone globally-extinct and
    #     #   removed from the system.
    #     lineage_death_rates = collections.defaultdict(list)
    #     lineage_death_habitat_localities = collections.defaultdict(list)
    #     lineage_counts = collections.Counter()
    #     for island in self.islands:
    #         for habitat in island.habitat_list:
    #             lineage_counts.update(habitat.lineages)
    #             if not habitat.lineages:
    #                 continue
    #             death_rate = self.diversification_model_e0 * (len(habitat.lineages) ** self.diversification_model_b)
    #             for lineage in habitat.lineages:
    #                 lineage_death_rates[lineage].append(death_rate)
    #                 lineage_death_habitat_localities[lineage].append(habitat)
    #     to_remove = []
    #     for lineage in lineage_death_rates:
    #         total_rate = sum(lineage_death_rates[lineage])
    #         if self.rng.uniform(0, 1) <= total_rate:
    #             # print(">>> {}: death: {}:{}".format(self.current_gen, lineage, lineage.label))
    #             selected_habitat_idx = weighted_index_choice(lineage_death_rates[lineage],
    #                     self.rng)
    #             to_remove.append( (lineage, lineage_death_habitat_localities[lineage][selected_habitat_idx]) )
    #     for lineage, habitat in to_remove:
    #         habitat.remove_lineage(lineage)
    #         lineage_counts.subtract([lineage])
    #     for lineage in lineage_counts:
    #         count = lineage_counts[lineage]
    #         if count == 0:
    #             if lineage is self.phylogeny.seed_node:
    #                 raise TotalExtinctionException()
    #             else:
    #                 self.phylogeny.prune_subtree(node=lineage,
    #                         update_splits=False, delete_outdegree_one=True)
    #                 if self.phylogeny.seed_node.num_child_nodes() == 0:
    #                     raise TotalExtinctionException()
    #         elif _DEBUG_MODE:
    #             ## sanity checking ...
    #             found = True
    #             for island in self.islands:
    #                 for habitat in island.habitat_list:
    #                     if lineage in habitat.lineages:
    #                         found = True
    #                         break
    #                 if found:
    #                     break
    #             assert found, lineage

    # def run_lineage_death_3(self):
    #     # - Each lineage has a probability of going globally-extinct given by
    #     #   the sum of the local death rates for the lineage across all
    #     #   habitats in which it occurs.
    #     # - A particular lineage's probability of going extinct in a
    #     #   particular habitat is independent of the any other lineage going
    #     #   extinct in the same or any other habitat.
    #     lineage_death_rates = collections.defaultdict(list)
    #     lineage_death_habitat_localities = collections.defaultdict(list)
    #     lineage_counts = collections.Counter()
    #     for island in self.islands:
    #         for habitat in island.habitat_list:
    #             lineage_counts.update(habitat.lineages)
    #             if not habitat.lineages:
    #                 continue
    #             death_rate = self.diversification_model_e0 * (len(habitat.lineages) ** self.diversification_model_b)
    #             for lineage in habitat.lineages:
    #                 lineage_death_rates[lineage].append(death_rate)
    #                 lineage_death_habitat_localities[lineage].append(habitat)
    #     for lineage in lineage_death_rates:
    #         total_rate = sum(lineage_death_rates[lineage])
    #         if self.rng.uniform(0, 1) <= total_rate:
    #             for island in self.islands:
    #                 for habitat in island.habitat_list:
    #                     if lineage in habitat.lineages:
    #                         habitat.remove_lineage(lineage)
    #             if lineage is self.phylogeny.seed_node:
    #                 raise TotalExtinctionException()
    #             else:
    #                 self.phylogeny.prune_subtree(node=lineage,
    #                         update_splits=False, delete_outdegree_one=True)
    #                 if self.phylogeny.seed_node.num_child_nodes() == 0:
    #                     raise TotalExtinctionException()

    # def run_simple_birth(self):
    #     tips = self.phylogeny.leaf_nodes()
    #     if not tips:
    #         raise TotalExtinctionException()
    #     if self.rng.uniform(0, 1) > self.global_per_lineage_birth_prob * len(tips):
    #         return
    #     lineage_habitat_localities = collections.defaultdict(list)
    #     for island in self.islands:
    #         for habitat in island.habitat_list:
    #             for lineage in habitat.lineages:
    #                 lineage_habitat_localities[lineage].append(habitat)
    #     splitting_lineage = self.rng.choice(list(lineage_habitat_localities.keys()))
    #     c0, c1 = splitting_lineage.diversify(finalize_distribution_label=True)
    #     if _DEBUG_MODE:
    #         try:
    #             self.phylogeny._debug_check_tree()
    #         except AttributeError:
    #             self.phylogeny.debug_check_tree()
    #     if self.rng.uniform(0, 1) <= self.global_lineage_niche_evolution_probability:
    #         c1.habitat_type = self.rng.choice([ h for h in self.habitat_types if h is not c1.habitat_type ])
    #     splitting_lineage_habitat_localities = lineage_habitat_localities[splitting_lineage]
    #     assert splitting_lineage_habitat_localities
    #     if len(splitting_lineage_habitat_localities) == 1:
    #         target = splitting_lineage_habitat_localities[0]
    #     else:
    #         target = self.rng.choice(splitting_lineage_habitat_localities)
    #     for habitat in splitting_lineage_habitat_localities:
    #         habitat.remove_lineage(splitting_lineage)
    #         if habitat is target:
    #             # sympatric speciation: "old" species retained in original habitat on island
    #             # new species added to new habitat on island
    #             habitat.island.add_lineage(lineage=c0, habitat_type=c0.habitat_type)
    #             habitat.island.add_lineage(lineage=c1, habitat_type=c1.habitat_type)
    #         else:
    #             habitat.island.add_lineage(lineage=c0, habitat_type=c0.habitat_type)

    # def run_density_dependent_birth(self):
    #     self.run_simple_birth() # TODO!!!

    # def run_simple_death(self):
    #     tips = self.phylogeny.leaf_nodes()
    #     if not tips:
    #         raise TotalExtinctionException()
    #     if self.rng.uniform(0, 1) > self.global_per_lineage_death_prob * len(tips):
    #         return
    #     extinguishing_lineage = self.rng.choice(tips)
    #     lineage_habitat_localities = []
    #     for island in self.islands:
    #         for habitat in island.habitat_list:
    #             if extinguishing_lineage in habitat.lineages:
    #                 lineage_habitat_localities.append(habitat)
    #     for loc in lineage_habitat_localities:
    #         loc.remove_lineage(extinguishing_lineage)
    #     if extinguishing_lineage is self.phylogeny.seed_node:
    #         raise TotalExtinctionException()
    #     else:
    #         self.phylogeny.prune_subtree(node=extinguishing_lineage,
    #                 update_splits=False, delete_outdegree_one=True)
    #         if self.phylogeny.seed_node.num_child_nodes() == 0:
    #             raise TotalExtinctionException()

    # def run_density_dependent_death(self):
    #     lineage_counts = collections.Counter()
    #     for island in self.islands:
    #         for habitat in island.habitat_list:
    #             lineage_counts.update(habitat.lineages)
    #             n = len(habitat.lineages)
    #             if n <= 0:
    #                 continue
    #             assert habitat.carrying_capacity is not None
    #             K = habitat.carrying_capacity
    #             if n > K:
    #                 while n > K:
    #                     lineage = self.rng.choice(list(habitat.lineages))
    #                     habitat.remove_lineage(lineage)
    #                     lineage_counts.subtract([lineage])
    #                     n -= 1
    #             else:
    #                 weight = 1.0 - ((K-n)//K)
    #                 prob = self.global_per_lineage_death_prob * weight
    #                 prob = self.diversification_model_s0 * (len(habitat.lineages) ** self.diversification_model_a)
    #                 # print("n={:2d}, K={:2d}, weight={:8.6f}, prob={:8.6f}".format(n, K, weight, prob))
    #                 if self.rng.uniform(0, 1) <= prob:
    #                     lineage = self.rng.choice(list(habitat.lineages))
    #                     habitat.remove_lineage(lineage)
    #                     lineage_counts.subtract([lineage])
    #     for lineage in lineage_counts:
    #         count = lineage_counts[lineage]
    #         if count == 0:
    #             if lineage is self.phylogeny.seed_node:
    #                 raise TotalExtinctionException()
    #             else:
    #                 self.phylogeny.prune_subtree(node=lineage,
    #                         update_splits=False, delete_outdegree_one=True)
    #                 if self.phylogeny.seed_node.num_child_nodes() == 0:
    #                     raise TotalExtinctionException()
    #         elif _DEBUG_MODE:
    #             ## sanity checking ...
    #             found = True
    #             for island in self.islands:
    #                 for habitat in island.habitat_list:
    #                     if lineage in habitat.lineages:
    #                         found = True
    #                         break
    #                 if found:
    #                     break
    #             assert found, lineage

