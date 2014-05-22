#! /usr/bin/env python

import os
import sys
import random
import unittest
from supertramp import utility
from supertramp import simulate

class HackedIsland(simulate.Island):

    def __init__(self, *args, **kwargs):
        self.dispersal_records = []
        simulate.Island.__init__(self, *args, **kwargs)

    def receive_migrant(self,
            lineage,
            habitat_type,
            from_island,
            from_habitat,
            ):
        self.dispersal_records.append(
                (lineage, habitat_type, from_island, from_habitat)
                )
        simulate.Island.receive_migrant(
                self,
                lineage=lineage,
                habitat_type=habitat_type,
                from_island=from_island,
                from_habitat=from_habitat,
                )

class HackedSupertrampSimulator(simulate.SupertrampSimulator):
    island_type = HackedIsland

    def clear_dispersal_records(self):
        for island in self.islands:
            island.dispersal_records = []

    def compile_dispersal_records(self):
        all_records = []
        for island in self.islands:
            all_records.extend(island.dispersal_records)
        return all_records

class DispersalRatesTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.logger = utility.RunLogger(
                name="supertramp-run",
                log_to_stderr=False,
                stderr_logging_level="info",
                log_to_file=False,
                log_stream=None,
                file_logging_level="debug")
        cls.random_seed = random.randint(0, sys.maxsize)
        cls.rng = random.Random(cls.random_seed)
        cls.logger.info("Random seed: {}".format(cls.random_seed))

    def get_simulator(self,
            num_islands,
            dispersal_rate,
            dispersal_model="unconstrained",
            num_habitat_types=1,
            niche_evolution_probability=0.0
            ):
        cls = self.__class__
        simulator = HackedSupertrampSimulator(
                num_islands=num_islands,
                num_habitat_types=num_habitat_types,
                diversification_model_s0=0.0,
                diversification_model_e0=0.0,
                diversification_model_a=-0.5,
                diversification_model_b=0.5,
                dispersal_model=dispersal_model,
                dispersal_rate=dispersal_rate,
                niche_evolution_probability=niche_evolution_probability,
                run_logger=cls.logger,
                tree_log=open(os.devnull, "w"),
                general_stats_log=open(os.devnull, "w"),
                rng=cls.rng,
                track_extinct_lineages=False,
                log_frequency=0,
                report_frequency=0,
            )
        return simulator

    def test_basic(self):
        for dispersal_rate in (1e-6, 1e-4, 1e-2):
            simulator = self.get_simulator(
                    num_islands=2,
                    dispersal_rate=dispersal_rate)
            for ngen in range(1000):
                simulator.clear_dispersal_records()
                simulator.execute_life_cycle()
                all_records = simulator.compile_dispersal_records()
                for record in all_records:
                    print(record)

if __name__ == "__main__":
    unittest.main()

