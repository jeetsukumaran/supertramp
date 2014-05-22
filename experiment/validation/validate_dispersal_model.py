#! /usr/bin/env python

import os
import sys
import random
import unittest
from supertramp import utility
from supertramp import simulate
from supertramp import monitor

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

    def __init__(self, *args, **kwargs):
        simulate.SupertrampSimulator.__init__(self, *args, **kwargs)

    def clear_dispersal_records(self):
        for island in self.islands:
            island.dispersal_records = []

    def _get_total_dispersals(self):
        n = 0
        for island in self.islands:
            n += len(island.dispersal_records)
        return n
    total_dispersals = property(_get_total_dispersals)

class RateTracker(object):

    def sample(self, simulator, record):
        count = record["delta_observed_total_disperals"]
        period = record["delta_gen"]
        rate = count/period
        record["observed_dispersal_rate"] = rate


class DispersalRatesValidator():

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

    def run(self):
        for dispersal_rate in (1e-6, 1e-4, 1e-2, 1e-1):
            simulator_monitor = monitor.SimulatorMonitor()
            simulator_monitor.add_attribute_tracker(
                    attr_name="current_gen",
                    field_name="gen",
                    sample_diffs=True,
                    )
            simulator_monitor.add_attribute_tracker(
                    attr_name="num_islands",
                    field_name="num_islands",
                    sample_diffs=False,
                    )
            simulator_monitor.add_attribute_tracker(
                    attr_name="global_dispersal_rate",
                    field_name="global_dispersal_rate",
                    sample_diffs=False,
                    )
            simulator_monitor.add_attribute_tracker(
                    attr_name="total_dispersals",
                    field_name="observed_total_disperals",
                    sample_diffs=True,
                    )
            simulator_monitor.trackers.append(RateTracker())
            # for dispersal_rate in (1e-6, 1e-4, 1e-2):
            for num_islands in (2,):
                simulator = self.get_simulator(
                        num_islands=num_islands,
                        dispersal_rate=dispersal_rate)
                for ngen in range(100):
                    simulator.run(100)
                    simulator_monitor.sample(simulator)
            df = simulator_monitor.as_data_frame()
            print("{: 10.4e}: {:6.4e}".format(dispersal_rate, df["observed_dispersal_rate"].mean()))

if __name__ == "__main__":
    d = DispersalRatesValidator()
    d.setUpClass()
    d.run()
