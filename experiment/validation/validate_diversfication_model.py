#! /usr/bin/env python

import os
import sys
import argparse
import itertools
import tempfile
import random
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from supertramp import utility
from supertramp import simulate

class DiversificationSubmodelValidator(object):

    def __init__(self,
            nreps=1,
            random_seed=None,
            auto_delete_simulation_log=True):
        self.test_logger = utility.RunLogger(
                name=self.__class__.__name__ + ".test",
                log_to_stderr=True,
                stderr_logging_level="info",
                log_to_file=True,
                file_logging_level="debug",
                )
        self.nreps = nreps
        self.random_seed = random_seed
        if self.random_seed is None:
            random_seed = random.randint(0, sys.maxsize)
        self.test_logger.info("Initializing with random seed: {}".format(random_seed))
        self.rng = random.Random(self.random_seed)
        self.delete_simulation_log_file = auto_delete_simulation_log

    def run(self):
        self.sim_log_stream = tempfile.NamedTemporaryFile(
                mode="w",
                dir=os.curdir,
                prefix=self.__class__.__name__ + ".simulation.",
                suffix=".log",
                delete=self.delete_simulation_log_file,
                )
        self.test_logger.info("Simulation log will be saved to: '{}' ({} on successful exit)".format(
            self.sim_log_stream.name,
            "auto-deleted" if self.delete_simulation_log_file else "not deleted",
            ))
        self.sim_logger = utility.RunLogger(
                name="supertramp-run",
                log_to_stderr=True,
                stderr_logging_level="info",
                log_to_file=True,
                log_stream=self.sim_log_stream,
                file_logging_level="debug",
                )
        s0e0_values = (1e-8, 1e-6, 1e-4, 1e-2)
        for s0e0_idx, (s0, e0) in enumerate(itertools.product(s0e0_values, s0e0_values)):
            for rep in range(self.nreps):
                self.test_logger.info("Starting replicate {rep} of {nreps} for diversification regime: s0={s0}, e0={e0}".format(
                    rep=rep+1,
                    nreps=self.nreps,
                    s0=s0,
                    e0=e0))

                model_params_d = self.get_model_params_dict()
                model_params_d["diversification_model_s0"] = s0
                model_params_d["diversification_model_e0"] = e0

                configd = {}
                configd.update(model_params_d)
                configd["run_logger"] = self.sim_logger
                configd["rng"] = self.rng
                configd["tree_log"] = tempfile.NamedTemporaryFile(delete=True)
                self.test_logger.info("Tree log: '{}'".format(configd["tree_log"].name))
                configd["general_stats_log"] = tempfile.NamedTemporaryFile(delete=True)
                self.test_logger.info("General stats log: '{}'".format(configd["general_stats_log"].name))
                configd["general_stats_log"].header_written = False

                supertramp_simulator = simulate.SupertrampSimulator(
                        name="supertramp", **configd)
                supertramp_simulator.run(1000)




    def get_model_params_dict(self):
        model_params_d = {}
        model_params_d["num_islands"] = 1
        model_params_d["num_habitat_types"] = 1
        model_params_d["diversification_model_a"] = -0.5
        model_params_d["diversification_model_b"] = 0.5
        model_params_d["dispersal_model"] = "unconstrained"
        model_params_d["dispersal_rate"] = 0.0
        model_params_d["niche_evolution_probability"] = 0.0
        return model_params_d

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--random-seed",
            default=None,
            help="Seed for random number generator engine.")
    parser.add_argument("-n", "--nreps",
            type=int,
            default=10,
            help="number of replicates (default = %(default)s).")
    parser.add_argument("--log-frequency",
            default=100,
            type=int,
            help="Frequency that background progress messages get written to the log (default = %(default)s).")
    parser.add_argument("--file-logging-level",
            default="debug",
            help="Message level threshold for file logs.")
    parser.add_argument("--screen-logging-level",
            default="info",
            help="Message level threshold for screen logs.")
    parser.add_argument("--debug-mode",
            action="store_true",
            default=False,
            help="Run in debugging mode.")
    args = parser.parse_args()

    validator = DiversificationSubmodelValidator(
            nreps=args.nreps,
            random_seed=args.random_seed)
    validator.run()


if __name__ == "__main__":
    main()
