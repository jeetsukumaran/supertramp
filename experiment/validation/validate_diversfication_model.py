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
from supertramp import monitor
from supertramp import utility
from supertramp import simulate

class DiversificationSubmodelValidator(object):

    def __init__(self,
            nreps=1,
            ngens=1e6,
            sample_frequency=1e3,
            random_seed=None,
            auto_delete_output_files=True):
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
        self.test_logger.info("||SUPERTRAMP-TEST|| Initializing with random seed: {}".format(random_seed))
        self.rng = random.Random(self.random_seed)
        self.auto_delete_output_files = auto_delete_output_files
        self.ngens = int(ngens)
        self.sample_frequency = int(sample_frequency)

        # create sampler
        self.simulator_monitor = monitor.SimulatorMonitor()

        # params
        self.simulator_monitor.add_attribute_tracker(
                attr_name="current_gen",
                field_name="gen",
                sample_diffs=False,
                )
        self.simulator_monitor.add_attribute_tracker(
                attr_name="diversification_model_s0",
                field_name="s0",
                sample_diffs=False,
                )
        self.simulator_monitor.add_attribute_tracker(
                attr_name="diversification_model_e0",
                field_name="e0",
                sample_diffs=False,
                )

        # data
        self.simulator_monitor.add_attribute_tracker(
                attr_name="num_extant_lineages",
                field_name="num_extant_lineages",
                sample_diffs=True,
                )
        self.simulator_monitor.add_attribute_tracker(
                attr_name="num_births",
                field_name="num_births",
                sample_diffs=True,
                )
        self.simulator_monitor.add_attribute_tracker(
                attr_name="num_extinctions",
                field_name="num_extinctions",
                sample_diffs=True,
                )
        self.simulator_monitor.add_attribute_tracker(
                attr_name="num_extirpations",
                field_name="num_extirpations",
                sample_diffs=True,
                )

    def get_model_params_dict(self):
        model_params_d = {}
        model_params_d["num_islands"] = 1
        model_params_d["num_habitat_types"] = 1
        model_params_d["diversification_model_a"] = -0.5
        model_params_d["diversification_model_b"] = 0.5
        model_params_d["dispersal_model"] = "unconstrained"
        model_params_d["dispersal_rate"] = 1e-3
        model_params_d["niche_evolution_probability"] = 0.0
        return model_params_d

    def run(self):
        self.sim_log_stream = tempfile.NamedTemporaryFile(
                mode="w",
                dir=os.curdir,
                prefix=self.__class__.__name__ + ".simulation.",
                suffix=".log",
                delete=self.auto_delete_output_files,
                )
        self.test_logger.info("||SUPERTRAMP-TEST|| Simulation log: '{}'".format( self.sim_log_stream.name,))
        self.sim_logger = utility.RunLogger(
                name="supertramp-run",
                log_to_stderr=True,
                stderr_logging_level="info",
                log_to_file=True,
                log_stream=self.sim_log_stream,
                file_logging_level="debug")
        # s0e0_values = (1e-8, 1e-6, 1e-4, 1e-2)
        # for s0e0_idx, (s0, e0) in enumerate(itertools.product(s0e0_values, s0e0_values)):
        s0e0_values = (
                (0.01, 0.0001),
                )
        for s0e0_idx, (s0, e0) in enumerate(s0e0_values):
            output_file_tag = "_s0={}_e0={}_".format(s0, e0)
            output_files_d = {}
            output_files_d["tree_log"] = tempfile.NamedTemporaryFile(
                    mode="w",
                    dir=os.curdir,
                    prefix=self.__class__.__name__ + output_file_tag,
                    suffix=".trees",
                    delete=self.auto_delete_output_files)
            output_files_d["general_stats_log"] = tempfile.NamedTemporaryFile(
                    mode="w",
                    dir=os.curdir,
                    prefix=self.__class__.__name__ + output_file_tag,
                    suffix=".stats.txt",
                    delete=self.auto_delete_output_files)
            self.test_logger.info("||SUPERTRAMP-TEST|| Tree log: '{}'".format(output_files_d["tree_log"].name))
            self.test_logger.info("||SUPERTRAMP-TEST|| General stats log: '{}'".format(output_files_d["general_stats_log"].name))
            current_rep = 0
            while current_rep < self.nreps:
                current_rep += 1
                num_attempts = 0
                while True:
                    num_attempts += 1
                    if num_attempts > 100:
                        self.test_logger.info("||SUPERTRAMP-TEST|| Maximum number of retries exceeded: abandoning replicate completely")
                        break
                    try:
                        self.test_logger.info("||SUPERTRAMP-TEST|| Starting replicate {rep} of {nreps} for diversification regime: s0={s0}, e0={e0}".format(
                            rep=current_rep,
                            nreps=self.nreps,
                            s0=s0,
                            e0=e0))
                        model_params_d = self.get_model_params_dict()
                        model_params_d["diversification_model_s0"] = s0
                        model_params_d["diversification_model_e0"] = e0
                        configd = {}
                        configd.update(model_params_d)
                        configd.update(output_files_d)
                        configd["run_logger"] = self.sim_logger
                        configd["rng"] = self.rng
                        configd["general_stats_log"].header_written = False
                        configd["log_frequency"] = 0
                        supertramp_simulator = simulate.SupertrampSimulator(
                                name="supertramp", **configd)
                        while supertramp_simulator.current_gen < self.ngens:
                            supertramp_simulator.run(self.sample_frequency)
                            self.simulator_monitor.sample(supertramp_simulator)
                        supertramp_simulator.report()
                    except simulate.TotalExtinctionException as e:
                        self.test_logger.info("||SUPERTRAMP-TEST|| Replicate {} of {}: [t={}] total extinction of all lineages before termination condition: {}".format(current_rep, self.nreps, supertramp_simulator.current_gen, e))
                        self.test_logger.info("||SUPERTRAMP-TEST|| Replicate {} of {}: restarting".format(current_rep, self.nreps))
                    else:
                        self.test_logger.info("||SUPERTRAMP-TEST|| Replicate {} of {}: completed to termination condition of {} generations".format(current_rep, self.nreps, self.ngens))
                        break
            df = self.simulator_monitor.as_data_frame()
            df = df[(df.s0 == s0) & (df.e0 == e0)]
            print("s0={s0}, e0={e0}: "
                  " num_extant_linages={num_extant_lineages};"
                  " num_births={num_births};"
                  " num_extirpations={num_extirpations};"
                  " num_extinctions={num_extinctions};"
                  " birth_rate={birth_rate};"
                  " extirpation_rate={extirpation_rate};"
                  " extinction_rate={extinction_rate};"
                  "".format(
                s0=s0,
                e0=e0,
                num_extant_lineages=df["num_extant_lineages"].mean(),
                num_births=df["num_births"].mean(),
                num_extirpations=df["num_extirpations"].mean(),
                num_extinctions=df["num_extinctions"].mean(),
                birth_rate=df["rate_of_change_num_births"].mean(),
                extirpation_rate=df["rate_of_change_num_extirpations"].mean(),
                extinction_rate=df["rate_of_change_num_extinctions"].mean(),
                ))

    def analyze(self):
        pass
        # df = self.simulator_monitor.as_data_frame()
        # print(df.describe())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--random-seed",
            default=None,
            help="Seed for random number generator engine.")
    parser.add_argument("-n", "--nreps",
            type=int,
            default=10,
            help="Number of replicates (default = %(default)s).")
    parser.add_argument("-g", "--ngens",
            type=int,
            default=10000,
            help="Number of generations to run in each replicate (default = %(default)s).")
    parser.add_argument("-s", "--sample-frequency",
            type=int,
            default=100,
            help="Frequency (number of generations) of samples to be taken during each replicate run (default = %(default)s).")
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
    parser.add_argument("--preserve-run-output",
            action="store_true",
            default=False,
            help="Do not clean up tree files, logs, etc.")
    parser.add_argument("--debug-mode",
            action="store_true",
            default=False,
            help="Run in debugging mode.")
    args = parser.parse_args()

    validator = DiversificationSubmodelValidator(
            nreps=args.nreps,
            ngens=args.ngens,
            sample_frequency=args.sample_frequency,
            random_seed=args.random_seed,
            auto_delete_output_files=not args.preserve_run_output)
    validator.run()
    validator.analyze()

if __name__ == "__main__":
    main()
