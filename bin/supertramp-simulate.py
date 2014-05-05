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

import os
import sys
import argparse
import random
from supertramp import supertramp

def main():
    parser = argparse.ArgumentParser(description="Biogeographical simulator")

    run_options = parser.add_argument_group("Run Options")
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
    landscape_options = parser.add_argument_group("Landscape Options")
    landscape_options.add_argument("--num-islands",
            type=int,
            default=4,
            help="number of islands (default = %(default)s).")
    landscape_options.add_argument("--num-habitat-types",
            type=int,
            default=4,
            help="number of habitat types per island (default = %(default)s).")
    diversification_param_options = parser.add_argument_group("Diversification Model Parameters")
    diversification_param_options.add_argument("-a", "--diversification-model-a",
            type=float,
            default=-0.5,
            help="'a' parameter of the diversfication model (default: %(default)s).")
    diversification_param_options.add_argument("-b", "--diversification-model-b",
            type=float,
            default=0.5,
            help="'b' parameter of the diversfication model (default: %(default)s).")
    diversification_param_options.add_argument("-s", "--diversification-model-s0",
            type=float,
            default=0.40,
            help="'s' parameter of the diversfication model (default: %(default)s).")
    diversification_param_options.add_argument("-e", "--diversification-model-e0",
            type=float,
            default=0.01,
            help="'e' parameter of the diversfication model (default: %(default)s).")
    taxon_cycle_param_options = parser.add_argument_group("Taxon Cycle (Sub-)Model Parameters")
    taxon_cycle_param_options.add_argument("--dispersal_model",
            type=str,
            default="unconstrained",
            choices=["constrained", "unconstrained"],
            help="Dispersal model: constrained or unconstrained by habitat")
    taxon_cycle_param_options.add_argument("-d", "--dispersal-rate",
            default=0.01,
            type=float,
            help="Dispersal rate (default = %(default)s).")
    taxon_cycle_param_options.add_argument("-y", "--niche-evolution-probability",
            default=0.01,
            type=float,
            help="Lineage (post-splitting) niche evolution probability (default = %(default)s).")
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
    output_options.add_argument("-r", "--report-frequency",
            default=None,
            type=int,
            help="Frequency that data is sampled from the simulation (default = None [final report only]).")
    args = parser.parse_args()
    argsd = vars(args)

    _DEBUG_MODE = argsd.pop("debug_mode")
    nreps = argsd.pop("nreps")
    ngens = argsd.pop("ngens")
    random_seed = argsd.pop("random_seed", None)
    if random_seed is None:
        random_seed = random.randint(0, sys.maxsize)

    configd = dict(argsd)
    configd["run_logger"] = supertramp.RunLogger(name="supertramp",
            log_path=args.output_prefix + ".log")
    configd["run_logger"].info("Initializing with random seed: {}".format(random_seed))
    configd["rng"] = random.Random(random_seed)
    configd["tree_log"] = open(configd["output_prefix"] + ".trees",
            "w")
    configd["general_stats_log"] = open(configd["output_prefix"] + ".general_stats.txt",
            "w")
    configd["general_stats_log"].header_written = False
    header_written = False
    rep = 0
    while rep < nreps:
        simulation_name="Run{}".format((rep+1))
        run_output_prefix = "{}.R{:04d}".format(configd["output_prefix"], rep+1)
        configd["run_logger"].info("Run {} of {}: starting".format(rep+1, nreps))
        supertramp_system = supertramp.System(
                name=simulation_name,
                **configd)
        supertramp_system.bootstrap()
        success = False
        while not success:
            try:
                success = supertramp_system.run(ngens)
            except supertramp.TotalExtinctionException:
                configd["run_logger"].info("Run {} of {}: [t={}] total extinction of all lineages before termination condition".format(rep+1, nreps, supertramp_system.current_gen))
                configd["run_logger"].info("Run {} of {}: restarting".format(rep+1, nreps))
                supertramp_system = supertramp.System(
                        name=simulation_name,
                        **configd)
                supertramp_system.bootstrap()
            else:
                configd["run_logger"].info("Run {} of {}: completed to termination condition of {} generations".format(rep+1, nreps, ngens))
                supertramp_system.report()
                break
        rep += 1

if __name__ == "__main__":
    main()
