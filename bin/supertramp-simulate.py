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
import logging
import supertramp
from supertramp import simulate
from supertramp import utility

def main():
    simulation_model_arg_parser = simulate.SupertrampSimulator.simulation_model_arg_parser()
    parser = argparse.ArgumentParser(
            parents=[simulation_model_arg_parser],
            description="{} Biogeographical simulator".format(supertramp.description())
            )

    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument("-z", "--random-seed",
            default=None,
            help="Seed for random number generator engine.")
    run_options.add_argument("-n", "--nreps",
            type=int,
            default=10,
            help="number of replicates (default = %(default)s).")
    run_options.add_argument("--log-frequency",
            default=100,
            type=int,
            help="Frequency that background progress messages get written to the log (default = %(default)s).")
    run_options.add_argument("--file-logging-level",
            default="debug",
            help="Message level threshold for file logs.")
    run_options.add_argument("--stderr-logging-level",
            default="info",
            help="Message level threshold for screen logs.")
    run_options.add_argument("--debug-mode",
            action="store_true",
            default=False,
            help="Run in debugging mode.")

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

    termination_condition_options = parser.add_argument_group("Termination Condition Options")
    termination_condition_options.add_argument("-g", "--ngens",
            type=int,
            default=10000,
            help="Number of generations to run (default = %(default)s).")

    args = parser.parse_args()

    argsd = vars(args)
    ngens = argsd.pop("ngens")
    nreps = argsd.pop("nreps")
    random_seed = argsd.pop("random_seed", None)
    output_prefix = argsd.pop("output_prefix", "supertramp_run")
    stderr_logging_level=argsd.pop("stderr_logging_level")
    file_logging_level=argsd.pop("file_logging_level")

    simulate.repeat_run_supertramp(
            model_params_d=argsd,
            ngens=ngens,
            nreps=nreps,
            output_prefix=output_prefix,
            random_seed=random_seed,
            stderr_logging_level=stderr_logging_level,
            file_logging_level=file_logging_level)

if __name__ == "__main__":
    main()

