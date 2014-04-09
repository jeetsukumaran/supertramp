#! /usr/bin/env python

import sys
import os
import json
import random
import argparse

kwyjibo_job_template = """\
#! /bin/bash
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -l h_vmem=8G
#$ -l virtual_free=8G
{commands}
"""

def main():
    """
    Assumptions
    -----------

    1 simulation generation = 100 years
    10000 simulation generations = 1e6 years
    Simulation run-time:        1000000 generations    = 1e8 years
    High speciation rate:       0.001   per generation = 0.1 per MY
    Med speciation rate:        0.0001  per generation = 0.01 per MY
    Low speciation rate:        0.00001 per generation = 0.001 per MY
    Dispersal rate:             0.01, 0.5, 1.0, 2.0, 10.0 x speciation rates
    Niche evolution prob:       0.001, 0.01, 0.10, 1.0
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--venv",
            default=None,
            help="Path to Python virtual environment.")
    parser.add_argument("-z", "--random-seed",
            default=None,
            help="Seed for random number generator engine.")
    parser.add_argument("--ngens",
            type=int,
            default=1000000,
            help="Number of generations to run (default = %(default)s).")
    parser.add_argument("--nreps",
            type=int,
            default=10,
            help="number of replicates (default = %(default)s).")
    args = parser.parse_args()

    if args.random_seed is None:
        args.random_seed = random.randint(0, sys.maxsize)
    rng = random.Random(args.random_seed)
    if args.venv is not None:
        venv_dir = os.path.expanduser(os.path.expandvars(args.venv))
        venv_activate = os.path.abspath(os.path.join(venv_dir, "bin", "activate"))
        if not os.path.exists(venv_activate):
            raise Exception("Virtual environment activation script not found: '{}'".format(venv_activate))
        source_venv = "source {}".format(venv_activate)
    else:
        source_venv = ""
    # python_path = "python3"
    # supertramp_path = os.path.abspath(os.path.join(
    #         os.path.dirname(__file__),
    #         "supertramp-simulate.py"))
    supertramp_path = "supertramp-simulate.py"

    dispersal_models = ["constrained", "unconstrained"]
    # birth_rates = [0.001, 0.0001, 0.00001]
    # dispersal_rate_factors = [0.01, 0.5, 1.0, 2.0, 10.0]
    # niche_evolution_probs = [0.001, 0.01, 0.10, 1.0]
    birth_rates = [0.01, 0.001, 0.0001]
    dispersal_rate_factors = [0.1, 1.0, 10.0]
    niche_evolution_probs = [0.01, 0.10, 1.0]
    run_manifest = {}
    for dm_idx, dispersal_model in enumerate(dispersal_models):
        for br_idx, birth_rate in enumerate(birth_rates):
            for drf_idx, dispersal_rate_factor in enumerate(dispersal_rate_factors):
                for nef_idx, niche_evolution_prob in enumerate(niche_evolution_probs):
                    dispersal_rate = dispersal_rate_factor * birth_rate
                    stem = "{dispersal_model}_b{birth_rate:6.5f}_r{dispersal_rate:6.5f}_n{niche_evolution_prob:6.5f}".format(
                            dispersal_model=dispersal_model,
                            birth_rate=birth_rate,
                            dispersal_rate=dispersal_rate,
                            niche_evolution_prob=niche_evolution_prob)
                    output_prefix = stem + ".run"
                    run_cmd = []
                    run_cmd.append(supertramp_path)
                    run_cmd.extend(["-z", str(rng.randint(0, sys.maxsize))])
                    run_cmd.extend(["--nreps", str(args.nreps)])
                    run_cmd.extend(["--log-frequency", "1000"])
                    run_cmd.extend(["--birth-probability", str(birth_rate)])
                    run_cmd.extend(["--death-probability", "0"])
                    run_cmd.extend(["--niche-evolution-probability", str(niche_evolution_prob)])
                    run_cmd.extend(["--dispersal-rate", str(dispersal_rate)])
                    run_cmd.extend(["--ngens", str(args.ngens)])
                    run_cmd.extend(["--output-prefix", output_prefix])
                    run_cmd.append(dispersal_model)
                    run_cmd = " ".join(run_cmd)
                    commands = []
                    if source_venv:
                        commands.append(source_venv)
                    commands.append(run_cmd)
                    job_filepath = stem + ".job"
                    with open(job_filepath, "w") as jobf:
                        template = kwyjibo_job_template
                        jobf.write(template.format(commands="\n".join(commands)))
                    run_manifest[output_prefix] = {
                            "dispersal_model"       : dispersal_model,
                            "birth_rate"            : birth_rate,
                            "dispersal_rate_factor" : dispersal_rate_factor,
                            "dispersal_rate"        : dispersal_rate,
                            "niche_evolution_prob"  : niche_evolution_prob,
                            "treefile"              : output_prefix + ".trees",
                            "logfile"              : output_prefix + ".log",
                            }
    with open("run-manifest.json", "w") as manifestf:
        json.dump(run_manifest, manifestf)
if __name__ == "__main__":
    main()




