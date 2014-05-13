#! /usr/bin/env python

import os
import sys
import fnmatch
import json
import argparse
import collections
import dendropy

def get_species_index_max_and_min(trees):
    species_idxs = set()
    for taxon in trees.taxon_namespace:
        label = taxon.label
        species_idx = label.split(".")[0][1:]
        species_idx = int(species_idx)
        species_idxs.add(species_idx)
    return min(species_idxs), max(species_idxs)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "source_paths",
            nargs="+",
            help="Path(s) to directory or directories with simulation generated files.")
    parser.add_argument(
            "-u", "--unmanaged-runs",
            action="store_true",
            default=False,
            help="Not a managed simulation (no 'run-manifest.json' available);"
            " simply summarize results of all tree files found in source path(s).")
    args = parser.parse_args()
    args.quiet = False

    param_keys = collections.OrderedDict()
    param_keys["disp.model"]      = "dispersal_model"
    param_keys["s0"]                   = "s0"
    param_keys["e0"]                   = "e0"
    param_keys["disp.rate"]       = "dispersal_rate"
    param_keys["niche.shift.prob"] = "niche_evolution_prob"

    summaries = []
    out = sys.stdout
    float_template = "{:1.1e}"
    if args.unmanaged_runs:
        tree_filepaths = []
        pattern = "*.trees"
        for source_path in args.source_paths:
            root_dir = os.path.expanduser(os.path.expandvars(source_path))
            for root, dirs, files in os.walk(root_dir):
                for filename in fnmatch.filter(files, pattern):
                    tree_filepaths.append(os.path.join(root, filename))
        if not tree_filepaths:
            sys.exit("No tree files found")
        col_width = max(len(tree_filepath) for tree_filepath in tree_filepaths)
        text_template = "{{:{}}}".format(col_width)
        row_template = text_template + ":   " + float_template
        for tree_filepath in tree_filepaths:
            trees = dendropy.TreeList.get_from_path(
                    tree_filepath,
                    "newick")
            value = get_species_index_max_and_min(trees)[1]
            out.write(row_template.format(tree_filepath, value))
            out.write("\n")
    else:
        for source_path in args.source_paths:
            source_dir = os.path.abspath(os.path.expanduser(os.path.expandvars(source_path)))
            run_manifest_path = os.path.join(source_dir, "run-manifest.json")
            if not os.path.exists(run_manifest_path):
                sys.exit("Manifest file not found: {}".format(run_manifest_path))
            with open(run_manifest_path, "r") as run_manifest_f:
                run_manifest = json.load(run_manifest_f)
            jobs = list(run_manifest.keys())
            header = list(param_keys.keys()) + ["num.lin"]
            col_width = max(len(h) for h in header)
            text_template = "{{:{}}}".format(col_width)
            header = [text_template.format(h) for h in header]
            out.write(" | ".join(header))
            out.write("\n")
            for job_idx, job in enumerate(jobs):
                params = {}
                for param_key in param_keys:
                    params[param_key] = run_manifest[job][param_keys[param_key]]
                run_data = run_manifest[job]
                tree_filepath = os.path.join(source_dir, run_data["treefile"])
                if not os.path.exists(tree_filepath):
                    sys.stderr.write("Skipping job {} of {} (missing): {}\n".format(job_idx+1, len(jobs), job))
                    continue
                if not args.quiet:
                    sys.stderr.write("Processing job {} of {}: {}\n".format(job_idx+1, len(jobs), job))
                values = [params["s0"], params["e0"], params["disp.rate"], params["niche.shift.prob"]]
                trees = dendropy.TreeList.get_from_path(
                        tree_filepath,
                        "newick")
                values.append(get_species_index_max_and_min(trees)[1])
                values = [float_template.format(v) for v in values]
                values.insert(0, params["disp.model"])
                values = [text_template.format(v) for v in values]
                out.write(" | ".join(values))
                out.write("\n")

if __name__ == "__main__":
    main()

