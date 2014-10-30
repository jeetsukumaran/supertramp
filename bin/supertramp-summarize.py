#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import csv
import dendropy

from supertramp import summarize

# library(adegenet)
# summary.df = read.table("processed/summary.txt", header=T)
# summary.df = na.omit(summary.df)
# groups = summary.df$dispersal.model
# cols.to.drop <- c(
#                   "dispersal.model",
#                   "birth.rate",
#                   "death.rate",
#                   "dispersal.rate",
#                   "niche.evolution.prob",
#                   "edges",
#                   "est.birth.rate",
#                   "length",
#                   "size"
#                   )
# predictors = summary.df[,!(names(summary.df) %in% cols.to.drop)]
# result = dapc(predictors, groups)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "source_paths",
            nargs="+",
            help="Path(s) to directory or directories with simulation generated files.")
    parser.add_argument(
            "-o", "--output-root-dir",
            default='processed',
            help="Output directory root (default: '%(default)s').")
    args = parser.parse_args()
    args.quiet = False
    args.group_processed_trees_by_model = False

    tree_processor = summarize.TreeProcessor()
    param_keys = collections.OrderedDict()
    param_keys["dispersal.model"]      = "dispersal_model"
    param_keys["birth.rate"]           = "birth_rate"
    param_keys["death.rate"]           = "death_rate"
    param_keys["dispersal.rate"]       = "dispersal_rate"
    param_keys["niche.evolution.prob"] = "niche_evolution_prob"

    summaries = []
    output_root_dir = args.output_root_dir
    output_dir = output_root_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    try:
        for source_path in args.source_paths:
            source_dir = os.path.abspath(os.path.expanduser(os.path.expandvars(source_path)))
            run_manifest_path = os.path.join(source_dir, "run-manifest.json")
            if not os.path.exists(run_manifest_path):
                sys.exit("Manifest file not found: {}".format(run_manifest_path))
            with open(run_manifest_path, "r") as run_manifest_f:
                run_manifest = json.load(run_manifest_f)
            jobs = list(run_manifest.keys())
            stats_fields = set()
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
                trees = dendropy.TreeList.get_from_path(
                        tree_filepath,
                        "newick")
                for tree in trees:
                    tree.treefile = tree_filepath
                summary_stat, sub_stats_fields = tree_processor.process_trees(
                        trees,
                        params=params,
                        summaries=summaries)
                stats_fields.update(sub_stats_fields)
                for color_scheme in ("by-island", "by-habitat"):
                    if args.group_processed_trees_by_model:
                        out_fname = job
                    else:
                        parts = job.rsplit("_", 1)
                        assert len(parts) == 2
                        model_name = parts[-1]
                        out_fname = parts[0]
                    colorized_trees_filepath = os.path.join(output_dir, "{}.processed.{}.{}.trees".format(out_fname, color_scheme, model_name))
                    with open(colorized_trees_filepath, "w") as trees_outf:
                        tree_processor.write_colorized_trees(trees_outf, trees, color_scheme)
    except KeyboardInterrupt:
        pass

    param_fields = list(param_keys.keys())
    stats_fields = sorted(list(stats_fields))
    all_fields = param_fields + stats_fields
    summary_stats_fpath = os.path.join(output_dir, "summary.txt")

    with open(summary_stats_fpath, "w") as summary_outf:
        writer = csv.DictWriter(summary_outf,
                fieldnames=all_fields,
                restval="NA",
                delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

if __name__ == "__main__":
    main()
