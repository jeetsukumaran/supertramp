#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import csv
import dendropy

from supertramp import summarize

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "source_paths",
            nargs="+",
            help="Path(s) to tree files.")
    parser.add_argument(
            "-o", "--output-root-dir",
            default='processed',
            help="Output directory root (default: '%(default)s').")
    parser.add_argument(
            "-x", "--exclude-first-island-as-continental-source-outside-of-analysis",
            action="store_true",
            default=False,
            help="Treat Island 0 as a continental source, and exclude it from analysis.")
    args = parser.parse_args()
    args.quiet = False
    args.group_processed_trees_by_model = False

    tree_processor = summarize.TreeProcessor()
    tree_processor.exclude_first_island_as_continental_source_outside_of_analysis = args.exclude_first_island_as_continental_source_outside_of_analysis
    summaries = []
    stats_fields = set()
    output_root_dir = args.output_root_dir
    output_dir = output_root_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    try:
        for source_idx, source_path in enumerate(args.source_paths):
            if not args.quiet:
                sys.stderr.write("Processing file {} of {}: {}\n".format(source_idx+1, len(args.source_paths), source_path))
            params = {}
            trees = dendropy.TreeList.get_from_path(source_path, "newick")
            for tree in trees:
                tree.treefile = source_path
                summary_stat, sub_stats_fields = tree_processor.process_trees(
                        trees,
                        params=params,
                        summaries=summaries)
                stats_fields.update(sub_stats_fields)
                for color_scheme in ("by-island", "by-habitat"):
                    colorized_trees_filepath = os.path.join(output_dir,
                            "processed.{}.trees".format(
                            color_scheme))
                    with open(colorized_trees_filepath, "w") as trees_outf:
                        tree_processor.write_colorized_trees(trees_outf, trees, color_scheme)
    except KeyboardInterrupt:
        pass

    stats_fields = sorted(list(stats_fields))
    summary_stats_fpath = os.path.join(output_dir, "summary.txt")

    with open(summary_stats_fpath, "w") as summary_outf:
        writer = csv.DictWriter(summary_outf,
                fieldnames=stats_fields,
                restval="NA",
                delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

if __name__ == "__main__":
    main()

