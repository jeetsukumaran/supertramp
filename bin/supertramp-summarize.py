#! /usr/bin/env python

import sys
import os
import argparse
import json
import collections
import dendropy

if dendropy.__version__.startswith("4"):
    _get_taxa = lambda x: x.taxon_namespace
else:
    _get_taxa = lambda x: x.taxon_set

class TreeColorizer(object):

    def __init__(self):

        self.island_colors = collections.defaultdict(lambda: "#666666")
        self.island_colors.update({
                "0001" : "#ff0000",
                "0010" : "#00cc00",
                "0100" : "#0000ff",
                "1000" : "#cccc00",
                })
        self.habitat_colors = {}
        self.habitat_colors.update({
                "001" : "#ff00ff",
                "010" : "#00ffff",
                "100" : "#ff8800",
                })

    def colorize_trees(self, trees, outf=None):
        self.encode_taxa(_get_taxa(trees))
        for tree in trees:
            for nd in tree:
                if nd.taxon is None and nd.label is None:
                    continue
                if nd.taxon and nd.taxon.habitat_color is not None:
                    nd.annotations.add_new("!color", nd.taxon.habitat_color)
                elif nd.label is not None:
                    self.encode_labeled_item(nd, annotate_habitat_color=True)
        if outf is not None:
            try:
                trees.write_to_stream(outf, "nexus")
            except AttributeError:
                self.write_nexus(trees, outf)

    def encode_labeled_item(self,
            t,
            annotate_island_color=False,
            annotate_habitat_color=False):
        if annotate_island_color and annotate_habitat_color:
            raise TypeError("Cannot simultaneously annotate island and habitat color")
        label_parts = t.label.split(".")
        t.island_code = label_parts[1]
        t.island_color = self.island_colors[t.island_code]
        if annotate_island_color and t.island_color is not None:
            t.annotations.add_new("!color", t.island_color)
        t.habitat_code = label_parts[2]
        t.habitat_color = self.habitat_colors.get(t.habitat_code, None)
        if annotate_habitat_color and t.habitat_color is not None:
            t.annotations.add_new("!color", t.habitat_color)

    def encode_taxa(self, taxa):
        for t in taxa:
            self.encode_labeled_item(t, annotate_island_color=True)

    def write_nexus(self, trees, outf):
        parts = []
        parts.append("#NEXUS\n")
        parts.append("BEGIN TAXA;")
        parts.append("    DIMENSIONS NTAX={};".format(len(_get_taxa(trees))))
        parts.append("    TAXLABELS")
        for t in _get_taxa(trees):
            if t.island_color is None:
                color = ""
            else:
                color = "[&!color={}]".format(t.island_color)
            parts.append("        '{}'{}".format(t.label, color))
        parts.append("    ;")
        parts.append("END;")
        parts.append("BEGIN TREES;")
        for tree_idx, tree in enumerate(trees):
            parts.append("    TREE T{} = [&R] {};".format(tree_idx, tree._as_newick_string()))
        parts.append("END;")
        outf.write("\n".join(parts))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "source_paths",
            nargs="+",
            help="Path to run output directories.")
    parser.add_argument('-o', '--output-prefix',
        action='store',
        dest='output_prefix',
        type=str,
        default='supertramp.summary',
        metavar='OUTPUT-FILE-PREFIX',
        help="Prefix for output files (default='%(default)s').")
    args = parser.parse_args()
    args.quiet = False

    tree_colorizer = TreeColorizer()
    for source_path in args.source_paths:
        source_dir = os.path.abspath(os.path.expanduser(os.path.expandvars(source_path)))
        run_manifest_path = os.path.join(source_dir, "run-manifest.json")
        if not os.path.exists(run_manifest_path):
            sys.exit("Manifest file not found: {}".format(run_manifest_path))
        with open(run_manifest_path, "r") as run_manifest_f:
            run_manifest = json.load(run_manifest_f)
        jobs = list(run_manifest.keys())
        for key_idx, key in enumerate(jobs):
            if not args.quiet:
                sys.stderr.write("Processing job {} of {}: {}\n".format(key_idx+1, len(jobs), key))
            run_data = run_manifest[key]
            tree_filepath = os.path.join(source_dir, run_data["treefile"])
            trees = dendropy.TreeList.get_from_path(
                    tree_filepath,
                    "newick")
            colorized_trees_filepath = args.output_prefix + ".{}.trees".format(key)
            with open(colorized_trees_filepath, "w") as trees_outf:
                tree_colorizer.colorize_trees(trees, trees_outf)
            # for tree_idx, tree in enumerate(trees):
            #     if not args.quiet:
            #         sys.stderr.write("Processing tree {} of {} in job {} of {}\n".format(tree_idx+1, len(trees), key_idx+1, len(jobs), key))


if __name__ == "__main__":
    main()
