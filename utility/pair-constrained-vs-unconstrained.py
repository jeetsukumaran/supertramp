#! /usr/bin/env python

import re
import os
import sys
import argparse
import collections
import shutil

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "source_paths",
            nargs="+",
            help="Path(s) to directory or directories with simulation generated files.")
    args = parser.parse_args()

    sim_file_pattern = re.compile(r"^(un)?constrained_(.*)")
    for source_path in args.source_paths:
        root_dir = os.path.expanduser(os.path.expandvars(source_path))
        dest_root = source_path + "_paired"
        file_sets = collections.defaultdict(set)
        for root, dirs, filenames in os.walk(root_dir):
            for filename in filenames:
                if sim_file_pattern.match(filename):
                    parts = sim_file_pattern.split(filename)
                    title = parts[2].rsplit(".", 1)[0]
                    file_sets[title].add((source_path, filename))
        if not file_sets:
            continue
        for title in file_sets:
            file_set = file_sets[title]
            dest_dir = os.path.join(dest_root, title)
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)
            for source_path, filename in file_set:
                shutil.copy2(os.path.join(source_path, filename),
                        os.path.join(dest_dir, filename))


if __name__ == "__main__":
    main()

