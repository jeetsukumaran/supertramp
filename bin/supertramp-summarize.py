#! /usr/bin/env python

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "source_path",
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


if __name__ == "__main__":
    main()
