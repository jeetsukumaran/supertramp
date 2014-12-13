#! /usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import sys
import os
import re
import collections
import dendropy

from supertramp import summarize

class TraitEvolutionRateEstimator(object):

    def __init__(self):
        self.tree_processor = summarize.TreeSummarizer()

    def write_bt_commands(
            self,
            filepath,
            state_symbols,
            ml=True,
            ):
        bt_commands = []
        bt_commands.append("1") # multstate
        if ml:
            bt_commands.append("1") # ml
        else:
            bt_commands.append("2") # mcmc
        if len(state_symbols) > 7:
            bt_commands.append("restrictall q{}{}".format(state_symbols[0], state_symbols[1]))
        bt_commands.append("run")
        run_file = open(filepath, "w")
        run_file.write("\n".join(bt_commands))
        run_file.write("\n")
        run_file.close()

    def create_analysis(
            self,
            prefix,
            taxon_keys,
            taxon_state_set_map,
            tree_file,
            ml=True):
        data_filepath =  prefix + ".txt"
        dataf = open(data_filepath, "w")
        symbols = analyses.NameToSymbolMap()
        for taxon in taxon_keys:
            row = [taxon.replace(" ", "_")]
            states = sorted([symbols[s] for s in taxon_state_set_map[taxon]])
            row.append("".join(states))
            dataf.write("{}\n".format("\t".join(row)))
        dataf.close()
        run_filepath = prefix + ".bayestraits"
        write_bt_commands(
                filepath=run_filepath,
                state_symbols=symbols.SYMBOLS,
                ml=ml)
        write_job_file(
                job_filepath=prefix + ".job",
                tree_filepath=tree_file,
                data_filepath=data_filepath,
                run_filepath=run_filepath)

    def estimate_niche_evolution_rate(self, trees):
        trees = self.tree_processor.encode_trees(trees)
        for tree, tree_idx in enumerate(trees):
            pass



