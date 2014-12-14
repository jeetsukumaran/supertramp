#! /usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import sys
import os
import re
import collections
import tempfile
import subprocess
import dendropy
from dendropy.utility import processio
from supertramp import postprocess

class TraitEvolutionRateEstimator(object):

    def __init__(self,
            exclude_first_island_as_continental_source_outside_of_analysis,
            drop_stunted_trees=True,
            ):
        self.tree_postprocessor = postprocess.TreePostProcessor(
                exclude_first_island_as_continental_source_outside_of_analysis=exclude_first_island_as_continental_source_outside_of_analysis,
                # drop_trees_not_occupying_all_islands=True,
                # drop_trees_not_occupying_all_habitats=True,
                drop_stunted_trees=drop_stunted_trees,
                )
        # self.tree_file = tempfile.NamedTemporaryFile()
        # self.data_file = tempfile.NamedTemporaryFile()
        # self.tree_file_name = self.tree_file.name
        # self.data_file_name = self.data_file.name
        self.tree_file_name = "x1.tre"
        self.data_file_name = "x1.txt"

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

    def write_job_file(
            self,
            job_filepath,
            tree_filepath,
            data_filepath,
            run_filepath):
        jobf = open(job_filepath, "w")
        jobf.write("""\
    #! /bin/bash
    #$ -cwd
    #$ -V
    #$ -S /bin/bash
    #$ -l h_vmem=8G
    #$ -l virtual_free=8G

    {app_path} {tree_path} {data_path} < {commands_path}
    """.format(
        app_path="BayesTraits",
        tree_path=os.path.abspath(tree_filepath),
        data_path=os.path.abspath(data_filepath),
        commands_path=os.path.abspath(run_filepath),
        ))

    def create_analysis(
            self,
            prefix,
            taxon_keys,
            taxon_state_set_map,
            tree_file,
            ml=True):
        data_filepath =  prefix + ".txt"
        dataf = open(data_filepath, "w")
        symbols = postprocess.NameToSymbolMap()
        for taxon in taxon_keys:
            row = [taxon.label.replace(" ", "_")]
            states = sorted([symbols[s] for s in taxon_state_set_map[taxon]])
            row.append("".join(states))
            dataf.write("{}\n".format("\t".join(row)))
        dataf.close()
        run_filepath = prefix + ".bayestraits"
        self.write_bt_commands(
                filepath=run_filepath,
                state_symbols=symbols.SYMBOLS,
                ml=ml)
        self.write_job_file(
                job_filepath=prefix + ".job",
                tree_filepath=tree_file,
                data_filepath=data_filepath,
                run_filepath=run_filepath)

    def estimate_niche_evolution_rate(self, trees):
        trees = self.tree_postprocessor.process_trees(trees)
        for tree_idx, tree in enumerate(trees):

            taxa = tree.poll_taxa()
            taxon_state_set_map = {}
            for taxon in taxa:
                taxon_state_set_map[taxon] = set()
                for idx, i in enumerate(taxon.habitat_code):
                    if i == "1":
                        taxon_state_set_map[taxon].add(str(idx+1))

            tree.taxon_namespace = dendropy.TaxonNamespace(taxa)
            for nd in tree:
                nd.label = None # BayesTraits gets confused with internal taxon labels, especially those with periods etc.
            tree.write_to_path(
                    self.tree_file_name,
                    "nexus",
                    translate_tree_taxa=True)

            name_to_symbol_map = postprocess.NameToSymbolMap()
            dataf = open(self.data_file_name, "w")
            for taxon in taxa:
                row = [taxon.label]
                states = sorted([name_to_symbol_map[s] for s in taxon_state_set_map[taxon]])
                row.append("".join(states))
                dataf.write("{}\n".format("\t".join(row)))
            dataf.close()

            bt_commands = []
            bt_commands.append("1") # multstate
            bt_commands.append("1") # ml; 2 == mcmc
            if True: #len(name_to_symbol_map.SYMBOLS) > 7:
                bt_commands.append("restrictall q{}{}".format(
                    name_to_symbol_map.SYMBOLS[0],
                    name_to_symbol_map.SYMBOLS[1]))
            bt_commands.append("run")
            bt_commands = "\n".join(bt_commands)
            p = subprocess.Popen(
                    ["BayesTraits", self.tree_file_name, self.data_file_name],
                    stdout=subprocess.PIPE,
                    stdin=subprocess.PIPE,
                    )
            stdout, stderr = processio.communicate(p, bt_commands)
            stdout = stdout.split("\n")
            result = dict(zip(stdout[-3].split("\t"), stdout[-2].split("\t")))
            del result['']
            print(result)
            # self.create_analysis(
            #         prefix="x{}".format(tree_idx),
            #         taxon_keys=taxa,
            #         taxon_state_set_map=tsm,
            #         tree_file="z{}.nexus".format(tree_idx),
            #         )



