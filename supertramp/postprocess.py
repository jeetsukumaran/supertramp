#! /usr/bin/env python

import sys
import os
import collections
import dendropy
from supertramp import BitVector


class OutOfRegionError(Exception):
    pass

class ColorAssigner(object):

    COLORS = (
        "#ff0000", "#00ff00", "#0000ff", "#ffff00", "#ff00ff", "#00ffff", "#000000",
        "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#808080",
        "#c00000", "#00c000", "#0000c0", "#c0c000", "#c000c0", "#00c0c0", "#c0c0c0",
        "#400000", "#004000", "#000040", "#404000", "#400040", "#004040", "#404040",
        "#200000", "#002000", "#000020", "#202000", "#200020", "#002020", "#202020",
        "#600000", "#006000", "#000060", "#606000", "#600060", "#006060", "#606060",
        "#a00000", "#00a000", "#0000a0", "#a0a000", "#a000a0", "#00a0a0", "#a0a0a0",
        "#e00000", "#00e000", "#0000e0", "#e0e000", "#e000e0", "#00e0e0", "#e0e0e0",
        "#666666",
        )
    COLOR_INDEXES = {}
    for idx, c in enumerate(COLORS):
        COLOR_INDEXES[c] = idx

    def __init__(self, offset=0):
        self.assigned_colors = {}
        self.offset = offset

    def __getitem__(self, code):
        try:
            return self.assigned_colors[code]
        except KeyError:
            on_bits = []
            for idx, i in enumerate(code):
                if i == "1":
                    on_bits.append(idx)
            if len(on_bits) > 1:
                self.assigned_colors[code] = "#666666"
            elif len(on_bits) == 0:
                raise OutOfRegionError
            else:
                self.assigned_colors[code] = self.COLORS[on_bits[0]+self.offset]
            return self.assigned_colors[code]

class NameToSymbolMap(object):

    SYMBOLS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    SYMBOL_INDEXES = {}
    for idx, s in enumerate(SYMBOLS):
        SYMBOL_INDEXES[s] = idx

    def __init__(self):
        self.assigned_symbols = {}

    def __getitem__(self, name):
        try:
            return self.assigned_symbols[name]
        except KeyError:
            sidx = len(self.assigned_symbols)
            assert sidx < len(NameToSymbolMap.SYMBOLS)
            self.assigned_symbols[name] = NameToSymbolMap.SYMBOLS[sidx]
            return self.assigned_symbols[name]

    def __iter__(self):
        for idx in range(len(self.assigned_symbols)):
            return NameToSymbolMap.SYMBOLS[idx]

    def __len__(self):
        return len(self.assigned_symbols)

class TreePostProcessor(object):

    def __init__(self,
            exclude_first_island_as_continental_source_outside_of_analysis,
            drop_stunted_trees,
            ):
        self.island_colors = ColorAssigner()
        self.habitat_colors = ColorAssigner()
        self.drop_stunted_trees = drop_stunted_trees
        self.exclude_first_island_as_continental_source_outside_of_analysis = exclude_first_island_as_continental_source_outside_of_analysis

    def decode_labeled_item_biogeography(self, t):
        if hasattr(t, "is_encoded") and t.is_encoded:
            return
        if self.exclude_first_island_as_continental_source_outside_of_analysis:
            p = t.label.split(".")
            t.label = ".".join([p[0], p[1][1:], p[2]])
        label_parts = t.label.split(".")
        t.island_code = label_parts[1]
        # if self.exclude_first_island_as_continental_source_outside_of_analysis:
        #     t.island_code = t.island_code[1:]
        try:
            t.island_color = self.island_colors[t.island_code]
            t.habitat_code = label_parts[2]
            t.habitat_color = self.habitat_colors[t.habitat_code]
            t.is_encoded = True
            t.out_of_region = False
        except OutOfRegionError:
            # taxon only occurs on first island, and needs to be removed
            t.out_of_region = True

    def decode_taxon_biogeography(self, taxa):
        if hasattr(taxa, "is_encoded") and taxa.is_encoded:
            return []
        taxa_to_exclude = []
        for t in taxa:
            self.decode_labeled_item_biogeography(t)
            if t.out_of_region:
                taxa_to_exclude.append(t)
        taxa.is_encoded = True
        return taxa_to_exclude

    def prune_trees(self, trees, taxa_to_exclude):
        if not taxa_to_exclude:
            return trees
        to_keep = dendropy.TreeList(taxon_namespace=trees.taxon_namespace)
        for tree in trees:
            try:
                tree.prune_taxa(taxa_to_exclude)
                if not self.drop_stunted_trees or tree.seed_node.num_child_nodes() > 1:
                    to_keep.append(tree)
            except AttributeError:
                # trying to prune root node
                pass
        for taxon in taxa_to_exclude:
            trees.taxon_namespace.remove_taxon(taxon)
        return to_keep

    def write_colorized_trees(self,
            outf,
            trees,
            schema,
            process_trees):
        if schema == "by-island":
            taxon_color_attr = "island_color"
            branch_color_attr = "habitat_color"
        elif schema == "by-habitat":
            taxon_color_attr = "habitat_color"
            branch_color_attr = "island_color"
        else:
            raise ValueError("Unrecognized schema: {}".format(schema))
        if process_trees:
            trees = self.process_trees(trees)
        for taxon in trees.taxon_namespace:
            taxon.annotations["!color"] = getattr(taxon, taxon_color_attr)
        trees.write_to_stream(outf, "nexus")

    def process_trees(self, trees):
        taxa_to_exclude = self.decode_taxon_biogeography(trees.taxon_namespace)
        trees = self.prune_trees(trees, taxa_to_exclude)
        return trees

