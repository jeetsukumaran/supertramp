#! /usr/bin/env python

import collections
import dendropy
from dendropy import treecalc
from supertramp import BitVector

if dendropy.__version__.startswith("4"):
    _get_taxa = lambda x: x.taxon_namespace
else:
    _get_taxa = lambda x: x.taxon_set

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
            else:
                self.assigned_colors[code] = self.COLORS[on_bits[0]+self.offset]
            return self.assigned_colors[code]

class TreeProcessor(object):

    def __init__(self):

        # self.island_colors = collections.defaultdict(lambda: "#666666")
        # self.island_colors.update({
        #         "000001" : "#ff0000",
        #         "000010" : "#00cc00",
        #         "000100" : "#0000ff",
        #         "001000" : "#cccc00",
        #         "010000" : "#cccc00",
        #         "100000" : "#cccc00",
        #         })
        # self.habitat_colors = {}
        # self.habitat_colors.update({
        #         "01" : "#ff00ff",
        #         "10" : "#00ffff",
        #         })
        self.island_colors = ColorAssigner()
        self.habitat_colors = ColorAssigner()
        self.drop_stunted_trees = True

    def write_colorized_trees(self, outf, trees, scheme):
        if scheme == "by-island":
            taxon_color_attr = "island_color"
            branch_color_attr = "habitat_color"
        elif scheme == "by-habitat":
            taxon_color_attr = "habitat_color"
            branch_color_attr = "island_color"
        else:
            raise ValueError("Unrecognized scheme: {}".format(scheem))
        taxa = _get_taxa(trees)
        self.encode_taxa(taxa)
        for taxon in taxa:
            taxon.annotations["!color"] = getattr(taxon, taxon_color_attr)
        # for tree in trees:
        #     for nd in tree:
        #         if nd.label is None:
        #             continue
        #         self.encode_labeled_item(nd)
        #         nd.annotations["!color"] = getattr(nd, branch_color_attr)
        trees.write_to_stream(outf, "nexus")

    def get_mean_patristic_distance(self, pdm, nodes):
        if len(nodes) == 1:
            return None, None
        weighted_dist = 0.0
        unweighted_dist = 0.0
        ncomps = 0
        for idx1, nd1 in enumerate(nodes[:-1]):
            for nd2 in nodes[idx1+1:]:
                assert nd1.taxon
                assert nd2.taxon
                weighted_dist += pdm(nd1.taxon, nd2.taxon)
                unweighted_dist += pdm.path_edge_count(nd1.taxon, nd2.taxon)
                ncomps += 1
        return weighted_dist/ncomps, unweighted_dist/ncomps

    def process_trees(self,
            trees,
            trees_outf=None,
            params=None,
            summaries=None):
        self.encode_taxa(_get_taxa(trees))
        stats_fields = set()
        for tree in list(trees):
            if self.drop_stunted_trees and tree.seed_node.num_child_nodes() <= 1:
                trees.remove(tree)
                continue
            num_tips = 0
            total_length = 0.0
            nodes_by_island = collections.defaultdict(list)
            nodes_by_habitat = collections.defaultdict(list)
            for nd in tree:
                # colorize
                if nd.taxon is None and nd.label is None:
                    continue
                if nd.label is not None:
                    self.encode_labeled_item(nd)
                # stats
                num_tips += 1
                total_length += nd.edge.length
                if nd.is_leaf():
                    island_code = nd.taxon.island_code
                    for idx, i in enumerate(island_code):
                        # island_idx = len(island_code) - idx
                        island_idx = idx
                        if i == "1":
                            nodes_by_island[island_idx].append(nd)
                    habitat_code = nd.taxon.habitat_code
                    for idx, i in enumerate(habitat_code):
                        # habitat_idx = len(habitat_code) - idx
                        habitat_idx = idx
                        if i == "1":
                            nodes_by_habitat[habitat_idx].append(nd)
            pdm = treecalc.PatristicDistanceMatrix(tree=tree)

            tree.stats = collections.defaultdict(lambda:"NA")
            if params is not None:
                tree.params = params.copy()
            tree.stats["size"] = num_tips
            tree.stats["length"] = total_length

            weighted_habitat_dist_total = 0.0
            weighted_habitat_dist_count = 0
            unweighted_habitat_dist_total = 0.0
            unweighted_habitat_dist_count = 0
            habitat_keys = sorted(nodes_by_habitat.keys())
            for habitat_key in habitat_keys:
                weighted_habitat_dist_key = "habitat.{}.pdist.weighted".format(habitat_key)
                unweighted_habitat_dist_key = "habitat.{}.pdist.unweighted".format(habitat_key)
                weighted, unweighted = self.get_mean_patristic_distance(pdm, nodes_by_habitat[habitat_key])
                if weighted is None:
                    assert unweighted is None
                    tree.stats[weighted_habitat_dist_key] = "NA"
                else:
                    weighted = weighted/total_length
                    weighted_habitat_dist_total += weighted
                    tree.stats[weighted_habitat_dist_key] = weighted
                    weighted_habitat_dist_count += 1
                if unweighted is None:
                    assert weighted is None
                    tree.stats[unweighted_habitat_dist_key] = "NA"
                else:
                    unweighted = unweighted/num_tips
                    unweighted_habitat_dist_total += unweighted
                    tree.stats[unweighted_habitat_dist_key] = weighted
                    unweighted_habitat_dist_count += 1
            assert weighted_habitat_dist_count == unweighted_habitat_dist_count
            try:
                tree.stats["habitat.mean.pdist.weighted"] = weighted_habitat_dist_total / weighted_habitat_dist_count
            except ZeroDivisionError:
                tree.stats["habitat.mean.pdist.weighted"] = "NA"
            try:
                tree.stats["habitat.mean.pdist.unweighted"] = unweighted_habitat_dist_total / unweighted_habitat_dist_count
            except ZeroDivisionError:
                tree.stats["habitat.mean.pdist.unweighted"] = "NA"

            weighted_island_dist_total = 0.0
            weighted_island_dist_count = 0
            unweighted_island_dist_total = 0.0
            unweighted_island_dist_count = 0
            island_keys = sorted(nodes_by_island.keys())
            for island_key in island_keys:
                weighted_island_dist_key = "island.{}.pdist.weighted".format(island_key)
                unweighted_island_dist_key = "island.{}.pdist.unweighted".format(island_key)
                weighted, unweighted = self.get_mean_patristic_distance(pdm, nodes_by_island[island_key])
                if weighted is None:
                    assert unweighted is None
                    tree.stats[weighted_island_dist_key] = "NA"
                else:
                    weighted = weighted/total_length
                    weighted_island_dist_total += weighted
                    tree.stats[weighted_island_dist_key] = weighted
                    weighted_island_dist_count += 1
                if unweighted is None:
                    assert weighted is None
                    tree.stats[unweighted_island_dist_key] = "NA"
                else:
                    unweighted = unweighted/num_tips
                    unweighted_island_dist_total += unweighted
                    tree.stats[unweighted_island_dist_key] = weighted
                    unweighted_island_dist_count += 1
            assert weighted_island_dist_count == unweighted_island_dist_count
            try:
                tree.stats["island.mean.pdist.weighted"] = weighted_island_dist_total / weighted_island_dist_count
            except ZeroDivisionError:
                tree.stats["island.mean.pdist.weighted"] = "NA"
            try:
                tree.stats["island.mean.pdist.unweighted"] = unweighted_island_dist_total / unweighted_island_dist_count
            except ZeroDivisionError:
                tree.stats["island.mean.pdist.unweighted"] = "NA"

            stats_fields.update(tree.stats.keys())

            if summaries is not None:
                sss = tree.stats.copy()
                sss.update(tree.params)
                summaries.append(sss)

        if trees_outf is not None:
            try:
                trees.write_to_stream(trees_outf, "nexus")
            except AttributeError:
                self.write_nexus(trees, trees_outf)
        return trees, stats_fields

    def encode_labeled_item(self, t):
        if hasattr(t, "is_encoded") and t.is_encoded:
            return
        label_parts = t.label.split(".")
        t.island_code = label_parts[1]
        t.island_color = self.island_colors[t.island_code]
        t.habitat_code = label_parts[2]
        t.habitat_color = self.habitat_colors[t.habitat_code]
        t.is_encoded = True

    def encode_taxa(self, taxa):
        if hasattr(taxa, "is_encoded") and taxa.is_encoded:
            return
        for t in taxa:
            self.encode_labeled_item(t)
        taxa.is_encoded = True


