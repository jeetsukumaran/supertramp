#! /usr/bin/env python

import collections
import dendropy
from dendropy import treecalc

if dendropy.__version__.startswith("4"):
    _get_taxa = lambda x: x.taxon_namespace
else:
    _get_taxa = lambda x: x.taxon_set

class TreeProcessor(object):

    def __init__(self):

        self.island_colors = collections.defaultdict(lambda: "#666666")
        self.island_colors.update({
                "000001" : "#ff0000",
                "000010" : "#00cc00",
                "000100" : "#0000ff",
                "001000" : "#cccc00",
                "010000" : "#cccc00",
                "100000" : "#cccc00",
                })
        self.habitat_colors = {}
        self.habitat_colors.update({
                "01" : "#ff00ff",
                "10" : "#00ffff",
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
        for tree in trees:
            num_tips = 0
            total_length = 0.0
            nodes_by_island = collections.defaultdict(list)
            nodes_by_habitat = collections.defaultdict(list)
            for nd in tree:
                # colorize
                if nd.taxon is None and nd.label is None:
                    continue
                if nd.taxon and nd.taxon.habitat_color is not None:
                    nd.annotations.add_new("!color", nd.taxon.habitat_color)
                elif nd.label is not None:
                    self.encode_labeled_item(nd, annotate_habitat_color=True)
                # stats
                num_tips += 1
                total_length += nd.edge.length
                if nd.is_leaf():
                    island_code = nd.taxon.island_code
                    for idx, i in enumerate(island_code):
                        island_idx = len(island_code) - idx
                        if i == "1":
                            nodes_by_island[island_idx].append(nd)
                    habitat_code = nd.taxon.habitat_code
                    for idx, i in enumerate(habitat_code):
                        habitat_idx = len(habitat_code) - idx
                        if i == "1":
                            nodes_by_habitat[habitat_idx].append(nd)
            pdm = treecalc.PatristicDistanceMatrix(tree=tree)

            tree.stats = collections.defaultdict(lambda:"NA")
            if params is not None:
                tree.stats.update(params)
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
                summaries.append(tree.stats.copy())

        if trees_outf is not None:
            try:
                trees.write_to_stream(trees_outf, "nexus")
            except AttributeError:
                self.write_nexus(trees, trees_outf)
        return trees, stats_fields

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

