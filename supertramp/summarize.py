#! /usr/bin/env python

import sys
import os
import collections
import dendropy
import subprocess
from dendropy.calculate import treemeasure
from dendropy.model import birthdeath
from dendropy.utility import processio
from supertramp import BitVector
from supertramp import postprocess

class Rcalculator(object):

    RESULT_FLAG_LEADER = "[!!!]"

    def __init__(self):
        pass

    def execute_rscript(self, script):
        cmd = []
        cmd.append("Rscript")
        cmd.append("--vanilla")
        cmd.append("-")
        p = subprocess.Popen(cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                )
        stdout, stderr = processio.communicate(p, script)
        if p.returncode != 0:
            print(script)
            for row in stderr.split("\n"):
                print("# {}".format(row))
            sys.exit(p.returncode)
        results = {}
        num_lines_with_results = 0
        for line in stdout.split("\n"):
            if not line.startswith(Rcalculator.RESULT_FLAG_LEADER):
                continue
            parts = line[len(Rcalculator.RESULT_FLAG_LEADER):].split("=")
            assert len(parts) == 2
            key = parts[0].strip()
            try:
                value = float(parts[1].strip())
            except ValueError as e:
                value = "NA"
            results[key] = value
            num_lines_with_results += 1
        return results

    def _compose_cophenetic_matrix(
            self,
            dists,
            taxon_names,
            byrow=True,
            ):
        return "matrix(c({data}), nrow={nrow}, byrow={byrow}, dimnames=list(c({names}),c({names})))".format(
            data=",".join("{}".format(d) for d in dists),
            nrow=len(taxon_names),
            byrow="T" if byrow else "F",
            names=",".join(taxon_names))

    def _compose_community_matrix(
            self,
            data,
            comm_names,
            taxon_names,
            byrow=True,
            ):
        return "matrix(c({data}), nrow={nrow}, byrow={byrow}, dimnames=list(c({comm_names}),c({taxon_names})))".format(
            data=",".join("{}".format(d) for d in data),
            nrow=len(comm_names),
            byrow="T" if byrow else "F",
            comm_names=",".join(comm_names),
            taxon_names=",".join(taxon_names)
            )

    def calc_ecological_stats(
            self,
            tree,
            patristic_distance_matrix,
            total_tree_length,
            total_tree_edges,
            nodes_by_island,
            nodes_by_habitat,
            disturbed_habitat_nodes,
            interior_habitat_nodes,
            ):
        pdm = patristic_distance_matrix
        leaf_taxa = set([leaf.taxon for leaf in tree.leaf_node_iter()])
        tree_taxa = [t for t in tree.taxon_namespace if t in leaf_taxa] # to maintain order
        taxon_names = []
        weighted_dists = []
        unweighted_dists = []
        normalized_weighted_dists = []
        normalized_unweighted_dists = []
        for taxon1 in tree_taxa:
            taxon_names.append("'{}'".format(taxon1.label))
            for taxon2 in tree_taxa:
                weighted_dist = pdm(taxon1, taxon2)
                unweighted_dist = pdm.path_edge_count(taxon1, taxon2)
                normalized_weighted_dist = weighted_dist / total_tree_length
                normalized_unweighted_dist = unweighted_dist / total_tree_edges
                weighted_dists.append(weighted_dist)
                unweighted_dists.append(unweighted_dist)
                normalized_weighted_dists.append(normalized_weighted_dist)
                normalized_unweighted_dists.append(normalized_unweighted_dist)
        rscript = []
        rscript.append("suppressMessages(library(picante))")
        for dists, dists_desc in (
                    # (weighted_dists, "weighted"),
                    # (unweighted_dists, "unweighted"),
                    (normalized_weighted_dists, "normalized_weighted"),
                    (normalized_unweighted_dists, "normalized_unweighted"),
                ):
            cophenetic_dist_matrix_str = self._compose_cophenetic_matrix(
                    dists=dists,
                    taxon_names=taxon_names)
            cophenetic_dist_matrix_name = "{}_cophenetic_dist_matrix".format(dists_desc)
            rscript.append("{} <- {}".format(cophenetic_dist_matrix_name, cophenetic_dist_matrix_str))
            for comm_prefix, comm_data, comm_desc in (
                ("I", nodes_by_island, "by_island"),
                ("H", nodes_by_habitat, "by_habitat"),
                    ):
                comm_names = []
                pa_data = []
                for idx in comm_data:
                    comm_taxa = [nd.taxon for nd in comm_data[idx]]
                    # print("# {}: {}: {} of {}: {}\n".format(comm_desc, idx, len(comm_taxa), len(tree_taxa), ", ".join(t.label for t in comm_taxa)))
                    comm_names.append("'{}{}'".format(comm_prefix, idx))
                    for taxon in tree_taxa:
                        if taxon in comm_taxa:
                            pa_data.append(1)
                        else:
                            pa_data.append(0)
                comm_pa_matrix_name = "community_{}".format(comm_desc)
                comm_pa_matrix_str = self._compose_community_matrix(
                        data=pa_data,
                        comm_names=comm_names,
                        taxon_names=["'{}'".format(t.label) for t in tree_taxa])
                rscript.append("{} <- {}".format(comm_pa_matrix_name, comm_pa_matrix_str))
                nruns = 100
                prefix = "{comm}.{dists}".format(comm=comm_pa_matrix_name,
                        dists=cophenetic_dist_matrix_name.replace("_cophenetic_dist_matrix", "")).replace("_", ".")
                out = "stdout()"
                # out = "'z.txt'"
                for stat_type in ("mpd", "mntd"):
                    rscript.append("result <- ses.{stat_type}({comm},{dists},null.model='taxa.labels',abundance.weighted=FALSE,runs={nruns})".format(
                        stat_type=stat_type,
                        comm=comm_pa_matrix_name,
                        dists=cophenetic_dist_matrix_name,
                        nruns=nruns))
                    rscript.append("result.df <- as.data.frame(result)")
                    if comm_desc == "by_habitat":
                        rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.Z.', rownames(result.df), ' = ', result.df${stat_type}.obs.z, '\\n', sep=''), {out})".format(
                            stat_type=stat_type,
                            result_flag=Rcalculator.RESULT_FLAG_LEADER,
                            prefix=prefix,
                            out=out,
                            ))
                        rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.p.', rownames(result.df), ' = ', result.df${stat_type}.obs.p, '\\n', sep=''), {out})".format(
                            stat_type=stat_type,
                            result_flag=Rcalculator.RESULT_FLAG_LEADER,
                            prefix=prefix,
                            out=out,
                            ))
                    rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.Z.mean', ' = ', mean(result.df${stat_type}.obs.z), '\\n', sep=''), {out})".format(
                        stat_type=stat_type,
                        result_flag=Rcalculator.RESULT_FLAG_LEADER,
                        prefix=prefix,
                        out=out,
                        ))
                    rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.p.mean', ' = ', mean(result.df${stat_type}.obs.p), '\\n', sep=''), {out})".format(
                        stat_type=stat_type,
                        result_flag=Rcalculator.RESULT_FLAG_LEADER,
                        prefix=prefix,
                        out=out,
                        ))
                    rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.Z.var', ' = ', var(result.df${stat_type}.obs.z), '\\n', sep=''), {out})".format(
                        stat_type=stat_type,
                        result_flag=Rcalculator.RESULT_FLAG_LEADER,
                        prefix=prefix,
                        out=out,
                        ))
                    rscript.append("write(paste('{result_flag}', '{prefix}.{stat_type}.obs.p.var', ' = ', var(result.df${stat_type}.obs.p), '\\n', sep=''), {out})".format(
                        stat_type=stat_type,
                        result_flag=Rcalculator.RESULT_FLAG_LEADER,
                        prefix=prefix,
                        out=out,
                        ))
        rscript = "\n".join(rscript)
        results = self.execute_rscript(rscript)
        tree.stats.update(results)
        return results

class TreeSummarizer(object):

    def __init__(self,
            exclude_first_island_as_continental_source_outside_of_analysis,
            drop_trees_not_occupying_all_habitats,
            drop_trees_not_occupying_all_islands,
            drop_stunted_trees,
            ):
        self.tree_postprocessor = postprocess.TreePostProcessor(
                exclude_first_island_as_continental_source_outside_of_analysis=exclude_first_island_as_continental_source_outside_of_analysis,
                drop_stunted_trees=drop_stunted_trees,
                )
        self.drop_trees_not_occupying_all_islands = drop_trees_not_occupying_all_islands
        self.drop_trees_not_occupying_all_habitats = drop_trees_not_occupying_all_habitats
        self.rcalc = Rcalculator()

    def get_mean_patristic_distance(self, pdm, nodes):
        if len(nodes) <= 1:
            return "NA", "NA"
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

    def summarize_trees(self,
            trees,
            trees_outf=None,
            params=None,
            summaries=None):
        trees = self.tree_postprocessor.process_trees(trees)
        stats_fields = set()

        # crucial assumption here is all trees from same landscape wrt to
        # number of islands and habitats
        representative_taxon = trees[0].taxon_namespace[0]
        community_by_disturbed_vs_interior_habitat = {}
        num_islands = len(representative_taxon.island_code)
        num_habitats = len(representative_taxon.habitat_code)
        # community_by_island = {}
        # community_by_habitat = {}
        # for i in num_islands:
        #     community_by_island[i] = {}
        # for i in num_habitats:
        #     community_by_habitat[i] = {}
        # community_by_disturbed_vs_interior_habitat[0] = {}
        # community_by_disturbed_vs_interior_habitat[1] = {}

        for tree in list(trees):
            num_tips = 0
            total_length = 0.0
            total_edges = 0
            nodes_by_island = collections.defaultdict(list)
            nodes_by_habitat = collections.defaultdict(list)
            disturbed_habitat_nodes = []
            interior_habitat_nodes = []

            all_tips = []
            for nd in tree:
                # colorize
                if nd.taxon is None and nd.label is None:
                    continue
                if nd.label is not None:
                    self.tree_postprocessor.decode_labeled_item_biogeography(nd)
                # stats
                total_edges += 1
                num_tips += 1
                total_length += nd.edge.length
                if nd.is_leaf():
                    all_tips.append(nd)
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
                            if habitat_idx == 0:
                                disturbed_habitat_nodes.append(nd)
                            else:
                                if nd not in interior_habitat_nodes:
                                    interior_habitat_nodes.append(nd)
            if len(nodes_by_island) < num_islands and self.drop_trees_not_occupying_all_islands:
                trees.remove(tree)
                continue
            if len(nodes_by_habitat) < num_habitats and self.drop_trees_not_occupying_all_habitats:
                trees.remove(tree)
                continue
            pdm = treemeasure.PatristicDistanceMatrix(tree=tree)
            tree.stats = collections.defaultdict(lambda:"NA")
            if params is not None:
                tree.params = params.copy()
            tree.stats["size"] = num_tips
            tree.stats["length"] = total_length
            tree.stats["edges"] = total_edges
            # node_ages = tree.internal_node_ages()
            # node_ages = [n/total_length for n in node_ages]
            # tree.stats["est.birth.rate"] = birthdeath.fit_pure_birth_model(internal_node_ages=node_ages)["birth_rate"]
            tree.stats["est.birth.rate"] = birthdeath.fit_pure_birth_model(tree=tree)["birth_rate"]

            weighted_disturbed, unweighted_disturbed = self.get_mean_patristic_distance(pdm, disturbed_habitat_nodes)
            weighted_interior, unweighted_interior = self.get_mean_patristic_distance(pdm, interior_habitat_nodes)
            tree.stats["weighted.disturbed.habitat.pd"] = weighted_disturbed
            tree.stats["unweighted.disturbed.habitat.pd"] = unweighted_disturbed
            tree.stats["weighted.interior.habitat.pd"] = weighted_interior
            tree.stats["unweighted.interior.habitat.pd"] = unweighted_interior
            try:
                tree.stats["weighted.disturbed.to.interior.habitat.pd"] = weighted_disturbed/weighted_interior
                tree.stats["unweighted.disturbed.to.interior.habitat.pd"] = unweighted_disturbed/unweighted_interior
            except (ZeroDivisionError, TypeError):
                tree.stats["weighted.disturbed.to.interior.habitat.pd"] = "NA"
                tree.stats["unweighted.disturbed.to.interior.habitat.pd"] = "NA"

            rstats = self.rcalc.calc_ecological_stats(
                    tree=tree,
                    patristic_distance_matrix=pdm,
                    total_tree_length=total_length,
                    total_tree_edges=total_edges,
                    nodes_by_island=nodes_by_island,
                    nodes_by_habitat=nodes_by_habitat,
                    disturbed_habitat_nodes=disturbed_habitat_nodes,
                    interior_habitat_nodes=interior_habitat_nodes,
                    )

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

    def write_colorized_trees(self, outf, trees, schema, is_trees_postprocessed):
        self.tree_postprocessor.write_colorized_trees(outf, trees, schema, is_trees_postprocessed)
