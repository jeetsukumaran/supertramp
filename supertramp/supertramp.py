#! /usr/bin/env python

##############################################################################
##
##  Copyright 2010-2014 Jeet Sukumaran.
##  All rights reserved.
##
##  Redistribution and use in source and binary forms, with or without
##  modification, are permitted provided that the following conditions are met:
##
##      * Redistributions of source code must retain the above copyright
##        notice, this list of conditions and the following disclaimer.
##      * Redistributions in binary form must reproduce the above copyright
##        notice, this list of conditions and the following disclaimer in the
##        documentation and/or other materials provided with the distribution.
##      * The names of its contributors may not be used to endorse or promote
##        products derived from this software without specific prior written
##        permission.
##
##  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
##  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
##  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
##  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN OR MARK T. HOLDER
##  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
##  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
##  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
##  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
##  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
##  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
##  POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

import io
import sys
import argparse
import random
import logging
import collections
import numpy

_LOGGING_LEVEL_ENVAR = "SUPERTRAMP_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "SUPERTRAMP_LOGGING_FORMAT"

class RunLogger(object):

    def __init__(self, **kwargs):
        self.name = kwargs.get("name", "RunLog")
        self._log = logging.getLogger(self.name)
        self._log.setLevel(logging.DEBUG)
        if kwargs.get("log_to_stderr", True):
            ch1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", logging.INFO))
            ch1.setLevel(stderr_logging_level)
            ch1.setFormatter(self.get_default_formatter())
            self._log.addHandler(ch1)
        if kwargs.get("log_to_file", True):
            log_stream = kwargs.get("log_stream", \
                open(kwargs.get("log_path", self.name + ".log"), "w"))
            ch2 = logging.StreamHandler(log_stream)
            file_logging_level = self.get_logging_level(kwargs.get("file_logging_level", logging.DEBUG))
            ch2.setLevel(file_logging_level)
            ch2.setFormatter(self.get_default_formatter())
            self._log.addHandler(ch2)

    def get_logging_level(self, level=None):
        if level in [logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING,
            logging.ERROR, logging.CRITICAL]:
            return level
        elif level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        if level_name == "NOTSET":
            level = logging.NOTSET
        elif level_name == "DEBUG":
            level = logging.DEBUG
        elif level_name == "INFO":
            level = logging.INFO
        elif level_name == "WARNING":
            level = logging.WARNING
        elif level_name == "ERROR":
            level = logging.ERROR
        elif level_name == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
        return level

    def get_default_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_rich_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(filename)s (%(lineno)d): %(levelname) 8s: %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_simple_formatter(self):
        return logging.Formatter("%(levelname) 8s: %(message)s")

    def get_raw_formatter(self):
        return logging.Formatter("%(message)s")

    def get_logging_formatter(self, format=None):
        if format is not None:
            format = format.upper()
        elif _LOGGING_FORMAT_ENVAR in os.environ:
            format = os.environ[_LOGGING_FORMAT_ENVAR].upper()
        if format == "RICH":
            logging_formatter = self.get_rich_formatter()
        elif format == "SIMPLE":
            logging_formatter = self.get_simple_formatter()
        elif format == "NONE":
            logging_formatter = self.get_raw_formatter()
        else:
            logging_formatter = self.get_default_formatter()
        if logging_formatter is not None:
            logging_formatter.datefmt='%H:%M:%S'

    def debug(self, msg, *args, **kwargs):
        self._log.debug(msg, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self._log.info(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self._log.warning(msg, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self._log.error(msg, *args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self._log.critical(msg, *args, **kwargs)


class HabitatType(object):

    def __init__(self, label):
        self.label = label

class Habitat(object):

    def __init__(self, habitat_type):
        self.habitat_type = habitat_type
        pass

class Island(object):

    def __init__(self,
            rng,
            habitat_types,
            dispersal_habitat_types=None):
        self.rng = rng
        self.habitat_types = habitat_types
        self.dispersal_habitat_types = dispersal_habitat_types
        self.habitats = {}
        self.lineages = collections.OrderedDict()
        for ht in self.habitat_types:
            h = Habitat(habitat_type=ht)
            self.habitats[ht] = h

class Lineage(object):

    counter = 0

    def __init__(self, parent=None):
        Lineage.counter += 1
        self.index = Lineage.counter
        self.label = "s{}".format(self.index)
        self.age = 0
        self.parent = parent
        self.child_nodes = []

    def add_age_to_tips(self, ngens=1):
        """
        Grows tree by adding ``ngens`` time unit(s) to all tips.
        """
        if self.child_nodes:
            for nd in self.leaf_iter():
                nd.age += 1
        else:
            self.age += 1

    def leaf_iter(self, filter_fn=None):
        """
        Returns an iterator over the leaf_nodes that are descendants of self
        (with leaves returned in same order as a post-order traversal of the
        tree).
        """
        if filter_fn:
            ff = lambda x: (not x.child_nodes) and filter_fn(x) or None
        else:
            ff = lambda x: (not x.child_nodes) and x or None
        for node in self.postorder_iter(ff):
            yield node

    def postorder_iter(self, filter_fn=None):
        """
        Postorder traversal of the self and its child_nodes.  Returns self
        and all descendants such that a node's child_nodes (and their
        child_nodes) are visited before node.  Filtered by filter_fn:
        node is only returned if no filter_fn is given or if filter_fn
        returns True.
        """
        stack = [(self, False)]
        while stack:
            node, state = stack.pop(0)
            if state:
                if filter_fn is None or filter_fn(node):
                    yield node
            else:
                stack.insert(0, (node, True))
                child_nodes = [(n, False) for n in node.child_nodes]
                child_nodes.extend(stack)
                stack = child_nodes

    def diversify(self):
        """
        Spawns two child lineages with self as parent.
        Returns tuple consisting of these two lineages.
        """
        if self.child_nodes:
            raise Exception("Trying to diversify internal node: {}: {}".format(self.label, ", ".join(c.label for c in self.child_nodes)))
        c1 = Lineage(parent=self)
        c2 = Lineage(parent=self)
        self.child_nodes.append(c1)
        self.child_nodes.append(c2)
        return (c1, c2)

    def __str__(self):
        return self.label

    def __repr__(self):
        return "<Lineage {}>".format(self.label)

    ###########################################################################
    ## Hacked-in NEWICK representation.

    def as_newick_string(self, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.
        For production purposes, use the the full-fledged 'as_string()'
        method of the object.
        """
        out = io.StringIO()
        self.write_newick(out, **kwargs)
        return out.getvalue()

    def write_newick(self, out, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.  For
        production purposes, use the the full-fledged 'write_to_stream()'
        method of the object.
        """
        child_nodes = self.child_nodes
        if child_nodes:
            out.write('(')
            f_child = child_nodes[0]
            for child in child_nodes:
                if child is not f_child:
                    out.write(',')
                child.write_newick(out, **kwargs)
            out.write(')')
        out.write(self.label)
        out.write(":{}".format(self.age))


class System(object):

    def __init__(self, random_seed=None):
        self.logger = RunLogger(name="supertramp")
        if random_seed is None:
            self.random_seed = random.randint(0, sys.maxsize)
        else:
            self.random_seed = random_seed
        self.log_frequency = 1
        self.current_gen = 0
        self.logger.info("Initializing with random seed {}".format(self.random_seed))
        self.rng = numpy.random.RandomState(seed=[self.random_seed])

        self.num_islands = 2
        self.habitat_type_labels = ["coastal", "interior", "deep"]
        self.dispersal_model = "unconstrained"
        self.habitat_types = []
        self.islands = []
        self.phylogeny = None
        self.bootstrap()

    def bootstrap(self):
        self.phylogeny = Lineage()
        for ht_label in self.habitat_type_labels:
            h = HabitatType(label=ht_label)
            self.habitat_types.append(h)
        if self.dispersal_model = "unconstrained":
            self.dispersal_habitat_types = list(self.habitat_types)
        else:
            self.dispersal_habitat_types = self.habitat_types[0]
        for isl_idx in range(self.num_islands):
            island = Island(
                    rng=self.rng,
                    habitat_types=self.habitat_types,
                    dispersal_habitat_types=self.dispersal_habitat_types)

    def execute_life_cycle(self):
        self.current_gen += 1
        self.phylogeny.add_age_to_tips(1)
            if self.current_gen % self.log_frequency == 0:
                self.logger.info("Executing life-cycle {}".format(self.current_gen))

