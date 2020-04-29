import os
import unittest
import numpy as np

import tskit
import pyslim
from graft import *


def run_slim_script(slimfile, args=''):
    print("running:: " + "slim -s 23 " + args + " " + slimfile)
    out = os.system("slim -s 23 " + args + " " + slimfile + ">/dev/null")
    return out


def get_examples(T1, T2):
    print(os.getcwd())
    # note: could make N and gens parameters here
    run_slim_script("tests/recipe.slim", args="-d N=100 -d gens=100 -d \"outfile='tests/data/root.trees'\"")
    run_slim_script("tests/recipe.slim", args="-d N=100 -d gens={} -d \"infile='tests/data/root.trees'\" -d 'outfile=\"tests/data/branch1.trees\"'".format(T1))
    run_slim_script("tests/recipe.slim", args="-d N=100 -d gens={} -d \"infile='tests/data/root.trees'\" -d 'outfile=\"tests/data/branch2.trees\"'".format(T2))
    ts1 = pyslim.load("tests/data/branch1.trees")
    ts2 = pyslim.load("tests/data/branch2.trees")
    return ts1, ts2


class TestFindSplitTime(unittest.TestCase):

    def test_simple_example(self):
        for (T1, T2) in [(100, 100), (100, 200), (200, 10)]:
            ts1, ts2 = get_examples(T1, T2)
            S1, S2 = find_split_time(ts1, ts2)
            self.assertEqual(S1, T1)
            self.assertEqual(S2, T2)


class TestMatchNodes(unittest.TestCase):

    def verify_match_nodes(self, ts1, ts2):
        T1, T2 = find_split_time(ts1, ts2)
        matched_nodes = match_nodes(ts1, ts2, T2)
        # note: should also check that we got all the nodes
        for (a, b) in matched_nodes:
            n1 = ts1.node(a)
            n2 = ts2.node(b)
            self.assertGreaterEqual(n2.time, T2)
            self.assertEqual(n1.metadata.slim_id, n2.metatdata.slim_id)

    def test_simple_example(self):
        ts1, ts2 = get_examples(100, 100)
        self.verify_match_nodes(ts1, ts2)


class TestGraft(unittest.TestCase):

    def test_simple_example(self):
        ts1, ts2 = get_examples(35, 12)
        T1, T2 = find_split_time(ts1, ts2)
        node_map21 = match_nodes(ts1, ts2, T2)

        # easier to deal with tskit treeseqs
        ts1 = ts1.tables.tree_sequence()
        ts2 = ts2.tables.tree_sequence()

        tsg, all_node_map2new, pop_map2new = graft(ts1, ts2, node_map21)

        # ts1 from grafted
        ts1g = tsg.simplify(ts1.samples())
        tables1g = ts1g.tables
        tables1 = ts1.simplify().tables
        self.assertEqual(tables1g.nodes, tables1.nodes)
        self.assertEqual(tables1g.edges, tables1.edges)
        self.assertEqual(tables1g.sites, tables1.sites)
        self.assertEqual(tables1g.mutations, tables1.mutations)

        # testing ts2
        samp = ts2.samples()
        new_samp = np.array([all_node_map2new[nid] for nid in list(samp)])
        tables2 = ts2.simplify().tables
        tables2g = tsg.simplify(new_samp).tables
        # the test for the nodes table will be more complicated
        # because of the changes pop id, indiv id, and time shift
        # self.assertEqual(tables2g.nodes, tables2.nodes)
        self.assertEqual(tables2g.edges, tables2.edges)
        self.assertEqual(tables2g.sites, tables2.sites)
        self.assertEqual(tables2g.mutations, tables2.mutations)
        # tests: graft one tree to another, simplify over just the nodes in original tree, test if the tables are the same
        # test: take a grafted tree, simplify over samples in first and second tree, graft them back together and test if the same
