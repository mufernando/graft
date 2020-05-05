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


def get_examples(T1, T2, gens=100, N=100):
    print(os.getcwd())
    # note: could make N and gens parameters here
    run_slim_script("tests/recipe.slim", args=f"-d N={N} -d gens={gens} -d \"outfile='tests/data/root.trees'\"")
    run_slim_script("tests/recipe.slim", args=f"-d N={N} -d gens={T1} -d \"infile='tests/data/root.trees'\" -d 'outfile=\"tests/data/branch1.trees\"'")
    run_slim_script("tests/recipe.slim", args=f"-d N={N} -d gens={T2} -d \"infile='tests/data/root.trees'\" -d 'outfile=\"tests/data/branch2.trees\"'")
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
        node_map21 = match_nodes(ts1, ts2, T2)
        # note: should also check that we got all the nodes
        for a,b in node_map21.items():
            n1 = ts1.node(b)
            n2 = ts2.node(a)
            self.assertGreaterEqual(n2.time, T2)
            self.assertEqual(n1.metadata.slim_id, n2.metadata.slim_id)

    def verify_simplification_nodes(self, ts1, ts2):
        T1, T2 = find_split_time(ts1, ts2)
        node_map21 = match_nodes(ts1, ts2, T2)
        nodes2 = list(node_map21.keys())
        nodes1 = list(node_map21.values())
        ts1, ts2 = reset_time(ts1.tables.tree_sequence(), ts2.tables.tree_sequence(), T1-T2)
        ts1s = ts1.simplify(nodes1)
        ts2s = ts2.simplify(nodes2)
        tables1s = ts1s.tables
        tables2s = ts2s.tables
        tables1s.provenances.clear()
        tables2s.provenances.clear()
        self.assertEqual(tables1s, tables2s)

    def test_simple_example(self):
        for (T1, T2) in [(100, 100), (100, 200), (200, 10)]:
            ts1, ts2 = get_examples(T1, T2)
            self.verify_match_nodes(ts1, ts2)
            self.verify_simplification_nodes(ts1, ts2)

class TestGraft(unittest.TestCase):

    def test_simple_example(self):
        ts1, ts2 = get_examples(35, 12, gens=4, N=4)
        T1, T2 = find_split_time(ts1, ts2)
        node_map21 = match_nodes(ts1, ts2, T2)

        # easier to deal with tskit treeseqs
        ts1 = ts1.tables.tree_sequence()
        ts2 = ts2.tables.tree_sequence()

        tsg, maps = graft(ts1, ts2, node_map21, T1, T2)
        node_map2new = maps[0]

        # resetting times so the trees are comparable
        ts1, ts2 = reset_time(ts1, ts2, T1-T2)

        # ts1 from grafted
        ts1g = tsg.simplify(ts1.samples())
        tables1g = ts1g.tables
        tables1 = ts1.simplify().tables

        # all tables but the provenance table should be the same
        tables1.provenances.clear()
        tables1g.provenances.clear()
        self.assertEqual(tables1g, tables1)

        # testing ts2
        samp = ts2.samples()
        new_samp = np.array([node_map2new[nid] for nid in list(samp)])
        tables2 = ts2.simplify(samp, filter_populations=False, filter_individuals=False).tables
        tables2g = tsg.simplify(new_samp, filter_populations=False, filter_individuals=False).tables
        # the test for the nodes table will be more complicated
        # because of the changes pop id, indiv ids
        # all tables but the provenance table should be the same
        tables2.provenances.clear()
        tables2g.provenances.clear()
        if tables2g.individuals != tables2.individuals:
            for i1, i2 in zip(tables2.individuals, tables2g.individuals):
                print(i1)
                print(i2)
        self.assertEqual(tables2g, tables2)

