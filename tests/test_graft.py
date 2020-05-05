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

    def verify_graft_simplification(self, ts, tsg, node_map):
        # check that if we simplify tsg = graft(ts, ts2) back to ts.samples,
        # we get the same trees: after simplification everything but provenance
        # should be the same
        nodes = ts.samples()
        nodesg = [node_map[n] for n in nodes]
        tables = ts.simplify(nodes).tables
        # remove this - sanity check
        for j, n in enumerate(nodes):
            ia = ts.node(n).individual
            ib = tables.nodes.individual[j]
            assert(ts.tables.individuals[ia].metadata == tables.individuals[ib].metadata)
        tablesg = tsg.simplify(nodesg).tables
        # remove this - sanity check
        for j, n in enumerate(nodesg):
            ia = tsg.node(n).individual
            ib = tablesg.nodes.individual[j]
            assert(tsg.tables.individuals[ia].metadata == tablesg.individuals[ib].metadata)
        tables.provenances.clear()
        tablesg.provenances.clear()
        # remove this assert stuff - sanity check
        assert(tables.nodes.num_rows == tablesg.nodes.num_rows)
        for j in range(len(nodes)):
            ia = tables.nodes.individual[j]
            ib = tablesg.nodes.individual[j]
            assert(tables.individuals[ia].metadata == tablesg.individuals[ib].metadata)
        if tablesg.individuals != tables.individuals:
            for j, (i1, i2) in enumerate(zip(tables.individuals, tablesg.individuals)):
                print(j, pyslim.decode_individual(i1.metadata))
                print(j, pyslim.decode_individual(i2.metadata))
        self.assertEqual(tables, tablesg)


    def test_simple_example(self):
        ts1, ts2 = get_examples(35, 12, gens=4, N=4)
        T1, T2 = find_split_time(ts1, ts2)
        node_map21 = match_nodes(ts1, ts2, T2)

        tsg, (node_map2new, pop_map2new, ind_map2new) = graft(ts1, ts2, node_map21, T1, T2)

        # resetting times so the trees are comparable
        ts1, ts2 = reset_time(ts1, ts2, T1-T2)

        # check that ts1 has not been changed by grafting
        self.verify_graft_simplification(ts1, tsg, node_map={n: n for n in ts1.samples()})

        print(ind_map2new)
        print('-----')

        # check that ts2 has not been changed by grafting
        full_sample_map = node_map2new.copy()
        # we'll need all the samples of ts2 in the node map:
        for n in ts2.samples():
            if n in node_map2new:
                pass
            else:
                full_sample_map[n] = n
        self.verify_graft_simplification(ts2, tsg, node_map=full_sample_map)
