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

def reset_time(ts1, ts2, dT):
    '''Resets one of the trees so depending on difference in
    maximum time ago between ts1 and ts2, dT'''
    if dT > 0:
        ts2 = add_time(ts2, dT)
    elif dT < 0:
        ts1 = add_time(ts1, abs(dT))
    return (ts1, ts2)

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
        # simplifying with filter_populations=False bc SLiM can start
        # with id 1 instead of zero
        tables = ts.simplify(nodes, filter_populations=False).tables
        tablesg = tsg.simplify(nodesg, filter_populations=False).tables
        tables.provenances.clear()
        tablesg.provenances.clear()
        self.assertEqual(len(tables.nodes), len(tablesg.nodes))
        # simplify does not put individuals in order
        # but it does for nodes. that is why we check
        # for equality of individuals as follows
        for j in range(len(tables.nodes)):
            na = node_asdict(tables.nodes[j])
            nb = node_asdict(tablesg.nodes[j])
            if not (na['individual'] == nb['individual'] == -1):
                self.assertFalse((na['individual'] == -1) or (nb['individual'] == -1))
                ia = tables.individuals[na['individual']]
                ib = tablesg.individuals[nb['individual']]
                self.assertEqual(ia.flags, ib.flags)
                self.assertEqual(ia.metadata, ib.metadata)
                self.assertTrue((ia.location==ib.location).all())
            if not (na['population'] == nb['population'] == -1):
                pa = tables.populations[na['population']]
                pb = tablesg.populations[nb['population']]
                self.assertEqual(pa, pb)
            na['population']=nb['population']
            na['individual']=nb['individual']
            self.assertEqual(na, nb)
        tables.individuals.clear()
        tablesg.individuals.clear()
        tables.populations.clear()
        tablesg.populations.clear()
        tables.nodes.clear()
        tablesg.nodes.clear()
        self.assertEqual(tables, tablesg)

    def verify_nodes_pop(self, ts1, ts2, tsg, new_nodes):
        popg = [n.population for n in tsg.nodes()]
        popnew = [p for i,p in enumerate(popg) if i in new_nodes]
        pop1 = [n.population for n in ts1.nodes()]
        pop2 = [n.population for n in ts2.nodes()]
        # test new nodes are in new pops
        self.assertTrue(set(pop1)-set(popnew) == set(pop1))
        # test there are npop1+npop2 in grafted ts
        self.assertTrue(len(set(pop1))+len(set(pop2))==len(set(popg)))

    def test_simple_example(self):
        ts1, ts2 = get_examples(35, 12, gens=4, N=4)
        T1, T2 = find_split_time(ts1, ts2)
        node_map21 = match_nodes(ts1, ts2, T2)

        tsg, (node_map2new, pop_map2new, ind_map2new) = graft(ts1, ts2, node_map21)

        # resetting times so the trees are comparable
        ts1, ts2 = reset_time(ts1, ts2, T1-T2)

        # check that nodes added to ts1 were assigned new pop
        added_nodes = list(set(node_map2new)-set(node_map21))
        new_nodes = np.array([node_map2new[n] for n in added_nodes])
        self.verify_nodes_pop(ts1, ts2, tsg, new_nodes)

        # check that ts1 has not been changed by grafting
        self.verify_graft_simplification(ts1, tsg, node_map={n: n for n in ts1.samples()})

        # check that ts2 has not been changed by grafting
        full_sample_map = node_map2new.copy()
        # we'll need all the samples of ts2 in the node map:
        for n in ts2.samples():
            if n in node_map2new:
                pass
            else:
                full_sample_map[n] = n
        self.verify_graft_simplification(ts2, tsg, node_map=full_sample_map)
