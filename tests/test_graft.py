import os
import unittest
import msprime
import pyslim
from graft import *


def run_slim_script(slimfile, args=''):
    print("running:: " + "slim -s 23 " + args + " " + slimfile)
    out = os.system("slim -s 23 " + args + " " + slimfile + ">/dev/null")
    return out


def get_slim_examples(
        T1,
        T2,
        gens=100,
        N=100,
        recipe_path="tests/recipe.slim"):
    print(os.getcwd())
    # note: could make N and gens parameters here
    run_slim_script(
        recipe_path,
        args=f"-d N={N} -d gens={gens} -d \"outfile='tests/data/root.trees'\"")
    run_slim_script(
        recipe_path,
        args=f"-d N={N} -d gens={T1} -d \"infile='tests/data/root.trees'\" -d 'outfile=\"tests/data/branch1.trees\"'")
    run_slim_script(
        recipe_path,
        args=f"-d N={N} -d gens={T2} -d \"infile='tests/data/root.trees'\" -d 'outfile=\"tests/data/branch2.trees\"'")
    ts1 = pyslim.load("tests/data/branch1.trees")
    ts2 = pyslim.load("tests/data/branch2.trees")
    return ts1, ts2


def get_msprime_mig_example(T=100, t=10, N=100, n=10):
    # we assume after the T the ts are completely independent
    # t is used to set the split within the two independent pops
    M = [
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0]
    ]
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=n),
        msprime.PopulationConfiguration(sample_size=n),
        msprime.PopulationConfiguration(sample_size=n),
        msprime.PopulationConfiguration(sample_size=n)
    ]
    demographic_events = [
        msprime.MassMigration(t, source=2, dest=0, proportion=1),
        msprime.MassMigration(t, source=3, dest=1, proportion=1),
        msprime.MassMigration(T, source=1, dest=0, proportion=1),
        msprime.CensusEvent(time=T)
    ]
    # dd = msprime.DemographyDebugger(
    #   population_configurations=population_configurations,
    #   migration_matrix=M,
    #   demographic_events=demographic_events)
    # dd.print_history()
    ts = msprime.simulate(
        Ne=N,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=M,
        length=2e4,
        recombination_rate=1e-8,
        mutation_rate=1e-8,
        record_migrations=True)
    return ts


def get_msprime_example(T=100, N=100, n=10):
    # we assume after the split the ts are completely independent
    M = [[0, 0], [0, 0]]
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=n),
        msprime.PopulationConfiguration(sample_size=n)
    ]
    demographic_events = [
        msprime.CensusEvent(time=T),
        msprime.MassMigration(T, source=1, dest=0, proportion=1)
    ]
    # dd = msprime.DemographyDebugger(
    #    population_configurations=population_configurations,
    #    migration_matrix=M,
    #    demographic_events=demographic_events)
    # dd.print_history()
    ts = msprime.simulate(
        Ne=N,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=M,
        length=2e4,
        recombination_rate=1e-8,
        mutation_rate=1e-8)
    return ts


def node_asdict(node):
    return {
        "time": node.time,
        "population": node.population,
        "individual": node.individual,
        "metadata": node.metadata,
        "flags": node.flags
    }


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
            ts1, ts2 = get_slim_examples(T1, T2)
            S1, S2 = find_split_time(ts1, ts2)
            self.assertEqual(S1, T1)
            self.assertEqual(S2, T2)


class TestMatchNodes(unittest.TestCase):

    def verify_match_nodes(self, ts1, ts2):
        T1, T2 = find_split_time(ts1, ts2)
        node_map21 = match_nodes(ts1, ts2, T2)
        # note: should also check that we got all the nodes
        for a, b in node_map21.items():
            n1 = ts1.node(b)
            n2 = ts2.node(a)
            self.assertGreaterEqual(n2.time, T2)
            self.assertEqual(n1.metadata.slim_id, n2.metadata.slim_id)

    def verify_simplification_nodes(self, ts1, ts2):
        T1, T2 = find_split_time(ts1, ts2)
        node_map21 = match_nodes(ts1, ts2, T2)
        nodes2 = list(node_map21.keys())
        nodes1 = list(node_map21.values())
        ts1, ts2 = reset_time(ts1.tables.tree_sequence(),
                              ts2.tables.tree_sequence(), T1 - T2)
        ts1s = ts1.simplify(nodes1)
        ts2s = ts2.simplify(nodes2)
        tables1s = ts1s.tables
        tables2s = ts2s.tables
        tables1s.provenances.clear()
        tables2s.provenances.clear()
        self.assertEqual(tables1s, tables2s)

    def test_simple_example(self):
        for (T1, T2) in [(100, 100), (100, 200), (200, 10)]:
            ts1, ts2 = get_slim_examples(T1, T2)
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
        # TODO: a test for provenances
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
                self.assertFalse(
                    (na['individual'] == -
                     1) or (
                        nb['individual'] == -
                        1))
                ia = tables.individuals[na['individual']]
                ib = tablesg.individuals[nb['individual']]
                self.assertEqual(ia.flags, ib.flags)
                self.assertEqual(ia.metadata, ib.metadata)
                self.assertTrue((ia.location == ib.location).all())
            if not (na['population'] == nb['population'] == -1):
                pa = tables.populations[na['population']]
                pb = tablesg.populations[nb['population']]
                self.assertEqual(pa, pb)
            na['population'] = nb['population']
            na['individual'] = nb['individual']
            self.assertEqual(na, nb)
        tables.individuals.clear()
        tablesg.individuals.clear()
        tables.populations.clear()
        tablesg.populations.clear()
        tables.nodes.clear()
        tablesg.nodes.clear()
        print(tables)
        print(tablesg)
        self.assertEqual(tables, tablesg)

    def verify_node_populations(self, ts1, ts2, node_map21):
        tsg, (node_map2new, pop_map2new, ind_map2new) = graft(
            ts1, ts2, node_map21)
        # check that ts1 pops remain unchanged
        for j, p in enumerate(ts1.populations()):
            self.assertEqual(p, tsg.population(j))
        # check that each ts2 pop maps to a unique tsg pop
        self.assertEqual(len(pop_map2new.keys()),
                         len(set(pop_map2new.values())))
        # check that ts2 pops are unchanged, and that newly added nodes have a
        # new population
        for n2 in range(ts2.num_nodes):
            p2 = ts2.node(n2).population
            ng = node_map2new[n2]
            pg = tsg.node(ng).population
            self.assertEqual(
                ts2.tables.populations[p2],
                tsg.tables.populations[pg])
            if n2 not in node_map21:
                self.assertGreaterEqual(pg, ts1.num_populations)

    def test_msprime_example(self):
        T = 100
        ts = get_msprime_example(T, 50, 2)
        #ts = get_msprime_mig_example(30, t=10, N=100, n=2)
        # simpifying to get things in the right order
        shared_nodes = [n.id for n in ts.nodes() if n.time >= T]
        pop1 = list(ts.samples(population=0))
        pop2 = list(ts.samples(population=1))
        ts1_samples = shared_nodes + pop1
        ts2_samples = shared_nodes + pop2
        assert len(ts1_samples) == len(ts2_samples)
        node_map21 = {i: i for i in range(len(shared_nodes))}
        ts1 = ts.simplify(ts1_samples)
        ts2 = ts.simplify(ts2_samples)
        tsg, (node_map2new, pop_map2new, ind_map2new) = graft(
            ts1, ts2, node_map21)

        # check that nodes added to ts1 were assigned new pop
        self.verify_node_populations(ts1, ts2, node_map21)

        # check that ts1 has not been changed by grafting
        self.verify_graft_simplification(
            ts1, tsg, node_map={
                n: n for n in ts1.samples()})

        # check that ts2 has not been changed by grafting
        full_sample_map = node_map2new.copy()
        # we'll need all the samples of ts2 in the node map:
        for n in ts2.samples():
            if n in node_map2new:
                pass
            else:
                full_sample_map[n] = n
        self.verify_graft_simplification(ts2, tsg, node_map=full_sample_map)

    def test_slim_nonwf_example(self):
        ts1, ts2 = get_slim_examples(
            10, 10, gens=100, N=100, recipe_path="tests/recipe_nonwf1.slim")
        T1, T2 = find_split_time(ts1, ts2)
        node_map21 = match_nodes(ts1, ts2, T2)

        tsg, (node_map2new, pop_map2new, ind_map2new) = graft(
            ts1, ts2, node_map21)

    def test_slim_example(self):
        ts1, ts2 = get_slim_examples(15, 30, gens=100, N=100)
        T1, T2 = find_split_time(ts1, ts2)
        node_map21 = match_nodes(ts1, ts2, T2)
        for k in node_map21:
            assert(k == node_map21[k])
        tsg, (node_map2new, pop_map2new, ind_map2new) = graft(
            ts1, ts2, node_map21)

        # resetting times so the trees are comparable
        ts1, ts2 = ts1.tables.tree_sequence(), ts2.tables.tree_sequence()
        ts1, ts2 = reset_time(ts1, ts2, T1 - T2)

        # check that nodes added to ts1 were assigned new pop
        self.verify_node_populations(ts1, ts2, node_map21)

        # check that ts1 has not been changed by grafting
        self.verify_graft_simplification(
            ts1, tsg, node_map={
                n: n for n in ts1.samples()})

        # check that ts2 has not been changed by grafting
        full_sample_map = node_map2new.copy()
        # we'll need all the samples of ts2 in the node map:
        for n in ts2.samples():
            if n in node_map2new:
                pass
            else:
                full_sample_map[n] = n
        self.verify_graft_simplification(ts2, tsg, node_map=full_sample_map)
