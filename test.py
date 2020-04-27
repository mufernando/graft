import numpy as np
import pyslim
import tskit
from init import *

ts1=pyslim.load("data/branch1.trees")
ts2=pyslim.load("data/branch2.trees")


# first, let's find the split times between the two trees
T1, T2 = find_split_time(ts1,ts2)

# now we can match the node ids
matched_nodes = match_nodes(ts1,ts2,T2)

# easier to deal with tskit treeseqs
ts1 = ts1.tables.tree_sequence()
ts2 = ts2.tables.tree_sequence()

tsg, all_nodes = graft(ts1,ts2, matched_nodes)

# ts1 from grafted
ts1g = tsg.simplify(ts1.samples())
tables1g = ts1g.tables
tables1 = ts1.simplify().tables
assert tables1g.nodes == tables1.nodes
assert tables1g.edges == tables1.edges

#testing ts2
samp = ts2.samples()
new_samp = np.array([all_nodes[nid==all_nodes[:,1]][0][0] for nid in list(samp)])
tables2 = ts2.simplify().tables
tables2g = tsg.simplify(new_samp).tables
assert tables2g.nodes == tables2.nodes
assert tables2g.edges == tables2.edges
# tests: graft one tree to another, simplify over just the nodes in original tree, test if the tables are the same
# test: take a grafted tree, simplify over samples in first and second tree, graft them back together and test if the same
