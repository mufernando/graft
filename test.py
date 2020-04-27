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

# tests: graft one tree to another, simplify over just the nodes in original tree, test if the tables are the same
# test: take a grafted tree, simplify over samples in first and second tree, graft them back together and test if the same
