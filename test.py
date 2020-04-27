import numpy as np
import pyslim
from init import *

ts1=pyslim.load("data/branch1.trees")
ts2=pyslim.load("data/branch2.trees")


# first, let's find the split times between the two trees
T1, T2 = find_split_time(ts1,ts2)

# now we can match the node ids
matched_nodes = match_nodes(ts1,ts2,T2)

# mapping parent and child nodes
map_parent = map2new[ts2.tables.edges.parent]
map_child = map2new[ts2.tables.edges.child]
# not sure what this is doing
map_edges = np.logical_and(map_parent>0, map_child>0)
new_tables.edges.append_columns(left=ts2.tables.edges.left[map_edges])
new_tables.edges.append_columns(right=ts2.tables.edges.right[map_edges])
new_tables.edges.append_columns(parent=map_parent[map_edges])
new_tables.edges.append_columns(child=map_child[map_edges])

# tests: graft one tree to another, simplify over just the nodes in original tree, test if the tables are the same
# test: take a grafted tree, simplify over samples in first and second tree, graft them back together and test if the same
