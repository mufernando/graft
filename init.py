import numpy as np
import pyslim

def get_slim_gens(ts):
    return np.array([p.slim_generation for p in ts.slim_provenances])

def node_asdict(node):
    return {
	"time" : node.time,
	"population" : node.population,
	"individual" : node.individual,
	"metadata" : node.metadata,
	"flags" : node.flags
    }

def find_split_time(ts1,ts2):
    """
    Given two SLiM tree sequences with shared history, this
    function returns the split times (in time ago) for each tree.
    """
    T1 = 0
    T2 = 0
    slim_gens1 = get_slim_gens(ts1)
    slim_gens2 = get_slim_gens(ts2)
    for i, (p1, p2) in enumerate(zip(ts1.provenances(),ts2.provenances())):
        if p1 != p2:
            assert i > 0, "No shared history between the two trees."
            T1 = abs(last_gen-slim_gens1[-1])
            T2 = abs(last_gen-slim_gens2[-1])
            break
        last_gen = ts1.slim_provenances[i].slim_generation
    return T1,T2

def match_nodes(ts1, ts2, T2=0):
    """
    Given two SLiM tree sequences, returns a list of pairs of node ids,
    [id1, id2], so that id1 and id2 are the node IDs in the two tree
    sequences that refer to the same node. If split time in ts3 (T2)
    is given, then only nodes before the split are considered
    """
    matched_nodes = []
    # getting slim ids for ts1
    sids_1 = np.array([n.metadata.slim_id for n in ts1.nodes()])
    # looping through nodes in ts2 and finding equiv in ts1
    for k, n in enumerate(ts2.nodes()):
        #print(k)
        # nodes should be equiv only at split or before it happened
        if n.time >= T2:
            # rn just checking slim ids, but would be good
            # to compare other features (like the subtree
            # simplifying on that node?)
            nid = np.where(n.metadata.slim_id == sids_1)[0]
            assert len(nid)== 1
            matched_nodes.append([nid[0],k])
    return np.array(matched_nodes)

def graft(ts1, ts2, matched_nodes):
    """
    Returns a tree sequence obtained by grafting together the
    two tree sequences along the nodes in ``matched_nodes``,
    which should be a list of pairs of node IDs...
    More precisely, ts2 is grafted onto ts1.
    """
    new_tables = ts1.tables
    new_nodes = []
    # need to loop throgh nodes in ts2, find the unmatched nodes
    # and add them to the ts1 node table
    for k, n in enumerate(ts2.nodes()):
        if not k in matched_nodes[:,1]:
            nid = new_tables.nodes.add_row(**node_asdict(n))
            new_nodes.append([nid,k])
    new_nodes=np.array(new_nodes)
    all_nodes=np.concatenate((new_nodes, matched_nodes), axis=0)
    # now we need to add the edges
    for i, e in enumerate(ts2.edges()):
        if (e.parent in new_nodes[:,1]) or e.child in new_nodes[:,1]:
            new_parent = all_nodes[e.parent==all_nodes[:,1]][0][0]
            new_child = all_nodes[e.child==all_nodes[:,1]][0][0]
            new_tables.edges.add_row(left = e.left, right = e.right, parent=new_parent, child=new_child)
    new_tables.sort()
    return new_tables.tree_sequence(), all_nodes
