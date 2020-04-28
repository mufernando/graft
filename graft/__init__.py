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
    return T1, T2

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

def add_time(ts, dt):
    '''
    This function returns a tskit.TreeSequence in which `dt`
    has been added to the times in all nodes.
    '''
    tables = ts.tables
    tables.nodes.set_columns(flags=tables.nodes.flags, time=tables.nodes.time+dt,population=tables.nodes.population,individual=tables.nodes.individual,metadata=tables.nodes.metadata,metadata_offset=tables.nodes.metadata_offset)
    return tables.tree_sequence()

def reset_time(ts1, ts2):
    '''
    Give two tskit.TreeSequences(), returns the respective tree
    sequences but now with correspondent times. That is, the times
    in the tree with the smallest `max_root_time` are shifted so
    that they are comparable to the other tree.
    '''
    time_diff = ts1.max_root_time - ts2.max_root_time
    if time_diff > 0:
        ts2 = add_time(ts2, time_diff)
    else:
        ts1 = add_time(ts1, abs(time_diff))
    return(ts1, ts2)

def graft(ts1, ts2, matched_nodes):
    """
    Returns a tree sequence obtained by grafting together the
    two tree sequences along the nodes in ``matched_nodes``,
    which should be a list of pairs of node IDs...
    More precisely, ts2 is grafted onto ts1. ``time_diff``
    indicates the time difference between ts1 and ts2. i.e.,
    the difference between the time of nodes that existed in
    both ts1 and ts2.
    """
    ts1, ts2 = reset_time(ts1, ts2)
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
            # translating the node ids from ts2 to new
            new_parent = all_nodes[e.parent==all_nodes[:,1]][0][0]
            new_child = all_nodes[e.child==all_nodes[:,1]][0][0]
            new_tables.edges.add_row(left = e.left, right = e.right, parent=new_parent, child=new_child)
    # grafting sites and muts
    new_muts = {}
    for k, m in enumerate(ts2.mutations()):
        if m.node in new_nodes[:,1]:
            # translating node id
            new_node = new_nodes[m.node==new_nodes[:,1]][0][0]
            # add site
            s = ts2.site(m.site)
            sid = new_tables.sites.add_row(position=s.position, ancestral_state="", metadata=s.metadata)
            # add mutation
            parent = m.parent
            # if parent was also new need to translate the id
            if m.parent in new_muts:
                parent = new_muts[m.parent]
            mid=new_tables.mutations.add_row(site=sid, node=new_node,derived_state=m.derived_state, parent=parent, metadata=m.metadata)
            new_muts[k] = mid
    new_tables.sort()
    new_tables.deduplicate_sites()
    return new_tables.tree_sequence(), all_nodes
