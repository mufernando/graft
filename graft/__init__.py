import numpy as np
import pyslim
import tskit
import json

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
    slim_gens1 = get_slim_gens(ts1)
    slim_gens2 = get_slim_gens(ts2)
    # counting SLiM prov before chains diff
    j = 0
    # finding the first diff between provenance chains
    for p1, p2 in zip(ts1.provenances(),ts2.provenances()):
        if p1 != p2:
            break
        record = json.loads(p1.record)
        if record["software"]["name"] == "SLiM":
            j += 1
    if j == 0:
        raise ValueError("No shared SLiM provenance entries.")
    last_slim_gen = ts1.slim_provenances[j-1].slim_generation
    T1 = abs(last_slim_gen - slim_gens1[-1])
    T2 = abs(last_slim_gen - slim_gens2[-1])
    return T1, T2

def match_nodes(ts1, ts2, T2=0):
    """
    Given two SLiM tree sequences, returns a dictionary relating
    the id in ts2 (key) to id in ts1 (item) for  node IDs in the
    two tree sequences that refer to the same node. If split time
    in ts2 (T2) is given, then only nodes before the split are
    considered. Note the only check of equivalency is the slim_id
    of the nodes.
    """
    node_map21 = {}
    # getting slim ids for ts1
    slim_ids1 = np.array([n.metadata.slim_id for n in ts1.nodes()])
    # looping through nodes in ts2 and finding equiv in ts1
    for k, n in enumerate(ts2.nodes()):
        # nodes should be equiv only at split or before it happened
        if n.time >= T2:
            # rn just checking slim ids, but would be good
            # to compare other features (like the subtree
            # simplifying on that node?)
            nid = np.where(n.metadata.slim_id == slim_ids1)[0]
            assert len(nid)== 1
            node_map21[k] = nid[0]
    return node_map21

def add_time(ts, dt):
    '''
    This function returns a tskit.TreeSequence in which `dt`
    has been added to the times in all nodes.
    '''
    tables = ts.tables
    nodes_dict = tables.nodes.asdict()
    nodes_dict['time'] = nodes_dict['time'] + dt
    tables.nodes.set_columns(**nodes_dict)
    return tables.tree_sequence()

def _check_shared_nodes(ts1, ts2, node_map21):
    '''
    Given two tree sequences with shared nodes as described in
    `node_map21`, test whether simplifying on those nodes gives
    you the same tree sequences.
    '''
    nodes2 = list(node_map21.keys())
    nodes1 = list(node_map21.values())
    ts1s = ts1.simplify(nodes1)
    ts2s = ts2.simplify(nodes2)
    tables1s = ts1s.tables
    tables2s = ts2s.tables
    tables1s.provenances.clear()
    tables2s.provenances.clear()
    assert tables1s == tables2s

def graft(ts1, ts2, node_map21):
    """
    Returns a tree sequence obtained by grafting together the
    two tree sequences along the nodes in ``node_map21``,
    which should be a dictionary mapping nodes in ts2 that are
    equivalent to nodes in ts1.
    More precisely, ts2 is grafted onto ts1.
    Populations of nodes new to ts1 are considered new in the
    grafted tree sequence. Map is returned.
    T1 and T2 are used to shift the time in the tree seqs.
    It is used in cases where the after split portion of
    are run for different number of generations in each of ts1
    and ts2. If this is not the case set T1=T2=0
    """
    # making sure the ts are tskit ts
    ts1, ts2 = ts1.tables.tree_sequence(), ts2.tables.tree_sequence()
    # checking shift in time ago between ts1 and ts2
    dt = [ts1.node(node_map21[a]).time - ts2.node(a).time for a in node_map21]
    if len(set(dt)) > 1:
        raise ValueError("Inconsistent time differences among the equivalent nodes.")
    dT = int(dt[0])
    if dT > 0:
        ts2 = add_time(ts2, dT)
    elif dT < 0:
        ts1 = add_time(ts1, abs(dT))
    # checking the trees are the same below the nodes_map21
    _check_shared_nodes(ts1, ts2, node_map21)
    # the grafted tree will bbe based off of ts1
    new_tables = ts1.tables
    # mapping nodes in ts2 to new nodes in the grafted tables
    node_map2new = {}
    node_map2new.update(node_map21)
    # mapping of individuals in ts2 to new
    ind_map2new = {}
    # adding the pops in ts2 to new
    pop_map2new = {}
    # need to loop throgh nodes in ts2, find the unmatched nodes
    # and add them to the ts1 node table
    for k, n in enumerate(ts2.nodes()):
        if not k in node_map21:
            if not n.population in pop_map2new:
                pop = ts2.tables.populations[n.population]
                pid = new_tables.populations.add_row(pop.metadata)
                pop_map2new[n.population]=pid
            # translating pop to new
            n.population = pop_map2new[n.population]
            # adding individual
            if n.individual >= 0:
                ind = ts2.individual(n.individual)
                iid = new_tables.individuals.add_row(flags=ind.flags, location=ind.location, metadata=ind.metadata)
                ind_map2new[n.individual] = iid
                n.individual=iid
            # addingn node
            nid = new_tables.nodes.add_row(**node_asdict(n))
            node_map2new[k] = nid
    # creating a set with nodes that are new to ts1
    new_nodes=set(node_map2new)-set(node_map21)
    # now we need to add the edges
    for i, e in enumerate(ts2.edges()):
        if (e.parent in new_nodes) or (e.child in new_nodes):
            if (e.parent in new_nodes) and (not e.child in new_nodes):
                raise ValueError("Cannot graft nodes above existing nodes.")
            # translating the node ids from ts2 to new
            new_parent = node_map2new[e.parent]
            new_child = node_map2new[e.child]
            new_tables.edges.add_row(left = e.left, right = e.right, parent=new_parent, child=new_child)
    # grafting sites and muts
    new_muts = {}
    for k, m in enumerate(ts2.mutations()):
        if m.node in new_nodes:
            # translating node id
            new_node = node_map2new[m.node]
            # add site
            s = ts2.site(m.site)
            sid = new_tables.sites.add_row(position=s.position, ancestral_state="", metadata=s.metadata)
            mid=new_tables.mutations.add_row(site=sid, node=new_node,derived_state=m.derived_state, parent=tskit.NULL, metadata=m.metadata)
            new_muts[k] = mid
    new_tables.sort()
    new_tables.deduplicate_sites()
    new_tables.build_index()
    new_tables.compute_mutation_parents()
    return new_tables.tree_sequence(), (node_map2new, pop_map2new, ind_map2new)
