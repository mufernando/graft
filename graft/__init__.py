import numpy as np
import pyslim
import tskit

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

def get_provenance(ts, only_last=True):
    '''
    Extracts model type, slim generation, and remembmered node count from either
    the last entry in the provenance table that is tagged with "program"="SLiM"
    (if ``only_last=True``) or a list of all of them (otherwise).

    :param SlimTreeSequence ts: The tree sequence.
    :param bool only_last: Whether to return only the last SLiM provenance entry,
        (otherwise, returns a list of all SLiM entries).
    :rtype ProvenanceMetadata:
    '''
    provenances = []
    for j, p in enumerate(ts.tables.provenances):
        this_record = json.loads(p.record)
        is_slim, this_file_version = _slim_provenance_version(this_record)
        if is_slim:
            record = this_record
            file_version = this_file_version
            if file_version == "0.1":
                out = ProvenanceMetadata(record['model_type'],
                                         record['generation'],
                                         file_version)
            else: # >= 0.2
                out = ProvenanceMetadata(record['parameters']['model_type'],
                                         record['slim']["generation"],
                                         file_version)
            provenances.append(out)

    if len(provenances) == 0:
        raise ValueError("Tree sequence contains no SLiM provenance entries"
                          "(or your pyslim is out of date).")
    if only_last:
        return provenances[-1]
    else:
        return provenances

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

def graft(ts1, ts2, node_map21):
    """
    Returns a tree sequence obtained by grafting together the
    two tree sequences along the nodes in ``node_map21``,
    which should be a dictionary mapping nodes in ts2 that are
    equivalent to nodes in ts1.
    More precisely, ts2 is grafted onto ts1.
    """
    ts1, ts2 = reset_time(ts1, ts2)
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
            else:
                n.population = pop_map2new[n.population]
            # adding individual
            if n.individual > 0:
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
            if (e.parent in new_nodes):
                assert e.child in new_nodes, "Parent is new but child is not."
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
