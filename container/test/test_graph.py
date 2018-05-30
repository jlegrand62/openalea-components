import random
import numpy as np
from nose import with_setup
from openalea.container import Graph

g = Graph()


def setup_func():
    for i in xrange(10):
        g.add_vertex(i)
    for i in xrange(9):
        g.add_edge(i, i + 1, i)


def teardown_func():
    g.clear()


# ##########################################################
#
# Graph concept
#
# ##########################################################
@with_setup(setup_func, teardown_func)
def test_source():
    for i in xrange(9):
        assert g.source(i) == i


@with_setup(setup_func, teardown_func)
def test_target():
    for i in xrange(9):
        assert g.target(i) == (i + 1)


@with_setup(setup_func, teardown_func)
def test_has_vertex():
    for i in xrange(10):
        assert g.has_vertex(i)


@with_setup(setup_func, teardown_func)
def test_has_edge():
    for i in xrange(9):
        assert g.has_edge(i)


@with_setup(setup_func, teardown_func)
def test_is_valid():
    assert g.is_valid()


# ##########################################################
#
# Vertex List Graph Concept
#
# ##########################################################
@with_setup(setup_func, teardown_func)
def test_vertices():
    assert list(g.vertices()) == range(10)


@with_setup(setup_func, teardown_func)
def test_nb_vertices():
    assert g.nb_vertices() == 10


@with_setup(setup_func, teardown_func)
def test_in_neighbors():
    for i in xrange(9):
        assert list(g.in_neighbors(i + 1)) == [i]


@with_setup(setup_func, teardown_func)
def test_out_neighbors():
    for i in xrange(9):
        assert list(g.out_neighbors(i)) == [i + 1]


@with_setup(setup_func, teardown_func)
def test_neighbors():
    for i in xrange(8):
        neis = list(g.neighbors(i + 1))
        assert i in neis
        assert i + 2 in neis


@with_setup(setup_func, teardown_func)
def test_nb_in_neighbors():
    for i in xrange(9):
        assert g.nb_in_neighbors(i + 1) == 1


@with_setup(setup_func, teardown_func)
def test_nb_out_neighbors():
    for i in xrange(9):
        assert g.nb_out_neighbors(i) == 1


@with_setup(setup_func, teardown_func)
def test_nb_neighbors():
    for i in xrange(8):
        assert g.nb_neighbors(i + 1) == 2


@with_setup(setup_func, teardown_func)
def test_edge():
    assert g.edge(0, 1) == 0
    assert g.edge(0, 2) is None


# ##########################################################
#
# Edge List Graph Concept
#
# ##########################################################
@with_setup(setup_func, teardown_func)
def test_edges():
    assert list(g.edges()) == range(9)


@with_setup(setup_func, teardown_func)
def test_nb_edges():
    assert g.nb_edges() == 9


@with_setup(setup_func, teardown_func)
def test_in_edges():
    for i in xrange(9):
        assert list(g.in_edges(i + 1)) == [i]


@with_setup(setup_func, teardown_func)
def test_out_edges():
    for i in xrange(9):
        assert list(g.out_edges(i)) == [i]


@with_setup(setup_func, teardown_func)
def test_vertex_edges():
    for i in xrange(8):
        neis = list(g.edges(i + 1))
        assert i in neis
        assert i + 1 in neis


@with_setup(setup_func, teardown_func)
def test_nb_in_edges():
    for i in xrange(9):
        assert g.nb_in_edges(i + 1) == 1


@with_setup(setup_func, teardown_func)
def test_nb_out_edges():
    for i in xrange(9):
        assert g.nb_out_edges(i) == 1


@with_setup(setup_func, teardown_func)
def test_nb_edges():
    for i in xrange(8):
        assert g.nb_edges(i + 1) == 2


# ##########################################################
#
# Mutable Vertex Graph concept
#
# ##########################################################
@with_setup(setup_func, teardown_func)
def test_add_vertex():
    assert g.add_vertex(100) == 100
    vid = g.add_vertex()
    assert g.has_vertex(vid)


@with_setup(setup_func, teardown_func)
def test_remove_vertex():
    g.remove_vertex(5)
    assert not g.has_vertex(5)
    assert not g.has_edge(4)
    assert not g.has_edge(5)
    assert 5 not in list(g.neighbors(6))
    assert 5 not in list(g.neighbors(4))


@with_setup(setup_func, teardown_func)
def test_clear():
    g.clear()
    assert g.nb_vertices() == 0
    assert g.nb_edges() == 0


# ##########################################################
#
# Mutable Edge Graph concept
#
# ##########################################################
@with_setup(setup_func, teardown_func)
def test_add_edge():
    assert g.add_edge(0, 9, 100) == 100
    eid = g.add_edge(2, 1)
    assert eid in list(g.in_edges(1))
    assert eid in list(g.out_edges(2))


@with_setup(setup_func, teardown_func)
def test_remove_edge():
    g.remove_edge(4)
    assert not g.has_edge(4)
    assert 4 not in list(g.neighbors(5))
    assert 5 not in list(g.neighbors(4))


@with_setup(setup_func, teardown_func)
def test_clear_edges():
    g.clear_edges()
    assert g.nb_vertices() == 10
    assert g.nb_edges() == 0


# ##########################################################
#
# Extend Graph concept
#
# ##########################################################
@with_setup(setup_func, teardown_func)
def test_extend():
    trans_vid, trans_eid = g.extend(g)
    assert len(trans_vid) == 10
    assert len(trans_eid) == 9


# ##########################################################
#
# Spatial distance
#
# ##########################################################

# TOPOLOGICAL DISTANCES between vids:
################################################################################
def create_graph(n_vids, n_eids):
    """
    Create a PropertyGraph with 'n_vids' vertex ids and 'n_eids' edge ids.
    Edges are linking vid_i with vid_i+1

    Parameters
    ----------
    n_vids : int
        number of vertex to add to the graph
    n_eids : int
        number of edge to add to the graph

    Returns
    -------
    the PropertyGraph
    """
    # Add vertices and save the list of ids
    vids = []
    for i in xrange(n_vids):
        vids.append(g.add_vertex(i))
    # Add edges and save the list of ids
    eids = []
    for i in xrange(n_eids):
        eids.append(g.add_edge(i, i + 1, i))

    return g


g = Graph()

n_vids = 10
n_eids = n_vids - 1


def setup_func():
    g = create_graph(n_vids, n_eids)


@with_setup(setup_func, teardown_func)
def test_topological_distance():
    topo = g.topological_distance(0)
    assert len(topo) == n_vids
    for i in xrange(n_vids):
        assert topo[i] == i


@with_setup(setup_func, teardown_func)
def test_topological_distance2():
    topo = g.topological_distance(0, lambda x, y: 2)
    assert len(topo) == n_vids
    for i in xrange(n_vids):
        assert topo[i] == i * 2


@with_setup(setup_func, teardown_func)
def test_topological_distance3():
    max_depth = 2
    topo = g.topological_distance(0, max_depth=max_depth)
    assert len(topo) == n_vids
    for i in xrange(max_depth):
        assert topo[i] == i
    for i in xrange(max_depth + 1, n_vids):
        assert topo[i] == float('inf')
    topo = g.topological_distance(0, max_depth=max_depth, return_inf=False)
    for i in xrange(max_depth + 1, n_vids):
        assert np.isnan(topo[i])


# ADJACENCY MATRIX between all vertices:
################################################################################
edge_dist = 1
no_edge_val = 0
reflexive_value = 0.01


@with_setup(setup_func, teardown_func)
def test_adjacency_matrix_non_oriented_non_reflexive():
    # init the TEST adjacency matrix
    if no_edge_val == 0:
        m = np.zeros([n_vids, n_vids])
    else:
        m = np.ones([n_vids, n_vids]) * no_edge_val
    # create the TEST adjacency matrix
    for e in g.edges():
        s, t = g.edge_vertices(e)
        m[s, t] = edge_dist
    # get the graph.adjacency_matrix
    adj = g.adjacency_matrix(edge_dist, no_edge_val, oriented=False,
                             reflexive=False, reflexive_value=0)
    assert np.alltrue(adj == m)


@with_setup(setup_func, teardown_func)
def test_adjacency_matrix_oriented_non_reflexive():
    # init the TEST adjacency matrix
    if no_edge_val == 0:
        m = np.zeros([n_vids, n_vids])
    else:
        m = np.ones([n_vids, n_vids]) * no_edge_val
    # create the TEST adjacency matrix
    for e in g.edges():
        s, t = g.edge_vertices(e)
        m[s, t] = m[t, s] = edge_dist
    # get the graph.adjacency_matrix
    adj = g.adjacency_matrix(edge_dist, no_edge_val, oriented=True,
                             reflexive=False, reflexive_value=0)
    assert np.alltrue(adj == m)


@with_setup(setup_func, teardown_func)
def test_adjacency_matrix_oriented_reflexive():
    # init the TEST adjacency matrix
    if no_edge_val == 0:
        m = np.zeros([n_vids, n_vids])
    else:
        m = np.ones([n_vids, n_vids]) * no_edge_val
    # create the TEST adjacency matrix
    for e in g.edges():
        s, t = g.edge_vertices(e)
        m[s, t] = m[t, s] = edge_dist
    for i in xrange(n_vids):
        m[i, i] = reflexive_value
    # get the graph.adjacency_matrix
    adj = g.adjacency_matrix(edge_dist, no_edge_val, oriented=True,
                             reflexive=True, reflexive_value=reflexive_value)
    assert np.alltrue(adj == m)


# FLOYD-WARSHALL DISTANCE MATRIX between all vertices:
################################################################################
@with_setup(setup_func, teardown_func)
def test_floyd_warshall_non_oriented_non_reflexive():
    m = np.zeros([n_vids, n_vids])

    for t in xrange(n_vids):
        for s in xrange(n_vids):
            m[s, t] = t - s
            m[t, s] = 'inf'
    for i in xrange(n_vids):
        m[i, i] = 'inf'

    fw = g.floyd_warshall(edge_dist, oriented=False, reflexive=False,
                          reflexive_value=0)
    assert np.alltrue(fw == m)
    # Add edge between two non-neighbors vertices, assert they are now "close":
    s, t = 0, 0
    while abs(t - s) <= 1:
        s, t = random.sample(range(n_vids), 2)
    g.add_edge(s, t)
    fw = g.floyd_warshall(edge_dist, oriented=False, reflexive=False,
                          reflexive_value=0)
    assert fw[s, t] == 1
    assert fw[t, s] == float('inf')


@with_setup(setup_func, teardown_func)
def test_floyd_warshall_oriented_non_reflexive():
    m = np.zeros([n_vids, n_vids])

    for t in xrange(n_vids):
        for s in xrange(n_vids):
            m[s, t] = t - s
            m[t, s] = t - s
    for i in xrange(n_vids):
        m[i, i] = 'inf'

    fw = g.floyd_warshall(edge_dist, oriented=True, reflexive=False,
                          reflexive_value=0)
    assert np.alltrue(fw == m)
    # Add edge between two non-neighbors vertices, assert they are now "close":
    s, t = 0, 0
    while abs(t - s) <= 1:
        s, t = random.sample(range(n_vids), 2)
    g.add_edge(s, t)
    fw = g.floyd_warshall(edge_dist, oriented=True, reflexive=False,
                          reflexive_value=0)
    assert fw[s, t] == 1
    assert fw[t, s] == t - s


@with_setup(setup_func, teardown_func)
def test_floyd_warshall_oriented_reflexive():
    m = np.zeros([n_vids, n_vids])

    for t in xrange(n_vids):
        for s in xrange(n_vids):
            m[s, t] = t - s
            m[t, s] = t - s
    for i in xrange(n_vids):
        m[i, i] = reflexive_value

    fw = g.floyd_warshall(edge_dist, oriented=True, reflexive=True,
                          reflexive_value=reflexive_value)
    assert np.alltrue(fw == m)

# TODO: test 'to_networkx' & 'from_networkx'
