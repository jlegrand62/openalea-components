__license__ = "Cecill-C"
__revision__ = " $Id$ "

# Test node module
# from openalea.core.graph.property_graph import PropertyGraph
from openalea.container import PropertyGraph

from nose import with_setup


g = PropertyGraph()
n_vids = 10
n_eids = n_vids - 1


def test_init():
    assert hasattr(g, 'metavidtypepropertyname')
    assert hasattr(g, 'metaeidtypepropertyname')
    assert g.metavidtypepropertyname == "valueproperty_as_vid"
    assert g.metaeidtypepropertyname == "valueproperty_as_eid"


def test_add_vtx():
    # Add vertices and save the list of ids
    vids = []
    for i in xrange(n_vids):
        vids.append(g.add_vertex(i))
    assert len(vids) == n_vids
    assert len(g._vertices) == n_vids


def test_add_edge():
    # Add edges and save the list of ids
    eids = []
    for i in xrange(n_eids):
        eids.append(g.add_edge(i, i + 1, i))
    assert len(eids) == n_eids
    assert len(g._edges) == n_eids


def test_add_empty_vtx_ppty():
    g.add_vertex_property('squared_id')
    assert 'squared_id' in g.vertex_property_names()
    g.remove_vertex_property('squared_id')
    assert 'squared_id' not in g.vertex_property_names()


def test_add_empty_edge_ppty():
    g.add_edge_property('sum_vids')
    assert 'sum_vids' in g.edge_property_names()
    g.remove_edge_property('sum_vids')
    assert 'sum_vids' not in g.edge_property_names()


def test_add_vtx_ppty():
    vids = g._vertices
    assert 'squared_id' not in g.vertex_property_names()
    # Add vertex property: the value is the squared value of the ids
    v = dict([(k, k ** 2) for k in vids])
    g.add_vertex_property('squared_id', v)
    assert 'squared_id' in g.vertex_property_names()
    assert len(g.vertex_property('squared_id')) == n_vids


def test_add_edge_ppty():
    eids = g._edges
    assert 'sum_vids' not in g.edge_property_names()
    # Add edge property: the value is the squared value of the ids
    v = dict([(k, sum(g.edge_vertices(k))) for k in eids])
    g.add_edge_property('sum_vids', v)
    assert 'sum_vids' in g.edge_property_names()
    assert len(g.edge_property('sum_vids')) == n_eids


def test_add_graph_ppty():
    # Add graph property:
    g.add_graph_property('name', "my_property_graph")
    g.add_graph_property('units', {'volume': 'real', 'barycenter': 'voxel'})
    assert len(list(g.graph_property_names())) == 2
    assert 'name' in g.graph_property_names()
    assert 'units' in g.graph_property_names()
    assert g.graph_property('name') == "my_property_graph"
    assert g.graph_property('units')['volume'] == "real"


def test_clear_graph():
    g.clear()
    assert len(g._vertices) == 0
    assert len(g._edges) == 0
    assert len(list(g.vertex_property_names())) == 0
    assert len(list(g.edge_property_names())) == 0
    assert len(list(g.graph_property_names())) == 0



# ##########################################################
#
# Nose setup and teardown methods test functions
#
# ##########################################################
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

    # Add vertex property: the value is the squared value of the ids
    v = dict([(k, k ** 2) for k in vids])
    g.add_vertex_property('squared_id', v)
    # Add edge property: the value is the squared value of the ids
    v = dict([(k, sum(g.edge_vertices(k))) for k in eids])
    g.add_edge_property('sum_vids', v)
    # Add a graph property:
    g.add_graph_property('units', {'squared_id': 'real', 'sum_vids': 'voxel'})
    return g


def setup_func():
    """
    Construct a PropertyGraph, from scratch.
    """
    g = create_graph(n_vids, n_eids)


def teardown_func():
    g.clear()


# ##########################################################
#
# PropertyGraph concept: EMPTY init
#
# ##########################################################
@with_setup(setup_func, teardown_func)
def test_vertex_ppty():
    assert len(list(g.vertex_property_names())) == 1
    assert "squared_id" in list(g.vertex_property_names())
    assert len(g.vertex_property("squared_id")) == n_vids


@with_setup(setup_func, teardown_func)
def test_edge_ppty():
    assert len(list(g.edge_property_names())) == 1
    assert "sum_vids" in list(g.edge_property_names())
    assert len(g.edge_property("sum_vids")) == n_eids


@with_setup(setup_func, teardown_func)
def test_vertex_ppty_values():
    for i in xrange(10):
        assert g.vertex_property('squared_id')[i] == i ** 2
    assert 11 not in g.vertex_property('squared_id')


@with_setup(setup_func, teardown_func)
def test_edge_ppty_values():
    for i in xrange(9):
        s_vid, t_vid = g.edge_vertices(i)
        assert g.edge_property('sum_vids')[i] == s_vid + t_vid
    assert n_eids + 1 not in g.edge_property('sum_vids')


@with_setup(setup_func, teardown_func)
def test_remove_vertex_ppty():
    g.remove_vertex_property("squared_id")
    assert "squared_id" not in list(g.vertex_property_names())


@with_setup(setup_func, teardown_func)
def test_remove_edge_ppty():
    g.remove_edge_property("sum_vids")
    assert "sum_vids" not in list(g.edge_property_names())


# ##########################################################
#
# PropertyGraph concept: graph init
#
# ##########################################################
@with_setup(setup_func, teardown_func)
def test_vertex_ppty_graph_init():
    gb = PropertyGraph(g)
    assert len(list(gb.vertex_property_names())) == 1
    assert "squared_id" in list(gb.vertex_property_names())
    assert len(gb.vertex_property("squared_id")) == n_vids


@with_setup(setup_func, teardown_func)
def test_edge_ppty():
    gb = PropertyGraph(g)
    assert len(list(gb.edge_property_names())) == 1
    assert "sum_vids" in list(gb.edge_property_names())
    assert len(gb.edge_property("sum_vids")) == n_eids


@with_setup(setup_func, teardown_func)
def test_vertex_ppty_values():
    gb = PropertyGraph(g)
    for i in xrange(10):
        assert gb.vertex_property('squared_id')[i] == i ** 2
    assert n_vids + 1 not in gb.vertex_property('squared_id')


@with_setup(setup_func, teardown_func)
def test_edge_ppty_values():
    gb = PropertyGraph(g)
    for i in xrange(9):
        s_vid, t_vid = gb.edge_vertices(i)
        assert gb.edge_property('sum_vids')[i] == s_vid + t_vid
    assert n_eids + 1 not in gb.edge_property('sum_vids')


@with_setup(setup_func, teardown_func)
def test_extend_property_graph():
    gb = PropertyGraph(g)
    gb.extend_property_graph(g)
    # test vertex property relabelling:
    v_ppty = gb.vertex_property('squared_id')
    for i in xrange(n_vids):
        assert v_ppty[i] == v_ppty[i + n_vids]
    # test edge property relabelling:
    e_ppty = gb.edge_property('sum_vids')
    for i in xrange(n_eids):
        assert e_ppty[i] == e_ppty[i + n_eids]



#TODO: test 'graph_property' relabelling, have to re-implement first!
#TODO: test 'domains' creation