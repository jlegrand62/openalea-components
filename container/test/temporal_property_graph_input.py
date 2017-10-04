# -*- python -*-
#
#       OpenAlea.image
#
#       Copyright 2012 INRIA - CIRAD - INRA
#
#       File author(s):  Leo Guignard <leo.guignard@inria.fr>
#                        Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#                        Frederic Boudon <frederic.boudon@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://openalea.gforge.inria.fr
#
################################################################################
"""Helps creating TemporalPropertyGraph"""

from random import sample, randint
from openalea.container import PropertyGraph, TemporalPropertyGraph

def create_random_PG(n_vids=10, n_eids_struct=10):
    g = PropertyGraph()
    # - Add vertices to the two PropertyGraph:
    vids = range(n_vids)
    for v in vids:
        g.add_vertex(v)
    # - Add edges to the two PropertyGraph:
    eids = []
    while len(eids) < n_eids_struct:
        vid_source, vid_target = sample(vids, 2)
        eids.append(g.add_edge(vid_source, vid_target))

    return g, vids, eids


def create_random_TPG(n_vids=10, n_eids_struct=10):
    """
    Create two identical random 'PropertyGraph', link them using random mapping.
    Temporal mapping is done for all vertices

    Parameters
    ----------
    n_vids
    n_eids_struct

    Returns
    -------
    a TemporalPropertyGraph
    """
    g = TemporalPropertyGraph()
    # - Create two 'PropertyGraph':
    p1, p2 = PropertyGraph(), PropertyGraph()
    # - Add vertices to the two PropertyGraph:
    vids = range(n_vids)
    for v in vids:
        p1.add_vertex(v)
        p2.add_vertex(v)
    # - Add STRUCTURAL edges to the two PropertyGraph:
    eids = []
    while len(eids) < n_eids_struct:
        vid_source, vid_target = sample(vids, 2)
        eids.append(p1.add_edge(vid_source, vid_target))
        p2.add_edge(vid_source, vid_target)

    # - Create a random mapping:
    mapping = dict([(vid, [sample(vids, randint(1, 2))]) for vid in vids])
    # - Create the 'TemporalPropertyGraph':
    g.temporal_extension([p1, p2], [mapping])
    return g, mapping


def create_TemporalGraph():
    """create a graph"""
    g = TemporalPropertyGraph()
    g1 = PropertyGraph()
    g2 = PropertyGraph()
    g3 = PropertyGraph()
    # -- Adding vertices to 'PropertyGraph'
    dict_0, dict_1, dict_2 = {}, {}, {}
    for i in range(2):
        g1.add_vertex()
    for i in range(5):
        g2.add_vertex()
    for i in range(13):
        g3.add_vertex()
    # -- Adding STRUCTURAL edges to 'PropertyGraph'
    g1.add_edge(0, 1)
    for i in range(4):
        g2.add_edge(i, i + 1)
    g2.add_edge(1, 3)
    for i in range(12):
        g3.add_edge(i, i + 1)
    g3.add_edge(0, 3)
    g3.add_edge(2, 6)
    g3.add_edge(3, 6)
    g3.add_edge(4, 6)
    g3.add_edge(8, 10)
    g3.add_edge(10, 12)
    # -- Creating lineage: TEMPORAL edges
    L12 = {
        0: [0, 1, 2],
        1: [3, 4]
    }
    L23 = {
        0: [0, 1],
        1: [2],
        2: [3, 4, 5, 6],
        3: [7, 8],
        4: [9, 10, 11, 12]
    }
    # -- Extending TemporalPropertyGraph with structural graph ('PropertyGraph') and linking them with lineage information.
    g.temporal_extension([g1, g2, g3], [L12, L23])
    return g


def create_TemporalGraph_partial_lineage():
    """
    Create a TemporalPropertyGraph with partial lineage.
    Simulate partial lost of lineage du to daughter cells being out of the acquisition frame.

    Here the is no daughters for TPG_vertex_id #6 (cell_id #4 in g2).
    """
    g = TemporalPropertyGraph()
    g1 = PropertyGraph()
    g2 = PropertyGraph()
    g3 = PropertyGraph()
    # -- Adding vertices to 'PropertyGraph'
    dict_0, dict_1, dict_2 = {}, {}, {}
    for i in range(2):
        g1.add_vertex()
    for i in range(5):
        g2.add_vertex()
    for i in range(13):
        g3.add_vertex()
    # -- Adding STRUCTURAL edges to 'PropertyGraph'
    g1.add_edge(0, 1)
    for i in range(4):
        g2.add_edge(i, i + 1)
    g2.add_edge(1, 3)
    for i in range(12):
        g3.add_edge(i, i + 1)
    g3.add_edge(0, 3)
    g3.add_edge(2, 6)
    g3.add_edge(3, 6)
    g3.add_edge(4, 6)
    g3.add_edge(8, 10)
    g3.add_edge(10, 12)
    # -- Creating lineage: TEMPORAL edges
    L12 = {
        0: [0, 1, 2],
        1: [3, 4]
    }
    L23 = {
        0: [0, 1],
        1: [2],
        2: [3, 4, 5, 6],
        3: [7, 8],
    }
    # -- Extending TemporalPropertyGraph with structural graph ('PropertyGraph') and linking them with lineage information.
    g.temporal_extension([g1, g2, g3], [L12, L23])
    return g


def create_TPG_VertexProperty():
    """ Test of laplacian function with TemporalPropertyGraph """
    g = create_TemporalGraph()
    prop = dict([(i, float(i)) for i in g.vertices()])
    g.add_vertex_property('property', prop)
    return g


def create_TemporalGraph_GraphProperty(set_as_ids=True):
    """create a TemporalGraph with extra 'graph_property' to test relabbeling of keys and labels as ids"""
    g = TemporalPropertyGraph()
    g1 = PropertyGraph()
    g2 = PropertyGraph()
    g3 = PropertyGraph()
    # -- Adding vertices to 'PropertyGraph'
    dict_1, dict_2, dict_3 = {}, {}, {}
    for i in range(2):
        g1.add_vertex()
        dict_1[i] = "t1_vid_" + str(i)
    g1.add_graph_property('time-point_vid', dict_1)
    for i in range(5):
        g2.add_vertex()
        dict_2[i] = "t2_vid_" + str(i)
    g2.add_graph_property('time-point_vid', dict_2)
    for i in range(13):
        g3.add_vertex()
        dict_3[i] = "t3_vid_" + str(i)
    g3.add_graph_property('time-point_vid', dict_3)
    # -- Setting keys to vid type: will be relabelled when creating the TPG 
    if set_as_ids:
        g1.set_graph_property_key_to_vid_type('time-point_vid')
        g2.set_graph_property_key_to_vid_type('time-point_vid')
        g3.set_graph_property_key_to_vid_type('time-point_vid')
    # -- Adding STRUCTURAL edges to 'PropertyGraph'
    g1.add_edge(0, 1)
    for i in range(4):
        g2.add_edge(i, i + 1)
    g2.add_edge(1, 3)
    for i in range(12):
        g3.add_edge(i, i + 1)
    g3.add_edge(0, 3)
    g3.add_edge(2, 6)
    g3.add_edge(3, 6)
    g3.add_edge(4, 6)
    g3.add_edge(8, 10)
    g3.add_edge(10, 12)
    # -- Creating lineage: TEMPORAL edges
    L12 = {
        0: [0, 1, 2],
        1: [3, 4]
    }
    L23 = {
        0: [0, 1],
        1: [2],
        2: [3, 4, 5, 6],
        3: [7, 8],
        4: [9, 10, 11, 12]
    }
    # -- Extending TemporalPropertyGraph with structural graph ('PropertyGraph') and linking them with lineage information.
    g.temporal_extension([g1, g2, g3], [L12, L23])
    return g

# ~ def set_graph_property_value_to_vid_type(self, propertyname, property_type = VertexProperty):
# ~ def set_graph_property_value_to_eid_type(self, propertyname,  property_type = VertexProperty):
# ~ def set_graph_property_key_to_vid_type(self, propertyname):
# ~ def set_graph_property_key_to_eid_type(self, propertyname):
