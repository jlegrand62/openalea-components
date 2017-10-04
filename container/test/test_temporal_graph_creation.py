# -*- python -*-
#
#       OpenAlea.image
#
#       Copyright 2012 INRIA - CIRAD - INRA
#
#       File author(s):  Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#                        Frederic Boudon <frederic.boudon@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://openalea.gforge.inria.fr
#
################################################################################
"""Test creation of TemporalPropertyGraph"""
from random import sample, randint

from nose import with_setup
from openalea.container import TemporalPropertyGraph
from openalea.container.temporal_property_graph import flatten

from temporal_property_graph_input import create_random_PG, create_random_TPG


def test_init_from_PGs():
    p1, vids1, eids1 = create_random_PG(10, 12)
    p2, vids2, eids2 = create_random_PG(20, 24)
    n_vids1, n_vids2 = len(vids1), len(vids2)
    n_eids1, n_eids2 = len(eids1), len(eids2)
    mapping = dict([(vid, [sample(vids2, randint(1, 2))]) for vid in vids1])
    n_temp_eids = len(list(flatten(mapping.values())))
    g = TemporalPropertyGraph(graph=[p1, p2], mappings=[mapping],
                              time_steps=[0, 1])
    assert 'edge_type' in g.edge_property_names()
    assert 'index' in g.vertex_property_names()
    assert hasattr(g, 'nb_time_points')
    assert hasattr(g, '_old_to_new_vids')
    assert len(g._old_to_new_vids) == 2
    assert len(list(flatten(g._old_to_new_vids))) == n_vids1 + n_vids2
    assert hasattr(g, '_old_to_new_eids')
    assert len(g._old_to_new_eids) == 2
    assert len(list(flatten(g._old_to_new_eids))) == n_eids1 + n_eids2
    assert len(list(g.vertices())) == n_vids1 + n_vids2
    assert len(g.vertex_property('index')) == n_vids1 + n_vids2
    assert len(list(g.edges())) == n_eids1 + n_eids2 + n_temp_eids
    assert len(g.edge_property('edge_type')) == n_eids1 + n_eids2 + n_temp_eids


n_vids = 20
n_eids = 30


def test_init():
    """ Test TPG creation from random 'PropertyGraph' and mapping."""
    g, mapping = create_random_TPG(n_vids, n_eids)
    n_temp_eids = len(list(flatten(mapping.values())))
    assert 'edge_type' in g.edge_property_names()
    assert 'index' in g.vertex_property_names()
    assert hasattr(g, 'nb_time_points')
    assert hasattr(g, '_old_to_new_vids')
    assert len(g._old_to_new_vids) == 2
    assert len(list(flatten(g._old_to_new_vids))) == 2 * n_vids
    assert hasattr(g, '_old_to_new_eids')
    assert len(g._old_to_new_eids) == 2
    assert len(list(flatten(g._old_to_new_eids))) == 2 * n_eids
    assert len(list(g.vertices())) == 2 * n_vids
    assert len(g.vertex_property('index')) == 2 * n_vids
    assert len(list(g.edges())) == 2 * n_eids + n_temp_eids
    assert len(g.edge_property('edge_type')) == 2 * n_eids + n_temp_eids


def test_edge_type():
    """ Test use of 'edge_type' parameter to retrieve edge ids."""
    g, mapping = create_random_TPG(n_vids, n_eids)
    n_temp_eids = len(list(flatten(mapping.values())))
    assert len(list(g.edges(edge_type='s'))) == n_eids * 2
    assert len(list(g.edges(edge_type='t'))) == n_temp_eids


def test_vertex_at_time():
    g, mapping = create_random_TPG(n_vids, n_eids)
    n_temp_eids = len(list(flatten(mapping.values())))
    assert len(g.vertex_at_time(0)) == n_vids
    assert len(g.vertex_at_time(1)) == n_vids


def test_export_TPG_2_networkx():
    """ Create a random graph and export it to networkx """
    g, mapping = create_random_TPG()
    nxg = g.to_networkx()


def test_import_TPG_from_networkx():
    """ Create a random graph and import it to networkx """
    g, mapping = create_random_TPG()
    nxg = g.to_networkx()
    gg = TemporalPropertyGraph().from_networkx(nxg)
    # ~ assert g==gg

# def test_TGP_display_by_networkx(display = False):
#     """ Test the display of a TPG  by networkx """
#     g = create_random_TPG()
#     nxg = g.to_networkx()
#     import matplotlib.pyplot as plt
#     nx.draw(nxg)
#     if display:
#         plt.show()
