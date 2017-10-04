# -*- coding: utf-8 -*-
#
#       Graph : graph package
#
#       Copyright or Copr. 2006 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       VPlants WebSite : https://gforge.inria.fr/projects/vplants/
#
"""
This module provide a simple pure python implementation
for a graph interface
do not implement copy concept
"""

__license__ = "Cecill-C"
__revision__ = " $Id$ "

from heapq import heappop, heappush

import networkx as nx
import numpy as np
from interface.graph import IEdgeListGraph
from interface.graph import IExtendGraph
from interface.graph import IGraph
from interface.graph import IMutableEdgeGraph
from interface.graph import IMutableVertexGraph
from interface.graph import IVertexListGraph
from interface.graph import InvalidEdge
from interface.graph import InvalidVertex
from utils import IdDict


class Graph(IGraph, IVertexListGraph, IEdgeListGraph, IMutableVertexGraph,
            IMutableEdgeGraph, IExtendGraph):
    """
    Directed graph with multiple links in this implementation:
      * vertices are tuple of edge_in, edge_out
      * edges are tuple of source, target
    """

    def __init__(self, graph=None, idgenerator="set"):
        """
        Graph object constructor.

        If graph is not None, make a copy of the topological structure of graph
        (i.e. don't use the same id)

        Parameters
        ----------
        graph : None|Graph
             if None (default) construct an  empty instance of the object, else
             extend this Graph object with 'graph'.
        """
        self._vertices = IdDict(idgenerator=idgenerator)
        self._edges = IdDict(idgenerator=idgenerator)
        if graph is None:
            print "Constructing EMPTY Graph object"
        else:
            dummy = self.extend(graph)

    def __str__(self):
        """
        Format printed instance type information.
        """
        s = "Object 'Graph' containing:"
        s += "\n  - {} vertices".format(len(self._vertices))
        s += "\n  - {} edges".format(len(self._edges))
        return s

    # ##########################################################
    #
    # DIRECTED Graph concept
    #
    # ##########################################################
    def source(self, eid):
        """
        Get the source vertex id of the directed edge 'eid'.

        Parameters
        ----------
        eid: int
            the id of the directed edge for which to get the source

        Returns
        -------
        vid : int, vid giving the edge source vertex
        """
        try:
            return self._edges[eid][0]
        except KeyError:
            raise InvalidEdge(eid)

    source.__doc__ = IGraph.source.__doc__

    def target(self, eid):
        """
        Get the target vertex id of the directed edge 'eid'.

        Parameters
        ----------
        eid: int
            the id of the directed edge for which to get the target

        Returns
        -------
        vid : int, vid giving the edge target vertex
        """
        try:
            return self._edges[eid][1]
        except KeyError:
            raise InvalidEdge(eid)

    target.__doc__ = IGraph.target.__doc__

    def edge_vertices(self, eid):
        """
        Get the source and target vertex ids of the directed edge 'eid'.

        Parameters
        ----------
        eid : int
            the id of the directed edge for which to get the target

        Returns
        -------
        (vid_source, vid_target) : tuple, vid pair giving the edge source and
        target vertex
        """
        try:
            return self._edges[eid]
        except KeyError:
            raise InvalidEdge(eid)

    edge_vertices.__doc__ = IGraph.edge_vertices.__doc__

    def edge(self, source, target):
        """
        Get the directed edge id using source and target vertex ids.

        Parameters
        ----------
        source : int
            the source vertex id
        target : int
            the target vertex id

        Returns
        -------
        eid : int, edge id corresponding to given source and target vertex ids
        """
        link_in, link_out = self._vertices[source]
        for eid in link_in:
            if self._edges[eid][0] == target:
                return eid
        for eid in link_out:
            if self._edges[eid][1] == target:
                return eid
        return None

    edge.__doc__ = IGraph.edge.__doc__

    def __contains__(self, vid):
        return self.has_vertex(vid)

    __contains__.__doc__ = IGraph.__contains__.__doc__

    def has_vertex(self, vid):
        """
        Test the presence of vertex id 'vid' in the object.

        Parameters
        ----------
        vid : int
            a vertex id to test

        Returns
        -------
        bool
        """
        return vid in self._vertices

    has_vertex.__doc__ = IGraph.has_vertex.__doc__

    def has_edge(self, eid):
        """
        Test the presence of edge id 'eid' in the object.

        Parameters
        ----------
        eid : int
            an edge id to test

        Returns
        -------
        bool
        """
        return eid in self._edges

    has_edge.__doc__ = IGraph.has_edge.__doc__

    def is_valid(self):
        """
        Returns
        -------
        True (always)
        """
        return True

    is_valid.__doc__ = IGraph.is_valid.__doc__

    # ##########################################################
    #
    # Vertex List Graph Concept
    #
    # ##########################################################
    def vertices(self):
        """
        Returns
        -------
        iterator on the list of vertex defined in the object.
        """
        return iter(self._vertices)

    vertices.__doc__ = IVertexListGraph.vertices.__doc__

    def __iter__(self):
        """
        Returns
        -------
        iterator on the list of vertex defined in the object.
        """
        return iter(self._vertices)

    __iter__.__doc__ = IVertexListGraph.__iter__.__doc__

    def nb_vertices(self):
        """
        Compute the number of vertices in the object.

        Returns
        -------
        nb_vertices : int
        """
        return len(self._vertices)

    nb_vertices.__doc__ = IVertexListGraph.nb_vertices.__doc__

    def __len__(self):
        """
        Returns
        -------
        Number of vertices in the object.
        """
        return self.nb_vertices()

    __len__.__doc__ = IVertexListGraph.__len__.__doc__

    # ##########################################################
    #
    # Edge List Graph Concept
    #
    # ##########################################################
    def _iter_edges(self, vid):
        """
        internal function that perform 'edges' with vid not None
        """
        link_in, link_out = self._vertices[vid]
        for eid in link_in:
            yield eid
        for eid in link_out:
            yield eid

    def edges(self, vid=None):
        """
        Get the edges list defined in the object.
        Can be filtered for edges linked to a given vertex 'vid'.

        Parameters
        ----------
        vid : int|None, optional
            if None (default), returns all edge ids, else those related to this
             vertex-id

        Returns
        -------
        iterator on the edge list.
        """
        if vid is None:
            return iter(self._edges)
        elif vid not in self:
            raise InvalidVertex(vid)
        else:
            return self._iter_edges(vid)

    edges.__doc__ = IEdgeListGraph.edges.__doc__

    def nb_edges(self, vid=None):
        """
        Compute the number of edge defined in the object.
        Can be filtered for edges in contact with a given vertex 'vid'.

        Parameters
        ----------
        vid : int|None, optional
            if None (default), compute for all edge ids, else for those related
            to this vertex-id

        Returns
        -------
        int : the number of edges found.
        """
        if vid is None:
            return len(self._edges)
        if vid not in self:
            raise InvalidVertex(vid)
        return len(self._vertices[vid][0]) + len(self._vertices[vid][1])

    nb_edges.__doc__ = IEdgeListGraph.nb_edges.__doc__

    def in_edges(self, vid):
        """
        Get the "in-edges" list defined in the object.
        Can be filtered for "in-edges" in contact with a given vertex 'vid'.

        Parameters
        ----------
        vid : int|None, optional
            if None (default), returns all "in-edges" ids, else those related to
             this vertex-id

        Returns
        -------
        iterator on the "in-edges" list.
        """
        if vid not in self:
            raise InvalidVertex(vid)
        for eid in self._vertices[vid][0]:
            yield eid

    in_edges.__doc__ = IEdgeListGraph.in_edges.__doc__

    def out_edges(self, vid):
        """
        Get the "out-edges" list defined in the object.
        Can be filtered for "out-edges" in contact with a given vertex 'vid'.

        Parameters
        ----------
        vid : int|None, optional
            if None (default), returns all "out-edges" ids, else those related
            to this vertex-id

        Returns
        -------
        iterator on the "out-edges" list.
        """
        if vid not in self:
            raise InvalidVertex(vid)
        for eid in self._vertices[vid][1]:
            yield eid

    out_edges.__doc__ = IEdgeListGraph.out_edges.__doc__

    def nb_in_edges(self, vid):
        """
        Compute the number of "in-edges" defined in the object.
        Can be filtered for "in-edges" in contact with a given vertex 'vid'.

        Parameters
        ----------
        vid : int|None, optional
            if None (default), compute for all "in-edges" ids, else for those
            related to this vertex-id

        Returns
        -------
        int : the number of "in-edges" found.
        """
        if vid not in self:
            raise InvalidVertex(vid)
        return len(self._vertices[vid][0])

    nb_in_edges.__doc__ = IEdgeListGraph.nb_in_edges.__doc__

    def nb_out_edges(self, vid):
        """
        Compute the number of "out-edges" defined in the object.
        Can be filtered for "out-edges" in contact with a given vertex 'vid'.

        Parameters
        ----------
        vid : int|None, optional
            if None (default), compute for all "out-edges" ids, else for those
            related to this vertex-id

        Returns
        -------
        int : the number of "out-edges" found.
        """
        if vid not in self:
            raise InvalidVertex(vid)
        return len(self._vertices[vid][1])

    nb_out_edges.__doc__ = IEdgeListGraph.nb_out_edges.__doc__

    def in_neighbors(self, vid):
        """
        In neighbors of 'vid' are vids with an edge directed toward 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id

        Returns
        -------
        iterator on "in-neighbors" list
        """
        if vid not in self:
            raise InvalidVertex(vid)
        neighbors_list = [self.source(eid) for eid in self._vertices[vid][0]]
        return set(neighbors_list)

    in_neighbors.__doc__ = IVertexListGraph.in_neighbors.__doc__

    def iter_in_neighbors(self, vid):
        """
        Return an iterator on the "in vertices" for vertex 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id

        Returns
        -------
        iterator : an iterator on the set of parent vertices of vertex 'vid'
        """
        return iter(self.in_neighbors(vid))

    def out_neighbors(self, vid):
        """
        Out neighbors of 'vid' are vids with an edge coming from 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id

        Returns
        -------
        iterator on "out-neighbors" list
        """
        if vid not in self:
            raise InvalidVertex(vid)
        neighbors_list = [self.target(eid) for eid in self._vertices[vid][1]]
        return set(neighbors_list)

    out_neighbors.__doc__ = IVertexListGraph.out_neighbors.__doc__

    def iter_out_neighbors(self, vid):
        """
        Return an iterator of the "out vertices" for vertex 'vid'

        Parameters
        ----------
        vid : int
            a vertex id

        Returns
        -------
        iterator : an iterator on the set of child vertices of vertex 'vid'
        """
        return iter(self.out_neighbors(vid))

    def neighbors(self, vid):
        """
        Neighbors of 'vid' are vids sharing a directed edge toward or outward
        'vid'.

        Parameters
        ----------
        vid : int
            a vertex id

        Returns
        -------
        "neighbors" list
        """
        return self.in_neighbors(vid) | self.out_neighbors(vid)

    neighbors.__doc__ = IVertexListGraph.neighbors.__doc__

    def iter_neighbors(self, vid):
        """
        Return the neighbors vertices of vertex 'vid'

        Parameters
        ----------
        vid : int
            a vertex id

        Returns
        -------
        iterartor : iterator on the set of neighbors vertices of vertex 'vid'
        """
        return iter(self.neighbors(vid))

    def nb_in_neighbors(self, vid):
        """
        Compute the number of "in-neighbors" for vertex 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id

        Returns
        -------
        int : the number of "in-neighbors"
        """
        return len(self.in_neighbors(vid))

    nb_in_neighbors.__doc__ = IVertexListGraph.nb_in_neighbors.__doc__

    def nb_out_neighbors(self, vid):
        """
        Compute the number of "out-neighbors" for vertex 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id

        Returns
        -------
        int : the number of "out-neighbors"
        """
        return len(self.out_neighbors(vid))

    nb_out_neighbors.__doc__ = IVertexListGraph.nb_out_neighbors.__doc__

    def nb_neighbors(self, vid):
        """
        Compute the number of "neighbors" for vertex 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id

        Returns
        -------
        int : the number of neighbors
        """
        return len(self.neighbors(vid))

    nb_neighbors.__doc__ = IVertexListGraph.nb_neighbors.__doc__

    def topological_distance(self, vid, edge_dist=lambda x, y: 1,
                             max_depth=float('inf'), full_dict=True,
                             return_inf=True):
        """
        Return vertex distances from 'vid' according to the 'edge_dist' cost
        function, up to a given 'max_depth'.

        Parameters
        ----------
        vid : int
            a vertex id
        edge_dist : function, optional
            the cost function to apply on each edge, by default 1
        max_depth : float, optional
            the maximum depth that we want to reach, by default use 'inf' (ie.
            not limited)
        full_dict : bool, optional
            if True (default), return the entire dictionary (with inf values)
        return_inf : bool, optional
            if True (default) return 'inf' values, else 'nan'.

        Returns
        -------
        dist_dict : dict of the vid_j distances to 'vid' for j in
        "set(vertices) - vid", {vid_j: edge_distance(vid, vid_j)}
        """
        dist = {}
        reduced_dist = {vid: 0}
        q = []

        inf_dist = float('inf')
        for k in self._vertices:
            dist[k] = inf_dist
            heappush(q, (dist[k], k))

        dist[vid] = 0
        heappush(q, (dist[vid], vid))
        treated = set()
        modif = True
        n = self.nb_vertices()
        while len(treated) != n and modif:
            modif = False
            actual_vid = heappop(q)
            while actual_vid[1] in treated and actual_vid[0] == inf_dist:
                actual_vid = heappop(q)

            if actual_vid[0] != inf_dist:
                actual_vid = actual_vid[1]
                treated.add(actual_vid)

                for nei in self.iter_neighbors(actual_vid):
                    d_nei = dist[nei]
                    cost = edge_dist(nei, actual_vid)
                    d_vid = dist[actual_vid]
                    if (((d_nei == inf_dist) or (d_nei > d_vid + cost)) and (
                                    d_vid + cost < max_depth + 1)):
                        d_nei = d_vid + cost
                        dist[nei] = d_nei
                        reduced_dist[nei] = d_nei
                        heappush(q, (d_nei, nei))
                    modif = True
        # ~ return (reduced_dist, dist)[full_dict], q
        if not return_inf:
            dist = dict(
                [(k, (v if v != inf_dist else np.nan)) for k, v in
                 dist.iteritems()])

        return (reduced_dist, dist)[full_dict]

    def adjacency_matrix(self, edge_dist=1, no_edge_val=0, oriented=False,
                         reflexive=False, reflexive_value=0):
        """
        Return the adjacency matrix of the graph, ie. the pairwise matrix
        presence or absence of adjacency (neighbors) between two vertex.
        Then :
         - adj_matrix[i,j] = 'edge_dist' account for an existing adjacency
           between vertex i and j;
         - adj_matrix[i,j] = 'no_edge_val' account for the absence of adjacency;
           between vertex i and j
         - adj_matrix[i,i] = 'reflexive_value' account for graph reflexivity;

        Parameters
        ----------
        edge_dist : int|float|function, optional
            cost function to apply between two edges, default : 1
        no_edge_val : int|float, optional
            cost to put if there is no edge between two vertices, default : 0
        oriented : bool, optional
            if False (default), the graph is considered non oriented
            (ie. edge[i, j] != edge[j, i]) and only the upper part of the
            symmetric matrix is returned,
            if True, the graph is considered oriented, ie. edge[i, j] != edge[j, i]
        reflexive : bool, optional
            if True (default is False), the graph is considered reflexive and we
            will put the cost or the cost function 'reflexive_value' on the
            diagonal of the adjacency matrix
        reflexive_value : int|float, optional
            value used by the cost function on the diagonal of the adjacency
            matrix, by default 0.

        Returns
        -------
        numpy.array : a NxN matrix
        """
        if isinstance(edge_dist, int) or isinstance(edge_dist, float):
            edge_dist_func = lambda g, x, y: edge_dist
        elif isinstance(edge_dist, type(lambda m: 1)):
            edge_dist_func = edge_dist
        else:
            err = "Got wrong type for 'edge_dist': {}".format(type(edge_dist))
            raise TypeError(err)

        n = self.nb_vertices()
        # adj_matrix = np.array(n * [n * [no_edge_val]])
        if no_edge_val == 0:
            adj_matrix = np.zeros([n, n])
        else:
            adj_matrix = np.ones([n, n]) * no_edge_val
        for edge in self.edges():
            v1, v2 = self.edge_vertices(edge)
            adj_matrix[v1, v2] = edge_dist_func(self, v1, v2)
            if oriented:
                adj_matrix[v2, v1] = edge_dist_func(self, v2, v1)

        if reflexive:
            if isinstance(reflexive_value, int) or isinstance(reflexive_value,
                                                              float):
                reflexive_func = lambda g, x, y: reflexive_value
            elif isinstance(reflexive_value, type(lambda m: 1)):
                reflexive_func = reflexive_value
            else:
                err = "Got wrong type for 'reflexive_value': {}".format(
                    type(reflexive_value))
                raise TypeError(err)
            for i in range(n):
                adj_matrix[i, i] = reflexive_func(self, i, i)

        return adj_matrix

    def floyd_warshall(self, edge_dist=1, oriented=False, reflexive=False,
                       reflexive_value=0):
        """
        Returns the pairwise topological distance matrix using the shortest path
        between vertices of the object.

        Parameters
        ----------
        edge_dist : int
            cost function to apply between two edges, default : 1
        oriented : bool
            if False (default), the graph is not oriented (edge_ij == edge_ji)
            and only the upper part of the matrix is returned, else the graph
            is oriented and every pairwise distance is computed.
        reflexive : bool, optional
            if True (default is False), the graph is considered reflexive and we
            will put the cost or the cost function 'reflexive_value' on the
            diagonal of the adjacency matrix
        reflexive_value : int|float, optional
            value used by the cost function on the diagonal of the adjacency
            matrix, by default 0.

        Returns
        -------
        NxN matrix, where N it the number of vids
        """
        adj_matrix = self.adjacency_matrix(edge_dist, float('inf'), oriented,
                                           reflexive, reflexive_value)
        n = self.nb_vertices()
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    adj_matrix[i, j] = min(adj_matrix[i, j],
                                           adj_matrix[i, k] + adj_matrix[k, j])
                    if not reflexive and (i == j):
                        adj_matrix[i, j] = float('inf')
        return adj_matrix

    def neighborhood(self, vid, depth=1):
        """
        Return the neighborhood of the vertex 'vid' at distance 'depth'.
        Beware, only the vids at distance 'depth' of 'vid' will be returned, ie. not
        those in between.

        Parameters
        ----------
        vid : int
            a vertex id
        depth : int, optional
            the distance at which we detect neighbors

        Returns
        -------
        neighbors_set : the set of the vertices at distance 'depth' of 'vid' (including vid)
        """
        # TODO: Check if this returns ONLY vids at distance 'depth'
        dist = self.topological_distance(vid, max_depth=depth, full_dict=False)
        return set(dist.keys())

    def iter_neighborhood(self, vid, depth):
        """
        Return an iterator on the neighborhood of the vertex 'vid' at distance 'depth'.
        Beware, only the vids at distance '' of 'vid' will be returned, ie. not
        those in between.

        Parameters
        ----------
        vid : int
            a vertex id
        depth : int
            the distance at which we detect neighbors

        Returns
        -------
        iterator : an iterator on the set of the vertices at distance n of the vertex 'vid'
        """
        return iter(self.neighborhood(vid, depth))

    # ##########################################################
    #
    # Mutable Vertex Graph concept
    #
    # ##########################################################
    def add_vertex(self, vid=None):
        """
        Add a vertex to the object, the id may be specified.

        Parameters
        ----------
        vid : int|None, optional
            if None (default) use the 'idgenerator' to compute next id, else use
            provided value

        Returns
        -------
        Nothing, add to the object.
        """
        return self._vertices.add((set(), set()), vid)

    add_vertex.__doc__ = IMutableVertexGraph.add_vertex.__doc__

    def remove_vertex(self, vid):
        """
        Remove a vertex from the object using its id, as well as the in/out
        edges linked to this vertex.

        Parameters
        ----------
        vid : int
            id of the vertex to remove

        Returns
        -------
        Nothing, remove from the object.
        """
        if vid not in self:
            raise InvalidVertex(vid)
        link_in, link_out = self._vertices[vid]
        for edge in list(link_in):
            self.remove_edge(edge)
        for edge in list(link_out):
            self.remove_edge(edge)
        del self._vertices[vid]

    remove_vertex.__doc__ = IMutableVertexGraph.remove_vertex.__doc__

    def clear(self):
        """
        Clear the object from its edges and vertices.

        Returns
        -------
        Nothing, empty the object.
        """
        self._edges.clear()
        self._vertices.clear()

    clear.__doc__ = IMutableVertexGraph.clear.__doc__

    # ##########################################################
    #
    # Mutable Edge Graph concept
    #
    # ##########################################################
    def add_edge(self, sid, tid, eid=None):
        # def add_edge(self, vid_pair, eid=None):
        """
        Add a directed edge between two vertex (source and target), the edge id 
        may be specified.

        Parameters
        ----------
        sid : int
            source vertex id
        tid : int
            target vertex id
        eid : int|None, optional
            id to assign to the edge, if None (default) use 'idgenerator' else
            use the given id

        Returns
        -------
        the id of the newly created edge
        """
        # sid, tid = vid_pair
        if sid not in self:
            raise InvalidVertex(sid)
        if tid not in self:
            raise InvalidVertex(tid)
        eid = self._edges.add((sid, tid), eid)
        self._vertices[sid][1].add(eid)
        self._vertices[tid][0].add(eid)
        return eid

    add_edge.__doc__ = IMutableEdgeGraph.add_edge.__doc__

    def remove_edge(self, eid):
        """
        Remove the edge from the object using its id 'eid'.

        Parameters
        ----------
        eid : int
            the edge id to remove

        Returns
        -------
        Nothing, remove edge from the object.
        """
        if not self.has_edge(eid):
            raise InvalidEdge(eid)
        sid, tid = self._edges[eid]
        self._vertices[sid][1].remove(eid)
        self._vertices[tid][0].remove(eid)
        del self._edges[eid]

    remove_edge.__doc__ = IMutableEdgeGraph.remove_edge.__doc__

    def clear_edges(self):
        """
        Removes all edges from the object.

        Returns
        -------
        Nothing, remove all edges from the object.
        """
        self._edges.clear()
        for vid, (in_edges, out_edges) in self._vertices.iteritems():
            in_edges.clear()
            out_edges.clear()

    clear_edges.__doc__ = IMutableEdgeGraph.clear_edges.__doc__

    # ##########################################################
    #
    # Extend Graph concept
    #
    # ##########################################################
    def extend(self, graph):
        """
        Extend the current object with a graph object.

        Parameters
        ----------
        graph : Graph|PropertyGraph
            a graph object to use to extend this one

        Returns
        -------
        trans_vid : dict
            translate 'graph' to object vertex ids
        trans_eid : dict
            translate 'graph' to object edge ids
        """
        # TODO: Shouldn't we check for existence of 'vertices', 'edges', 'source' & 'target' attribute in 'graph'?
        # TODO: or use 'isinstance(graph, Graph) or 'isinstance(graph, PropertyGraph)' ?
        # vertex adding
        trans_vid = {}
        for vid in list(graph.vertices()):
            trans_vid[vid] = self.add_vertex()

        # edge adding
        trans_eid = {}
        for eid in list(graph.edges()):
            sid = trans_vid[graph.source(eid)]
            tid = trans_vid[graph.target(eid)]
            trans_eid[eid] = self.add_edge(sid, tid)

        return trans_vid, trans_eid

    extend.__doc__ = IExtendGraph.extend.__doc__

    def sub_graph(self, vids):
        """
        Create a sub-graph using a list of vertex ids 'vids'. All edges found
        between given 'vids' will also be returned.

        Parameters
        ----------
        vids : list
            a graph object to use to extend this one

        Returns
        -------
        a sub-graph of the same type than the object.
        """
        # Create a copy of the object:
        from copy import deepcopy
        vids = set(vids)
        # Clear the copy:
        result = deepcopy(self)
        result._vertices.clear()
        result._edges.clear()
        # Re-build the copy, with the necessary vertex and edge ids:
        s = self.source
        t = self.target
        for key, edges in self._vertices.items():
            if key in vids:
                inedges, outedges = edges
                sortedinedges = set([e for e in inedges if s(e) in vids])
                sortedoutedges = set([e for e in outedges if t(e) in vids])
                result._vertices.add((sortedinedges, sortedoutedges), key)
                for eid in sortedoutedges:
                    result._edges.add(self._edges[eid], eid)

        return result

    def to_networkx(self):
        """
        Return a NetworkX Graph object from current object.

        Returns
        -------
        a NetworkX graph object (nx.Graph).
        """
        nx_graph = nx.Graph()
        nx_graph.add_nodes_from(self.vertices())
        nx_graph.add_edges_from(
            ((self.source(eid), self.target(eid)) for eid in self.edges()))

        return nx_graph

    def from_networkx(self, nx_graph):
        """
        Return a Graph from a NetworkX Directed graph 'nx_graph'.

        Parameters
        ----------
        graph : networkx.graph
            a networkx directed graph
        """
        # Clear the object:
        self.clear()
        # Transform NetworkX graph to "directed graph" if not:
        if not nx_graph.is_directed():
            nx_graph = nx_graph.to_directed()
        # Add vertices:
        for vid in nx_graph.nodes_iter():
            self.add_vertex(vid)
        # Add edges:
        for source, target in nx_graph.edges_iter():
            d = nx_graph[source][target]
            eid = d.get('eid')
            self.add_edge(source, target, eid)

        return self.__str__()
