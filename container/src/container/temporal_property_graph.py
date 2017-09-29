# -*- python -*-
#
#       OpenAlea.Container
#
#       Copyright 2011 - 2013 INRIA - CIRAD - INRA
#
#       File author(s): Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#                       Christophe Pradal
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://openalea.gforge.inria.fr
#
################################################################################
"""
This module provide a class that extends the PropertyGraph with new types of
edges: temporal edges.
"""

__license__ = "Cecill-C"
__revision__ = " $Id$ "

import collections

from property_graph import *

try:
    from networkx import greedy_color
except ImportError:
    raise ImportError(
        "NetworkX library cannot be found! Please install package 'python-networkx'.")


def iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False


# TODO: check the use of edge_property 'old_label', not sure it is necessary
class TemporalPropertyGraph(PropertyGraph):
    """
    Extension of PropertyGraph to temporal case, by defining two types of edges,
    spatial ('s') and temporal ('t').

    Uses four dictionaries to maintain these properties:
     * edge_property 'edge_type': define type of edge, ie. spatial or temporal;
     * edge_property 'old_label': define 'old_label' edge ppty: {eid: (label_i, label_j)}
     * vertex_property 'old_label' : define the translation dict between the id of the ;
     * vertex_property 'index' : define the temporal index of each vid;
    """
    STRUCTURAL = 's'
    TEMPORAL = 't'

    def __init__(self, graph=None, **kwds):
        """
        TemporalPropertyGraph constructor.
        """

        PropertyGraph.__init__(self, graph, idgenerator='max', **kwds)
        if graph is None:
            # Introduce 'edge_type' ppty: {eid: STRUCTURAL|TEMPORAL}
            self.add_edge_property('edge_type')
            # Introduce 'old_label' edge ppty: {eid: (label_i, label_j)}
            self.add_edge_property('old_label')
            # Introduce 'old_label' vertex ppty: {vid: label_i}
            self.add_vertex_property('old_label')
            # Introduce temporal 'index' vertex ppty: {vid: time_point}
            self.add_vertex_property('index')
            # Save number of time point added to the object
            self.nb_time_points = 0
            # list of (N_tp) tuples of two dict vertex & edges:
            # (self._old_to_new_ids[tp][0] = old_to_new_vids; self._old_to_new_ids[tp][1] = old_to_new_eids).
            self._old_to_new_ids = []
        else:
            try:
                self._edge_property['edge_type'] = graph._edge_property[
                    'edge_type']
            except KeyError:
                # Introduce 'edge_type' ppty: {eid: STRUCTURAL|TEMPORAL}
                self.add_edge_property('edge_type')
            try:
                self._edge_property['old_label'] = self._edge_property[
                    'old_label']
            except KeyError:
                # Introduce 'old_label' edge ppty: {eid: (label_i, lable_j)}
                self.add_edge_property('old_label')
            try:
                self._vertex_property['old_label'] = self._vertex_property[
                    'old_label']
            except KeyError:
                # Introduce 'old_label' vertex ppty: {vid: (label_i, lable_j)}
                self.add_vertex_property('old_label')
            try:
                self._vertex_property['index'] = self._vertex_property['index']
            except KeyError:
                # Introduce temporal 'index' vertex ppty: {vid: time_point}
                self.add_vertex_property('index')
            try:
                self.nb_time_points = getattr(graph, 'nb_time_points')
            except AttributeError:
                # Save number of time point added to the object
                self.nb_time_points = 0
            try:
                self._old_to_new_ids = getattr(graph, '_old_to_new_ids')
            except AttributeError:
                # list of (N_tp) tuples of two dict vertex & edges:
                # (self._old_to_new_ids[tp][0] = old_to_new_vids; self._old_to_new_ids[tp][1] = old_to_new_eids).
                self._old_to_new_ids = []

    def __str__(self):
        """
        Format returned instance type informations.
        """
        s = "Object 'TemporalPropertyGraph' containing:"
        s += "\n  - {} time-points".format(len(self.nb_time_points))
        s += "\n  - {} vertices".format(len(self._vertices))
        s += "\n  - {} edges".format(len(self._edges))
        s += "\n  - {} vertex properties".format(len(self._vertex_property))
        s += "\n  - {} edge properties".format(len(self._edge_property))
        s += "\n  - {} edge properties".format(len(self._graph_property))
        return s

    def temporal_extension(self, graphs, mappings, time_steps=None,
                           disable_lineage_checking=True):
        # def extend(self, graphs, mappings, time_steps = None, disable_lineage_checking=True):
        """
        Extend the TemporalPropertyGraph with PropertyGraph `graphs` and `mappings`.
        Each graph in the `graphs` list contains :structural edges:.
        Each lineage in the `mappings` list define :temporal edges: between two `graphs`.

        Parameters
        ----------
        graphs : list
            a list of PropertyGraph
        mappings : list
            a list of dictionaries linking the vids of consecutive graphs
        time_steps : list, optional
            absolute acquisition time starting at t=0
        disable_lineage_checking : bool, optional
            if True (default), the lineage is NOT checked before creation, see notes

        Notes
        -----
        * `mapping[0]` relate `graph[0]` vertex-ids to `graph[1]` vertex-ids
        * hence: `len(graphs) == len(mappings)-1` should be true
        * `time_steps` will be used to compute time-derivative scalars
        * 'lineage checking' make sure a vid can not have two mothers
        """
        # - Usual paranoia (avoid useless computation):
        assert len(graphs) == len(mappings) + 1
        if time_steps is not None:
            assert len(graphs) == len(time_steps)
        # - First append a spatial graph (PropertyGraph):
        self.append(graphs[0])
        self.nb_time_points += 1
        # - Now loop over PropertyGraph & mapping (temporal relation between graphs):
        for g, m in zip(graphs[1:], mappings):
            self.append(g, m, disable_lineage_checking)
            self.nb_time_points += 1
        # - Finally save the first property: 'time_steps'
        if time_steps is not None:
            self.add_graph_property('time_steps', time_steps)

        return self._old_to_new_ids

    # TODO: rename to 'append_time_point' ?
    def append(self, graph, mapping=None, disable_lineage_checking=False):
        """
        Append a PropertyGraph `graph` to the TemporalPropertyGraph at next time-point
        using temporal relations within `mapping`.
        This require to first add a PropertyGraph without `mapping`.

        Parameters
        ----------
        graph : PropertyGraph
            a PropertyGraph to append to the previous one with `mapping`
        mapping : list, optional
            if None (default) add the `graph` as first PropertyGraph (the TPG should be
            empty), else a dictionaries linking the vids of consecutive graphs
        disable_lineage_checking : bool, optional
            if True (default), the lineage is checked before creation, see notes

        Notes
        -----
        * 'graph[t]' the last PropertyGraph (at time-point `t`) of the TPG
        * 'graph[t+1]' the new PropertyGraph to add (`graph`)
        * then `mapping` relate 'graph[t]' vertex-ids to 'graph[t+1]' vertex-ids
        * 'lineage checking' make sure a vid can not have two mothers
        """
        # Usual paranoia: check the existence of a first graph when trying to link it to a second one!
        if mapping:
            assert len(
                self._old_to_new_ids) >= 1, 'To create temporal edges between two graphs, first add a graph without mapping.'
        # Get our current index (position) in time using the length of the 'old_to_new_ids' dictionary:
        current_index = len(self._old_to_new_ids)
        # Get useful values in shorter variable names:
        edge_types = self.edge_property('edge_type')
        old_edge_labels = self.edge_property('old_label')
        old_vertex_labels = self.vertex_property('old_label')
        indices = self.vertex_property('index')

        # Add and translate the vertex and edge ids of the second graph:
        relabel_ids = Graph.extend(self, graph)
        old_to_new_vids, old_to_new_eids = relabel_ids
        self._old_to_new_ids.append(
            relabel_ids)  # was put later, here should be good too!
        # Relabel the edge and vertex property:
        self._relabel_and_add_vertex_edge_properties(graph, old_to_new_vids,
                                                     old_to_new_eids)

        # Update properties on graph
        temporalgproperties = self.graph_properties()

        # while on a property graph, graph_property are just dict of dict,
        # on a temporal property graph, graph_property are dict of list of dict
        # to keep the different values for each time point.
        for gname in graph.graph_property_names():
            if gname in [self.metavidtypepropertyname,
                         self.metavidtypepropertyname]:
                temporalgproperties[gname] = graph.graph_property(gname)
            else:
                newgproperty = graph._relabel_and_add_graph_property(gname,
                                                                     old_to_new_vids,
                                                                     old_to_new_eids)
                temporalgproperties[gname] = temporalgproperties.get(gname,
                                                                     []) + [
                                                 newgproperty]

        # ~ self._old_to_new_ids.append(relabel_ids)

        # set edge_type property for structural edges
        for old_eid, eid in old_to_new_eids.iteritems():
            edge_types[eid] = self.STRUCTURAL
            old_edge_labels[eid] = old_eid

        for old_vid, vid in old_to_new_vids.iteritems():
            old_vertex_labels[vid] = old_vid
            indices[vid] = current_index

        def use_sub_lineage(mother, daughters, on_ids_source, on_ids_target):
            found_sub_lineage = False
            tmp_daughters = []
            for d in daughters:
                if iterable(d):
                    found_sub_lineage = True
                    tmp_d = []
                    if "sub_lineage" not in self.graph_properties():
                        self.add_graph_property("sub_lineage")
                    for sub_d in d:
                        if iterable(sub_d):
                            use_sub_lineage(mother, sub_d, on_ids_source,
                                            on_ids_target)
                        else:
                            tmp_d.append(on_ids_target[0][sub_d])
                    tmp_daughters.append(tmp_d)
                else:
                    tmp_daughters.append(on_ids_target[0][d])
            if found_sub_lineage:
                self.graph_property("sub_lineage").update(
                    {on_ids_source[0][mother]: tmp_daughters})

        def flatten(l):
            """
            Flatten everything that's a list!

            Examples
            --------
            >>> list(flatten([1,2,3,4]))
            >>> [1, 2, 3, 4]

            >>> list(flatten([[1,2],[3,4]]))
            >>> [1, 2, 3, 4]

            >>> list(flatten([[1,[2,3]],4]))
            >>> [1, 2, 3, 4]
            """
            for el in l:
                if isinstance(el, collections.Iterable) and not isinstance(el,
                                                                           basestring):
                    for sub in flatten(el):
                        yield sub
                else:
                    yield el

        if mapping:
            if disable_lineage_checking:
                print "Lineage checking is disabled!"
            unused_lineage, no_ancestor, no_all_desc = {}, {}, {}
            on_ids_source, on_ids_target = self._old_to_new_ids[
                                           -2:]  # get the last two image2graph dict (translate_ids_Image2Graph @ t_n-1 & t_n)
            for k, l in mapping.iteritems():
                l_flat = list(flatten(
                    l))  # flatten the lineage (i.e. [[1,2],3] -> [1,2,3])
                if disable_lineage_checking:
                    for v in l_flat:
                        try:
                            eid = self.add_edge(on_ids_source[0][k],
                                                on_ids_target[0][v])
                            edge_types[eid] = self.TEMPORAL
                        except KeyError:
                            unused_lineage.update({k: v})
                # Check if the mother cell and ALL daugthers are present in their respective topological graph : WE DON'T WANT TO CREATE A PARTIAL LINEAGE !!!!
                elif k in on_ids_source[0] and (
                    sum([v in on_ids_target[0].keys() for v in l_flat]) == len(
                        l_flat)):
                    use_sub_lineage(k, l, on_ids_source, on_ids_target)
                    for v in l_flat:
                        eid = self.add_edge(on_ids_source[0][k],
                                            on_ids_target[0][v])
                        edge_types[eid] = self.TEMPORAL
                else:
                    unused_lineage.update({k: l})
                if k not in on_ids_source[0]:
                    no_ancestor.update({k: l})
                if not (
                    sum([v in on_ids_target[0].keys() for v in l_flat]) == len(
                        l_flat)):
                    no_all_desc.update({k: l})
            if unused_lineage != {} and not disable_lineage_checking:
                print "Detected partial lineage info from t{} to t{} !!".format(
                    current_index - 1, current_index)
                print "   - {} lineage infos could not be used...".format(
                    len(unused_lineage))
                print "   - this represent {}% ({}/{}) of the initially provided mapping...".format(
                    round(float(len(unused_lineage)) / len(mapping), 3) * 100,
                    len(unused_lineage), len(mapping))
            elif unused_lineage != {}:
                print "Some lineages were not used, most likely because the cells they were related to have been deleted before topological-graph computation."
            if no_ancestor != {}:
                print "   - {} have missing ancestors in their topological graph (at time {}).".format(
                    len(no_ancestor), current_index - 1)
            if no_all_desc != {}:
                print "   - {} have missing descendants in their topological graph (at time {}).".format(
                    len(no_all_desc), current_index)

        return relabel_ids

    def clear(self):
        """
        Clear the object of all vetex, edge and graph properties.

        Returns
        -------
        Nothing, edit object
        """
        PropertyGraph.clear(self)
        self._old_to_new_ids = []

    # TODO: 'vertexpair2edge_map' should exists in Graph (only params are different, 'edge_type' & 'time_point')
    def vertexpair2edge_map(self, vids=None, time_point=None, edge_type=None):
        """
        Compute a dictionary that map pairs of vertex id to edge id.
        It requires the existence of a 'label' property

        Returns
        -------
        dictionary {(vid_i, vid_j): eid_ij}
        """
        # TODO: time-point is not taken into account!!! Do it, code 'edge_at_time'
        return dict([((self.source(eid), self.target(eid)), eid) for eid in
                     self.edges(vids, edge_type)])

    # TODO: 'vertexpair2edge_map' should exists in Graph (only params are different, 'edge_type' & 'time_point')
    def edge2vertexpair_map(self, vids=None, time_point=None, edge_type=None):
        """
        Compute a dictionary that map pairs of vertex id to edge id.
        It requires the existence of a 'label' property

        Parameters
        ----------
        time_point : int, optional
            provide it to restrict returned dictionary to this time-point

        Returns
        -------
        dictionary {eid_ij: (vid_i, vid_j)}
        """
        v2e = self.vertexpair2edge_map(vids, time_point, edge_type)
        return dict([(j, i) for i, j in v2e.iteritems()])

    def in_neighbors(self, vid, edge_type=None):
        """
        Return the list of "in vertices" for vertex 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses all types

        Returns
        -------
        neighbors_set: the set of parent vertices of vertex 'vid'
        """

        if vid not in self.vertices():
            raise InvalidVertex(vid)

        s = self.source
        vids = self._vertices
        if edge_type is None:
            neighbors_set = set([s(eid) for eid in self._vertices[vid][0]])
        else:
            edge_type = self.__to_set(edge_type)
            e_type = self._edge_property['edge_type']
            neighbors_set = set(
                [s(eid) for eid in vids[vid][0] if e_type[eid] in edge_type])
        return neighbors_set

    def iter_in_neighbors(self, vid, edge_type=None):
        """
        Return an iterator on the "in vertices" for vertex 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses all types

        Returns
        -------
        iterator : an iterator on the set of parent vertices of vertex 'vid'
        """
        return iter(self.in_neighbors(vid, edge_type))

    def out_neighbors(self, vid, edge_type=None):
        """
        Return the list of "out vertices" for vertex 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses all types

        Returns
        -------
        neighbors_set : the set of child vertices of vertex 'vid'
        """
        if vid not in self:
            raise InvalidVertex(vid)

        if edge_type is None:
            neighbors_set = set(
                [self.target(eid) for eid in self._vertices[vid][1]])
        else:
            edge_type = self.__to_set(edge_type)
            edge_type_property = self._edge_property['edge_type']
            neighbors_set = set(
                [self.target(eid) for eid in self._vertices[vid][1] if
                 edge_type_property[eid] in edge_type])
        return neighbors_set

    def iter_out_neighbors(self, vid, edge_type=None):
        """
        Return an iterator of the "out vertices" for vertex 'vid'

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses all types

        Returns
        -------
        iterator : an iterator on the set of child vertices of vertex 'vid'
        """
        return iter(self.out_neighbors(vid, edge_type))

    def neighbors(self, vid, edge_type=None):
        """
        Return the neighbors vertices of vertex 'vid'

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses all types

        Returns
        -------
        neighbors_set : the set of neighbors vertices of vertex 'vid'
        """
        return self.in_neighbors(vid, edge_type) | self.out_neighbors(vid,
                                                                      edge_type)

    def iter_neighbors(self, vid, edge_type=None):
        """
        Return the neighbors vertices of vertex 'vid'

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses all types

        Returns
        -------
        iterartor : iterator on the set of neighbors vertices of vertex 'vid'
        """
        return iter(self.neighbors(vid, edge_type))

    def in_edges(self, vid, edge_type=None):
        """
        Return in edges of the vertex 'vid'

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses all types

        Returns
        -------
        edge_list : the set of the in edges of the vertex 'vid'
        """
        if vid not in self:
            raise InvalidVertex(vid)

        if not edge_type:
            edge_list = set([eid for eid in self._vertices[vid][0]])
        else:
            edge_type = self.__to_set(edge_type)
            e_ppty = self._edge_property['edge_type']
            edge_list = set([eid for eid in self._vertices[vid][0] if
                             e_ppty[eid] in edge_type])
        return edge_list

    def iter_in_edges(self, vid, edge_type=None):
        """
        Return an iterator on "in edges" for vertex 'vid'

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses all types

        Returns
        -------
        iterator : an iterator on the set of the in edges of the vertex 'vid'
        """
        return iter(self.in_edges(vid, edge_type))

    def out_edges(self, vid, edge_type=None):
        """
        Return "out edges" ids of vertex 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses all types

        Returns
        -------
        edge_list : the set of the out edges of the vertex 'vid'
        """
        if vid not in self:
            raise InvalidVertex(vid)

        if edge_type is None:
            edge_list = set([eid for eid in self._vertices[vid][1]])
        else:
            edge_type = self.__to_set(edge_type)
            e_ppty = self._edge_property['edge_type']
            edge_list = set([eid for eid in self._vertices[vid][1] if
                             e_ppty[eid] in edge_type])
        return edge_list

    def iter_out_edges(self, vid, edge_type=None):
        """
        Return iterator on "out edge" ids of vertex 'vid'.

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses
            all types

        Returns
        -------
        iterator : an iterator on the set of the in edges of the vertex 'vid'
        """
        return iter(self.out_edges(vid, edge_type))

    def edges(self, vid=None, edge_type=None):
        """
        Return edges of the vertex 'vid', can be filtered by type using 'edge_type'.
        If vid is None (default), returned edges will not be filtered according
        to vertex they link.

        Parameters
        ----------
        vid : int, optional
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses
            all types

        Returns
        -------
        edge_list : the set of the edges of the vertex 'vid'
        """
        if vid is None:
            if edge_type is not None:
                return set([eid for eid in self._edges.keys() if
                            self.edge_property('edge_type')[eid] == edge_type])
            else:
                return set([eid for eid in self._edges.keys()])
        return self.out_edges(vid, edge_type) | self.in_edges(vid, edge_type)

    def iter_edges(self, vid, edge_type=None):
        """
        Return in edges of the vertex 'vid'
        If vid=None, return all edges of the graph

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses
            all types

        Returns
        -------
        iterator : an iterator on the set of the edges of the vertex 'vid'
        """
        return iter(self.edges(vid, edge_type))

    def neighborhood(self, vid, depth=1, edge_type=None):
        """
        Return the neighborhood of the vertex 'vid' at distance 'depth'.
        Beware, only the vids at distance 'depth' of 'vid' will be returned, ie.
        not those in between.

        Parameters
        ----------
        vid : int
            a vertex id
        depth : int, optional
            the distance at which we detect neighbors
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses
            all types

        Returns
        -------
        neighbors_set : the set of the vertices at distance 'depth' of 'vid'
        (including vid)
        """
        # TODO: Check if this returns ONLY vids at distance 'depth'
        dist = self.topological_distance(vid, edge_type=edge_type,
                                         max_depth=depth, full_dict=False)
        return set(dist.keys())

    def iter_neighborhood(self, vid, depth, edge_type=None):
        """
        Return an iterator on the neighborhood of the vertex 'vid' at distance
        'depth'. Beware, only the vids at distance '' of 'vid' will be returned,
         ie. not those in between.

        Parameters
        ----------
        vid : int
            a vertex id
        depth : int
            the distance at which we detect neighbors
        edge_type : str | set | None, optional
            edge type or set of edge types to consider, if None (default) uses
            all types

        Returns
        -------
        iterator : an iterator on the set of the vertices at distance n of the
        vertex 'vid'
        """
        return iter(self.neighborhood(vid, depth, edge_type))

    def topological_distance(self, vid, edge_type=None,
                             edge_dist=lambda x, y: 1, max_depth=float('inf'),
                             full_dict=True, return_inf=True):
        """
        Return each vertices distance from the vertex 'vid' according to
        'edge_dist' a cost function.

        Parameters
        ----------
        vid : int
            a vertex id
        edge_type : str|set|None, optional
            edge type or set of edge types to consider, if None (default) uses
            all types
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

                for nei in self.iter_neighbors(actual_vid, edge_type):
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

    def adjacency_matrix(self, edge_type=None, edge_dist=1, no_edge_val=0,
                         oriented=False, reflexive=False, reflexive_value=0):
        """
        Return the adjacency matrix of the graph.

        Parameters
        ----------
        edge_type : str|set|None, optional
            edge type or set of edge types to consider, if None (default) uses all types
        edge_dist : int, optional
            cost function to apply between two edges, default : 1
        no_edge_val : int, optional
            cost to put if there is no edge between two vertices, default : 0
        oriented : bool, optional
            if True, the graph is considered oriented and we always add an
            edge j -> i if i -> j exists
        reflexive : bool, optional
            if True (default is False), the graph is considered reflexive and we
            will put the cost or the cost function 'reflexive_value' on the
            diagonal of the adjacency matrix
        reflexive_value : int, optional
            value used by the cost function on the diagonal of the adjacency matrix

        Returns
        -------
        numpy.array : a NxN matrix
        """
        if not isinstance(edge_dist, type(lambda m: 1)):
            edge_dist_func = lambda g, x, y: edge_dist
        else:
            edge_dist_func = edge_dist

        n = self.nb_vertices()
        adj_matrix = np.array(n * [n * [no_edge_val]])
        for edge in self.edges(edge_type=edge_type):
            v1, v2 = self.edge_vertices(edge)
            adj_matrix[v1, v2] = edge_dist_func(self, v1, v2)
            if not oriented:
                adj_matrix[v2, v1] = edge_dist_func(self, v2, v1)

        if reflexive:
            if not isinstance(reflexive_value, type(lambda m: 1)):
                reflexive_func = lambda g, x, y: reflexive_value
            else:
                reflexive_func = reflexive_value
            for i in range(n):
                adj_matrix[i, i] = reflexive_func(self, i, i)

        return adj_matrix

    def floyd_warshall(self, edge_type=None, edge_dist=1, oriented=False):
        """
        Returns an adjacency matrix (distance matrix) using the shortest path between
        vertices of the object.
        Used edges can be filtered by 'edge_type'

        Parameters
        ----------
        edge_type : str | set | None
            edge type or set of edge types to consider, if None (default) uses all types
        edge_dist : int
            cost function to apply between two edges, default : 1
        oriented : bool
            define if the graph is oriented, default : False.
            If True, means i->j != j->i
        """
        adj_matrix = self.adjacency_matrix(edge_type, edge_dist,
                                                 float('inf'), oriented)
        n = self.nb_vertices()
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    adj_matrix[i, j] = min(adj_matrix[i, j],
                                                 adj_matrix[i, k] +
                                                 adj_matrix[k, j])
        return adj_matrix

    def vertex_temporal_index(self, vid):
        """
        Return the temporal index of a given 'vid'.
        
        Parameters
        ----------
        vid : int
            id of the vertex
        
        Returns
        -------
        integer `self.vertex_property('index')[vid]`
        """
        if isinstance(vid, list):
            return [self.vertex_temporal_index(v) for v in vid]
        else:
            return self.vertex_property('index')[vid]

    def children(self, vid):
        """
        Return the children list of the vertex 'vid', if any.

        Parameters
        ----------
        vid : int
            id of the vertex

        Returns
        -------
        children_list: the set of the children of the vertex vid
        """
        return self.out_neighbors(vid, 't')

    def iter_children(self, vid):
        """
        Return an iterator on the children list for the vertex 'vid'.

        Parameters
        ----------
        vid : int
            id of the vertex

        Returns
        -------
        children_iterator: an iterator on the set of 'vid' children
        """
        return iter(self.children(vid))

    def has_children(self, vid):
        """
        Return True if the vertex 'vid' has a child or children, else False.

        Parameters
        ----------
        vid : int
            id of the vertex
        """
        if self.children(vid) != set():
            return True
        else:
            return False

    def parent(self, vid):
        """
        Return parents of the vertex vid

        Parameters
        ----------
        vid : int
            id of the vertex

        Returns
        -------
        parents_list: the set of the parents of the vertex vid
        """
        return self.in_neighbors(vid, 't')

    def iter_parent(self, vid):
        """
        Return an iterator on the parents of the vertex 'vid'.

        Parameters
        ----------
        vid : int
            id of the vertex

        Returns
        -------
        parent_iterator:  an iterator on the set of 'vid' parent
        """
        return iter(self.parent(vid))

    def has_parent(self, vid):
        """
        Return True if the vertex 'vid' has a parent, else False.

        Parameters
        ----------
        vid : int
            id of the vertex
        """
        if self.parent(vid) != set():
            return True
        else:
            return False

    def sibling(self, vid):
        """
        Return sibling(s) list of the vertex 'vid'.

        Parameters
        ----------
        vid : int
            id of the vertex

        Returns
        -------
        sibling_list: the set of sibling of the vertex vid
        """
        if self.parent(vid):
            return self.children(self.parent(vid).pop()) - {vid}
        else:
            return None

    def iter_sibling(self, vid):
        """
        Return an iterator on the siblings of the vertex vid

        Parameters
        ----------
        vid : int
            id of the vertex

        Returns
        -------
        sibling_iterator:  an iterator on the set of 'vid' siblings
        """
        return iter(self.sibling(vid))

    # TODO: use 'vid' instead of 'vids'...
    def descendants(self, vids, rank=None):
        """
        Return the 0, 1, ..., nth descendants of the vertex 'vids'.
        This set contains 'vids' its/their descendants at 'rank' and also all other
        vertex ids in between.

        Parameters
        ----------
        vids : int|list|set
            ids of the vertex to use
        rank : int|None, optional
            if None (default) returns all descendant up to the last time-point,
            else is the temporal distance to use when searching descendants

        Returns
        -------
        descendant_set: the set of the 0, 1, ..., nth descendants of 'vids'
        """
        edge_type = 't'
        neighbs = set()
        vids = self.__to_set(vids)
        # If temporal distance is null, returns 'vids'
        if rank == 0:
            return vids
        # If temporal distance is one, use 'children'
        # TODO: use 'self.children()'
        elif rank == 1:
            for vid in vids:
                neighbs |= (self.out_neighbors(vid, edge_type) | {vid})
            return neighbs
        # Else
        else:
            if rank is None:
                rank = self.nb_time_points - 1
            for vid in vids:
                neighbs |= (self.descendants(self.out_neighbors(vid, edge_type),
                                             rank - 1) | {vid})
                if list(neighbs) == self._vertices.keys():
                    return neighbs
        return neighbs

    # TODO: use 'vid' instead of 'vids'...
    def iter_descendants(self, vids, rank=None):
        """
        Return the 0, 1, ..., nth descendants of the vertex 'vids'.
        This set contains 'vids' its/their descendants at 'rank' and also all other
        vertex ids in between.

        Parameters
        ----------
        vids : int|list|set
            ids of the vertex to use
        rank : int|None, optional
            if None (default) returns all descendant up to the last time-point,
            else is temporal distance to use when searching descendants

        Returns
        -------
        iterator: an iterator on the set of the 0, 1, ..., nth descendants of 'vid'
        """
        return iter(self.descendants(vids, rank))

    # TODO: use 'vid' instead of 'vids'...
    def rank_descendants(self, vids, rank=1):
        """
        Return the descendants of the vertex vid only at a given rank.

        Parameters
        ----------
        vid : int|list|set
            id or list of ids of the vertex
        rank : int|None, optional
            if None (default) returns all descendant up to the last time-point,
            else is temporal distance to look at for 'vid' descendants.

        Returns
        -------
        descendant_list: the set of the rank-descendant of vid or a list of set
        """
        if isinstance(vids, list):
            return [self.rank_descendants(v, rank) for v in vids]
        else:
            return self.descendants(vids, rank) - self.descendants(vids,
                                                                   rank - 1)

    def has_descendants(self, vid, rank=1):
        """
        Return True if 'vid' has at least a descendant at `rank`.
        
        Parameters
        ----------
        vid : int
            id or list of ids of the vertex
        rank : int
            temporal distance to look at for 'vid' descendants.
        """
        return self.rank_descendants(vid, rank) != set()

    # TODO: use 'vid' instead of 'vids'...
    def ancestors(self, vids, rank=None):
        """
        Return the 0, 1, ..., nth ancestors of the vertex 'vids'.
        This set contains 'vids' its/their ancestors at 'rank' and also all other
        vertex ids in between.

        Parameters
        ----------
        vids : int|list|set
            ids of the vertex to use
        rank : int|None, optional
            if None (default) returns all ancestors up to the last time-point,
            else is temporal distance to use when searching ancestors

        Returns
        -------
        ancestors_set: the set of the 0, 1, ..., nth ancestors of vertex 'vids'
        """
        edge_type = 't'
        neighbs = set()
        vids = self.__to_set(vids)
        if rank == 0:
            return vids
        elif rank == 1:
            for vid in vids:
                neighbs |= (self.in_neighbors(vid, edge_type) | {vid})
            return neighbs
        else:
            if rank is None:
                rank = self.nb_time_points - 1
            for vid in vids:
                neighbs |= (
                self.ancestors(self.in_neighbors(vid, edge_type), rank - 1) | {
                vid})
                if list(neighbs) == self._vertices.keys():
                    return neighbs
        return neighbs

    # TODO: use 'vid' instead of 'vids'...
    def iter_ancestors(self, vids, rank):
        """
        Return an iterator on the 0, 1, ..., nth ancestors of the vertex 'vid'.

        Parameters
        ----------
        vid : int
            id or list of ids of the vertex
        rank : int
            if None (default) returns all ancestors up to the last time-point,
            else is temporal distance to look at for 'vid' ancestors.

        Returns
        -------
        iterator: an iterator on the set of the 0, 1, ..., nth ancestors of the vertex vid
        """
        return iter(self.ancestors(vids, rank))

    # TODO: use 'vid' instead of 'vids'...
    def rank_ancestors(self, vid, rank=1):
        """
        Return the ancestors of the vertex 'vid' only at given 'rank'.

        Parameters
        ----------
        vid : int
            id or list of ids of the vertex
        rank : int
            temporal distance to look at for 'vid' ancestors.

        Returns
        -------
        descendant_list: the set of the rank-ancestor of the vertex vid or a list of set
        """
        if isinstance(vid, list):
            return [self.rank_ancestors(v, rank) for v in vid]
        else:
            return self.ancestors(vid, rank) - self.ancestors(vid, rank - 1)

    def has_ancestors(self, vid, rank=1):
        """
        Return True if the vertex 'vid' has at least one ancestor at 'rank'.
        
        Parameters
        ----------
        vid : int
            id or list of ids of the vertex
        n : int
            temporal distance to look at for 'vid' ancestors.
        """
        return self.rank_ancestors(vid, rank) != set()

    def exist_relative_at_rank(self, vid, rank, verbose=False):
        """
        Check if there is at least a temporal relative of 'vid' at given `rank`.
        Starts from the vid's temporal index.
        Positive and negative rank, respectively check for descendants or ancestors.

        Parameters
        ----------
        vids : list
            vertex id of origin
        rank : int
            temporal distance used on edges to gather neighbors
        verbose : bool, optional
            if True (default False) enable verbosity

        Note
        ----
        * `rank` >= 1 check for descendants
        * `rank` <= 1 check for ancestors
        """
        ## Obvious FALSE cases:
        # rank > max_self_temporal_depth - vid_index:
        if rank > self.nb_time_points - 1 - self.vertex_property('index')[vid]:
            if verbose: print "WARNING: rank > max_self_temporal_depth - vid_index !"
            return False
        # Null rank:
        if rank == 0:
            if verbose: print "WARNING: Null rank !"
            return False
        ## Non-obvious cases:
        if (rank > 0):
            try:
                return self.descendants(vid, rank) - self.descendants(vid,
                                                                      rank - 1) != set()
            except:
                return False
        if (rank < 0):
            try:
                return self.ancestors(vid, abs(rank)) - self.ancestors(vid, abs(
                    rank) - 1) != set()
            except:
                return False

    def exist_all_relative_at_rank(self, vid, rank, verbose=False):
        """
        Check if ALL temporal relative of 'vid' exists at given `rank`.
        Starts from the vid's temporal index.
        Positive rank check if 'vid' have all its descendants lineaged up to `rank`.
        Negative rank check if 'vid' rank-ancestor (from self.vertex_temporal_index(vid)-abs(rank)) is lineaged up to 'vid' time-point.
        It is the same as checking that all rank-siblings of vid are lineaged.

        Parameters
        ----------
        vids : list
            vertex id of origin
        rank : int
            temporal distance used on edges to gather neighbors
        verbose : bool, optional
            if True (default False) enable verbosity

        Notes
        -----
        * THIS ASSUME THE ABSENCE OF APOPTOSIS IN THE TISSUE
        """
        if rank == 0 or abs(rank) > self.nb_time_points - 1:
            return False

        if rank == 1 or rank == -1:
            return self.exist_relative_at_rank(vid, rank)

        if (rank > 0):
            descendants_at_rank = {1: self.descendants(vid)}
            if descendants_at_rank[1] == set():
                return False
            for r in xrange(2, rank + 1):
                for v in descendants_at_rank[r - 1]:
                    if self.descendants(v) == set():
                        return False
                    if r in descendants_at_rank:
                        descendants_at_rank[r].update(self.descendants(v))
                    else:
                        descendants_at_rank[r] = self.descendants(v)
            return True

        if (rank < 0):
            vid_rank_ancestor = list(self.rank_ancestors(vid, abs(rank)))
            if len(vid_rank_ancestor) == 1:
                vid_rank_ancestor = vid_rank_ancestor[0]
                return self.exist_all_relative_at_rank(vid_rank_ancestor,
                                                       abs(rank))
            elif len(vid_rank_ancestor) > 1:
                print "More than ONE ancestor?!! Sorry but... WEIRD!!!!"
            else:
                return False

    def _lineaged_as_ancestor(self, time_point=None, rank=1):
        """
        Return a list of vertex, lineaged as ancestors at a given 'rank' (ie. has
        descendants up to this rank).
        These vertex can be filtered for a given 'time_point'.
        
        Parameters
        ----------
        time_point : int, optional
            returned vids will belong to this temporal index 
        rank : int, optional
            temporal distance at which the returned vids should have descendants

        Returns
        -------
        a list of vids possessing descendants up to a given rank, and defined at
        'time_point' (if given).
        """
        if time_point is None:
            return [k for k in self.vertices() if self.has_descendants(k, rank)]
        else:
            return [k for k, v in self.vertex_property('index').iteritems() if
                    v == time_point and self.has_descendants(k, rank)]

    def _lineaged_as_descendant(self, time_point=None, rank=1):
        """
        Return a list of vertex, lineaged as descendants at a given 'rank' (ie. has
        ancestors up to this rank).
        These vertex can be filtered for a given 'time_point'.
        
        Parameters
        ----------
        time_point : int, optional
            returned vids will belong to this temporal index 
        rank : int, optional
            temporal distance at which the returned vids should have ancestors

        Returns
        -------
        a list of vids possessing ancestors up to a given rank, and defined at
        'time_point' (if given).
        """
        if time_point is None:
            return [k for k in self.vertices() if self.has_ancestors(k, rank)]
        else:
            return [k for k, v in self.vertex_property('index').iteritems() if
                    v == time_point and self.has_ancestors(k, rank)]

    def _fully_lineaged_vertex(self, time_point=None):
        """
        Return a list of fully lineaged vertex (from a given `time_point` if not None), i.e. lineaged from start to end.

        Parameters
        ----------
        time_point : int, optional
            returned vids will belong to this temporal index 
        """
        rank = self.nb_time_points - 1
        flv = self.descendants([k for k in self.vertex_at_time(0) if
                                self.exist_all_relative_at_rank(k, rank)], rank)
        if time_point is None:
            return flv
        else:
            return [vid for vid in flv if
                    self.vertex_temporal_index(vid) == time_point]

    def lineaged_vertex(self, fully_lineaged=False, as_ancestor=False,
                        as_descendant=False, lineage_rank=1):
        """
        Return ids of lineaged vertices, with differents type of lineage possible:
         - a full lineage, i.e. only vids with a lineage from the first to the last time-point (fully_lineaged=True);
         - a lineage over several ranks, i.e. only vids with a lineage from the vid to the vid+lineage_rank time-point (fully_lineaged=False, lineage_rank=int);
         - an 'ancestors' based lineage (as_ancestor = True), i.e. only vids lineaged as ancestor (over lineage_rank if not None);
         - an 'descendants' based lineage (as_ancestor = True), i.e. only vids lineaged as descendants (over lineage_rank if not None).

        Parameters
        ----------
        fully_lineaged : bool, optional
            if True, returned vids are lineaged from their time-point to the
            last one, else return vids have at least a parent or a child(ren);
        as_ancestor : bool, optional
            if True, return vertices lineaged as ancestor (has a descendant);
        as_descendant : bool, optional
            if True, return vertices lineaged as descendant (has an ancestor);
        lineage_rank : int, optional
            usefull if you want to check the lineage for a different rank than 
            the rank-1 temporal neighborhood.
        """
        if as_ancestor:
            vids_anc = self._lineaged_as_ancestor(time_point=None,
                                                  rank=lineage_rank)
        else:
            vids_anc = self.vertices()
        if as_descendant:
            vids_desc = self._lineaged_as_descendant(time_point=None,
                                                     rank=lineage_rank)
        else:
            vids_desc = self.vertices()

        if fully_lineaged:
            vids = self._fully_lineaged_vertex(time_point=None)
        else:
            all_rel = self.exist_all_relative_at_rank
            vids = [k for k in self.vertices() if
                    (all_rel(k, lineage_rank) or all_rel(k, -lineage_rank))]
        return list(set(vids) & set(vids_anc) & set(vids_desc))

    def _all_vertex_at_time(self, time_point):
        """
        Return a list containing all vertex assigned to a given `time_point`

        Parameters
        ----------
        time_point : int, optional
            returned vids will belong to this temporal index 
        """
        index = self.vertex_property('index')
        return [k for k, v in index.iteritems() if v == time_point]

    def vertex_at_time(self, time_point, lineaged=False, fully_lineaged=False,
                       as_ancestor=False, as_descendant=False, lineage_rank=1):
        """
        Return vertices ids corresponding to a given time point in the TPG.
        Optional parameters can be used to filter the list of vertices ids returned.

        Parameters
        ----------
        time_point : int
            integer giving which time point should be considered;
        lineaged : bool
            if True, return vertices having at least a parent or a child(ren);
        fully_lineaged : bool, optional
            if True, returned vids are lineaged from their time-point to the
            last one, else return vids have at least a parent or a child(ren);
        as_ancestor : bool, optional
            if True, return vertices lineaged as ancestor (has a descendant);
        as_descendant : bool, optional
            if True, return vertices lineaged as descendant (has an ancestor);
        lineage_rank : int, optional
            usefull if you want to check the lineage for a different rank than 
            the rank-1 temporal neighborhood.
        """
        if lineaged or fully_lineaged:
            vids = self.lineaged_vertex(fully_lineaged, as_ancestor,
                                        as_descendant, lineage_rank)
            return list(set(vids) & set(self._all_vertex_at_time(time_point)))
        else:
            return self._all_vertex_at_time(time_point)

    def vertex_property_at_time(self, vertex_property, time_point,
                                lineaged=False, fully_lineaged=False,
                                as_ancestor=False, as_descendant=False,
                                lineage_rank=1):
        """
        Return the `vertex_property``for a given `time_point`.
        May be conditionned by extra temporal property: `lineaged`, `fully_lineaged`, `as_parent`, `as_children`.

        Parameters
        ----------
        vertex_property : str
            the name of an existing 'vertex_property' to extract;
        time_point : int
            integer giving which time point should be considered;
        lineaged : bool
            if True, return vertices having at least a parent or a child(ren);
        fully_lineaged : bool, optional
            if True, returned vids are lineaged from their time-point to the
            last one, else return vids have at least a parent or a child(ren);
        as_ancestor : bool, optional
            if True, return vertices lineaged as ancestor (has a descendant);
        as_descendant : bool, optional
            if True, return vertices lineaged as descendant (has an ancestor);
        lineage_rank : int, optional
            usefull if you want to check the lineage for a different rank than 
            the rank-1 temporal neighborhood.
        """
        vids = self.vertex_at_time(time_point, lineaged, fully_lineaged,
                                   as_ancestor, as_descendant)
        ppty = self.vertex_property(vertex_property)
        return dict([(k, ppty[k]) for k in vids if k in ppty])

    def edge_at_time(self, time_point, edge_type=None, vids=None):
        """
        Return the edges linking vertices defined at given time-point.
        Both source and target vids should belong to the given time-point.

        Parameters
        ----------
        time_point : int
            integer giving which time point should be considered;
        edge_type : str, optional
            filter the type of edge returned
        vids : list, optional
            if not None, filter the list of eids returned;

        Returns
        -------
        list of edge ids defined at given time-point.

        Notes
        -----
        * parameter edge_type='t' will return nothing since none will be defined
        only at the given time-point
        * filtering the returned edges using 'vids' impose edges to be defined
        between these vids, ie. both vid source and target of each edge should
        belong to this list.
        """
        s = self.source
        t = self.target
        vids_tp = self.vertex_at_time(time_point)
        if vids is not None:
            vids_tp = list(set(vids_tp) & set(vids))

        edges = self.edges(vids, edge_type)
        return [eid for eid in edges if s(eid) in vids_tp and t(eid) in vids_tp]

    def edge_property_at_time(self, edge_property, time_point, edge_type=None,
                              vids=None):
        """
        Return the edge property 'edge_property' for a given time-point.

        Parameters
        ----------
        edge_property : str
            a string refering to an existing 'graph.edge_property' to extract;
        time_point : int
            time-point for which to return the `edge_property:;
        edge_type : str, optional
            filter the type of edge returned
        vids : list, optional
            if not None, filter the list of eids returned;

        Notes
        -----
        * parameter edge_type='t' will return nothing since none will be defined
        only at the given time-point
        * filtering the returned edges using 'vids' impose edges to be defined
        between these vids, ie. both vid source and target of each edge should
        belong to this list.
        """
        edges = self.edge_at_time(time_point, edge_type, vids)
        e_ppty = self.edge_property(edge_property)
        return dict([(eid, e_ppty[eid]) for eid in edges if eid in e_ppty])

    def spatial_graph_at_time(self, time_point, vids=None, vtx_ppty=None,
                              edge_ppty=None, graph_ppty=None):
        """
        Returns a spatial graph (PropertyGraph) from a given time-point.

        Parameters
        ----------
        time_point : int
            time-point for which to return the `edge_property:;
        vids : list
            if not None, filter the list of vids returned;
        vtx_ppty : str|list(str)
            a vertex property to return in the PropertyGraph;
        edge_ppty : str|list(str)
            an edge property to return in the PropertyGraph;
        graph_ppty : str|list(str)
            a graph property to return in the PropertyGraph;

        Returns
        -------
        PropertyGraph from given `time_point`. Can contain properties.
        """
        # TODO: Should also return `edge_ppty` and/or `graph_ppty` when required.
        # Check the parameters:
        vids_tp = self.vertex_at_time(time_point)
        if vids is not None:
            vids_tp = list(set(vids_tp) & set(vids))
            if vids_tp == []:
                raise ValueError("No vids left after filtering!")

        result = PropertyGraph()
        # Extract the spatial graph:
        for key, edges in self._vertices.items():
            if key in vids_tp:
                inedges, outedges = edges
                sortedinedges = set([eid for eid in inedges if (
                    (self.source(eid) in vids_tp) and (
                        self.edge_property('edge_type')[eid] == 's'))])
                sortedoutedges = set([eid for eid in outedges if (
                    (self.target(eid) in vids_tp) and (
                        self.edge_property('edge_type')[eid] == 's'))])
                result._vertices.add((sortedinedges, sortedoutedges), key)
                for eid in sortedoutedges:
                    result._edges.add(self._edges[eid], eid)
        # Extract the required vertex_property
        if vtx_ppty is not None:
            if isinstance(vtx_ppty, str):
                vtx_ppty = [vtx_ppty]
            for ppty in vtx_ppty:
                try:
                    ppty_values = [p for k, p in
                                   self.vertex_property(ppty).iteritems() if
                                   k in vids_tp]
                except:
                    print "Could not access to '{}' vertex_property...".format(
                        ppty)
                    pass
                else:
                    result.add_vertex_property(ppty, ppty_values)
        else:
            print 'No vertex property will be returned in the PropertyGraph'
        # Extract the required edge_property
        if edge_ppty is not None:
            if isinstance(edge_ppty, str):
                edge_ppty = [edge_ppty]
            for ppty in edge_ppty:
                try:
                    ppty_values = [p for k, p in
                                   self.edge_property(ppty).iteritems() if
                                   k in vids_tp]
                except:
                    print "Could not access to '{}' edge_property...".format(
                        ppty)
                    pass
                else:
                    result.add_edge_property(ppty, ppty_values)
        else:
            print 'No edge property will be returned in the PropertyGraph'

        return result

    def greedy_spatial_colormap(self, time_point, vids=None):
        """
        Uses `networkx.greedy_color` to extract a neighbor diverging colormap.
        """
        g = self.spatial_graph_at_time(time_point, vids)
        nxg = g.to_networkx()
        tp_color_dict = greedy_color(nxg)
        return tp_color_dict

    # TODO: '[add|extend|remove]_[vertex|edge|graph]_property' is defined in PropertyGraph, add 'time_point' params
    # TODO: create missing func to override function from PropertyGraph

    # TODO: Rename properties to override function from PropertyGraph
    def extend_vertex_ppty_from_vid_dict(self, name, vid_dict):
        """
        Extend a vertex property in the graph (dictionary) using SpatialImage
        labels as keys.

        Parameters
        ----------
        name : str
            the name of the vertex property to add
        vid_dict : dict
            dictionary {label: label_value} to add

        Note
        ----
        if you don't provide the 'mlabel2vertex', you need to give 'time_point'
        """
        # TODO: should not exists... uses 'self.extend_vertex_property' from PropertyGraph
        missing_vertex = list(set(vid_dict.keys()) - set(self.vertices()))
        if missing_vertex != []:
            print "The vid-dictionary '{}' contains {} unknown vertices:\n{}".format(
                name, len(missing_vertex), missing_vertex)
            vid_dict = dict([(k, v) for k, v in vid_dict.iteritems() if
                             k in self.vertices()])

        if name not in self.vertex_properties():
            self.add_vertex_property(name, vid_dict)
        else:
            self.vertex_property(name).update(vid_dict)
        return

    # TODO: Rename properties to override function from PropertyGraph
    def add_edge_ppty_from_eid_dict(self, name, eid_dict,
                                    try_erase_first=False):
        """
        Add an edge property in the graph using an eid dict {eid: value}.

        Parameters
        ----------
        name : str
            the name of the vertex property to add
        eid_dict : dict
            dictionary {eid: ppty_value} to add
        try_erase_first : bool
            try to erase the property before adding it

        Note
        ----
        Use it when transforming 'eid data' during TissueGraph build-up.
        """
        # TODO: should not exists... uses 'self.add_edge_property' from PropertyGraph
        if name in self.edge_properties() and not try_erase_first:
            raise ValueError("Existing edge property '{}'".format(name))
        if try_erase_first:
            print "You asked to try to erase property '{}'!".format(name)
            self.remove_edge_property(name)

        self.add_edge_property(name, eid_dict)
        return

    # TODO: Rename properties to override function from PropertyGraph
    def update_graph_ppty_from_dict(self, name, dictionary):
        """
        Update a graph property in the graph using using a dictionary.

        Parameters
        ----------
        name : str
            the name of the vertex property to update
        dictionary : dict
            dictionary containing values to update the self property
        """
        if name not in self.graph_properties():
            print "Adding self_property '{}'".format(name)
            self.add_graph_property(name, dictionary)
        else:
            self.graph_property(name).update(dictionary)
        return
