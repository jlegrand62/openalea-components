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
"""This module provide a class that extends the PropertyGraph with different types of edges"""

__license__ = "Cecill-C"
__revision__ = " $Id$ "

import warnings
import collections
import numpy as np
from property_graph import *

try:
    from networkx import greedy_color
except ImportError:
    raise ImportError("NetworkX library cannot be found! Please install package 'python-networkx'.")


def iterable(obj):
    try :
        iter(obj)
        return True
    except TypeError,te:
        return False


class TemporalPropertyGraph(PropertyGraph):
    """
    Simple implementation of PropertyGraph using
    dict as properties and two dictionaries to
    maintain these properties
    """
    STRUCTURAL = 's'
    TEMPORAL = 't'

    def __init__(self, graph=None, **kwds):
        """
        TemporalPropertyGraph constructor.
        """
        PropertyGraph.__init__(self, graph, idgenerator='max',**kwds)
        if graph is None:
            # Introduce 'edge_type' ppty: {eid: STRUCTURAL|TEMPORAL}
            self.add_edge_property('edge_type')
            # Introduce 'old_label' edge ppty: {eid: (label_i, lable_j)}
            self.add_edge_property('old_label')
            # Introduce 'old_label' vertex ppty: {vid: (label_i, lable_j)}
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
                self._edge_property['edge_type'] = graph._edge_property['edge_type']
            except:
                # Introduce 'edge_type' ppty: {eid: STRUCTURAL|TEMPORAL}
                self.add_edge_property('edge_type')
            try:
                self._edge_property['old_label'] = self._edge_property['old_label']
            except:
                # Introduce 'old_label' edge ppty: {eid: (label_i, lable_j)}
                self.add_edge_property('old_label')
            try:
                self._vertex_property['old_label'] = self._vertex_property['old_label']
            except:
                # Introduce 'old_label' vertex ppty: {vid: (label_i, lable_j)}
                self.add_vertex_property('old_label')
            try:
                self._vertex_property['index'] = self._vertex_property['index']
            except:
                # Introduce temporal 'index' vertex ppty: {vid: time_point}
                self.add_vertex_property('index')
            try:
                self.nb_time_points = getattr(graph, 'nb_time_points')
            except:
                # Save number of time point added to the object
                self.nb_time_points = 0
            try:
                self._old_to_new_ids = getattr(graph, '_old_to_new_ids')
            except:
                # list of (N_tp) tuples of two dict vertex & edges:
                # (self._old_to_new_ids[tp][0] = old_to_new_vids; self._old_to_new_ids[tp][1] = old_to_new_eids).
                self._old_to_new_ids = []

    def __str__(self):
        """
        Format returned instance type informations.
        """
        print "Object TemporalPropertyGraph:"
        print "  - {} time-points".format(len(self.nb_time_points))
        print "  - {} vertices".format(len(self._vertices))
        print "  - {} edges".format(len(self._edges))
        print "  - {} vertex properties".format(len(self._vertex_property))
        print "  - {} edge properties".format(len(self._edge_property))
        print "  - {} edge properties".format(len(self._graph_property))

    def temporal_extension(self, graphs, mappings, time_steps = None, disable_lineage_checking=True):
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
        lineage_checking : bool, optional
            if True (default), the lineage is checked before creation, see notes

        Notes
        -----
        * `mapping[0]` relate `graph[0]` vertex-ids to `graph[1]` vertex-ids
        * hence: `len(graphs) == len(mappings)-1` should be true
        * `time_steps` will be used to compute time-derivative scalars
        * 'lineage checking' make sure a vid can not have two mothers
        """
        # - Usual paranoia (avoid useless computation):
        assert len(graphs) == len(mappings)+1
        if time_steps is not None:
            assert len(graphs) == len(time_steps)
        # - First append a spatial graph (PropertyGraph):
        self.append(graphs[0])
        self.nb_time_points += 1
        # - Now loop over PropertyGraph & mapping (temporal relation between graphs):
        for g, m in zip(graphs[1:],mappings):
            self.append(g,m, disable_lineage_checking)
            self.nb_time_points += 1
        # - Finally save the first property: 'time_steps'
        if time_steps is not None:
            self.add_graph_property('time_steps',time_steps)

        return self._old_to_new_ids


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
        lineage_checking : bool, optional
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
            assert len(self._old_to_new_ids) >= 1, 'To create temporal edges between two graphs, first add a graph without mapping.'
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
        self._old_to_new_ids.append(relabel_ids) # was put later, here should be good too!
        # Relabel the edge and vertex property:
        self._relabel_and_add_vertex_edge_properties(graph, old_to_new_vids, old_to_new_eids)

        # Update properties on graph
        temporalgproperties = self.graph_properties()

        # while on a property graph, graph_property are just dict of dict,
        # on a temporal property graph, graph_property are dict of list of dict
        # to keep the different values for each time point.
        for gname in graph.graph_property_names():
            if gname in [self.metavidtypepropertyname,self.metavidtypepropertyname]:
                temporalgproperties[gname] = graph.graph_property(gname)
            else:
                newgproperty = graph.translate_graph_property(gname, old_to_new_vids, old_to_new_eids)
                temporalgproperties[gname] = temporalgproperties.get(gname,[])+[newgproperty]

        #~ self._old_to_new_ids.append(relabel_ids)

        # set edge_type property for structural edges
        for old_eid, eid in old_to_new_eids.iteritems():
            edge_types[eid] = self.STRUCTURAL
            old_edge_labels[eid] = old_eid

        for old_vid, vid in old_to_new_vids.iteritems():
            old_vertex_labels[vid] = old_vid
            indices[vid] = current_index

        def use_sub_lineage(mother, daughters, on_ids_source, on_ids_target):
            found_sub_lineage=False; tmp_daughters = []
            for d in daughters:
                if iterable(d):
                    found_sub_lineage=True; tmp_d = []
                    if "sub_lineage" not in self.graph_properties():
                        self.add_graph_property("sub_lineage")
                    for sub_d in d:
                        if iterable(sub_d):
                            use_sub_lineage(mother, sub_d, on_ids_source, on_ids_target)
                        else:
                            tmp_d.append(on_ids_target[0][sub_d])
                    tmp_daughters.append(tmp_d)
                else:
                    tmp_daughters.append(on_ids_target[0][d])
            if found_sub_lineage:
                self.graph_property("sub_lineage").update({on_ids_source[0][mother]:tmp_daughters})

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
                if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
                    for sub in flatten(el):
                        yield sub
                else:
                    yield el

        if mapping:
            if disable_lineage_checking:
                print "Lineage checking is disabled!"
            unused_lineage, no_ancestor, no_all_desc = {}, {}, {}
            on_ids_source, on_ids_target = self._old_to_new_ids[-2:] #get the last two image2graph dict (translate_ids_Image2Graph @ t_n-1 & t_n)
            for k, l in mapping.iteritems():
                l_flat = list(flatten(l)) # flatten the lineage (i.e. [[1,2],3] -> [1,2,3])
                if disable_lineage_checking:
                    for v in l_flat:
                        try:
                            eid = self.add_edge(on_ids_source[0][k], on_ids_target[0][v])
                            edge_types[eid] = self.TEMPORAL
                        except KeyError:
                            unused_lineage.update({k:v})
                # Check if the mother cell and ALL daugthers are present in their respective topological graph : WE DON'T WANT TO CREATE A PARTIAL LINEAGE !!!!
                elif on_ids_source[0].has_key(k) and ( sum([v in on_ids_target[0].keys() for v in l_flat]) == len(l_flat) ):
                    use_sub_lineage(k, l, on_ids_source, on_ids_target)
                    for v in l_flat:
                        eid = self.add_edge(on_ids_source[0][k], on_ids_target[0][v])
                        edge_types[eid] = self.TEMPORAL
                else:
                    unused_lineage.update({k:l})
                if not on_ids_source[0].has_key(k):
                    no_ancestor.update({k:l})
                if not ( sum([v in on_ids_target[0].keys() for v in l_flat]) == len(l_flat) ):
                    no_all_desc.update({k:l})
            if unused_lineage != {} and not disable_lineage_checking:
                print "Detected partial lineage info from t{} to t{} !!".format(current_index-1,current_index)
                print "   - {} lineage infos could not be used...".format(len(unused_lineage))
                print "   - this represent {}% ({}/{}) of the initially provided mapping...".format(round(float(len(unused_lineage))/len(mapping),3)*100,len(unused_lineage),len(mapping))
                #~ print "It is most likely that you are trying to add lineage between non-existant vertex in your spatial graphs!"
                #~ print "Check if you are not using an out-dated graph and erase temporary files (TPG creation...)."
                #~ print "Or maybe the lignage is detected as 'incomplete' because some cells have been removed before topological-graph computation."
            elif unused_lineage != {}:
                print "Some lineages were not used, most likely because the cells they were related to have been deleted before topological-graph computation."
            if no_ancestor != {}:
                print "   - {} have missing ancestors in their topological graph (at time {}).".format(len(no_ancestor), current_index-1)
            if no_all_desc != {}:
                print "   - {} have missing descendants in their topological graph (at time {}).".format(len(no_all_desc), current_index)

        return relabel_ids

    def clear(self):
        PropertyGraph.clear(self)
        self._old_to_new_ids = []

    def __to_set(self, s):
        if not isinstance(s, set):
            if isinstance(s, list):
                s=set(s)
            else:
                s=set([s])
        return s

    def vertex_temporal_index(self, vid):
        """Return the temporal index of a vid `self.vertex_property('index')[vid]`."""
        if isinstance(vid, list):
            return [self.vertex_temporal_index(v) for v in vid]
        else:
            return self.vertex_property('index')[vid]

    def children(self, vid):
        """ Return children of the vertex vid

        Args:
          vid: a vertex id

        Returns:
          children_list: the set of the children of the vertex vid
        """
        return self.out_neighbors(vid, 't')

    def iter_children(self, vid):
        """ Return children of the vertex vid

        Args:
          vid: a vertex id

        Returns:
          iterator: an iterator on the set of the children of the vertex vid
        """
        return iter(self.children(vid))

    def has_children(self, vid):
        """
        Return True if the vid `vid` has a child or children.
        """
        if self.children(vid) != set():
            return True
        else:
            return False

    def parent(self, vid):
        """ Return parents of the vertex vid

        Args:
          vid: a vertex id

        Returns:
          parents_list: the set of the parents of the vertex vid
        """
        return self.in_neighbors(vid, 't')

    def iter_parent(self, vid):
        """ Return parents of the vertex vid

        Args:
          vid: a vertex id

        Returns:
          iterator: the set of the children of the vertex vid
        """
        return iter(self.parent(vid))

    def has_parent(self, vid):
        """
        Return True if the vid `vid` has a parent.
        """
        if self.parent(vid) != set():
            return True
        else:
            return False

    def sibling(self, vid):
        """ Return sibling of the vertex vid

        Args:
          vid: a vertex id

        Returns:
          sibling_list: the set of sibling of the vertex vid
        """
        if self.parent(vid):
            return self.children(self.parent(vid).pop())-set([vid])
        else:
            return None

    def iter_sibling(self, vid):
        """ Return of the vertex vid

        Args:
          vid: a vertex id

        Returns:
          iterator: an iterator on the set of sibling of the vertex vid
        """
        return iter(self.sibling(vid))

    def descendants(self, vids, n = None):
        """ Return the 0, 1, ..., nth descendants of the vertex vid

        Args:
          vids: a set of vertex id

        Returns:
          descendant_list: the set of the 0, 1, ..., nth descendant of the vertex vid
        """
        edge_type='t'
        neighbs=set()
        vids=self.__to_set(vids)
        if n==0 :
            return vids
        elif n==1 :
            for vid in vids:
                neighbs |= (self.out_neighbors(vid, edge_type) | set([vid]))
            return neighbs
        else :
            if n is None :
                n = self.nb_time_points-1
            for vid in vids :
                neighbs |= (self.descendants(self.out_neighbors(vid, edge_type), n-1) | set([vid]))
                if list(neighbs)==self._vertices.keys():
                    return neighbs
        return neighbs

    def iter_descendants(self, vids, n = None):
        """ Return the 0, 1, ..., nth descendants of the vertex vid

        Args:
          vids: a set of vertex id

        Returns:
          iterator: an iterator on the set of the 0, 1, ..., nth descendants of the vertex vid
        """
        return iter(self.descendants(vids, n))

    def rank_descendants(self, vid, rank=1):
        """ Return the descendants of the vertex vid only at a given rank
        Args:
          vid: a vertex id or a set of vertex id
        Returns:
          descendant_list: the set of the rank-descendant of the vertex vid or a list of set
        """
        if isinstance(vid,list):
            return [self.rank_descendants(v,rank) for v in vid]
        else:
            return self.descendants(vid,rank)- self.descendants(vid,rank-1)

    def has_descendants(self, vid, rank=1):
        """
        Return True if the vid `vid` has at least a descendant at `rank`.
        """
        return self.rank_descendants(vid,rank) != set()

    def rank_descendants(self, vid, rank=1):
        """ Return the descendants of the vertex vid only at a given rank
        Args:
          vid: a vertex id or a set of vertex id
        Returns:
          descendant_list: the set of the rank-descendant of the vertex vid or a list of set
        """
        if isinstance(vid,list):
            return [self.rank_descendants(v,rank) for v in vid]
        else:
            return self.descendants(vid,rank)- self.descendants(vid,rank-1)

    def has_descendants(self, vid, rank=1):
        """
        Return True if the vid `vid` has at least a descendant at `rank`.
        """
        return self.rank_descendants(vid,rank) != set()

    def ancestors(self, vids, n = None):
        """Return the 0, 1, ..., nth ancestors of the vertex vid

        Args:
          vids: a set of vertex id

        Returns:
          anestors_list: the set of the 0, 1, ..., nth ancestors of the vertex vid
        """
        edge_type='t'
        neighbs=set()
        vids=self.__to_set(vids)
        if n==0 :
            return vids
        elif n==1 :
            for vid in vids:
                neighbs |= (self.in_neighbors(vid, edge_type) | set([vid]))
            return neighbs
        else :
            if n is None:
                n = self.nb_time_points-1
            for vid in vids :
                neighbs |= (self.ancestors(self.in_neighbors(vid, edge_type), n-1) | set([vid]))
                if list(neighbs)==self._vertices.keys():
                    return neighbs
        return neighbs

    def iter_ancestors(self, vids, n):
        """ Return the 0, 1, ..., nth ancestors of the vertex vid

        Args:
          vids: a set of vertex id

        Returns:
          iterator: an iterator on the set of the 0, 1, ..., nth ancestors of the vertex vid
        """
        return iter(self.ancestors(vids, n))

    def rank_ancestors(self, vid, rank=1):
        """ Return the ancestor of the vertex vid only at a given rank
        Args:
          vid: a vertex id or a set of vertex id
        Returns:
          descendant_list: the set of the rank-ancestor of the vertex vid or a list of set
        """
        if isinstance(vid,list):
            return [self.rank_ancestors(v,rank) for v in vid]
        else:
            return self.ancestors(vid,rank)- self.ancestors(vid,rank-1)

    def has_ancestors(self, vid, rank=1):
        """
        Return True if the vid `vid` has at least an ancestor at `rank`.
        """
        return self.rank_ancestors(vid,rank) != set()


    def _lineaged_as_ancestor(self, time_point=None, rank=1):
        """ Return a list of vertex lineaged as ancestors."""
        if time_point is None:
            return [k for k in self.vertices() if self.has_descendants(k,rank)]
        else:
            return [k for k,v in self.vertex_property('index').iteritems() if v==time_point and self.has_descendants(k,rank)]

    def _lineaged_as_descendant(self, time_point=None, rank=1):
        """ Return a list of vertex lineaged as descendants."""
        if time_point is None:
            return [k for k in self.vertices() if self.has_ancestors(k,rank)]
        else:
            return [k for k,v in self.vertex_property('index').iteritems() if v==time_point and self.has_ancestors(k,rank)]

    def _fully_lineaged_vertex(self, time_point=None):
        """
        Return a list of fully lineaged vertex (from a given `time_point` if not None), i.e. lineaged from start to end.
        """
        rank = self.nb_time_points-1
        flv = self.descendants([k for k in self.vertex_at_time(0) if self.exist_all_relative_at_rank(k, rank)], rank)
        if time_point is None:
            return flv
        else:
            return [vid for vid in flv if self.vertex_temporal_index(vid)==time_point]

    def lineaged_vertex(self, fully_lineaged=False, as_ancestor=False, as_descendant=False, lineage_rank=1):
        """
        Return ids of lineaged vertices, with differents type of lineage possible:
         - a full lineage, i.e. only vids with a lineage from the first to the last time-point (fully_lineaged=True);
         - a lineage over several ranks, i.e. only vids with a lineage from the vid to the vid+lineage_rank time-point (fully_lineaged=False, lineage_rank=int);
         - an 'ancestors' based lineage (as_ancestor = True), i.e. only vids lineaged as ancestor (over lineage_rank if not None);
         - an 'descendants' based lineage (as_ancestor = True), i.e. only vids lineaged as descendants (over lineage_rank if not None).

        :Parameter:
           fully_lineaged: (bool) : if True (and lineage_rank is None), return vertices lineaged from the first to the last time-point (or from vid_time_point to vid_time_point + lineage_rank), otherwise return vertices having at least a parent or a child(ren);
           as_parent: (bool) : if True, return vertices lineaged as parents;
           as_children: (bool) : if True, return vertices lineaged as children;
         - 'lineage_rank' (int): usefull if you want to check the lineage for a different rank than the rank-1 temporal neighborhood.
        """
        if as_ancestor:
            vids_anc = self._lineaged_as_ancestor(time_point=None, rank=lineage_rank)
        else:
            vids_anc = self.vertices()
        if as_descendant:
            vids_desc = self._lineaged_as_descendant(time_point=None, rank=lineage_rank)
        else:
            vids_desc = self.vertices()

        if fully_lineaged:
            vids = self._fully_lineaged_vertex(time_point=None)
        else:
            vids = [k for k in self.vertices() if (self.exist_all_relative_at_rank(k, lineage_rank) or self.exist_all_relative_at_rank(k, -lineage_rank))]
        return list( set(vids) & set(vids_anc) & set(vids_desc) )


    def _all_vertex_at_time(self, time_point):
        """ Return a list containing all vertex assigned to a given `time_point`."""
        return [k for k,v in self.vertex_property('index').iteritems() if v==time_point]

    def vertex_at_time(self, time_point, lineaged=False, fully_lineaged=False, as_ancestor=False, as_descendant=False, lineage_rank=1):
        """
        Return vertices ids corresponding to a given time point in the TPG.
        Optional parameters can be used to filter the list of vertices ids returned.

        Args:
           time_point: (int) : integer giving which time point should be considered;
           lineaged: (bool) : if True, return vertices having at least a parent or a child(ren);
           lineage_rank (int): usefull if you want to check the lineage for a different rank than the rank-1 temporal neighborhood;
           fully_lineaged: (bool) : if True, return vertices linked from the beginning to the end of the graph;
           as_parent: (bool) : if True, return vertices lineaged as parents;
           as_children: (bool) : if True, return vertices lineaged as children.
        """
        if lineaged or fully_lineaged:
            vids = self.lineaged_vertex(fully_lineaged, as_ancestor, as_descendant, lineage_rank)
            return list(set(vids) & set(self._all_vertex_at_time(time_point)))
        else:
            return self._all_vertex_at_time(time_point)

    def vertex_property_at_time(self, vertex_property, time_point, lineaged=False, fully_lineaged=False, as_ancestor=False, as_descendant=False, lineage_rank=1):
        """
        Return the `vertex_property``for a given `time_point`.
        May be conditionned by extra temporal property: `lineaged`, `fully_lineaged`, `as_parent`, `as_children`.

        Args:
           vertex_property: (str): a string refering to an existing 'graph.vertex_property' to extract;
           time_point: (int) : integer giving which time point should be considered;
           lineaged: (bool) : if True, return vertices having at least a parent or a child(ren);
           lineage_rank (int): usefull if you want to check the lineage for a different rank than the rank-1 temporal neighborhood;
           fully_lineaged: (bool) : if True, return vertices linked from the beginning to the end of the graph;
           as_parent: (bool) : if True, return vertices lineaged as parents;
           as_children: (bool) : if True, return vertices lineaged as children.
        """
        vids = self.vertex_at_time(time_point, lineaged, fully_lineaged, as_ancestor, as_descendant)
        return dict([(k,self.vertex_property(vertex_property)[k]) for k in vids if self.vertex_property(vertex_property).has_key(k)])

    def edge_property_at_time(self, edge_property, time_point):
        """
        Return a subpart of graph.edge_property(`edge_property`) at a given time-point.

        Args:
           edge_property: (str): a string refering to an existing 'graph.edge_property' to extract;
           time_point` (int): time-point for which to return the `edge_property:;

        :Examples:
         graph.edge_property_with_image_labelpairs('wall_area', 0)
        """
        eid2vidpair = self.edge2vertexpair_map(time_point)
        tmp_dict = {}
        for eid,(label1,label2) in eid2vidpair.iteritems():
            if self.edge_property(edge_property).has_key(eid):
                tmp_dict[eid] = self.edge_property(edge_property)[eid]
        return tmp_dict

    def spatial_graph_at_time(self, timepoint, vids=None, vtx_ppty=None, egde_ppty=None, graph_ppty=None):
        """
        Returns a spatial graph from a given time-point.

        Args:
          timepoint (int) - time point for which we will return the PropertyGraph;
          vids (list) - list of vids to keep, other vids from this time-point will not be returned;
          vtx_ppty (str|list(str)) - a vertex property to return in the PropertyGraph;
          egde_ppty (str|list(str)) - an edge property to return in the PropertyGraph;
          graph_ppty (str|list(str)) - a graph property to return in the PropertyGraph;
        :TODO:
            Should also return `egde_ppty` and/or `graph_ppty` when required.

        :Returns: a PropertyGraph from given `time_point`. Can contain properties.
        """
        result = PropertyGraph()
        vids_tp = self.vertex_at_time(timepoint)
        if vids is not None:
            vids_tp = list(set(vids_tp) & set(vids))
            if vids_tp == []:
                raise ValueError("No vids left after filtering!")

        for key, edges in self._vertices.items():
            if key in vids_tp:
                inedges, outedges = edges
                sortedinedges = set([eid for eid in inedges if ((self.source(eid) in vids_tp) and (self.edge_property('edge_type')[eid]=='s'))])
                sortedoutedges = set([eid for eid in outedges if ((self.target(eid) in vids_tp) and (self.edge_property('edge_type')[eid]=='s'))])
                result._vertices.add((sortedinedges,sortedoutedges), key)
                for eid in sortedoutedges:
                    result._edges.add(self._edges[eid], eid)

        if vtx_ppty is not None:
            if isinstance(vtx_ppty, str):
                vtx_ppty = [vtx_ppty]
            for ppty in vtx_ppty:
                try:
                    ppty_values = [p for k,p in self.vertex_property(ppty).iteritems() if k in vids_tp]
                except:
                    print "Could not access to '{}' vertex_property...".format(ppty)
                    pass
                else:
                    add_vertex_ppty_from_label_dict(result, ppty, ppty_values)

        return result

    def greedy_spatial_colormap(self, time_point, vids=None):
        """
        Uses `networkx.greedy_color` to extract a neighbor diverging colormap.
        """
        g = self.spatial_graph_at_time(time_point, vids)
        nxg = g.to_networkx()
        tp_color_dict = greedy_color(nxg)
        return tp_color_dict


    def vertexpair2edge_map(self, vids=None, time_point=None, edge_type=None):
        """
        Compute a dictionary that map pairs of vertex id to edge id.
        It requires the existence of a 'label' property
    
        Returns
        -------
        dictionary {(vid_i, vid_j): eid_ij}
        """
        return dict([((self.source(eid),self.target(eid)),eid) for eid in self.edges(vids, edge_type)])

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
        return dict([(j,i) for i,j in v2e.iteritems()])

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
        missing_vertex = list(set(vid_dict.keys())-set(self.vertices()))
        if missing_vertex != []:
            print "The vid-dictionary '{}' contains {} unknown vertices:\n{}".format(name, len(missing_vertex), missing_vertex)
            vid_dict = dict([(k,v) for k,v in vid_dict.iteritems() if k in self.vertices()])
    
        if name not in self.vertex_properties():
            self.add_vertex_property(name, vid_dict)
        else:
            self.vertex_property(name).update(vid_dict)
        return

    def add_edge_ppty_from_eid_dict(self, name, eid_dict, try_erase_first=False):
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
        if name in self.edge_properties() and not try_erase_first:
            raise ValueError("Existing edge property '{}'".format(name))
        if try_erase_first:
            print "You asked to try to erase property '{}'!".format(name)
            self.remove_edge_property(name)
    
        self.add_edge_property(name, eid_dict)
        return

    def extend_edge_ppty_from_eid_dict(self, name, eid_dict):
        """
        Entend an edge property in the graph using an eid dict {eid: value}.
    
        Parameters
        ----------
        name : str
            the name of the vertex property to add
        eid_dict : dict
            dictionary {eid: ppty_value} to add
    
        Note
        ----
        Use it when transforming 'eid data' during TissueGraph build-up.
        """
        if name not in list(self.edge_property_names()):
            self.add_edge_ppty_from_eid_dict(name, eid_dict)
        else:
            self.edge_property(name).update(eid_dict)
        return

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

    def exist_relative_at_rank(self, vid, rank, verbose=False):
        """
        Check if there is at least a temporal relative of `vid` at given `rank`.
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
        if rank > self.nb_time_points-1 - self.vertex_property('index')[vid]:
            if verbose: print "WARNING: rank > max_self_temporal_depth - vid_index !"
            return False
        # Null rank:
        if rank == 0:
            if verbose: print "WARNING: Null rank !"
            return False
        ## Non-obvious cases:
        if (rank > 0):
            try :
                return self.descendants(vid,rank)-self.descendants(vid,rank-1) != set()
            except:
                return False
        if (rank < 0):
            try :
                return self.ancestors(vid,abs(rank))-self.ancestors(vid,abs(rank)-1) != set()
            except:
                return False

    def exist_all_relative_at_rank(self, vid, rank):
        """
        Check if ALL temporal relative of `vid` exists at given `rank`.
        Starts from the vid's temporal index.
        Positive rank check if `vid` have all its descendants lineaged up to `rank`.
        Negative rank check if `vid` rank-ancestor (from self.vertex_temporal_index(vid)-abs(rank)) is lineaged up to `vid` time-point.
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
        if rank == 0 or abs(rank) > self.nb_time_points-1:
            return False
    
        if rank == 1 or rank == -1:
            return self.exist_relative_at_rank(vid, rank)
    
        if (rank > 0):
            descendants_at_rank = {}
            descendants_at_rank[1] = self.descendant(vid)
            if descendants_at_rank[1] == set():
                return False
            for r in xrange(2,rank+1):
                for v in descendants_at_rank[r-1]:
                    if self.descendant(v) == set():
                        return False
                    if descendants_at_rank.has_key(r):
                        descendants_at_rank[r].update(self.descendant(v))
                    else:
                        descendants_at_rank[r] = self.descendant(v)
            return True
    
        if (rank < 0):
            vid_rank_ancestor = list(self.rank_ancestors(vid,abs(rank)))
            if len(vid_rank_ancestor) == 1:
                vid_rank_ancestor = vid_rank_ancestor[0]
                return self.exist_all_relative_at_rank(vid_rank_ancestor, abs(rank))
            elif len(vid_rank_ancestor) > 1:
                print "More than ONE ancestor?!! Sorry but... WEIRD!!!!"
            else:
                return False
