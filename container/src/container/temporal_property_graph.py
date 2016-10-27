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

import warnings, numpy as np
from property_graph import *

from vplants.tissue_analysis.temporal_graph_analysis import translate_keys_Graph2Image, exist_all_relative_at_rank

class TemporalPropertyGraph(PropertyGraph):
    """
    Simple implementation of PropertyGraph using
    dict as properties and two dictionaries to
    maintain these properties
    """
    STRUCTURAL = 's'
    TEMPORAL = 't'

    def __init__(self, graph=None, **kwds):
        PropertyGraph.__init__(self, graph, idgenerator='max',**kwds)
        self.add_edge_property('edge_type')

        # list of dict
        # each dict define the mapping between the new and the old vid.
        # old label define both graph index and local id.
        self.add_edge_property('old_label')
        self.add_vertex_property('old_label')
        self.add_vertex_property('index')
        self.nb_time_points = 0
        # list of (N_tp) tuples of two dict vertex & edges:
        # (self._old_to_new_ids[tp][0] = old_to_new_vids; self._old_to_new_ids[tp][1] = old_to_new_eids).
        self._old_to_new_ids = []


    def extend(self, graphs, mappings, time_steps = None, disable_lineage_checking=True):
        """
        Extend the structure with graphs and mappings.
        Each graph contains structural edges. 
        Mapping define dynamic edges between two graphs.

        Args:
            - graphs - a list of PropertyGraph
            - mappings - a list defining the dynamic or temporal edges between two graphs.
            - time_steps - time corresponding to each graph

        :warning:: len(graphs) == len(mappings)-1
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
        Append a (spatial) graph to the tpg structure with respect to a given temporal mapping.
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
        relabel_ids = Graph.extend(self,graph)
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
            :Examples:
            >>> list(flatten([1,2,3,4]))
            >>> [1, 2, 3, 4]

            >>> list(flatten([[1,2],[3,4]]))
            >>> [1, 2, 3, 4]

            >>> list(flatten([[1,[2,3]],4]))
            >>> [1, 2, 3, 4]
            """
            import collections
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
            on_ids_source, on_ids_target = self._old_to_new_ids[-2:] #get the last two image2graph dict (label2vertex @ t_n-1 & t_n)
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
        """ Return the temporal index of a vid `self.vertex_property('index')[vid]`."""
        if isinstance(vid, list):
            return [self.vertex_temporal_index(v) for v in vid]
        else:
            return self.vertex_property('index')[vid]

    def children(self, vid):
        """ Return children of the vertex vid

        Args:
          vid: : a vertex id

        Returns:
          children_list: : the set of the children of the vertex vid
        """
        return self.out_neighbors(vid, 't')

    def iter_children(self, vid):
        """ Return children of the vertex vid

        Args:
          vid: : a vertex id

        Returns:
          iterator: : an iterator on the set of the children of the vertex vid
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
          vid: : a vertex id

        Returns:
          parents_list: : the set of the parents of the vertex vid
        """
        return self.in_neighbors(vid, 't')

    def iter_parent(self, vid):
        """ Return parents of the vertex vid

        Args:
          vid: : a vertex id

        Returns:
          iterator: : the set of the children of the vertex vid
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
          vid: : a vertex id

        Returns:
          sibling_list: : the set of sibling of the vertex vid
        """
        if self.parent(vid):
            return self.children(self.parent(vid).pop())-set([vid])
        else:
            return None

    def iter_sibling(self, vid):
        """ Return of the vertex vid

        Args:
          vid: : a vertex id

        Returns:
          iterator: : an iterator on the set of sibling of the vertex vid
        """
        return iter(self.sibling(vid))

    def descendants(self, vids, n = None):
        """ Return the 0, 1, ..., nth descendants of the vertex vid

        Args:
          vids: : a set of vertex id

        Returns:
          descendant_list: : the set of the 0, 1, ..., nth descendant of the vertex vid
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
          vids: : a set of vertex id

        Returns:
          iterator: : an iterator on the set of the 0, 1, ..., nth descendants of the vertex vid
        """
        return iter(self.descendants(vids, n))

    def rank_descendants(self, vid, rank=1):
        """ Return the descendants of the vertex vid only at a given rank
        Args:
          vid: : a vertex id or a set of vertex id
        Returns:
          descendant_list: : the set of the rank-descendant of the vertex vid or a list of set
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
          vid: : a vertex id or a set of vertex id
        Returns:
          descendant_list: : the set of the rank-descendant of the vertex vid or a list of set
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
          vids: : a set of vertex id

        Returns:
          anestors_list: : the set of the 0, 1, ..., nth ancestors of the vertex vid
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
          vids: : a set of vertex id

        Returns:
          iterator: : an iterator on the set of the 0, 1, ..., nth ancestors of the vertex vid
        """
        return iter(self.ancestors(vids, n))

    def rank_ancestors(self, vid, rank=1):
        """ Return the ancestor of the vertex vid only at a given rank
        Args:
          vid: : a vertex id or a set of vertex id
        Returns:
          descendant_list: : the set of the rank-ancestor of the vertex vid or a list of set
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
        flv = self.descendants([k for k in self.vertex_at_time(0) if exist_all_relative_at_rank(self, k, rank)], rank)
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
            vids = [k for k in self.vertices() if (exist_all_relative_at_rank(self, k, lineage_rank) or exist_all_relative_at_rank(self, k, -lineage_rank))]
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
         - 'lineage_rank' (int): usefull if you want to check the lineage for a different rank than the rank-1 temporal neighborhood;
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
         - 'lineage_rank' (int): usefull if you want to check the lineage for a different rank than the rank-1 temporal neighborhood;
           fully_lineaged: (bool) : if True, return vertices linked from the beginning to the end of the graph;
           as_parent: (bool) : if True, return vertices lineaged as parents;
           as_children: (bool) : if True, return vertices lineaged as children.
        """
        vids = self.vertex_at_time(time_point, lineaged, fully_lineaged, as_ancestor, as_descendant)
        return dict([(k,self.vertex_property(vertex_property)[k]) for k in vids if self.vertex_property(vertex_property).has_key(k)])

    def vertex_property_with_image_labels(self, vertex_property, time_point, lineaged=False, fully_lineaged=False, as_ancestor=False, as_descendant=False, lineage_rank=1):
        """
        Return a subpart of graph.vertex_property(`vertex_property`) with relabelled keys into "images labels" thanks to the dictionary graph.vertex_property('old_labels').
        Since "images labels" can be similar (not unique), it is mandatory to give a `time_point`.
        Additional parameters can be given and are related to the 'self.vertex_property_at_time' parameters.

        Args:
           vertex_property: (str): a string refering to an existing 'graph.vertex_property' to extract;
           time_point` (int): time-point for which to return the `vertex_property:;
           lineaged: (bool) : if True, return vertices having at least a parent or a child(ren);
         - 'lineage_rank' (int): usefull if you want to check the lineage for a different rank than the rank-1 temporal neighborhood;
           fully_lineaged: (bool) : if True, return vertices linked from the beginning to the end of the graph;
           as_parent: (bool) : if True, return vertices lineaged as parents;
           as_children: (bool) : if True, return vertices lineaged as children.
        Returns:
         - *key_ vertex/cell label, *values_ `vertex_property`

        :Examples:
         graph.vertex_property_with_image_labels('volume', 0)

        """
        return translate_keys_Graph2Image(self, self.vertex_property_at_time(vertex_property, time_point, lineaged, fully_lineaged, as_ancestor, as_descendant), time_point)

    def edge_property_at_time(self, edge_property, time_point):
        """
        Return a subpart of graph.edge_property(`edge_property`) at a given time-point.

        Args:
           edge_property: (str): a string refering to an existing 'graph.edge_property' to extract;
           time_point` (int): time-point for which to return the `edge_property:;

        :Examples:
         graph.edge_property_with_image_labelpairs('wall_area', 0)
        """
        eid2vidpair = edge2vertexpair_map(self, time_point)
        tmp_dict = {}
        for eid,(label1,label2) in eid2vidpair.iteritems():
            if self.edge_property(edge_property).has_key(eid):
                tmp_dict[eid] = self.edge_property(edge_property)[eid]
        return tmp_dict

    def edge_property_with_vids_labelpairs(self, edge_property, time_point):
        """
        Return a subpart of graph.edge_property(`edge_property`) with relabelled key-pair into "vids labelpairs" thanks to the dictionary from `edge2vertexpair_map`.

        Args:
           edge_property: (str): a string refering to an existing 'graph.edge_property' to extract;
           time_point` (int): time-point for which to return the `edge_property:;

        :Examples:
         graph.edge_property_with_image_labelpairs('wall_area', 0)
        """
        eid2vidpair = edge2vertexpair_map(self, time_point)
        tmp_dict = {}
        for eid,(label1,label2) in eid2vidpair.iteritems():
            if self.edge_property(edge_property).has_key(eid):
                tmp_dict[(label1,label2)] = self.edge_property(edge_property)[eid]
        return tmp_dict

    def edge_property_with_image_labelpairs(self, edge_property, time_point):
        """
        Return a subpart of graph.edge_property(`edge_property`) with relabelled key-pair into "images labelpairs" thanks to the dictionary graph.vertex_property('old_labels').

        Args:
           edge_property: (str): a string refering to an existing 'graph.edge_property' to extract;
           time_point` (int): time-point for which to return the `edge_property:;

        :Examples:
         graph.edge_property_with_image_labelpairs('wall_area', 0)
        """
        eid2labelpair = edge2labelpair_map(self, time_point)
        return dict([(tuple(sorted([label1,label2])),self.edge_property(edge_property)[eid]) for eid,(label1,label2) in eid2labelpair.iteritems() if self.edge_property(edge_property).has_key(eid)])

    def domain_vids(self, domain_name):
        """
        Return a list of vids (TPG vertex id type) that belong to the domain `domain_name` according to graph.
        Args:
           domain_name` (str) : the name of a domain (e.g. from `self.add_vertex_to_domain():)
        """
        if 'domains' in list(self.vertex_properties()):
            return sorted(list(set([k for k,v in self.vertex_property('domains').iteritems() for r in v if r==domain_name])))
        else:
            print "No property 'domains' added to the graph yet!"
            return None



def label2vertex_map(graph, time_point = None):
    """
        Compute a dictionary that map label to vertex id.
        It requires the existence of a 'label' vertex property

        :rtype: dict
    """
    if isinstance(graph, TemporalPropertyGraph):
        assert time_point is not None
        return dict([(j,i) for i,j in graph.vertex_property('label').iteritems() if graph.vertex_property('index')[i]==time_point])
    else:
        return dict([(j,i) for i,j in graph.vertex_property('label').iteritems()])

def vertex2label_map(graph, time_point = None):
    """
        Compute a dictionary that map label to vertex id.
        It requires the existence of a 'label' vertex property

        :rtype: dict
    """
    if isinstance(graph, TemporalPropertyGraph):
        assert time_point is not None
        return dict([(i,j) for i,j in graph.vertex_property('label').iteritems() if graph.vertex_property('index')[i]==time_point])
    else:
        return dict([(i,j) for i,j in graph.vertex_property('label').iteritems()])

def label2vertex(graph, labels, time_point = None):
    """
        Translate label as vertex id.
        It requires the existence of a 'label' vertex property

        :rtype: dict
    """
    label2vertexmap = label2vertex_map(graph, time_point)
    if isinstance(labels,list):
        return [label2vertexmap[label] for label in labels]
    else:
        return label2vertexmap[labels]

def labelpair2edge_map(graph, time_point = None):
    """
        Compute a dictionary that map pair of old_labels to edge id.
        It requires the existence of a 'label' property

        :rtype: dict
    """
    mvertex2label = vertex2label_map(graph, time_point)

    return dict([((mvertex2label[graph.source(eid)],mvertex2label[graph.target(eid)]),eid) for eid in graph.edges()
     if (mvertex2label.has_key(graph.source(eid)) and mvertex2label.has_key(graph.target(eid)))] )

def edge2labelpair_map(graph, time_point):
    """
        Compute a dictionary that map pair of old_labels to edge id.
        It requires the existence of a 'label' property

        :rtype: dict
    """
    mvertex2label = vertex2label_map(graph, time_point)

    return dict([(eid, (mvertex2label[graph.source(eid)],mvertex2label[graph.target(eid)])) for eid in graph.edges()
     if (mvertex2label.has_key(graph.source(eid)) and mvertex2label.has_key(graph.target(eid)))] )

def vertexpair2edge_map(graph):
    """
        Compute a dictionary that map pair of vertex id to edge id.
        It requires the existence of a 'label' property

        :rtype: dict
    """
    return dict([((graph.source(eid),graph.target(eid)),eid) for eid in graph.edges()])

def edge2vertexpair_map(graph, time_point = None):
    """
        Compute a dictionary that map pair of vertex id to edge id.
        It requires the existence of a 'label' property

        :rtype: dict
    """
    if time_point is None:
        return dict([(eid,(graph.source(eid),graph.target(eid))) for eid in graph.edges()])
    else:
        e2v = {}
        for eid in graph.edges():
            vid1 = graph.source(eid)
            vid2 = graph.target(eid)
            if (graph.vertex_temporal_index(vid1) == time_point) and (graph.vertex_temporal_index(vid2) == time_point):
                e2v[eid] = sorted((vid1, vid2))
        return e2v

def add_vertex_property_from_dictionary(graph, name, dictionary, mlabel2vertex = None, time_point = None, overwrite = False):
    """
        Add a vertex property with name 'name' to the graph build from an image.
        The values of the property are given as by a dictionary where keys are vertex labels.
        If overwrite is true, the property will be removed before adding the new values!
    """
    if isinstance(graph, TemporalPropertyGraph):
        assert time_point is not None
    if mlabel2vertex is None:
        mlabel2vertex = label2vertex_map(graph, time_point)
    if name in graph.vertex_properties() and not overwrite:
        raise ValueError("Existing vertex property '{}'".format(name))
    if overwrite:
        print "You asked to overwrite property '{}', it will be removed first!".format(name)
        graph.remove_vertex_property(name)

    graph.add_vertex_property(name)
    graph.vertex_property(name).update( dict([(mlabel2vertex[k], dictionary[k]) for k in dictionary]) )
    return "Done."

def add_vertex_property_from_label_and_value(graph, name, labels, property_values, mlabel2vertex = None, overwrite = False):
    """
        Add a vertex property with name 'name' to the graph build from an image.
        The values of the property are given as two lists.
        First one gives the label in the image and second gives the value of the property.
        Labels are first translated in id of the graph and values are assigned to these ids in the graph
        If overwrite is true, the property will be removed before adding the new values!
    """
    if mlabel2vertex is None:
        mlabel2vertex = label2vertex_map(graph)
    if name in graph.vertex_properties() and not overwrite:
        raise ValueError("Existing vertex property '{}'".format(name))
    if overwrite:
        print "You asked to overwrite property '{}', it will be removed first!".format(name)
        graph.remove_vertex_property(name)

    graph.add_vertex_property(name)
    graph.vertex_property(name).update(dict([(mlabel2vertex[i], v) for i,v in zip(labels,property_values)]))
    return "Done."

def add_vertex_property_from_label_property(graph, name, label_property, mlabel2vertex = None, overwrite = False):
    """
        Add a vertex property with name 'name' to the graph build from an image.
        The values of the property are given as a dictionnary associating a label and a value.
        Labels are first translated in id of the graph and values are assigned to these ids in the graph
        If overwrite is true, the property will be removed before adding the new values!
    """
    if mlabel2vertex is None:
        mlabel2vertex = label2vertex_map(graph)
    if name in graph.vertex_properties() and not overwrite:
        raise ValueError("Existing vertex property '{}'".format(name))
    if overwrite:
        print "You asked to overwrite property '{}', it will be removed first!".format(name)
        graph.remove_vertex_property(name)

    graph.add_vertex_property(name)
    graph.vertex_property(name).update(dict([(mlabel2vertex[i], v) for i,v in label_property.iteritems()]))
    return "Done."

def add_edge_property_from_dictionary(graph, name, dictionary, mlabelpair2edge = None, overwrite = False):
    """
        Add an edge property with name 'name' to the graph build from an image.
        The values of the property are given as by a dictionary where keys are vertex labels.
        If overwrite is true, the property will be removed before adding the new values!
    """
    if mlabelpair2edge is None:
        mlabelpair2edge = labelpair2edge_map(graph)
    if name in graph.edge_properties() and not overwrite:
        raise ValueError("Existing edge property '{}'".format(name))
    if overwrite:
        print "You asked to overwrite property '{}', it will be removed first!".format(name)
        graph.remove_edge_property(name)

    graph.add_edge_property(name)
    graph.edge_property(name).update( dict([(mlabelpair2edge[k], dictionary[k]) for k in dictionary]) )
    return "Done."

def add_edge_property_from_eid_dictionary(graph, name, dictionary, overwrite = False):
    """
        Add an edge property with name 'name' to the graph build from an image.
        The values of the property are given as by a dictionary where keys are vertex labels.
        If overwrite is true, the property will be removed before adding the new values!
    """
    if name in graph.edge_properties() and not overwrite:
        raise ValueError("Existing edge property '{}'".format(name))
    if overwrite:
        print "You asked to overwrite property '{}', it will be removed first!".format(name)
        graph.remove_edge_property(name)

    graph.add_edge_property(name)
    graph.edge_property(name).update(dictionary)
    return "Done."

def add_edge_property_from_label_and_value(graph, name, label_pairs, property_values, mlabelpair2edge = None, overwrite = False):
    """
        Add an edge property with name 'name' to the graph build from an image.
        The values of the property are given as two lists.
        First one gives the pair of labels in the image that are connected and the second list gives the value of the property.
        Pairs of labels are first translated in edge ids of the graph and values are assigned to these ids in the graph
        If overwrite is true, the property will be removed before adding the new values!
    """
    if mlabelpair2edge is None:
        mlabelpair2edge = labelpair2edge_map(graph)
    if name in graph.edge_properties() and not overwrite:
        raise ValueError("Existing edge property '{}'".format(name))
    if overwrite:
        print "You asked to overwrite property '{}', it will be removed first!".format(name)
        graph.remove_edge_property(name)

    graph.add_edge_property(name)
    graph.edge_property(name).update(dict([(mlabelpair2edge[labelpair], value) for labelpair,value in zip(label_pairs,property_values)]))
    return "Done."

def add_edge_property_from_label_property(graph, name, labelpair_property, mlabelpair2edge = None, overwrite = False):
    """
        Add an edge property with name 'name' to the graph build from an image.
        The values of the property are given as a dictionnary associating a pair of label and a value.
        Pairs of labels are first translated in edge ids of the graph and values are assigned to these ids in the graph
        If overwrite is true, the property will be removed before adding the new values!
    """
    if mlabelpair2edge is None:
        mlabelpair2edge = labelpair2edge_map(graph)
    if name in graph.edge_properties() and not overwrite:
        raise ValueError("Existing edge property '{}'".format(name))
    if overwrite:
        print "You asked to overwrite property '{}', it will be removed first!".format(name)
        graph.remove_edge_property(name)

    graph.add_edge_property(name)
    graph.edge_property(name).update(dict([(mlabelpair2edge[labelpair], value) for labelpair,value in labelpair_property.iteritems()]))
    return "Done."

def extend_edge_property_from_dictionary(graph, name, dictionary, time_point = None, mlabelpair2edge = None):
    """
        Add an edge property with name 'name' to the graph build from an image.
        The values of the property are given as by a dictionary where keys are labelpairs.
    """
    if isinstance(graph, TemporalPropertyGraph):
        assert time_point is not None
    if mlabelpair2edge is None:
        mlabelpair2edge = labelpair2edge_map(graph, time_point)
    if name not in graph.edge_properties():
        graph.add_edge_property(name)

    print "The dictionary extending edge property '{}' contains {} tuples, against {} edges for time-point #{} !".format(name, len(dictionary.keys()), len(mlabelpair2edge.keys()), time_point)
    edges_not_found = []
    for k in dictionary:
        sk = tuple(sorted(k)) # make sure the keys of the `dictionary` are sorted tuples
        if sk in mlabelpair2edge.keys():
            graph.edge_property(name).update({mlabelpair2edge[sk]: dictionary[k]})
        else:
            edges_not_found.append(sk)
        
    if len(edges_not_found) != 0:
        print "{} vid-pairs have no edges between them in the graph: {}".format(len(edges_not_found), edges_not_found)
    return "Done."

def extend_vertex_property_from_dictionary(graph, name, dictionary, mlabel2vertex = None, time_point = None):
    """
        Add a vertex property with name 'name' to the graph build from an image.
        The values of the property are given as by a dictionary where keys are labels.
    """
    if isinstance(graph, TemporalPropertyGraph):
        assert time_point is not None
    if mlabel2vertex is None:
        mlabel2vertex = label2vertex_map(graph, time_point)
    if name not in graph.vertex_properties():
        graph.add_vertex_property(name)

    missing_vertex = list(set(dictionary.keys())-set(mlabel2vertex.keys()))
    if missing_vertex != []:
        print "The dictionary extending vertex property '{}' contains {} labels, against {} vertices for time-point #{} !".format(name, len(dictionary.keys()), len(mlabel2vertex.keys()), time_point)
    graph.vertex_property(name).update( dict([(mlabel2vertex[k], dictionary[k]) for k in dictionary if k in mlabel2vertex.keys()]) )
    return "Done."

def extend_graph_property_from_dictionary(graph, name, dictionary):
    """
    """
    if name not in graph.graph_properties():
        print "Adding graph_property '{}'".format(name)
        graph.add_graph_property(name)

    graph.graph_property(name).update(dictionary)
    return "Done."

def add_bool_atlas_data(graph, time_sorted_exp_patterns_fnames, exp_patterns_time_steps, exp_pattern_as_cid=True, patterns_path=None):
    """
    Add ATLAS data, "boolean" data giving all cells (values;cids/vids) expressing a given gene (keys).
    
    Args:
       graph: (TemporalPropertyGraph) - graph to complete
       time_sorted_exp_patterns_fnames: (list) - list of strings, defining filenames of expression patterns;
       exp_patterns_time_steps: (list) - list of int, defining to which time-step the filenames of expression patterns relates to;
       exp_pattern_as_cid: (bool) - says weither the cell ids are given as cids (image) or vids (tpg);
       patterns_path: (str) - path where to find the filenames of expression patterns;
    """
    from vplants.tissue_analysis.misc import load_pickled_object
    from vplants.tissue_analysis.temporal_graph_analysis import translate_ids_Image2Graph
    # We list the time-steps to manually associate the correct time-point ID to each ATLAS:
    time_steps = graph.graph_property('time_steps')
    atlas_tp = [time_steps.index(ts) for ts in exp_patterns_time_steps]
    print "Detected 'time-steps' association with 'patterns files':\n{}".format(zip([time_steps[ts] for ts in atlas_tp],time_sorted_exp_patterns_fnames))

    ATLAS_genes = []
    # Loop over the ATLAS to transfert genes expression levels to the TPG:
    for tp, patterns_name in zip(atlas_tp, time_sorted_exp_patterns_fnames):
        print "\n"
        gene_dict = load_pickled_object(patterns_path+patterns_name)
        print "Detected the following list of genes: {}.".format(gene_dict.keys())
        ATLAS_genes += gene_dict.keys()
        for gene_name, cids in gene_dict.iteritems():
            gene_name = str(gene_name)
            print "Exporting '{}' pattern at time point {}...".format(gene_name, tp)
            print "Found {} initial cell-ids".format(len(cids)),
            if exp_pattern_as_cid:
                vids = translate_ids_Image2Graph(graph, cids, tp)
                ncids = len(cids); nvids=len(vids)
                print "and {}% ({}/{}) could be translated into vertex-ids!".format(round(float(nvids)/ncids,3)*100, nvids, ncids)
            else:
                vids = cids
            graph.add_vertex_to_domain(vids, gene_name)

    return graph


def iterable(obj):
    try :
        iter(obj)
        return True
    except TypeError,te:
        return False
