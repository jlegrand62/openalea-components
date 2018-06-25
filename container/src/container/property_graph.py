# -*- python -*-
#
#       OpenAlea.Core
#
#       Copyright 2006-2009 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#                       Fred Theveny <frederic.theveny@cirad.fr>
#                       Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://openalea.gforge.inria.fr
#
################################################################################
"""This module provide a set of concepts to add properties to graph elements"""

__license__ = "Cecill-C"
__revision__ = " $Id$ "

import warnings

from graph import Graph
from interface.property_graph import IPropertyGraph, PropertyError

try:
    import networkx as nx
except ImportError:
    raise ImportError(
        "NetworkX library cannot be found! Please install package 'python-networkx'.")

VertexProperty, EdgeProperty, GraphProperty = range(3)
VertexIdType, EdgeIdType, ValueType = range(3)

MISSING_PPTY = "Property '{}' is NOT defined as {} property!"
EXISTING_PPTY = "Property '{}' is already defined as {} property!"
EMPTY_PPTY = "Creating EMPTY {} property named '{}'!"
EMPTY_EXTEND = "Trying to extend {} property '{}' with empty {}!"
NULL_VALUES = "Trying to update {} property '{}' with NULL value!"

class PropertyGraph(IPropertyGraph, Graph):
    """
    Simple implementation of 'IPropertyGraph'.
    Extend class 'Graph' to add vertex, edge and graph properties.
    
    Usage
    -----
    The different "properties" are stored as dictionary under:
      * self._vertex_property save vertex properties with vertex-ids as keys;
      * self._edge_property save edge properties with edge-ids as keys;
      * self._graph_property save graph properties with a string as key.

    To define some keys or values of "graph properties" as vertex-ids or 
    edge-ids, use one of the following:
      - 'self.set_graph_property_key_to_vid_type()'
      - 'self.set_graph_property_key_to_eid_type()'
      - 'self.set_graph_property_value_to_vid_type()'
      - 'self.set_graph_property_value_to_eid_type()'
    If keys or values are defined as vid or eid type for a graph property, they 
    will be translated in case of object "extension", ie. when calling
    'self.extend_property_graph()'.
    
    Notes
    -----
    Defining keys or values as vid or eid type (for a graph property) is also 
    used for "temporal extension" by TemporalPropertyGraph, allowing to re-label
    them with the new vids ('old_to_new_vids') and eids ('old_to_new_eids').
    """
    meta_vid_key_type = "key_as_vid"
    meta_eid_key_type = "key_as_eid"
    meta_vid_value_type = "value_as_vid"
    meta_eid_value_type = "value_as_eid"

    def __init__(self, graph=None, **kwargs):
        """
        Parameters
        ----------
        graph : Graph|PropertyGraph, optional
            by default (None), initialize an empty PropertyGraph, else try to
            re-use the '_vertex_property', '_edge_property' & '_graph_property'
            attributes of the given object!

        kwargs
        ------
        idgenerator: 'list'|'max'|'set'
            indicate how the vertex and edge ids will be generated
            by the object, by default use 'set'
        verbose: bool
            control the verbosity of the object initialization, by default use
            silent mode (ie. False)
        """
        # - Parsing keyword-arguments:
        idgenerator = kwargs.get('idgenerator', "set")
        verbose = kwargs.get('verbose', False)
        # - Call 'Graph' object init:
        Graph.__init__(self, graph, idgenerator)
        # - Add the new object attributes:
        # '_vertex_property', '_edge_property' & '_graph_property'
        # TODO: also create '_vertex_property_unit', '_edge_property_unit' & '_graph_property_unit' ??
        if graph is None:
            self._vertex_property = {}
            self._edge_property = {}
            self._graph_property = {}
            if verbose:
                print "Initialized empty PropertyGraph object!"
        else:
            try:
                self._vertex_property = graph._vertex_property
                self._edge_property = graph._edge_property
                self._graph_property = graph._graph_property
            except AttributeError:
                raise TypeError("Unknown type of graph '{}', missing vertex, edge or graph property attribute!".format(type(graph)))
            else:
                if verbose:
                    print "Initialized ", self.__str__()

    def __str__(self):
        """
        Format returned object information.
        """
        s = "'PropertyGraph' object containing:\n"
        s += "  - {} vertices\n".format(len(self._vertices))
        s += "  - {} edges\n".format(len(self._edges))
        s += "  - {} vertex properties\n".format(len(self._vertex_property))
        s += "  - {} edge properties\n".format(len(self._edge_property))
        s += "  - {} graph properties\n".format(len(self._graph_property))
        return s

    ############################################################################
    # Overridden methods of Graph object:
    ############################################################################

    def remove_vertex(self, vid):
        """
        Remove vertex if 'vid' from the object.
        It is also removed from any vertex property dictionary!

        Parameters
        ----------
        vid : int
            id of the vertex to remove

        Returns
        -------
        Nothing, edit object
        """
        # Remove the vertex 'vid' from the object using 'Graph' function
        Graph.remove_vertex(self, vid)
        # Remove the key 'vid', and its associated value, for all properties:
        for prop in self._vertex_property.itervalues():
            prop.pop(vid, None)
        return
    # remove_vertex.__doc__ = Graph.remove_vertex.__doc__

    def remove_edge(self, eid):
        """
        Remove the edge 'eid' from the object.
        It is also removed from any edge property dictionary!

        Parameters
        ----------
        eid : int
            id of the edge to remove

        Returns
        -------
        Nothing, edit object
        """
        # Remove the edge 'eid' from the object using 'Graph' function
        Graph.remove_edge(self, eid)
        # Remove the key 'eid', and its associated value, for all properties:
        for e_ppty in self._edge_property.itervalues():
            e_ppty.pop(eid, None)
        return
    # remove_edge.__doc__ = Graph.remove_edge.__doc__

    def clear(self):
        """
        Clear the object of all vetex, edge and graph properties.

        Returns
        -------
        Nothing, edit object
        """
        Graph.clear(self)
        # Clear vertex properties:
        vid_ppty = list(self.vertex_property_names())
        for ppty in vid_ppty:
            self.remove_vertex_property(ppty)
        # Clear edge properties:
        eid_ppty = list(self.edge_property_names())
        for ppty in eid_ppty:
            self.remove_edge_property(ppty)
        # Clear graph properties:
        graph_ppty = list(self.graph_property_names())
        for ppty in graph_ppty:
            self.remove_graph_property(ppty)
        return
    # clear.__doc__ = Graph.clear.__doc__

    def clear_edges(self):
        """
        Remove all edges from the object.

        Returns
        -------
        Nothing, edit object
        """
        Graph.clear_edges(self)
        for prop in self._edge_property.itervalues():
            prop.clear()
        return
    # clear_edges.__doc__ = Graph.clear_edges.__doc__

    ############################################################################
    # VERTEX associated properties:
    ############################################################################

    def vertex_property_names(self):
        """
        Return a key-iterator of vertex property names in the object.

        Returns
        -------
        iterator of self._vertex_property.keys()
        """
        return self._vertex_property.iterkeys()
    # vertex_property_names.__doc__ = IPropertyGraph.vertex_property_names.__doc__

    def vertex_properties(self):
        """
        Returns all the vertex properties in the object as a dictionary
        with as keys, the property names and as values, a dictionary of corresponding
        vid values

        Returns
        -------
        a dictionary {'ppty_name': {vid: vid_ppty_value}}
        """
        return self._vertex_property
    # vertex_properties.__doc__ = IPropertyGraph.vertex_properties.__doc__

    def vertex_property(self, ppty_name, vids=None):
        """
        Return a dictionary of vertex values for given property name, can be 
        filtered by vertex-ids (`vids`) if not None.
        If an integer is given as `vids`, only its value is returned!

        Parameters
        ----------
        ppty_name: str
            the name of an existing vertex property
        vids: int|list|set|None, optional
            by default (None), return all values associated to this property;
            list/set of vertex-ids used to filter the keys of the returned 
            dictionary;
            an integer returns a value;

        Returns
        -------
        the associated values if `vids` is an integer, else the vertex-ids 
        dictionary for given `vids`
        """
        # - Assert the property exists:
        if not self._vertex_property.has_key(ppty_name):
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'vertex'))
        # - Return it for desired vid keys:
        ppty = self._vertex_property[ppty_name]
        if isinstance(vids, int):
            return ppty[vids]
        elif isinstance(vids, list) or isinstance(vids, set):
            return {k: v for k, v in ppty.items() if k in vids}
        else:
            return ppty
    # vertex_property.__doc__=IPropertyGraph.vertex_property.__doc__

    ############################################################################
    # EDGE associated properties:
    ############################################################################

    def edge_property_names(self):
        """
        Return a key-iterator of edges property names in the object.

        Returns
        -------
        iterator of self._edges_property.keys()
        """
        return self._edge_property.iterkeys()
    # edge_property_names.__doc__ = IPropertyGraph.edge_property_names.__doc__

    def edge_properties(self):
        """
        Returns all the edges properties in the object as a dictionary
        with as keys, the property names and as values, a dictionary of corresponding
        eid values

        Returns
        -------
        a dictionary {'ppty_name': {eid: eid_ppty_value}}, for all edge properties
        """
        return self._edge_property
    #  edge_properties.__doc__ = IPropertyGraph. edge_properties.__doc__

    def edge_property(self, ppty_name, eids=None):
        """
        Return a dictionary of edge values for given property name, can be 
        filtered by edge-ids (`eids`) if not None.
        If an integer is given as `eids`, only its value is returned!

        Parameters
        ----------
        ppty_name: str
            the name of an existing edge property;
        eids: int|list|set|None, optional
            by default (None), return all values associated to this property
            list/set of edge-ids used to filter the keys of the returned 
            dictionary;
            an integer returns a value;

        Returns
        -------
        the associated values if `eids` is an integer, else the edge-ids 
        dictionary for given `eids`
        """
        # - Assert the property exists:
        if not self._edge_property.has_key(ppty_name):
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'edge'))
        # - Return it for desired eid keys:
        ppty = self._edge_property[ppty_name]
        if isinstance(eids, int):
            return ppty[eids]
        elif isinstance(eids, list) or isinstance(eids, set):
            return {k: v for k, v in ppty.items() if k in eids}
        else:
            return ppty
    # edge_property.__doc__ = IPropertyGraph.edge_property.__doc__

    ############################################################################
    # GRAPH associated properties:
    ############################################################################

    def graph_property_names(self):
        """
        Return a key-iterator of graph property names in the object.

        Returns
        -------
        iterator of self._graph_property.keys()
        """
        return self._graph_property.iterkeys()

    def graph_properties(self):
        """
        Returns all the graph properties in the object as a dictionary
        with as keys, the graph property names and as values, a dictionary of
        corresponding graph values

        Returns
        -------
        a dictionary {'ppty_name': graph_property} for all graph properties
        """
        return self._graph_property

    def graph_property(self, ppty_name):
        """
        Return the corresponding graph property value.

        Parameters
        ----------
        ppty_name: str
            the name of an existing graph property

        Returns
        -------
        self._graph_property[ppty_name]
        """
        try:
            return self._graph_property[ppty_name]
        except KeyError:
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'graph'))
    # graph_property.__doc__ = IPropertyGraph.graph_property.__doc__

    ############################################################################
    # Hidden functions:
    ############################################################################

    def _ids_type(self, ids_type):
        """
        Hidden function using 'edge' or 'vertex' string to return list of
        corresponding vertex or edge ids, as well as the plural form.

        Parameters
        ----------
        ids_type: str
            type of dictionary, can be either 'edge' or 'vertex'

        Returns
        -------
        ids : list
            the list of corresponding vertex or edge ids
        t : str
            the plural form of the given id type
        """
        if ids_type == "vertex":
            ids = self.vertices()
            t = "vertices"
        elif ids_type == "edge":
            ids = self.edges()
            t = "edges"
        else:
            raise ValueError("Unknown 'ids_type': {}!".format(ids_type))
        return ids, t

    def _missing_ids(self, ppty_dict, ids_type):
        """
        Hidden function used to check vertex/edge-ids definition in given 
        property dictionary `ppty_dict`.

        Parameters
        ----------
        ppty_dict: dict
            an eid/vid based dictionary
        ids_type: str
            type of dictionary, can be either 'edge' or 'vertex'

        Returns
        -------

        """
        ids, t = self._ids_type(ids_type)
        
        missing_vertex = list(set(ppty_dict.keys()) - set(ids))
        if missing_vertex:
            n = len(missing_vertex)
            print "Given dictionary contains {} unknown {}!".format(n, t)
            return {k: v for k, v in ppty_dict.items() if k in ids}
        else:
            return ppty_dict

    def _input_dict(self, ppty_name, ppty_dict, ids_type, none_ok):
        """
        Hidden function to tests given `ppty_dict`.
        
        Parameters
        ----------
        ppty_name: str
            the name of the edge/vertex property
        ppty_dict: dict
            the dictionary to test
        ids_type: str
            type of dictionary, can be either 'edge' or 'vertex'
        none_ok: bool
            define if given `ppty_dict` can be None and converted to empty 
            dictionary

        Returns
        -------
        
        """
        ids, t = self._ids_type(ids_type)

        if ppty_dict is None:
            if none_ok:
                print EMPTY_PPTY.format(t, ppty_name)
                return {}
            else:
                raise TypeError(EMPTY_PPTY.format(t, ppty_name))
        
        if not isinstance(ppty_dict, dict):
            raise TypeError("Input '{}' is not a {} dictionary!".format(ppty_name, ids_type))
        
        # - Return a dictionary with existing eids/vids only:
        return self._missing_ids(ppty_dict, ids_type)

    def _extend_ppty(self, ppty_name, ppty_dict, id_dict, ids_type):
        """
        Performs the extension of an edge or vertex property.

        Parameters
        ----------
        ppty_name : str
            name of the property, used for printing
        ppty_dict : dict
            edge or vertex id dictionary to extend
        id_dict : dict
            the edge or vertex id dictionary to use for extension
        ids_type : str
            type of property, to choose among 'egde' and 'vertex'

        Returns
        -------
        ppty_dict : dict
            the extended dictionary
        """
        duplicated_ids = [k for k in id_dict.keys() if k in ppty_dict.keys()]
        new_ids = list(set(id_dict.keys()) - set(duplicated_ids))
        if new_ids:
            ppty_dict.update({k: id_dict[k] for k in new_ids})
            print "Extended {} property '{}' with {} values!".format(ids_type,
                                                                     ppty_name,
                                                                     len(
                                                                         new_ids))
        else:
            print "No new values found to extend {} property '{}'!".format(
                ids_type, ppty_name)
        # - Print if exist some duplicated keys:
        if duplicated_ids:
            print "Found {} {}-ids (out of {}) with a defined value for property '{}'.".format(
                len(duplicated_ids), ids_type, len(id_dict), ppty_name)

        return ppty_dict

    ############################################################################
    # VERTEX property creation and edition functions :
    ############################################################################

    def add_vertex_property(self, ppty_name, vid_dict=None):
        """
        Add a new vertex property, with vid_dict if it is a vid dictionary, or leaves
        it empty if None.

        Parameters
        ----------
        ppty_name : str
            the name of the property to add
        vid_dict : dict|None
            a dictionary with vid as keys, None will add an empty property

        Returns
        -------
        Nothing, edit object
        """
        # - Assert the property does NOT exist:
        if ppty_name in self._vertex_property:
            raise PropertyError(EXISTING_PPTY.format(ppty_name, 'vertex'))
        # - Check the input 'vid_dict', and that all keys are in the graph:
        vid_dict = self._input_dict(ppty_name, vid_dict, 'vertex', True)
        # - Finally add the vertex property to the graph:
        self._vertex_property[ppty_name] = vid_dict
        return
    # add_vertex_property.__doc__ = IPropertyGraph.add_vertex_property.__doc__

    def extend_vertex_property(self, ppty_name, vid_dict):
        """
        Extend an existing vertex property 'ppty_name' using the vid 
        dictionary 'vid_dict'.

        Parameters
        ----------
        ppty_name : str
            a string matching an existing property
        vid_dict : dict
            a dictionary {vid: vid_vid_dict}

        Returns
        -------
        Nothing, edit object

        Notes
        ----
        * 'ppty_name' should exist
        * 'vid_dict' cannot be an empty dictionary
        * if a vid already has a value associated fo vertex property
        'ppty_name' we do not update this value
        """
        # - Assert the property does exist:
        if ppty_name not in self._vertex_property:
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'vertex'))
        # - Check the input 'vid_dict', and that all keys are in the graph:
        vid_dict = self._input_dict(ppty_name, vid_dict, 'vertex', False)
        # - Extend the vertex property 'ppty_name':
        ppty_dict = self._vertex_property[ppty_name]
        ppty_dict = self._extend_ppty(ppty_name, ppty_dict, vid_dict, 'vertex')
        self._vertex_property[ppty_name] = ppty_dict
        return

    def remove_vertex_property(self, ppty_name):
        """
        Remove the vertex property 'ppty_name' of the object.

        Parameters
        ----------
        ppty_name : str
            name of the vertex property to remove

        Returns
        -------
        Nothing, edit object
        """
        # - Try to remove any existing associated graph property 'unit':
        try:
            del self._graph_property['units'][ppty_name]
        except KeyError:
            pass
        # - Try to remove the vertex property, and if so, print about it:
        try:
            n = len(self.vertex_property(ppty_name))
            del self._vertex_property[ppty_name]
        except KeyError:
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'vertex'))
        else:
            print "Removed vertex property '{}' (n={})".format(ppty_name, n)
        return
    # remove_vertex_property.__doc__ = IPropertyGraph.remove_vertex_property.__doc__

    def update_vertex_property(self, ppty_name, values, add_missing=True):
        """
        Update an existing vertex property, with non-null values.
        Can add missing property (instead of failing) if 'add_missing' is True.

        Parameters
        ----------
        ppty_name: str
            the name of the property to add
        values: any non-null
            any type of values or object, except null values like empty list,
            set, dictionary or 'None'
        add_missing: bool, optional
            if True (default), add the property `ppty_name` if missing,
            else may fail if not defined

        Returns
        -------
        Nothing, edit object
        """
        # - Assert the property does exist:
        if ppty_name not in self._vertex_property:
            if add_missing:
                self.add_vertex_property(ppty_name, values)
            else:
                raise PropertyError(MISSING_PPTY.format(ppty_name, 'vertex'))
        # - Raise an error when no 'values' are given:
        if not values:
            raise PropertyError(NULL_VALUES.format('vertex', ppty_name))
        # - Update the vertex property dictionary:
        self._vertex_property[ppty_name].update(values)
        return

    def clear_vertex_property(self, ppty_name):
        """
        Clear any values associated with the vertex property 'ppty_name'.

        Parameters
        ----------
        ppty_name : str
            name of the vertex property to clear

        Returns
        -------
        Nothing, edit object
        """
        # - Try to clear any existing associated graph property 'unit':
        try:
            self._graph_property['units'][ppty_name] = None
        except KeyError:
            pass
        # - Try to clear the vertex property, and if so, print about it:
        try:
            n = len(self.vertex_property(ppty_name))
            self._vertex_property[ppty_name] = {}
        except KeyError:
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'vertex'))
        else:
            print "Cleared vertex property '{}' of n={} values!".format(ppty_name, n)
        return

    ############################################################################
    # EDGE property creation and edition functions :
    ############################################################################

    def add_edge_property(self, ppty_name, eid_dict=None):
        """
        Add a new edge property, with eid_dict if it is an eid dictionary, or leaves
        it empty if None.

        Parameters
        ----------
        ppty_name : str
            the name of the property to add
        eid_dict : dict|None
            a dictionary with eid as keys, None will add an empty property

        Returns
        -------
        Nothing, edit object
        """
        # - Assert the property does NOT exist:
        if ppty_name in self._edge_property:
            raise PropertyError(EXISTING_PPTY.format(ppty_name, 'edge'))
        # - Check the input 'eid_dict', and that all keys are in the graph:
        eid_dict = self._input_dict(ppty_name, eid_dict, 'edge', True)
        # - Finally add the edge property to the graph:
        self._edge_property[ppty_name] = eid_dict
        return
    # add_edge_property.__doc__ = IPropertyGraph.add_edge_property.__doc__

    def extend_edge_property(self, ppty_name, eid_dict):
        """
        Extend an existing edge property 'ppty_name' with 'eid_dict', an eid dictionary.

        Parameters
        ----------
        ppty_name : str
            a string matching an exiting property
        eid_dict : dict
            a dictionary {eid: eid_values}

        Returns
        -------
        Nothing, edit object

        Notes
        ----
        * 'ppty_name' should exist
        * 'eid_dict' cannot be an empty dictionary
        """
        # - Assert the property does exist:
        if ppty_name not in self._edge_property:
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'edge'))

        # - Check the input 'eid_dict', and that all keys are in the graph:
        eid_dict = self._input_dict(ppty_name, eid_dict, 'edge', False)

        # - Update edge property 'ppty_name':
        ppty_dict = self._edge_property[ppty_name]
        ppty_dict = self._extend_ppty(ppty_name, ppty_dict, eid_dict, 'edge')
        self._vertex_property[ppty_name] = ppty_dict
        return

    def remove_edge_property(self, ppty_name):
        """
        Remove the edge property 'ppty_name' of the object.

        Parameters
        ----------
        ppty_name : str
            name of the edge property to remove

        Returns
        -------
        Nothing, edit object
        """
        # - Try to remove any existing associated graph property 'unit':
        try:
            del self._graph_property['units'][ppty_name]
        except KeyError:
            pass
        # - Try to remove the edge property, and if so, print about it:
        try:
            n = len(self.edge_property(ppty_name))
            del self._edge_property[ppty_name]
        except KeyError:
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'edge'))
        else:
            print "Removed edge property '{}' (n={})".format(ppty_name, n)
        return
    # remove_edge_property.__doc__ = IPropertyGraph.remove_edge_property.__doc__

    def update_edge_property(self, ppty_name, values, add_missing=True):
        """
        Update an existing edge property, with non-null values.
        Can add missing property (instead of failing) if 'add_missing' is True.

        Parameters
        ----------
        ppty_name: str
            the name of the property to add
        values: any non-null
            any type of values or object, except null values like empty list,
            set, dictionary or 'None'
        add_missing: bool, optional
            if True (default), add the property `ppty_name` if missing,
            else may fail if not defined

        Returns
        -------
        Nothing, edit object
        """
        # - Assert the property does exist:
        if ppty_name not in self._edge_property:
            if add_missing:
                self.add_edge_property(ppty_name, values)
            else:
                raise PropertyError(MISSING_PPTY.format(ppty_name, 'edge'))
        # - Raise an error when no 'values' are given:
        if not values:
            raise PropertyError(NULL_VALUES.format('edge', ppty_name))
        # - Update the edge property dictionary:
        self._edge_property[ppty_name].update(values)
        return

    def clear_edge_property(self, ppty_name):
        """
        Clear any values associated with the edge property 'ppty_name'.

        Parameters
        ----------
        ppty_name : str
            name of the edge property to clear

        Returns
        -------
        Nothing, edit object
        """
        # - Try to clear any existing associated graph property 'unit':
        try:
            self._graph_property['units'][ppty_name] = None
        except KeyError:
            pass
        # - Try to clear the edge property, and if so, print about it:
        try:
            n = len(self.edge_property(ppty_name))
            self._edge_property[ppty_name] = {}
        except KeyError:
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'edge'))
        else:
            print "Cleared edge property '{}' of n={} values!".format(ppty_name, n)
        return

    ############################################################################
    # GRAPH property creation and edition functions :
    ############################################################################

    def add_graph_property(self, ppty_name, values=None):
        """
        Add a new graph property, empty if values is None else add it to the object.

        Parameters
        ----------
        ppty_name : str
            the name of the property to add
        values : Any|None
            any type or object

        Returns
        -------
        Nothing, edit object
        """
        # - Assert the property does NOT exist:
        if ppty_name in self._graph_property:
            raise PropertyError(EXISTING_PPTY.format(ppty_name, 'graph'))
        # - Print when no 'values' are given:
        if values is None:
            print EMPTY_PPTY.format('graph', ppty_name)
        # - Add the graph property:
        self._graph_property[ppty_name] = values
        return

    def extend_graph_property(self, ppty_name, values):
        """
        Extend an existing graph property 'ppty_name' with 'values', the defined
        values should be an iterable type like 'list', 'set' or 'dict'.

        Parameters
        ----------
        ppty_name : str
            a string matching an exiting property
        values : Any
            any type or object

        Returns
        -------
        Nothing, edit object

        Notes
        ----
        * 'ppty_name' should exist
        * 'values' cannot be an empty dictionary
        """
        # - Assert the property does exist:
        if ppty_name not in self._graph_property:
            raise PropertyError(EXISTING_PPTY.format(ppty_name, 'graph'))
        # - Raise an error when no 'values' are given:
        if values is None:
            raise PropertyError("Trying to extend graph property '{}' with None value!".format(ppty_name))

        if isinstance(self.graph_property(ppty_name), list):
            if values:
                self._graph_property[ppty_name].extend(values)
            else:
                raise PropertyError(EMPTY_EXTEND.format('graph', ppty_name, 'list'))
        elif isinstance(self.graph_property(ppty_name), set):
            if values:
                self._graph_property[ppty_name].update(values)
            else:
                raise PropertyError(EMPTY_EXTEND.format('graph', ppty_name, 'set'))
        elif isinstance(self.graph_property(ppty_name), dict):
            try:
                assert values != {}
            except AssertionError:
                raise PropertyError(EMPTY_EXTEND.format('graph', ppty_name, 'dict'))
            else:
                ppty_dict = self.graph_property(ppty_name)
                ppty = {k: v for k, v in values.items() if k not in ppty_dict}
                self._graph_property[ppty_name].update(ppty)
        else:
            raise TypeError("Unable to extend graph property '{}' (type:{}) with this type of values: {}".format(ppty_name, type(self._graph_property[ppty_name]), type(values)))
        return

    def remove_graph_property(self, ppty_name):
        """
        Remove the graph property 'ppty_name' of the object.

        Parameters
        ----------
        ppty_name : str
            name of the graph property to remove

        Returns
        -------
        Nothing, edit object
        """
        # - Try to remove any existing associated graph property 'unit':
        try:
            del self._graph_property['units'][ppty_name]
        except KeyError:
            pass
        # - Try to remove the graph property, and if so, print about it:
        try:
            n = len(self.graph_property(ppty_name))
            del self._graph_property[ppty_name]
        except KeyError:
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'graph'))
        else:
            print "Removed graph property '{}' (n={})".format(ppty_name, n)
        return

    def update_graph_property(self, ppty_name, values, add_missing=True):
        """
        Update an existing graph property, with non-null values.
        Can add missing property (instead of failing) if 'add_missing' is True.

        Parameters
        ----------
        ppty_name: str
            the name of the property to add
        values: any non-null
            any type of values or object, except null values like empty list,
            set, dictionary or 'None'
        add_missing: bool, optional
            if True (default), add the property `ppty_name` if missing,
            else may fail if not defined

        Returns
        -------
        Nothing, edit object
        """
        # - Assert the property does exist:
        if ppty_name not in self._graph_property:
            if add_missing:
                self.add_graph_property(ppty_name, values)
            else:
                raise PropertyError(MISSING_PPTY.format(ppty_name, 'graph'))
        # - Raise an error when no 'values' are given:
        if not values:
            raise PropertyError(NULL_VALUES.format('graph', ppty_name))
        # - Update the graph property:
        ppty = self._graph_property[ppty_name]
        if isinstance(ppty, dict) or isinstance(ppty, set):
            self._graph_property[ppty_name].update(values)
        else:
            try:
                n_old = len(ppty)
                assert not isinstance(ppty, str)
            except:
                n_old = 1
            try:
                n_new = len(values)
                assert not isinstance(values, str)
            except:
                n_new = 1
            s = "Replaced the graph property '{}' ".format(ppty_name)
            s += "(n={}) ".format(n_old) if n_old > 1 else "('{}') ".format(ppty)
            s += "by {} values.".format(n_new) if n_new > 1 else "by '{}'.".format(ppty)
            print s
            self._graph_property[ppty_name] = values
        return

    def clear_graph_property(self, ppty_name):
        """
        Clear any values associated with the graph property 'ppty_name'.

        Parameters
        ----------
        ppty_name : str
            name of the graph property to clear

        Returns
        -------
        Nothing, edit object
        """
        # - Try to clear any existing associated graph property 'unit':
        try:
            self._graph_property['units'][ppty_name] = None
        except KeyError:
            pass
        # - Try to clear the graph property, and if so, print about it:
        try:
            n = len(self.graph_property(ppty_name))
            self._graph_property[ppty_name] = {}
        except KeyError:
            raise PropertyError(MISSING_PPTY.format(ppty_name, 'graph'))
        else:
            print "Cleared graph property '{}' of n={} values!".format(ppty_name, n)
        return

    ############################################################################
    # GRAPH extension and relabelling functions :
    ############################################################################

    @staticmethod
    def _translate_property(dict_values, trans_vid, trans_eid,
                            key_translation, value_translation):
        """
        Translate edge and vertex properties (dictionary) using 'trans_vid' &
        'trans_eid' dict, and "type"
        Typical use is with 'PropertyGraph.extend_property_graph'

        Parameters
        ----------
        dict_values : dict
            dictionary of values to translate
        trans_vid : dict
            vertex "translation" dictionary, linking current object vertex ids
            to 'graph' vertex ids, ie. {self_vid: 'graph'_vid}
        trans_eid : dict
            edge "translation" dictionary, linking current object edge ids to
            'graph' edge ids, ie. {self_eid: 'graph'_eid}
        key_translation : int
            indicate the "type" of the keys to translate: VertexIdType (0),
            EdgeIdType (1) or ValueType (2)
        value_translation : int
            indicate the "type" of the values to translate: VertexIdType (0),
            EdgeIdType (1) or ValueType (2)

        Returns
        -------
        trans_dict : dict
            translated dictionary
        """
        # - Define identity translation function, used with 'ValueType'
        id_value = lambda value: value

        # - Define vertex-ids translation function, used with 'VertexIdType'
        def translate_vid(vid):
            if isinstance(vid, list):
                return [trans_vid[i] for i in vid]
            elif isinstance(vid, tuple):
                return tuple([trans_vid[i] for i in vid])
            elif isinstance(vid, int):
                return trans_vid[vid]
            else:
                t = type(vid)
                err_message = "Got wrong 'vid' type: {}... ".format(t)
                err_message += "Choose among: int, list or tuple."
                raise TypeError(err_message)

        # - Define edge-ids translation function, used with 'EdgeIdType'
        def translate_eid(eid):
            if isinstance(eid, list):
                return [trans_eid[i] for i in eid]
            elif isinstance(eid, tuple):
                return tuple([trans_eid[i] for i in eid])
            elif isinstance(eid, int):
                return trans_eid[eid]
            else:
                t = type(eid)
                err_message = "Got wrong 'eid' type: {}... ".format(t)
                err_message += "Choose among: int, list or tuple."
                raise TypeError(err_message)

        # - Define a dict of translation function to use by type of ids:
        translator = {ValueType: id_value, VertexIdType: translate_vid,
                      EdgeIdType: translate_eid}
        # - Get the type of translation for keys & values:
        k_trans = translator[key_translation]
        v_trans = translator[value_translation]
        return {k_trans(k): v_trans(v) for k, v in dict_values.iteritems()}

    def _relabel_and_add_vertex_properties(self, graph, trans_vid, verbose):
        """
        Renumber and add vertex properties contained in 'graph'.
        The "translation" of graph.vertices() has already been performed and is
        given, respectively, by 'trans_vid'

        Parameters
        ----------
        graph : PropertyGraph
            a graph with edge and vertex property to renumber and add to the
            object
        trans_vid : dict
            vertex "translation" dictionary, linking current object vertex ids
            to 'graph' vertex ids, ie. {self_vid: 'graph'_vid}

        Returns
        -------
        Nothing, add 'graph' vertex properties to current object.
        """
        # update properties on vertex
        for ppty_name in graph.vertex_property_names():
            if verbose:
                print "Relabelling vertex property '{}'...".format(ppty_name)
            if ppty_name not in self._vertex_property:
                self.add_vertex_property(ppty_name, {})
            # Define translation type for 'ppty_name' values
            value_translator = graph.get_graph_property_value_type(ppty_name)
            # translate accordingly, using 'VertexIdType' for keys:
            new_ppty = self._translate_property(
                graph.vertex_property(ppty_name),
                trans_vid, {}, VertexIdType,
                value_translator)
            # update object graph property 'ppty_name':
            self._vertex_property[ppty_name].update(new_ppty)
            if verbose:
                n_given =  len(graph.vertex_property(ppty_name))
                n_obtained = len(new_ppty)
                pc = round(n_obtained/float(n_given) * 100, 1)
                print "{}% ({}/{}) translated!".format(pc, n_given, n_obtained)
        return

    def _relabel_and_add_edge_properties(self, graph, trans_eid, verbose):
        """
        Renumber and add edge properties contained in 'graph'.
        The "translation" of graph.edges() has already been performed and is
        given by 'trans_eid'

        Parameters
        ----------
        graph : PropertyGraph
            a graph with edge and vertex property to renumber and add to the
            object
        trans_eid : dict
            edge "translation" dictionary, linking current object edge ids to
            'graph' edge ids, ie. {self_eid: 'graph'_eid}

        Returns
        -------
        Nothing, add 'graph' edge properties to current object.
        """
        # update properties on edges
        for ppty_name in graph.edge_property_names():
            if verbose:
                print "Relabelling edge property '{}'...".format(ppty_name)
            if ppty_name not in self._edge_property:
                self.add_edge_property(ppty_name, {})
            # Define translation type for 'ppty_name' values
            value_translator = graph.get_graph_property_value_type(ppty_name)
            # translate accordingly, using 'EdgeIdType' for keys:
            new_ppty = self._translate_property(graph.edge_property(ppty_name),
                                                {}, trans_eid, EdgeIdType,
                                                value_translator)
            # update object graph property 'ppty_name':
            self._edge_property[ppty_name].update(new_ppty)
            if verbose:
                n_given =  len(graph.edge_property(ppty_name))
                n_obtained = len(new_ppty)
                pc = round(n_obtained/float(n_given) * 100, 1)
                print "{}% ({}/{}) translated!".format(pc, n_given, n_obtained)
        return

    def _relabel_and_add_graph_properties(self, graph, trans_vid, trans_eid, verbose):
        """
        Translate a graph property according to meta info

        Parameters
        ----------
        graph : PropertyGraph
            a graph with edge and vertex property to renumber and add to the
            object
        trans_vid : dict
            vertex "translation" dictionary, linking current object vertex ids
            to 'graph' vertex ids, ie. {self_vid: 'graph'_vid}
        trans_eid : dict
            edge "translation" dictionary, linking current object edge ids to
            'graph' edge ids, ie. {self_eid: 'graph'_eid}

        Returns
        -------
        Nothing, add 'graph' graph properties to current object.
        """
        for ppty_name in graph.graph_property_names():
            old_ppty = graph.graph_property(ppty_name)
            # - Case where 'old_ppty' is a dictionary:
            if isinstance(old_ppty, dict):
                if ppty_name not in self._graph_property:
                    self.add_graph_property(ppty_name, {})
                # - Define translation type for 'ppty_name' keys
                key_translator = graph.get_graph_property_key_type(ppty_name)
                # - Define translation type for 'ppty_name' values
                value_translator = graph.get_graph_property_value_type(ppty_name)
                new_ppty = self._translate_property(old_ppty, trans_vid,
                                                    trans_eid,
                                                    key_translator,
                                                    value_translator)
                # update object graph property 'ppty_name':
                self._graph_property[ppty_name].update(new_ppty)
            # - Case where 'old_ppty' is a list:
            elif isinstance(old_ppty, list):
                if ppty_name not in self._graph_property:
                    self.add_graph_property(ppty_name, [])
                # - Define translation type for 'ppty_name' values
                value_translator = graph.get_graph_property_value_type(ppty_name)
                # create a temprorary dict with "fake" keys:
                old_ppty = {n: v for n, v in enumerate(old_ppty)}
                new_ppty = self._translate_property(old_ppty, trans_vid,
                                                    trans_eid,
                                                    ValueType,
                                                    value_translator)
                # update object graph property 'ppty_name':
                self._graph_property[ppty_name].append(new_ppty.values())
            else:
                if ppty_name not in self._graph_property:
                    self.add_graph_property(ppty_name, old_ppty)
                elif isinstance(old_ppty, int) or isinstance(old_ppty, bool):
                    try:
                        assert old_ppty == self.graph_property(ppty_name)
                    except AssertionError:
                        err = "Found different values for graph property"
                        err += " '{}':".format(ppty_name)
                        err += " {} (old value)".format(self.graph_property(ppty_name))
                        err += " {} (new value)".format(old_ppty)
                        raise ValueError(err)

        return

    def extend_property_graph(self, graph, verbose=True):
        # def extend(self, graph):
        """
        Extend the current object ('self') with another 'graph' of type Graph or
        PropertyGraph.

        Parameters
        ----------
        graph : Graph|PropertyGraph
            the object used to extent the current one.

        Returns
        -------
        two dictionaries ({vid: vid_g}, {eid: eid_g}) specifying correspondences
        between "'graph' id" (vid_g, eid_g) and "'self' id"

        Notes
        -----
        * differ from Graph.extend() since relabelling of 'graph' vids and eids
        should be applied to any vertex or edge properties from this graph before .
        """
        try:
            assert isinstance(graph, Graph) or isinstance(graph, PropertyGraph)
        except AssertionError:
            err = "Input 'graph' should be a 'Graph' or a 'PropertyGraph' object"
            err += ", got '{}!".format(type(graph))
            raise TypeError(err)
        # Renumber and add the vertex and edge ids of 'graph':
        trans_vid, trans_eid = Graph.extend(self, graph)
        # relabel the vertex, edge and graph properties:
        self._relabel_and_add_vertex_properties(graph, trans_vid, verbose)
        self._relabel_and_add_edge_properties(graph, trans_eid, verbose)
        self._relabel_and_add_graph_properties(graph, trans_vid, trans_eid, verbose)
        return trans_vid, trans_eid

    # extend_property_graph.__doc__ = Graph.extend.__doc__
    # extend.__doc__ = Graph.extend.__doc__

    def _set_graph_property_value_to_type(self, ppty_name, value_type):
        """
        Associate self._graph_property['ppty_name'] value type to one of these
        types: VertexIdtype, EdgeIdType.

        Parameters
        ----------
        ppty_name : str
           graph property name for which to set the type of values
        key_type : int
           define values 'VertexIdtype' (0) or 'EdgeIdType' (1)

        Returns
        -------
        Nothing, add meta info to the object.
        """
        if value_type == VertexIdType:
            value_type = self.meta_vid_value_type
        elif value_type == EdgeIdType:
            value_type = self.meta_eid_value_type
        else:
            raise ValueError(
                "Unknown type '{}' for 'value_type'".format(type(value_type)))
        # Add the graph property 'ppty_name' to the corresponding 'value_type':
        if value_type not in self._graph_property:
            self.add_graph_property(value_type, [ppty_name])
        else:
            self.extend_graph_property(value_type, ppty_name)
        return

    def set_graph_property_value_to_vid_type(self, ppty_name):
        """
        Associate self._graph_property['ppty_name'] value type to VertexIdtype.

        Parameters
        ----------
        ppty_name : str
            graph property name for which to declare the values as vertex-ids

        Returns
        -------
        Nothing, add meta info to the object.
        """
        return self._set_graph_property_value_to_type(ppty_name, VertexIdType)

    def set_graph_property_value_to_eid_type(self, ppty_name):
        """
        Associate self._graph_property['ppty_name'] value type to EdgeIdType.


        Parameters
        ----------
        ppty_name : str
            graph property name for which to declare the values as edge-ids

        Returns
        -------
        Nothing, add meta info to the object.
        """
        return self._set_graph_property_value_to_type(ppty_name, EdgeIdType)

    def _set_graph_property_key_to_type(self, ppty_name, key_type):
        """
        Associate self._graph_property['ppty_name'] key type to one of these 
        types: VertexIdtype, EdgeIdType.
        
        Parameters
        ----------
        ppty_name : str
            graph property name for which to set the type of values
        key_type : int
            define values 'VertexIdtype' (0) or 'EdgeIdType' (1)

        Returns
        -------
        Nothing, add meta info to the object.
        """
        if key_type == VertexIdType:
            key_type = self.meta_vid_key_type
        elif key_type == EdgeIdType:
            key_type = self.meta_eid_key_type
        else:
            raise ValueError(
                "Unknown type '{}' for 'key_type'".format(type(key_type)))
        # Add the graph property 'ppty_name' to the corresponding 'key_type':
        if key_type not in self._graph_property:
            self.add_graph_property(key_type, [ppty_name])
        else:
            self.extend_graph_property(key_type, ppty_name)
        return

    def set_graph_property_key_to_vid_type(self, ppty_name):
        """
        Associate self._graph_property['ppty_name'] key type to one of these 
        types: VertexIdtype, EdgeIdType.

        Parameters
        ----------
        ppty_name : str
            graph property name for which to declare the keys as vertex-ids

        Returns
        -------
        Nothing, add meta info to the object.
        """
        return self._set_graph_property_key_to_type(ppty_name, VertexIdType)

    def set_graph_property_key_to_eid_type(self, ppty_name):
        """
        Associate self._graph_property['ppty_name'] key type to one of these
        types: VertexIdtype, EdgeIdType.

        Parameters
        ----------
        ppty_name : str
            graph property name for which to declare the keys as edge-ids

        Returns
        -------
        Nothing, add meta info to the object.
        """
        return self._set_graph_property_key_to_type(ppty_name, EdgeIdType)

    def get_graph_property_value_type(self, ppty_name):
        """
        Return meta info on property 'ppty_name' value type.

        Parameters
        ----------
        ppty_name : str
            graph property name to use

        Returns
        -------
        VertexIdType (0), EdgeIdType (1) or ValueType (2)
        """
        # - Test if is a 'VertexIdType' value:
        try:
            prop = self.graph_property(self.meta_vid_value_type)
        except PropertyError:
            pass
        else:
            if ppty_name in prop:
                return VertexIdType
        # - Test if is an 'EdgeIdType' value:
        try:
            prop = self.graph_property(self.meta_eid_value_type)
        except PropertyError:
            pass
        else:
            if ppty_name in prop:
                return EdgeIdType
        # - If neither an EdgeIdType nor a VertexIdType, return ValueType
        return ValueType

    def get_graph_property_key_type(self, ppty_name):
        """
        Return meta info on graph property 'ppty_name' key type.
        
        Parameters
        ----------
        ppty_name : str
            graph property name to use

        Returns
        -------
        VertexIdType (0), EdgeIdType (1) or ValueType (2)
        """
        # - Test if is a 'VertexIdType' value:
        try:
            prop = self.graph_property(self.meta_vid_key_type)
        except PropertyError:
            pass
        else:
            if ppty_name in prop:
                return VertexIdType
        # - Test if is an 'EdgeIdType' value:
        try:
            prop = self.graph_property(self.meta_eid_key_type)
        except PropertyError:
            pass
        else:
            if ppty_name in prop:
                return EdgeIdType
        # - If neither an EdgeIdType nor a VertexIdType, return ValueType
        return ValueType

    ############################################################################
    # EXPORT/IMPORT with NetworkX library
    ############################################################################

    def to_networkx(self):
        """
        Return a NetworkX Graph object from current object.

        Returns
        -------
        a NetworkX graph object (nx.Graph).
        """
        s = self.source
        t = self.target
        # - Initialize a NetworkX graph:
        graph = nx.Graph()
        # - Add vertices and edges:
        graph.add_nodes_from(self.vertices())
        graph.add_edges_from(((s(eid), t(eid)) for eid in self.edges()))

        # - Add graph, vertex and edge properties
        # -- Copy 'graph_properties', if any:
        for k, v in self.graph_properties().iteritems():
            graph.graph[k] = v

        # -- Copy 'vertex_properties', if any:
        vp = self._vertex_property
        for prop in vp:
            for vid, value in vp[prop].iteritems():
                graph.node[vid][prop] = value

        # -- Copy 'edge_properties', if any:
        ep = self._edge_property
        # --- Starts with the 'eid' values (not defined in NetworkX)
        for eid in self.edges():
            graph.edges[(s(eid), t(eid))]['eid'] = eid
        # --- Get the defined 'edge_properties':
        for prop in ep:
            for eid, value in ep[prop].iteritems():
                graph.edges[(s(eid), t(eid))][prop] = value

        return graph

    def from_networkx(self, nx_graph):
        """
        Return a PropertyGraph from a NetworkX Directed graph 'graph'.

        Parameters
        ----------
        nx_graph : networkx.graph
            a networkx directed graph

        Returns
        -------
        a PropertyGraph.
        """
        self.clear()
        g = self

        if not nx_graph.is_directed():
            nx_graph = nx_graph.to_directed()

        vp = self._vertex_property
        for vid in nx_graph.nodes_iter():
            g.add_vertex(vid)
            d = nx_graph.node[vid]
            for k, v in d.iteritems():
                vp.setdefault(k, {})[vid] = v

        ep = self._edge_property
        for source, target in nx_graph.edges_iter():
            d = nx_graph[source][target]
            eid = d.get('eid')
            eid = g.add_edge(source, target, eid)
            for k, v in d.iteritems():
                if k != 'eid':
                    ep.setdefault(k, {})[eid] = v

        gp = self._graph_property
        gp.update(nx_graph.graph)

        return g

    ############################################################################
    # EXPORT graph as a CSV
    ############################################################################
    # TODO: rewrite using Pandas library?
    def to_csv(self, graph_name, ppty2export=None, out_fname=None,
               datetime2fname=True):
        """
        Export vertex properties 'ppty2export' to a csv named 'out_fname'.
        Cells are given by row and properties by columns.
        For length-D vectors properties, like barycenters, each D value will be
        outputed in a separate column.
        By default export all properties.

        Parameters
        ----------
        graph_name : str
            name or id of the 'PropertyGraph'.
        ppty2export : list
            list of vid associated properties to export. None, by default,
            export them all.
        out_fname : str
            the name of the csv file!
        datetime2fname : bool
            if True (default) add the date to the csv filename.

        Note
        ----
        Based on the original work of V.Mirabet.
        """
        # Init the CSV header with ['graph_name'; 'vid';]:
        csv_head = "graph_name" + ";" + "vid" + ";"
        # Get the possible ppty to export crossing required 'ppty_sublist' &
        # available ones (g.vertex_property_names):
        if ppty2export is not None:
            ppty2export = list(
                set(self.vertex_property_names()) & set(ppty2export))
        else:
            ppty2export = sorted(self.vertex_property_names())

        # Add these ppty names to export to the CSV header:
        for ppty in ppty2export:
            if ppty == "barycenter":
                csv_head += "bary_x;bary_y;bary_z;"
            elif ppty == "barycenter_voxel":
                csv_head += "bary_vox_x;bary_vox_y;bary_vox_z;"
            elif ppty == "inertia_values_3D":
                csv_head += "inertia_val_1;inertia_val_2;inertia_val_3;"
            else:
                csv_head += ppty + ";"
        csv_head += "\n"  # CSV header endline!

        # Now loop the vids to get their feature values:
        # thus creating a line for each cell with their associated features.
        csv = csv_head
        for vid in sorted(self.vertices()):
            # Add 'flower_id':
            csv += graph_name + ";"
            # Add 'vid':
            csv += str(vid) + ";"
            # Add cell feature values:
            for ppty in ppty2export:
                # In case of length-3 vectors:
                if ppty in ["barycenter", "barycenter_voxel",
                            "inertia_values_3D"]:
                    if vid in self.vertex_property(ppty):
                        for i in range(3):
                            csv += str(self.vertex_property(ppty)[vid][i]) + ";"
                    else:
                        csv += ";;;"
                else:
                    if vid in self.vertex_property(ppty):
                        pc = self.vertex_property(ppty)[vid]
                        if type(pc) == str:
                            csv += pc + ";"
                        else:
                            try:
                                csv += str(float(pc)) + ";"
                            except Exception as e:
                                # print pc, e
                                csv += ";"
                    else:
                        csv += ";"
            # Remove the last ';' and put an endline '\n' instead:
            csv = csv[:-1] + "\n"

        if out_fname is None:
            out_fname = str(graph_name)
            datetime2fname = True
        # Add date time if requested by 'date=True' argument:
        if datetime2fname:
            import time
            m = time.localtime()
            out_fname = out_fname + "_" + str(m.tm_year) + "_" + str(
                m.tm_mon) + "_" + str(m.tm_mday) + "_" + str(
                m.tm_hour) + "_" + str(m.tm_min)

        # Finally write the WHOLE 's' string
        f = open(out_fname + ".csv", "w")
        f.write(csv)
        f.close()
        print "Done witting '{}' file!".format(out_fname + ".csv")
        return

    ############################################################################
    # Misc: Should probably be elsewhere !
    ############################################################################

    @staticmethod
    def _to_set(s):
        """
        Hidden method to transforms 's' into a list.

        Parameters
        ----------
        s : int|list|set
            a variable to transform into a set.

        Returns
        -------
        a set from 's'
        """
        if not isinstance(s, set):
            if isinstance(s, list):
                s = set(s)
            elif isinstance(s, int):
                s = {s}
            elif isinstance(s, str) and len(s)==1:
                s = {s}
            else:
                raise TypeError(
                    "Unable to transform input type {} into a set!".format(
                        type(s)))
        return s

    def get_vertex_ppty_type(self, vtx_ppty):
        """
        Return the type (scalar|vector|tensor) of a given vertex property.

        Parameters
        ----------
        vtx_ppty : str
            name of an existing vertex property
        """
        # TODO: move to TissueGraph? class.

        # TODO: Not really working except for scalars...
        # TODO: ...should be more strict when defining variables types!!!

        # Get every 'vtx_ppty' values:
        values_list = [v for v in self.vertex_property(vtx_ppty).values() if
                       v is not None]
        # Now get their types and simplify this list to its 'unique values' using 'set()'.
        types_set = set([type(v) for v in values_list])

        if len(types_set) != 1:
            raise warnings.warn(
                "More than ONE type detected for vertex property '{}'! Please check it!".format(
                    vtx_ppty))

        first_val = values_list[0]
        print first_val
        first_type = type(first_val)

        try:
            sh = first_val.shape
        except TypeError:
            pass
        else:
            if sh != (1, 1) and sh != ():
                return "tensor", sh
            else:
                return "scalar"

        try:
            le = first_val.len
        except TypeError:
            return "scalar"
        else:
            if le != 1:
                return "vector", le
            else:
                return "scalar"
