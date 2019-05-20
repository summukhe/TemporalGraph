import logging
import numpy as np
from .graph import *
from copy import deepcopy
from temporal_graph.spatial_ds import *

__version__ = "1.0"
__all__ = ['PathFilter', 'GeometricPathFilter', 'AllPathIterator']


class PathFilter:
    def __init__(self):
        pass

    def __call__(self, g, n, n_1, n_2=None):
        return True


class GeometricPathFilter(PathFilter):
    def __init__(self,
                 forward=True,
                 weight_cutoff=None,
                 distance_cutoff=None):
        self.__logger = logging.getLogger(name='temporal_graph.GeometricPathFilter')
        PathFilter.__init__(self)
        self.__forward_only = forward
        if weight_cutoff is not None:
            self.__logger.debug('Setting weight cutoff (%f)!!' % weight_cutoff)
            self.__weight_cutoff = weight_cutoff
        else:
            self.__weight_cutoff = None

        if distance_cutoff is not None:
            self.__logger.debug('Setting distance cutoff!!')
            if isinstance(distance_cutoff, DistanceCutoff):
                self._distance_cutoff = distance_cutoff
            elif isinstance(distance_cutoff, np.float) or \
                    isinstance(distance_cutoff, np.double) or \
                    isinstance(distance_cutoff, np.int):
                self._distance_cutoff = DistanceCutoff(def_cutoff=distance_cutoff)
        else:
            self._distance_cutoff = None

    @property
    def is_forward_only(self):
        return self.__forward_only

    @property
    def has_weight_cutoff(self):
        return self.__weight_cutoff is not None

    @property
    def weight_cutoff(self):
        return self.__weight_cutoff

    @property
    def has_distance_cutoff(self):
        return self._distance_cutoff is not None

    @property
    def distance_cutoff(self):
        return self._distance_cutoff

    def __call__(self, g, n, n_1, n_2=None):
        if not isinstance(g, GeometricGraph3d):
            self.__logger.debug("Not a geometric graph!")
            return False
        if not g.is_vertex(n) or not g.is_vertex(n_1) or not g.has_edge(n_1, n):
            self.__logger.debug("improper edge defined (%s : %s)" % (n_1, n))
            return False
        if (n_2 is not None) and (not g.is_vertex(n_2) or not g.has_edge(n_2, n_1)):
            self.__logger.debug("improper edge defined (%s : %s)" % (n_2, n_1))
            return False
        if self.has_weight_cutoff and g.weight(n_1, n) < self.weight_cutoff:
            self.__logger.debug("failed at weight cutoff (%f < %f)" % (g.weight(n_1, n), self.weight_cutoff))
            return False
        if self.has_distance_cutoff:
            d = distance(g.attribute(n), g.attribute(n_1))
            if d > self._distance_cutoff():
                self.__logger.debug("failed at distance cutoff (%f)" % d)
                return False
        if self.is_forward_only and n_2 is not None:
            v1 = connecting_vector(g.attribute(n_1), g.attribute(n)).unit_vector
            v2 = connecting_vector(g.attribute(n_2), g.attribute(n_1)).unit_vector
            if dotp(v1, v2) > 0:
                return False
        return True


class AllPathIterator:
    def __init__(self,
                 g,
                 visitor):
        assert isinstance(g, WeightedGraph)
        assert isinstance(visitor, PathFilter)
        self.__g = deepcopy(g)
        self.__visitor = visitor
        self.__nodes = {v: False for v in g.vertices}
        self.__stack = []
        self.__paths = []
        self.__stop_vertex = {}

    def is_vertex(self, v):
        return self.__g.is_vertex(v)

    def all_path(self,
                 start_vertex,
                 stop_vertex,
                 min_path_length=2,
                 max_path_length=None):
        assert self.__g.is_vertex(start_vertex)
        if isinstance(stop_vertex, list):
            self.__stop_vertex = set(stop_vertex)
        else:
            self.__stop_vertex = {stop_vertex}
        assert start_vertex not in self.__stop_vertex
        for s in self.__stop_vertex:
            assert self.__g.is_vertex(s)
        if max_path_length is None:
            max_path_length = self.__g.order // 2
        self.__all_paths(start_vertex,
                         min_path_length=min_path_length,
                         max_path_length=max_path_length)
        paths = deepcopy(self.__paths)
        self.__paths = []
        return paths

    def __all_paths(self,
                    start_vertex,
                    min_path_length=2,
                    max_path_length=10):
        assert self.__g.is_vertex(start_vertex)
        if not self.__nodes[start_vertex]:
            append = False
            if len(self.__stack) > 0:
                v = start_vertex
                v_1 = self.__stack[-1]
                v_2 = self.__stack[-2] if len(self.__stack) > 1 else None
                if self.__visitor(self.__g, v, v_1, v_2):
                    append = True
            else:
                append = True
            if append:
                self.__stack.append(start_vertex)
                self.__nodes[start_vertex] = True
                if start_vertex in self.__stop_vertex:
                    if len(self.__stack) >= min_path_length:
                        path_dict = {'path': deepcopy(self.__stack),
                                     'weights': []}
                        for i in range(1, len(self.__stack)):
                            path_dict['weights'].append(self.__g.weight(self.__stack[i-1], self.__stack[i]))
                        self.__paths.append(path_dict)
                elif len(self.__stack) < max_path_length:
                    for v in self.__g.out_neighbors(start_vertex):
                        self.__all_paths(v, min_path_length=min_path_length)
                self.__stack.pop(-1)
                self.__nodes[start_vertex] = False


