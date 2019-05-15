"""
Basic utility for graph representation
"""

import numpy as np
import pandas as pd

__version__ = "1.0"
__all__ = ['SimpleGraph', 'WeightedGraph']


class SimpleGraph(object):
    def __init__(self, directed=False):
        self.__directed = directed
        self.__vertex = set()
        self.__in_edges = dict()
        self.__out_edges = dict()
        self.__edge_count = 0
        self.__edge_id = 0

    def add_vertex(self, vertex):
        if vertex not in self.__vertex:
            self.__vertex.add(vertex)
            self.__in_edges[vertex] = dict()
            self.__out_edges[vertex] = dict()

    def has_edge(self, v1, v2):
        return v1 in self.__vertex and \
               v2 in self.__vertex and \
               v2 in self.__out_edges[v1]

    def edge_id(self, v1, v2):
        if self.has_edge(v1, v2):
            return self.__out_edges[v1][v2]
        return None

    def add_edge(self, src, dst):
        self.add_vertex(src)
        self.add_vertex(dst)
        if not self.has_edge(src, dst):
            self.__edge_id += 1
            self.__edge_count += 1
            self.__out_edges[src][dst] = self.__edge_id
            self.__in_edges[dst][src] = self.__edge_id
            if self.__directed is False:
                self.__out_edges[dst][src] = self.__edge_id
                self.__in_edges[src][dst] = self.__edge_id
            return True
        return False

    def del_edge(self, src, dst):
        if self.has_edge(src, dst):
            self.__edge_count -= 1
            del self.__out_edges[src][dst]
            del self.__in_edges[dst][src]
            if self.__directed is False:
                del self.__out_edges[dst][src]
                del self.__in_edges[src][dst]
            return True
        return False

    def del_vertex(self, v):
        if v in self.__vertex:
            for u in self.__out_edges[v].keys():
                del self.__in_edges[u][v]
                self.__edge_count -= 1

            for u in self.__in_edges[v].keys():
                del self.__out_edges[u][v]
                if self.__directed is True:
                    self.__edge_count -= 1
            del self.__out_edges[v]
            del self.__in_edges[v]
            self.__vertex.remove(v)
            return True
        return False

    @property
    def order(self):
        return len(self.__vertex)

    @property
    def size(self):
        return self.__edge_count

    @property
    def edges(self):
        E = list()
        for v1 in self.__out_edges.keys():
            for v2 in self.__out_edges[v1].keys():
                if (self.__directed is False and v1 < v2) or (self.__directed is True):
                    E.append((v1, v2, self.__out_edges[v1][v2]))
        return [(u, v) for u, v, _ in sorted(E, key=lambda x: x[2])]

    @property
    def vertices(self):
        return sorted(list(self.__vertex))

    @property
    def adjacency(self):
        vs = self.vertices
        adj = pd.DataFrame(0, index=vs, columns=vs, dtype=np.int)
        for v1 in vs:
            for v2 in self.__out_edges[v1].keys():
                adj.loc[v1, v2] = 1
        return adj

    def in_degree(self, v):
        if v in self.__vertex:
            return len(self.__in_edges[v])
        return 0

    def out_degree(self, v):
        if v in self.__vertex:
            return len(self.__out_edges[v])
        return 0

    def neighbors(self, v):
        nbr = []
        if v in self.__vertex:
            for u in self.__out_edges[v].keys():
                nbr.append(u)
            if self.__directed is True:
                for u in self.__in_edges[v].keys():
                    nbr.append(u)
        return sorted(list(set(nbr)))

    def out_neigbors(self, v):
        nbr = []
        if v in self.__vertex:
            for u in self.__out_edges[v].keys():
                nbr.append(u)
        return sorted(nbr)

    def in_neighbors(self, v):
        nbr = []
        if v in self.__vertex:
            for u in self.__in_edges[v].keys():
                nbr.append(u)
        return sorted(nbr)

    def is_vertex(self, v):
        return v in self.__vertex

    @property
    def is_directed(self):
        return self.__directed

    def components(self):
        def find_root(node, root):
            while node != root[node][0]:
                node = root[node][0]
            return node, root[node][1]

        root = {}
        for node in self.__out_edges.keys():
            root[node] = (node, 0)
        for vi in self.__out_edges.keys():
            for vj in self.__out_edges[vi]:
                (root_vi, depth_i) = find_root(vi, root)
                (root_vj, depth_j) = find_root(vj, root)
                if root_vi != root_vj:
                    mMin = root_vi
                    mMax = root_vj
                    if depth_i > depth_j:
                        mMin = root_vj
                        mMax = root_vi
                    root[mMax] = (mMax, max(root[mMin][1] + 1, root[mMax][1]))
                    root[mMin] = (root[mMax][0], -1)
        ret = {}
        for vi in self.__out_edges.keys():
            if root[vi][0] == vi:
                ret[vi] = []
        for vi in self.__out_edges.keys():
            ret[find_root(vi, root)[0]].append(vi)
        return list(ret.values())

    def is_connected(self, vi, vj):
        if vi in self.__vertex and vj in self.__vertex:
            c = self.components()
            for l in c:
                if vi in l and vj in l:
                    return True
        return False


class WeightedGraph:
    def __init__(self, directed=False, def_weight=1.0):
        self.__G = SimpleGraph(directed=directed)
        self.__def_wt = def_weight
        self.__vattrib_lookup = dict()
        self.__eweight_lookup = dict()

    def add_vertex(self, v, attribute=None):
        self.__G.add_vertex(v)
        if (v not in self.__vattrib_lookup) or (attribute is not None):
            self.__vattrib_lookup[v] = attribute

    def add_edge(self, src, dst, weight=None):
        self.add_vertex(src)
        self.add_vertex(dst)
        self.__G.add_edge(src, dst)
        e_id = self.__G.edge_id(src, dst)
        if (e_id not in self.__eweight_lookup) or (weight is not None):
            e_wt = self.__def_wt if weight is None else weight
            self.__eweight_lookup[e_id] = e_wt

    def del_edge(self, src, dst):
        e_id = self.__G.edge_id(src, dst)
        if e_id is not None:
            del self.__eweight_lookup[e_id]
        self.__G.del_edge(src, dst)

    def del_vertex(self, v):
        if self.__G.is_vertex(v):
            v_out = self.__G.out_neigbors(v)
            v_in = self.__G.in_neighbors(v)
            for u in v_out:
                e_id = self.__G.edge_id(v, u)
                del self.__eweight_lookup[e_id]
            for u in v_in:
                e_id = self.__G.edge_id(u, v)
                if e_id is not None and e_id in self.__eweight_lookup:
                    del self.__eweight_lookup[e_id]
            del self.__vattrib_lookup[v]
            return self.__G.del_vertex(v)
        return False

    @property
    def size(self):
        return self.__G.size

    @property
    def order(self):
        return self.__G.order

    @property
    def is_directed(self):
        return self.__G.is_directed

    def weight(self, src, dst):
        e_id = self.__G.edge_id(src, dst)
        if e_id is None:
            return None
        return self.__eweight_lookup[e_id]

    def attribute(self, v):
        if v in self.__vattrib_lookup:
            return self.__vattrib_lookup[v]
        return None

    @property
    def weight_matrix(self):
        vs = self.vertices
        wts = pd.DataFrame(0, index=vs, columns=vs)
        for vi in vs:
            for vj, eid in self.__out_edges[vi].items():
                wts.loc[vi, vj] = self.__eweight_lookup[eid]
        return wts

    @property
    def edges(self):
        return self.__G.edges

    @property
    def vertices(self):
        return self.__G.vertices

    def in_degree(self, v):
        return self.__G.in_degree(v)

    def in_neighbors(self, v):
        return self.__G.in_neighbors(v)

    def out_degree(self, v):
        return self.__G.out_degree(v)

    def out_neighbors(self, v):
        return self.__G.out_neigbors(v)

    def is_vertex(self, v):
        return self.__G.is_vertex(v)

    def has_edge(self, src, tgt):
        return self.__G.has_edge(src, tgt)

    def components(self):
        return self.__G.components()

    def is_connected(self, vi, vj):
        return self.__G.is_connected(vi, vj)

    @property
    def adjacency(self):
        return self.__G.adjacency


