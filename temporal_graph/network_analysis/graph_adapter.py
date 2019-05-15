import igraph
import logging
from .graph import *

__all__ = ['to_igraph', 'from_igraph']


def to_igraph(g):
    logger = logging.getLogger('network_analysis.to_igraph')
    if isinstance(g, SimpleGraph):
        logger.debug('Found SimpleGraph instance!!')
        return simple_graph_to_igraph(g)
    elif isinstance(g, WeightedGraph):
        logger.debug('Found WeightedGraph instance!!')
        return weighted_graph_to_igraph(g)
    logger.error('Does not know how to handle the graph instance')
    raise Exception('Invalid Graph Object!!')


def simple_graph_to_igraph(g):
    assert isinstance(g, SimpleGraph)
    vertices = g.vertices
    edges = g.edges
    ig = igraph.Graph(directed=g.is_directed)
    for v in vertices:
        ig.add_vertex(v)
    for u, v in edges:
        ig.add_edge(u, v)
    return ig


def weighted_graph_to_igraph(g):
    assert isinstance(g, WeightedGraph)
    vertices = g.vertices
    edges = g.edges
    wts = [g.weight(u,v) for u, v in edges]
    attr = [g.attribute(v) for v in vertices]
    assert len(edges) == len(wts)

    ig = igraph.Graph(directed=g.is_directed)
    for v in vertices:
        ig.add_vertex(v)
    for u, v in edges:
        ig.add_edge(u, v)
    ig.vs['attribute'] = attr
    ig.es['weight'] = wts
    return ig


def from_igraph(g):
    assert isinstance(g, igraph.Graph)
    weighted = g.is_weighted()
    directed = g.is_directed()

    if weighted:
        cg = WeightedGraph(directed=directed)
    else:
        cg = SimpleGraph(directed=directed)

    if g.is_named():
        vlist = g.vs['name']
    else:
        vlist = list(range(g.vcount()))

    if 'attribute' in g.vs.attribute_names():
        attr = g.vs['attribute']
    else:
        attr = [None for i in range(g.vcount())]

    for i,v in enumerate(vlist):
        if weighted:
            cg.add_vertex(v, attribute=attr[i])
        else:
            cg.add_vertex(v)

    weights = None
    if weighted:
        weights = g.es['weight']

    for i, edge in enumerate(g.get_edgelist()):
        if weighted:
            cg.add_edge(vlist[edge[0]], vlist[edge[1]], weight=weights[i])
        else:
            cg.add_edge(vlist[edge[0]], vlist[edge[1]])

    return cg

