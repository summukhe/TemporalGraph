import logging
import numpy as np
import pandas as pd
from .graph import *
from .graph_adapter import *


__all__ = ['diameter', 'betweenness', 'edge_betweenness',
           'shortest_path', 'get_all_shortest_paths',
           'between_groups_centrality', 'gomory_hu_cuts',
           'maxflow', 'weight_inversion']


def diameter(g, weight=True):
    ig = to_igraph(g)
    weighted = ig.is_weighted() and weight
    weights = 'weight' if weighted else None
    return ig.diameter(directed=ig.is_directed(),
                       unconn=ig.is_connected(),
                       weights=weights)


def betweenness(g, vertex_set=None, weight=True, max_length=None):
    logger = logging.getLogger('network_analysis.betweenness')
    ig = to_igraph(g)
    if (max_length is not None) and (max_length > g.order//2):
        logger.warning('Long max_length, cutoff provided.')

    if vertex_set is None:
        logger.info('No vertex set defined defaulting to all vertex evaluation')
        vertex_set = g.vertices
    assert isinstance(vertex_set, list)
    vertex_set = [v for v in vertex_set if g.is_vertex(v)]
    weighted = ig.is_weighted() and weight
    weights = 'weight' if weighted else None
    vlist = ig.vs['name'] if ig.is_named() else list(range(ig.vcount()))
    centrality = ig.betweenness(vertices=vertex_set,
                                directed=g.is_directed,
                                cutoff=max_length,
                                weights=weights)
    return {vlist[i]: c for i, c in enumerate(centrality)}


def edge_betweenness(g, weight=True, max_length=None):
    logger = logging.getLogger('network_analysis.edge_betweenness')
    ig = to_igraph(g)
    weighted = ig.is_weighted() and weight
    vlist = ig.vs['name'] if ig.is_named() else list(range(ig.vcount()))
    weights = 'weight' if weighted else None
    edge_list = [(vlist[u], vlist[v]) for u, v in ig.get_edgelist()]
    ebet = ig.edge_betweenness(directed=ig.is_directed(), cutoff=max_length, weights=weights)
    assert len(edge_list) == len(ebet)
    return [(edge_list[i][0], edge_list[i][1], c) for i, c in enumerate(ebet)]


def shortest_path(g, src=None, tgt=None, weight=True):
    ig = to_igraph(g)
    if src is None:
        src = g.vertices
    if tgt is None:
        tgt = g.vertices

    if not isinstance(src, list):
        src = [src]
    if not isinstance(tgt, list):
        tgt = [tgt]

    for v in src:
        assert g.is_vertex(v)

    for v in tgt:
        assert g.is_vertex(v)

    weighted = ig.is_weighted() and weight
    if weighted:
        res = pd.DataFrame(ig.shortest_paths(source=src,
                                             target=tgt,
                                             weights='weight'),
                           index=src,
                           columns=tgt)
    else:
        res = pd.DataFrame(ig.shortest_paths(source=src,
                                             target=tgt),
                           index=src,
                           columns=tgt)
    return res


def get_all_shortest_paths(g, src=None, tgt=None, weight=True, min_path_length=2):
    logger = logging.getLogger('network_analysis.get_all_shortest_paths')
    ig = to_igraph(g)
    if ig.is_named():
        vlist = ig.vs['name']
    else:
        vlist = list(range(ig.vcount()))

    if min_path_length < 2:
        logger.info('Expect minimum shortest path distance >= 2')
        min_path_length = 2

    if src is None:
        src = g.vertices
    if tgt is None:
        tgt = g.vertices

    if not isinstance(src, list):
        src = [src]
    if not isinstance(tgt, list):
        tgt = [tgt]

    for v in src:
        assert g.is_vertex(v)

    for v in tgt:
        assert g.is_vertex(v)

    weighted = ('weight' in ig.es.attribute_names()) and weight

    paths = dict()
    for v in src:
        if weighted:
            path_list = ig.get_all_shortest_paths(v=v, to=tgt, weights='weight')
        else:
            path_list = ig.get_all_shortest_paths(v=v, to=tgt)
        for path in path_list:
            assert isinstance(path, list)
            if len(path) >= min_path_length:
                assert vlist[path[0]] == v
                end_vertex = vlist[path[-1]]
                if v not in paths:
                    paths[v] = dict()
                if end_vertex not in paths[v]:
                    paths[v][end_vertex] = list()
                paths[v][end_vertex].append([vlist[n] for n in path])
    return paths


def between_groups_centrality(g,
                              group1,
                              group2,
                              weight=True,
                              min_path_length=3,
                              scale=1.):
    logger = logging.getLogger('get_between_groups_centrality')
    assert isinstance(group1, list) and isinstance(group2, list)
    paths = get_all_shortest_paths(g,
                                   src=group1,
                                   tgt=group2,
                                   weight=weight,
                                   min_path_length=min_path_length)
    node_stat, path_count = dict(), 0
    for v in paths:
        for u in paths[v]:
            npaths = len(paths[v][u])
            for i in range(npaths):
                path_length = len(paths[v][u][i])
                if path_length > 2:
                    path_count += 1
                    for n in paths[v][u][i][1:-1]:
                        if n not in node_stat:
                            node_stat[n] = 0
                        node_stat[n] += 1
    logger.debug('Total number of shortest paths are (%d)' % path_count)
    for n in node_stat:
        node_stat[n] = (node_stat[n] * scale / path_count)
    return node_stat


def gomory_hu_cuts(g, weight=True):
    logger = logging.getLogger(name="network_analysis.gomory_hu_cuts")
    ig = to_igraph(g)
    weighted = ig.is_weighted() and weight

    if weighted:
        logger.debug('Given graph is weighed!!')
    if ig.is_directed():
        logger.debug('Given graph is directed in nature!!')

    if weighted:
        gtree = ig.gomory_hu_tree(capacity='weight')
    else:
        gtree = ig.gomory_hu_tree()

    if ig.is_named():
        vlist = ig.vs['name']
    else:
        vlist = list(range(ig.vcount()))

    edge_list = gtree.get_edgelist()
    flows = gtree.es['flow']
    assert len(flows) == len(edge_list)

    gomory_cuts = []
    for i, edge in enumerate(edge_list):
        gomory_cuts.append((vlist[edge[0]], vlist[edge[1]], flows[i]))
    return gomory_cuts


def maxflow(g, src=None, tgt=None, weight=True):
    logger = logging.getLogger(name="network_analysis.maxflow")
    ig = to_igraph(g)

    if src is None:
        logger.info('No source list specified defaulting to all vertices')
        src = g.vertices
    elif not isinstance(src, list) and g.is_vertex(src):
        logger.info('Only single source specified')
        src = [src]
    assert isinstance(src, list)
    src = [v for v in src if g.is_vertex(v)]
    if tgt is None:
        logger.info('No target list specified defaulting to all vertices!')
        tgt = g.vertices
    elif not isinstance(tgt, list) and g.is_vertex(tgt):
        logger.info('Only single target specified!')
        tgt = [tgt]
    assert isinstance(tgt, list)
    tgt = [v for v in src if g.is_vertex(v)]

    assert len(src) > 0 and len(tgt) > 0
    weighted = ig.is_weighted() and weight

    capacities = dict()
    vlist = ig.vs['name']
    elist = [(vlist[u], vlist[v]) for u, v in ig.get_edgelist()]
    vmap = {v: i for i, v in enumerate(vlist)}
    capacity = 'weight' if weighted else None
    for v in src:
        capacities[v] = dict()
        for u in tgt:
            if u != v:
                capacities[v][u] = dict()
                g_cut = ig.maxflow(source=vmap[v], target=vmap[u], capacity=capacity)
                capacities[v][u]['flow'] = g_cut.value
                capacities[v][u]['cut'] = [elist[c] for c in g_cut.cut]
                capacities[v][u]['cover'] = [[vlist[n] for n in cover] for cover in list(g_cut.as_cover())]
    return capacities


def weight_inversion(g, buffer=1.):
    assert isinstance(g, WeightedGraph)
    assert buffer >= 0.

    vlist = g.vertices
    attr = [g.attribute(v) for v in vlist]

    elist = g.edges
    wlist = [g.weight(u,v) for u,v in elist]
    max_wt, min_wt = np.max(wlist), np.minimum(np.min(wlist),0)
    inv_wlist = (max_wt - min_wt - np.array(wlist) + buffer).tolist()
    assert len(inv_wlist) == len(elist)
    gcopy = WeightedGraph(directed=g.is_directed)

    for i,v in enumerate(vlist):
        gcopy.add_vertex(v, attribute=attr[i])

    for i, edge in enumerate(elist):
        gcopy.add_edge(src=edge[0], dst=edge[1], weight=inv_wlist[i])

    return gcopy