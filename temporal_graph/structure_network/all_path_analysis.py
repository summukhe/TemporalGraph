import numpy as np
from temporal_graph.network_analysis import *


__version__ = "1.0"
__all__ = ['all_path_centrality', 'all_path_maxflow', 'all_path_maxflow_score']


def all_path_centrality(g,
                        source,
                        target,
                        forward_path=True,
                        minimum_weight=0.5,
                        maximum_distance=8,
                        ntrails=25,
                        min_path_length=3,
                        max_path_length=12,
                        allowed_path_overlap=0.5,
                        backtrack_fraction=0.5):
    assert isinstance(g, GeometricGraph3d)
    assert isinstance(source, list)
    assert isinstance(target, list)
    assert isinstance(allowed_path_overlap, float)
    assert (allowed_path_overlap >= 0) and (allowed_path_overlap <= 1)
    assert isinstance(ntrails, int) and ntrails > 0
    for n in source + target:
        assert g.is_vertex(n)
    assert minimum_weight >= 0
    assert maximum_distance > 0
    path_filter = GeometricPathFilter(forward=forward_path,
                                      weight_cutoff=minimum_weight,
                                      distance_cutoff=maximum_distance)
    path_iterator = AllPathIterator(g,
                                    visitor=path_filter,
                                    path_separation=allowed_path_overlap)

    all_paths, scores, total = [], [], 0
    for v in source:
        result = path_iterator.all_path(v,
                                        stop_vertex=target,
                                        max_paths=ntrails,
                                        min_path_length=min_path_length,
                                        max_path_length=max_path_length,
                                        retrack_fraction=backtrack_fraction)
        if len(result) > 0:
            for res in result:
                assert isinstance(res, dict)
                assert ('path' in res) and ('weights' in res)
                lim_weight = np.min(res['weights'])
                all_paths.append(res['path'])
                scores.append(lim_weight)
                total += np.sum(res['weights'])
    assert len(all_paths) == len(scores)
    cen = {}
    for i, path in enumerate(all_paths):
        for node in path[1:-1]:
            if node not in cen:
                cen[node] = 0
            cen[node] += scores[i]
    return total, cen


def all_path_maxflow(g,
                     source,
                     target,
                     forward_path=True,
                     ntrails=25,
                     minimum_weight=0.5,
                     maximum_distance=8,
                     min_path_length=3,
                     max_path_length=12,
                     allowed_path_overlap=0.5,
                     backtrack_fraction=0.5):
    assert isinstance(g, GeometricGraph3d)
    assert isinstance(source, list)
    assert isinstance(target, list)
    for n in source + target:
        assert g.is_vertex(n)
    assert minimum_weight >= 0
    assert maximum_distance > 0
    assert isinstance(ntrails,int) and ntrails > 0
    assert isinstance(allowed_path_overlap, float)
    assert (allowed_path_overlap >= 0.0) and (allowed_path_overlap <= 1.0)
    path_filter = GeometricPathFilter(forward=forward_path,
                                      weight_cutoff=minimum_weight,
                                      distance_cutoff=maximum_distance)
    path_iterator = AllPathIterator(g, visitor=path_filter, path_separation=allowed_path_overlap)

    all_edges, vertices = dict(), set()
    for v in source:
        result = path_iterator.all_path(v,
                                        stop_vertex=target,
                                        ntrails=ntrails,
                                        min_path_length=min_path_length,
                                        max_path_length=max_path_length,
                                        retrack_fraction=backtrack_fraction)
        if len(result) > 0:
            for res in result:
                assert isinstance(res, dict)
                assert ('path' in res) and ('weights' in res)
                for i in range(1, len(res['path'])):
                    u = res['path'][i-1]
                    v = res['path'][i]
                    w = res['weights'][i-1]
                    vertices.add(u)
                    vertices.add(v)
                    if u not in all_edges:
                        all_edges[u] = {}
                    if v not in all_edges[u]:
                        all_edges[u][v] = w
    if len(vertices) > 0:
        ig = WeightedGraph(directed=g.is_directed)
        for u in all_edges:
            for v in all_edges[u]:
                ig.add_edge(u, v, weight=all_edges[u][v])
        return maxflow(ig, src=source, tgt=target, weight=True)
    return None


def all_path_maxflow_score(g,
                           source,
                           target,
                           forward_path=True,
                           ntrails=25,
                           minimum_weight=0.5,
                           maximum_distance=8,
                           min_path_length=3,
                           max_path_length=12,
                           allowed_path_overlap=0.5,
                           backtrack_fraction=0.5):
    cuts = all_path_maxflow(g,
                            source,
                            target,
                            forward_path=forward_path,
                            ntrails=ntrails,
                            minimum_weight=minimum_weight,
                            maximum_distance=maximum_distance,
                            min_path_length=min_path_length,
                            max_path_length=max_path_length,
                            allowed_path_overlap=allowed_path_overlap,
                            backtrack_fraction=backtrack_fraction)
    scores, total = {}, 0
    if (cuts is not None) and isinstance(cuts, dict):
        for u in cuts.keys():
            for v in cuts[u].keys():
                if ('flow' in cuts[u][v]) and ('cut' in cuts[u][v]):
                    for x, y in cuts[u][v]['cut']:
                        if x not in scores:
                            scores[x] = 0
                        if y not in scores:
                            scores[y] = 0
                        scores[x] += cuts[u][v]['flow']
                        scores[y] += cuts[u][v]['flow']
                        total += 1
        for s in scores.keys():
            scores[s] = scores[s] * 1.0 / total
    return total, scores
