import igraph
import numpy as np
from copy import deepcopy
from temporal_graph.pdb_processor import *
from temporal_graph.network_analysis import *
from temporal_graph.structure_network.structure_graph import *

__all__ = ['between_site_residues_by_stpath',
           'between_site_residues_by_mincut']


def between_site_residues_by_stpath(pdb_structure,
                                    site1,
                                    site2,
                                    potential='energy',
                                    contact_radius=12):
    assert potential in {'energy', 'mj', 'charmm'}
    assert isinstance(pdb_structure, PDBStructure) or isinstance(pdb_structure, CaTrace)
    assert isinstance(site1, list) and isinstance(site2, list)
    residue_key = [pdb_structure.key(r) for r in pdb_structure.residue_ids]
    for s in site1 + site2:
        assert s in residue_key
    g = contact_graph(pdb_structure,
                      cutoff=contact_radius,
                      potential=potential)
    g_inv = weight_inversion(g)
    node_stats = between_groups_centrality(g_inv,
                                           group1=site1,
                                           group2=site2,
                                           weight=True,
                                           scale=100)
    return node_stats


def robust_residue_importance_between_sites_by_stpath(pdb_structure,
                                                      site1,
                                                      site2,
                                                      min_distance=4.5,
                                                      max_distance=12,
                                                      nsample=15,
                                                      potential='charmm'):
    assert min_distance < max_distance
    assert nsample > 1
    assert energy in ['energy', 'mj', 'charmm']
    assert isinstance(pdb_structure, CaTrace) or isinstance(pdb_structure, PDBStructure)
    residue_ids = pdb_structure.residue_ids
    residue_keys = {pdb_structure.key(r) for r in residue_ids}
    for s in site1 + site2:
        assert s in residue_keys
    for distance in np.linspace(start=min_distance,
                                stop=max_distance,
                                num=nsample,
                                endpoint=True):
        value = distance


def between_site_residues_by_mincut(pdb_structure,
                                    site1,
                                    site2,
                                    potential='charmm',
                                    contact_radius=12):
    assert isinstance(pdb_structure, PDBStructure) or \
           isinstance(pdb_structure, CaTrace)
    assert isinstance(site1, list) and isinstance(site2, list)
    residue_key = [pdb_structure.key(r) for r in pdb_structure.residue_ids]
    for s in site1 + site2:
        assert s in residue_key
    g = contact_graph(pdb_structure,
                      cutoff=contact_radius,
                      potential=potential)
    cuts = maxflow(g, src=site1, tgt=site2, weight=True)
    assert isinstance(cuts, dict)
    residue_marking, valid_edges = dict(), 0
    all_sites = set(site1 + site2)
    for u in cuts.keys():
        for v in cuts[u].keys():
            assert 'cut' in cuts[u][v]
            if len(cuts[u][v]['cut']) > 0:
                for x, y in cuts[u][v]['cut']:
                    valid = False
                    if x not in all_sites:
                        if x not in residue_marking:
                            residue_marking[x] = 0
                        residue_marking[x] += 1
                        valid = True
                    if y not in all_sites:
                        if y not in residue_marking:
                            residue_marking[y] = 0
                        residue_marking[y] += 1
                        valid = True
                    if valid:
                        valid_edges += 1

    return valid_edges, residue_marking
