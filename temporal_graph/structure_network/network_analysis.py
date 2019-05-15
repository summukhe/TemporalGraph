import igraph
import numpy as np
from temporal_graph.pdb_processor import *
from temporal_graph.network_analysis import *
from temporal_graph.structure_network.structure_graph import *

__all__ = ['between_site_residues_by_stpath',
           'between_site_residues_by_mincut']


def between_site_residues_by_stpath(pdb_structure,
                                    site1,
                                    site2,
                                    type='energy',
                                    contact_radius=12):
    assert isinstance(pdb_structure, PDBStructure)
    assert isinstance(site1, list) and isinstance(site2, list)
    residue_key = [pdb_structure.key(r) for r in pdb_structure.residue_ids]
    for s in site1 + site2:
        assert s in residue_key
    if type == 'energy':
        g = contact_energy_graph(pdb_structure,
                                 contact_radius=contact_radius)
    elif type == 'mj':
        g = contact_graph(pdb_to_catrace(pdb_structure),
                          cutoff=DistanceCutoff(def_cutoff=contact_radius),
                          potential='mj')
    g_inv = weight_inversion(g)
    node_stats = between_groups_centrality(g_inv,
                                           group1=site1,
                                           group2=site2,
                                           weight=True,
                                           scale=100)
    return node_stats


def between_site_residues_by_mincut(pdb_structure,
                                    site1,
                                    site2,
                                    type='energy',
                                    contact_radius=12):
    assert isinstance(pdb_structure, PDBStructure)
    assert isinstance(site1, list) and isinstance(site2, list)
    residue_key = [pdb_structure.key(r) for r in pdb_structure.residue_ids]
    for s in site1 + site2:
        assert s in residue_key
    if type == 'energy':
        g = contact_energy_graph(pdb_structure,
                                 contact_radius=contact_radius)
    elif type == 'mj':
        g = contact_graph(pdb_to_catrace(pdb_structure),
                          cutoff=DistanceCutoff(def_cutoff=contact_radius),
                          potential='mj')
    g_inv = weight_inversion(g)
    cuts = maxflow(g_inv, src=site1, tgt=site2, weight=True)
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
