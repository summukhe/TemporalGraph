import logging
from copy import deepcopy
from .structure_graph import *
from temporal_graph.pdb_processor import *
from temporal_graph.network_analysis import *

__all__ = ['mutation_evaluation']


def mutation_evaluation(pdb_structure,
                        residue_id,
                        site1,
                        site2,
                        method='mincut',
                        contact_radius=12,
                        potential='mj'):
    assert method in ['mincut', 'centrality']
    assert isinstance(pdb_structure, PDBStructure) or isinstance(pdb_structure, CaTrace)
    current_residue = pdb_structure.get_amino(residue_id).name(one_letter_code=True)
    if method == 'mincut':
        result = mutation_evaluation_by_mincut(pdb_structure,
                                               residue_id,
                                               site1,
                                               site2,
                                               contact_radius=contact_radius,
                                               potential=potential)
    else:
        result = mutation_evaluation_by_centrality(pdb_structure,
                                                   residue_id,
                                                   site1,
                                                   site2,
                                                   contact_radius=contact_radius,
                                                   potential=potential)
    assert current_residue in result
    return result


def mutation_evaluation_by_centrality(pdb_structure,
                                      resid,
                                      site1,
                                      site2,
                                      contact_radius=12,
                                      potential='mj'):
    logger = logging.getLogger(name="structure_network.mutation_evaluation_by_centrality")
    assert potential in ['mj']
    if isinstance(pdb_structure, PDBStructure):
        logger.debug('Extracting CA trace from the PDB structure')
        structure = pdb_to_catrace(pdb_structure)
    else:
        structure = deepcopy(pdb_structure)
    assert isinstance(structure, CaTrace)
    assert isinstance(site1, list) and isinstance(site2, list)
    n = structure.size
    assert n > 1
    residue_ids = structure.residue_ids
    residue_key = [structure.key(r) for r in residue_ids]
    assert resid in residue_ids
    for s in site1 + site2:
        assert s in residue_key
    all_aminos = valid_amino_acids(one_letter=True)

    result = {}
    for aa in all_aminos:
        logger.debug('Mutating residue %d to amino acid %s' % (resid, aa))
        structure.set_amino(resid, aa_type=aa)
        logger.debug('Current residue key: %s' % structure.key(resid))
        cg = contact_graph(structure,
                           cutoff=DistanceCutoff(def_cutoff=contact_radius),
                           potential=potential)
        g_inv = weight_inversion(cg)
        node_stats = between_groups_centrality(g_inv,
                                               group1=site1,
                                               group2=site2,
                                               weight=True,
                                               scale=100)
        result[aa] = node_stats
    return result


def mutation_evaluation_by_mincut(pdb_structure,
                                  resid,
                                  site1,
                                  site2,
                                  contact_radius=12,
                                  potential="mj"):
    logger = logging.getLogger(name="structure_network.mutation_evaluation_by_mincut")
    assert potential in ['mj']
    if isinstance(pdb_structure, PDBStructure):
        logger.info('Extracting CA trace from the PDB structure')
        structure = pdb_to_catrace(pdb_structure)
    else:
        structure = deepcopy(pdb_structure)
    assert isinstance(structure, CaTrace)
    assert isinstance(site1, list) and isinstance(site2, list)
    n = structure.size
    assert n > 1
    residue_ids = structure.residue_ids
    residue_key = [structure.key(r) for r in residue_ids]
    assert resid in residue_ids
    for s in site1 + site2:
        assert s in residue_key
    all_aminos = valid_amino_acids(one_letter=True)
    all_sites = set(site1 + site2)

    result = dict()
    for aa in all_aminos:
        logger.debug('Mutating residue %d to amino acid %s' % (resid, aa))
        structure.set_amino(resid, aa_type=aa)
        logger.debug('Current residue key: %s' % structure.key(resid))
        cg = contact_graph(structure,
                           cutoff=DistanceCutoff(def_cutoff=contact_radius),
                           potential=potential)
        g_inv = weight_inversion(cg)
        flow = maxflow(g_inv, src=site1, tgt=site2, weight=True)
        assert isinstance(flow, dict)
        marked_residues = dict()
        for u in flow.keys():
            for v in flow[u].keys():
                assert 'cut' in flow[u][v]
                for x, y in flow[u][v]['cut']:
                    if x not in all_sites:
                        if x not in marked_residues:
                            marked_residues[x] = 0
                        marked_residues[x] += 1
                    if y not in all_sites:
                        if y not in marked_residues:
                            marked_residues[y] = 0
                        marked_residues[y] += 1
        result[aa] = deepcopy(marked_residues)
    return result

