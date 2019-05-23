import logging
import numpy as np
from copy import deepcopy
from .structure_graph import *
from .all_path_analysis import *
from temporal_graph.spatial_ds import *
from temporal_graph.pdb_processor import *
from temporal_graph.network_analysis import *

__all__ = ['mutation_evaluation', 'MutationEffectScore', 'all_path_mutation_evaluation']


class MutationEffectScore:
    def __init__(self, residue_id, base_amino, true_energy=False, npath=1):
        assert isinstance(residue_id, np.int)
        assert isinstance(npath, np.int) and (npath > 0)
        self.__residue_id = residue_id
        self.__ref_amino = base_amino
        self.__scores = dict()
        self.__energy = true_energy
        self.__npath = npath
        self.__logger = logging.getLogger(name="temporal_graph.MutationEffectScore")

    @property
    def normalized_score(self):
        ref_aa = self.ref_amino
        all_aa = self.keys
        assert ref_aa in self.__scores
        effect = dict()
        if self.__energy:
            ref_score = deepcopy(self.__scores[ref_aa])
            ref_moment = np.sum(list(ref_score.values()))
            ref_size = len(ref_score)
            for aa in self.__scores.keys():
                tgt_score = deepcopy(self.__scores[aa])
                tgt_moment = np.sum(list(tgt_score.values()))
                tgt_size = len(tgt_score)
                effect[aa] = (tgt_moment * (ref_size/tgt_size) - ref_moment)/self.__npath
        else:
            def solve_rank_order(x):
                assert isinstance(x, dict)
                order = list(reversed(sorted([(k, v) for k, v in x.items()], key=lambda p: p[1])))
                rank, rank_lookup, last_val = 1, {}, None
                for i, pair in enumerate(order):
                    assert (last_val is None) or (last_val >= pair[1])
                    rank_lookup[pair[0]] = rank
                    if last_val == pair[1]:
                        rank += 1
                    last_val = pair[1]
                return rank, rank_lookup

            scores_ref = self.__scores[ref_aa]
            max_rank, rank_ref = solve_rank_order(scores_ref)
            for aa in all_aa:
                scores_tgt = self.__scores[aa]
                mr, rank_tgt = solve_rank_order(scores_tgt)
                f = mr / max_rank
                common_set = set(list(rank_ref.keys())).intersection(set(list(rank_tgt.keys())))
                s = 0
                for i in common_set:
                    p = (scores_tgt[i] - scores_ref[i])**2
                    s += p * (rank_tgt[i] - rank_ref[i]*f)
                effect[aa] = s
        return effect

    @property
    def is_complete(self):
        return len(self.__scores) == len(valid_amino_acids(one_letter=True))

    @property
    def base_residue(self):
        return self.__residue_id

    @property
    def ref_amino(self):
        return self.__ref_amino

    @property
    def keys(self):
        return list(self.__scores.keys())

    @property
    def to_dict(self):
        return self.__scores

    def __getitem__(self, item):
        assert item in self.__scores
        return self.__scores

    def __setitem__(self, key, value):
        assert key in valid_amino_acids(one_letter=True)
        assert isinstance(value, dict)
        self.__scores[key] = value

    def __contains__(self, key):
        return key in self.__scores.keys()


def mutation_evaluation(pdb_structure,
                        residue_id,
                        site1,
                        site2,
                        method='mincut',
                        contact_radius=12,
                        potential='mj'):
    assert method in ['mincut', 'centrality']
    assert isinstance(pdb_structure, PDBStructure) or isinstance(pdb_structure, CaTrace)
    if isinstance(pdb_structure, PDBStructure):
        structure = pdb_to_catrace(pdb_structure)
    else:
        structure = deepcopy(pdb_structure)
    current_residue = structure.get_amino(residue_id).name(one_letter_code=True)
    if method == 'mincut':
        result = mutation_evaluation_by_mincut(structure,
                                               residue_id,
                                               site1,
                                               site2,
                                               contact_radius=contact_radius,
                                               potential=potential)
    else:
        result = mutation_evaluation_by_centrality(structure,
                                                   residue_id,
                                                   site1,
                                                   site2,
                                                   contact_radius=contact_radius,
                                                   potential=potential)
    assert current_residue in result
    return current_residue, result


def mutation_evaluation_by_centrality(pdb_structure,
                                      resid,
                                      site1,
                                      site2,
                                      contact_radius=12,
                                      potential='mj'):
    logger = logging.getLogger(name="structure_network.mutation_evaluation_by_centrality")
    assert potential in ['mj', 'charmm']
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
    assert resid in residue_ids
    for s in site1 + site2:
        assert s in residue_ids
    all_aminos = valid_amino_acids(one_letter=True)
    score = MutationEffectScore(resid, structure.get_amino(resid).name(one_letter_code=True))
    for aa in all_aminos:
        logger.debug('Mutating residue %d to amino acid %s' % (resid, aa))
        structure.set_amino(resid, aa_type=aa)
        logger.debug('Current residue key: %s' % structure.key(resid))
        cg = contact_graph(structure,
                           cutoff=DistanceCutoff(def_cutoff=contact_radius),
                           potential=potential)
        g_inv = weight_inversion(cg)
        site1_keys = [structure.key(r) for r in site1]
        site2_keys = [structure.key(r) for r in site2]
        node_stats = between_groups_centrality(g_inv,
                                               group1=site1_keys,
                                               group2=site2_keys,
                                               weight=True,
                                               scale=100)
        score[aa] = node_stats
    return score


def mutation_evaluation_by_mincut(pdb_structure,
                                  resid,
                                  site1,
                                  site2,
                                  contact_radius=12,
                                  potential="mj"):
    logger = logging.getLogger(name="structure_network.mutation_evaluation_by_mincut")
    assert potential in ['mj', 'charmm']
    if isinstance(pdb_structure, PDBStructure):
        logger.info('Extracting CA trace from the PDB structure')
        structure = pdb_to_catrace(pdb_structure)
    else:
        structure = deepcopy(pdb_structure)
    assert isinstance(structure, CaTrace)
    assert isinstance(site1, list) and isinstance(site2, list)
    assert structure.size > 1
    residue_ids = structure.residue_ids
    assert resid in residue_ids
    for s in site1 + site2:
        assert s in residue_ids
    all_aminos = valid_amino_acids(one_letter=True)
    all_sites = set(site1 + site2)
    result = MutationEffectScore(resid, structure.get_amino(resid).name(one_letter_code=True))
    for aa in all_aminos:
        logger.debug('Mutating residue %d to amino acid %s' % (resid, aa))
        structure.set_amino(resid, aa_type=aa)
        logger.debug('Current residue key: %s' % structure.key(resid))
        cg = contact_graph(structure,
                           cutoff=DistanceCutoff(def_cutoff=contact_radius),
                           potential=potential)
        site1_key = [structure.key(r) for r in site1]
        site2_key = [structure.key(r) for r in site2]
        flow = maxflow(cg, src=site1_key, tgt=site2_key, weight=True)
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


def all_path_mutation_evaluation(pdb_structure,
                                 residue_id,
                                 site1,
                                 site2,
                                 method='centrality',
                                 flow_forward=True,
                                 ntrails=25,
                                 minimum_energy=0.5,
                                 maximum_distance=12,
                                 min_path_length=3,
                                 max_path_length=10,
                                 allowed_path_overlap=0.5,
                                 backtrack_fraction=0.5):
    logger = logging.getLogger(name="temporal_network.all_path_mutation_evaluation")
    if isinstance(pdb_structure, PDBStructure):
        logger.info('Extracting CA trace from the PDB structure')
        structure = pdb_to_catrace(pdb_structure)
    else:
        structure = deepcopy(pdb_structure)
    assert isinstance(structure, CaTrace)
    assert method in {'mincut', 'centrality'}
    assert isinstance(site1, list) and isinstance(site2, list)
    assert isinstance(backtrack_fraction, float) and isinstance(allowed_path_overlap, float)
    assert (backtrack_fraction >= 0.0) and (backtrack_fraction <= 1.0)
    assert structure.size > 1
    residue_ids = structure.residue_ids
    assert residue_id in residue_ids
    for s in site1 + site2:
        assert s in residue_ids
    all_aminos = valid_amino_acids(one_letter=True)
    curr_residue = structure.get_amino(residue_id)
    result = MutationEffectScore(residue_id,
                                 structure.get_amino(residue_id).name(one_letter_code=True),
                                 true_energy=True,
                                 npath=ntrails)
    for aa in all_aminos:
        logger.debug('Mutating residue %d to amino acid %s' % (residue_id, aa))
        structure.set_amino(residue_id, aa_type=aa)
        logger.debug('Current residue key: %s' % structure.key(residue_id))
        cg = contact_graph(structure,
                           cutoff=maximum_distance + 1.0,
                           potential='charmm')
        site1_key = [structure.key(r) for r in site1]
        site2_key = [structure.key(r) for r in site2]
        if method == "centrality":
            logger.debug('Calculating all path trace based centrality measure!!')
            __, scores = all_path_centrality(cg,
                                             source=site1_key,
                                             target=site2_key,
                                             forward_path=flow_forward,
                                             ntrails=ntrails,
                                             maximum_distance=maximum_distance,
                                             minimum_weight=minimum_energy,
                                             min_path_length=min_path_length,
                                             max_path_length=max_path_length,
                                             allowed_path_overlap=allowed_path_overlap,
                                             backtrack_fraction=backtrack_fraction)
        elif method == "mincut":
            logger.debug('Calculating all path trace based mincut method!!')
            __, scores = all_path_maxflow_score(cg,
                                                source=site1_key,
                                                target=site2_key,
                                                forward_path=flow_forward,
                                                ntrails=ntrails,
                                                maximum_distance=maximum_distance,
                                                minimum_weight=minimum_energy,
                                                min_path_length=min_path_length,
                                                max_path_length=max_path_length,
                                                allowed_path_overlap=allowed_path_overlap,
                                                backtrack_fraction=backtrack_fraction)
        result[aa] = deepcopy(scores)
    return curr_residue, result






