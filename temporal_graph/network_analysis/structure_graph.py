import numpy as np
from .graph import WeightedGraph
from temporal_graph.pdb_processor import *


__author__ = "Sumanta Mukherjee"
__all__ = ['DistanceCutoff', 'contact_graph']


class DistanceCutoff:
    def __init__(self, def_cutoff=8.0):
        self.__cutoff = def_cutoff

    def __call__(self, amino1, amino2):
        if isinstance(amino1, str):
            amino1 = get_amino(amino1)
        if isinstance(amino2, str):
            amino2 = get_amino(amino2)
        assert isinstance(amino1,AminoAcid) and isinstance(amino2,AminoAcid)
        return self.__cutoff


def contact_graph(catrace, cutoff=DistanceCutoff(), potential='mj'):
    assert isinstance(catrace, CaTrace)
    assert isinstance(cutoff, DistanceCutoff)
    resids = catrace.residue_ids()
    cgraph = WeightedGraph(directed=False)
    for ri in resids:
        amino_i = catrace.get_amino(ri)
        x_i, y_i, z_i = catrace.xyz(ri)
        for rj in resids:
            if ri < rj:
                amino_j = catrace.get_amino(rj)
                x_j, y_j, z_j = catrace.xyz(rj)
                c = cutoff(amino_i, amino_j)
                d = np.sqrt( (x_i-x_j)**2 + (y_i-y_j)**2 + (z_i-z_j)**2)
                if d <= c:
                    p = get_pair_potential(amino_i, amino_j, pot_type=potential)
                    cgraph.add_edge( '%s%d' % (amino_i, ri), '%s%d' % (amino_j, rj), weight=p)
    return cgraph