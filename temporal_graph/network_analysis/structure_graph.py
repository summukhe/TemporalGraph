import numpy as np
from copy import deepcopy
from .graph import WeightedGraph
from temporal_graph.pdb_processor import *
from temporal_graph.network_process import *


__author__ = "Sumanta Mukherjee"
__all__ = ['DistanceCutoff', 'contact_graph', 'mark_diffusion_on_structure']


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


def contact_graph(ca_trace, cutoff=DistanceCutoff(), potential='mj'):
    assert isinstance(ca_trace, CaTrace)
    assert isinstance(cutoff, DistanceCutoff)
    res_ids = ca_trace.residue_ids()
    c_graph = WeightedGraph(directed=False)
    for ri in res_ids:
        amino_i = ca_trace.get_amino(ri)
        x_i, y_i, z_i = ca_trace.xyz(ri)
        for rj in res_ids:
            if ri < rj:
                amino_j = ca_trace.get_amino(rj)
                x_j, y_j, z_j = ca_trace.xyz(rj)
                c = cutoff(amino_i, amino_j)
                d = np.sqrt((x_i-x_j)**2 + (y_i-y_j)**2 + (z_i-z_j)**2)
                if d <= c:
                    p = get_pair_potential(amino_i, amino_j, pot_type=potential)
                    c_graph.add_edge('%s%d' % (amino_i, ri),
                                     '%s%d' % (amino_j, rj),
                                     weight=p)
    return c_graph


def mark_diffusion_on_structure(ca_trace,
                                pulse,
                                signal_nodes,
                                blocked_nodes=[],
                                diffusion_step=None,
                                report_step=10,
                                dissipation_rate=0.0,
                                cutoff=DistanceCutoff(),
                                potential='mj'):
    assert isinstance(ca_trace, CaTrace)
    assert isinstance(pulse, SignalInput)
    assert isinstance(signal_nodes, list)
    if diffusion_step is None:
        diffusion_step = ca_trace.size()
    assert (report_step > 0) and (report_step < diffusion_step)

    diffusion_process = DiffusionOnStaticGraph(dissipation_rate=dissipation_rate)
    c_graph = contact_graph(ca_trace, cutoff=cutoff, potential=potential)
    diffusion_process.start(g=c_graph,
                            pulse=pulse,
                            start_nodes=signal_nodes,
                            blocked=blocked_nodes)
    diffusion_process.n_steps(diffusion_step)
    res_ids = ca_trace.residue_ids()
    nodes = ["%s%d" % (ca_trace.get_amino(r),r) for r in res_ids]
    n_size = diffusion_step // report_step
    process_trj = [deepcopy(ca_trace) for i in range(n_size)]
    diffusion_values = np.array([diffusion_process.get_history(n)['value'] for n in nodes])
    cmax = np.max(diffusion_values)
    diffusion_values = diffusion_values * (100./cmax)
    for i in range(n_size):
        for j, r in enumerate(res_ids):
            process_trj[i].b_factor(r, diffusion_values[j, i*report_step])
    return process_trj
