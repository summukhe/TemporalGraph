import numpy as np
from copy import deepcopy
from temporal_graph.network_process import *
from temporal_graph.structure_network import *
from temporal_graph.pdb_processor import *

__version__ = "1.0"
__all__ = ["mark_diffusion_on_structure"]


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
        diffusion_step = ca_trace.size
    assert (report_step > 0) and (report_step < diffusion_step)

    diffusion_process = DiffusionOnStaticGraph(dissipation_rate=dissipation_rate)
    c_graph = contact_graph(ca_trace, cutoff=cutoff, potential=potential)
    diffusion_process.start(g=c_graph,
                            pulse=pulse,
                            start_nodes=signal_nodes,
                            blocked=blocked_nodes)
    diffusion_process.n_steps(diffusion_step)
    res_ids = ca_trace.residue_ids
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


