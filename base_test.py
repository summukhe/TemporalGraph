import os
import json
import logging
import numpy as np
from temporal_graph import *
from matplotlib import pyplot as plt


if __name__ == "__main__":
    pdbfile="/home/sumanta/Project/structural_dynamics/second_round_analysis/" + \
            "dynamics/electrostatics_analysis/5MFX/mode_catraj_7.pdb"
    atpsitefile="/home/sumanta/Project/structural_dynamics/second_round_analysis/" + \
            "ATPSites/5MFX.txt"
    rnasitefile="/home/sumanta/Project/structural_dynamics/second_round_analysis/" + \
            "RNASites/5MFX.txt"

    assert os.path.isfile(pdbfile)

    trajectory = read_trajectory_catrace(pdbfile)
    assert len(trajectory) > 0
    # idx = np.random.randint(len(trajectory))
    idx = 0
    residue_ids = trajectory[idx]['A'].residue_ids()

    print('Running for snapshot: %d' % idx)

    assert os.path.isfile(atpsitefile)
    with open(atpsitefile, "r") as f:
        atpsite_residues = [int(r) for r in f.readlines()]
    atpsite_residues = ['%s%d' % (trajectory[idx]['A'].get_amino(r), r) for r in atpsite_residues if r in residue_ids]

    assert os.path.isfile(rnasitefile)
    with open(rnasitefile, "r") as f:
        rnasite_residues = [int(r) for r in f.readlines()]
    rnasite_residues = ['%s%d' % (trajectory[idx]['A'].get_amino(r),r) for r in rnasite_residues if r in residue_ids]

    g = contact_graph(trajectory[idx]['A'])
    #pulse = PulseSignal(ts=0, te=4)
    pulse = PeriodicPulse(ts=0, span=5, periods=25, strength=10)
    n = np.minimum(1000, g.order())
    process = DiffusionOnStaticGraph(dissipation_rate=0.1)
    process.start(g, atpsite_residues, pulse)
    process.n_steps(n=n)

    p = np.minimum(10, len(rnasite_residues))
    for r in rnasite_residues[:p]:
        hist = process.get_smooth_history(r)
        plt.plot(hist['time'], hist['value'])
    plt.legend(loc='best', labels=rnasite_residues[:p])
    plt.show()

    v_info = {v: process.information_carried(v) for v in process.vertices()}
    v_ordered = reversed([k for k, v in sorted(v_info.items(), key=lambda kv: (kv[1], kv[0]))])
    for i, v in enumerate(v_ordered):
        print("%d: %s" % (i, v))

    diffusion_trj = mark_diffusion_on_structure(ca_trace=trajectory[idx]['A'],
                                                pulse=pulse,
                                                signal_nodes=atpsite_residues,
                                                dissipation_rate=0.1)
    with open('/tmp/u.pdb', 'w') as fh:
        for i, trj in enumerate(diffusion_trj):
            fh.write('MODEL     %4d\n' % (i+1))
            trj.write(fh)
            fh.write('ENDMDL\n')
    exit(0)

    # logging.basicConfig(level=logging.DEBUG)
    cnodes = g.vertices()
    np.random.shuffle(cnodes)

    scores = analyze_vertex_contribution(g,
                                         start_nodes=atpsite_residues,
                                         pulse=pulse,
                                         target_nodes=rnasite_residues,
                                         cycles=n,
                                         contrib_nodes=cnodes[:4],
                                         min_variation_score=1e-6,
                                         dissipation_rate=0.0)
    print(json.dumps(scores))

