import os
import numpy as np
from temporal_graph import *
from matplotlib import pyplot as plt


if __name__ == "__main__":
    pdbfile="/home/sumanta/Project/structural_dynamics/second_round_analysis/" + \
            "dynamics/electrostatics_analysis/1A1V/mode_catraj_7.pdb"
    atpsitefile="/home/sumanta/Project/structural_dynamics/second_round_analysis/" + \
            "ATPSites/1A1V.txt"
    rnasitefile="/home/sumanta/Project/structural_dynamics/second_round_analysis/" + \
            "RNASites/1A1V.txt"

    assert os.path.isfile(pdbfile)

    trajectory = read_trajectory_catrace(pdbfile)
    assert len(trajectory) > 0
    idx = np.random.randint(len(trajectory))
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
    process = DiffusionOnStaticGraph(dissipation_rate=0.3)
    #pulse = PulseSignal(ts=0, te=5)
    pulse = PeriodicPulse(ts=0, span=5, periods=25, strength=100.)
    process.start(g, atpsite_residues, pulse)
    n = np.minimum(1000, g.order())
    for i in range(n):
        process.step()

    p = np.minimum(10, len(rnasite_residues))
    for r in rnasite_residues[:p]:
        hist = process.get_smooth_history(r)
        plt.plot(hist['time'], hist['value'])
    plt.legend(loc='best', labels=rnasite_residues[:p])
    plt.show()

