import os
import json
import numpy as np
from temporal_graph import *
from matplotlib import pyplot as plt


if __name__ == "__main__":
    pdb_file= "data/example.pdb"
    atp_site_file= "data/site1.txt"
    rna_site_file= "data/site2.txt"
    rtf_file = "data/charmm_top_prot.rtf"
    prm_file = "data/nonbonded_prot.prm"
    out_pdb_file = "scratch/marked.pdb"

    assert os.path.isfile(pdb_file)
    trajectory = read_pdb(pdb_file)
    idx = 0
    residue_ids = trajectory[idx]['A'].residue_ids()

    print('Running for snapshot: %d' % idx)

    assert os.path.isfile(atp_site_file)
    with open(atp_site_file, "r") as f:
        atpsite_residues = [int(r) for r in f.readlines()]
    atp_site_residues = [trajectory[idx]['A'].key(r) for r in atpsite_residues if r in residue_ids]

    assert os.path.isfile(rna_site_file)
    with open(rna_site_file, "r") as f:
        rnasite_residues = [int(r) for r in f.readlines()]
    rna_site_residues = [trajectory[idx]['A'].key(r) for r in rnasite_residues if r in residue_ids]

    ec, cuts = between_site_residues_by_mincut(trajectory[0]['A'],
                                               site1=atp_site_residues,
                                               site2=rna_site_residues)
    print("Total valid edges: %d" % ec)
    print(cuts)

    g = contact_energy_graph(trajectory[0]['A'], contact_radius=8)
    pulse = PeriodicPulse(ts=0, span=5, periods=25, strength=10)
    n = np.minimum(1000, g.order)
    process = DiffusionOnStaticGraph(dissipation_rate=0.1)
    process.start(g, atp_site_residues, pulse)
    process.n_steps(n=n)

    p = np.minimum(10, len(rna_site_residues))
    for r in rna_site_residues[:p]:
        hist = process.get_history(r)
        plt.plot(hist['time'], hist['value'])
    plt.legend(loc='best', labels=rna_site_residues[:p])
    plt.show()

    diffusion_trj = mark_diffusion_on_structure(ca_trace=pdb_to_catrace(trajectory[idx]['A']),
                                                pulse=pulse,
                                                signal_nodes=atp_site_residues,
                                                dissipation_rate=0.1)
    with open(out_pdb_file, 'w') as fh:
        for i, trj in enumerate(diffusion_trj):
            fh.write('MODEL     %4d\n' % (i+1))
            trj.write(fh)
            fh.write('ENDMDL\n')






