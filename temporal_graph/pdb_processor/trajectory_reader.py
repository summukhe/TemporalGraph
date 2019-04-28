"""
This code file includes the base code for reading CA coordinate
trace from trajectory pdb files.
"""

import os
from .pdb_structure import CaTrace

__author__  = "Sumanta Mukherjee"
__version__ = "1.0"

__all__ = ['read_catrace_from_trajectory']


def read_catrace_from_trajectory(pdbfile):
    assert os.path.isfile(pdbfile)
    chain, modelname = None, None
    pdbname = os.path.basename(pdbfile).split('.')[0]
    structures = dict()
    trajectory = list()

    with open(pdbfile, "r") as f:
        for line in f.readlines():
            if line.startswith('MODEL'):
                assert len(line.split()) == 2
                modelname = 'Model_%d' % int(line.split()[1])
            if line.startswith('ATOM'):
                if line[13:15] == "CA":
                    if chain is None:
                        chain = line[21]

                    if chain != line[21]:
                        chain = line[21].strip()
                    chain = 'A' if chain == '' else chain
                    if chain not in structures:
                        structures[chain] = list()
                    alt = line[16]
                    if alt in [' ', 'A']:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        resId = int(line[22:26].strip())
                        resName = line[17:20]

                        structures[chain].append({'x': x,
                                                  'y': y,
                                                  'z': z,
                                                  'resid': resId,
                                                  'resname': resName})
            if (line.startswith('TER') or line.startswith('END')) and len(structures) > 0:
                pair = dict()
                for chain, seq in structures.items():
                    name = modelname if modelname is not None else pdbname
                    pair[chain] = CaTrace(name, chain, seq)
                trajectory.append(pair)
                structures.clear()
    if len(structures) > 0:
        pair = dict()
        for chain, seq in structures.items():
            name = modelname if modelname is not None else pdbname
            pair[chain] = CaTrace(name, chain, seq)
        trajectory.append(pair)
        structures.clear()
    return trajectory


"""
if __name__ == "__main__":
    pdbfile = '/home/sumanta/Project/structural_dynamics/second_round_analysis/dynamics/electrostatics_analysis/1A1V/mode_catraj_10.pdb'
    trajectory = read_catrace_from_trajectory(pdbfile)
    print( 'Number of snapshots: %d' % len(trajectory))
    for i in range(len(trajectory)):
        print('%s' % trajectory[i]['A'])
"""
