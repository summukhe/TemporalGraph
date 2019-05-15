"""
This code file includes the base code for reading CA coordinate
trace from trajectory pdb files.
"""

import os
from temporal_graph.pdb_processor import CaTrace, PDBStructure, pdb_to_catrace

__version__ = "1.0"

__all__ = ['read_trajectory_catrace', 'read_pdb']


def read_trajectory_catrace(pdbfile):
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
            if line.startswith('ENDMDL') and len(structures) > 0:
                pair = dict()
                for chain, seq in structures.items():
                    name = modelname if modelname is not None else pdbname
                    pair[chain] = CaTrace(name, seq, chain)
                trajectory.append(pair)
                structures.clear()
    if len(structures) > 0:
        pair = dict()
        for chain, seq in structures.items():
            name = modelname if modelname is not None else pdbname
            pair[chain] = CaTrace(name, seq, chain)
        trajectory.append(pair)
        structures.clear()
    return trajectory


def read_pdb(pdbfile):
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
                    resName = line[17:20].strip()
                    atom_id = int(line[6:11].strip())
                    atom_name=line[12:16].strip()
                    bfactor = float(line[60:66].strip())
                    if atom_name == "OXT":
                        atom_name = "O"
                    structures[chain].append({'x': x,
                                              'y': y,
                                              'z': z,
                                              'resid': resId,
                                              'resname': resName,
                                              'atomname': atom_name,
                                              'atomid': atom_id,
                                              'bfactor': bfactor})
            if (line.startswith('TER') or line.startswith('END')) and len(structures) > 0:
                pair = dict()
                for chain, seq in structures.items():
                    name = modelname if modelname is not None else pdbname
                    pair[chain] = PDBStructure(name, seq, chain)
                trajectory.append(pair)
                structures.clear()
    if len(structures) > 0:
        pair = dict()
        for chain, seq in structures.items():
            name = modelname if modelname is not None else pdbname
            pair[chain] = PDBStructure(name, seq, chain)
        trajectory.append(pair)
        structures.clear()
    return trajectory


if __name__ == "__main__":
    src_file_path = os.path.dirname(os.path.abspath(__file__))
    pdbfile = os.path.join(src_file_path, "..", "..", "data", "example.pdb")
    trajectory = read_pdb(pdbfile)
    print('Number of snapshots: %d' % len(trajectory))
    protein_rep1 = '%s' % trajectory[0]['A']
    protein_rep2 = '%s' % pdb_to_catrace(trajectory[0]['A'])
    assert protein_rep1 == protein_rep2
