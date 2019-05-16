import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import temporal_graph

__version__ = "1.0"
__all__ = ['is_valid_file', 'is_integer', 'read_site_residue']


def is_valid_file(parser, x):
    if not os.path.isfile(x):
        parser.error("File (%s) does not exists" % x)
    return x


def is_integer(parser, x):
    if isinstance(x, np.int) or int(x) != 0:
        return int(x)
    parser.error('Not an integer value (%s)' % x)


def read_site_residue(site_file, structure, key=True):
    assert os.path.isfile(site_file)
    assert isinstance(structure, temporal_graph.PDBStructure) or is_valid_file(structure, temporal_graph.CaTrace)
    with open(site_file, "r") as f:
        site_residues = [int(r) for r in f.readlines()]
    residue_ids = structure.residue_ids
    if key:
        site_reskey = [structure.key(r) for r in site_residues if r in residue_ids]
    else:
        site_reskey = [r for r in site_residues if r in residue_ids]
    return site_reskey

