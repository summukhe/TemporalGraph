"""
This program file defines the basic class object
for storing pdb 3D structure. The basic class structure
only allows to read and access the coordinates and related
information in a systematic way.
"""

import os
from .amino_acids import get_amino


__author__  = "Sumanta Mukherjee"
__version__ = "1.0"

__all__ = ['CaTrace', 'do_mutation']


class CaTrace:
    """
    This is a basic pdb structure class. This class reads only
    CA trace atom trace from a given pdb file. This class relies
    on the ATOM tags of the record. In case of modified residues
    with HETATM tag, this class may results in missing residues
    in the trace.
    """
    def __init__(self, name, chainId, entry):
        self.__name = name
        self.__bfactor = dict()
        self.__chain_id = chainId
        self.__structure = dict()
        self.__residues = dict()

        assert isinstance(entry, list)
        keys = ['resid', 'resname', 'x', 'y', 'z']
        for item in entry:
            assert isinstance(item, dict)
            for k in keys:
                assert k in item
            self.__structure[item['resid']] = {'x': item['x'],
                                               'y': item['y'],
                                               'z': item['z']}
            self.__residues[item['resid']] = get_amino(item['resname'])
            self.__bfactor[item['resid']] = item['bfactor'] if 'bfactor' in item else 100

    def size(self):
        return len(self.__structure)

    def sequence(self):
        s = list()
        for k in sorted([int(r) for r in self.__residues.keys()]):
            s.append(self.__residues[k])
        return s

    def residue_ids(self):
        return sorted([int(k) for k in self.__structure.keys()])

    def xyz(self, resid):
        residue_id = int(resid)
        if residue_id not in self.__structure:
            raise KeyError('Invalid residue id: %d' % residue_id)
        return self.__structure[residue_id]['x'], self.__structure[residue_id]['y'], self.__structure[residue_id]['z']

    def b_factor(self, resid, value=None):
        residue_id = int(resid)
        if residue_id not in self.__bfactor:
            raise KeyError('Invalid residue id: %d' % residue_id)
        if value is None:
            return self.__bfactor[residue_id]
        else:
            self.__bfactor[residue_id] = value
            return value

    def get_amino(self, resid):
        residue_id = int(resid)
        if residue_id not in self.__residues:
            raise KeyError('Invalid residue id: %d' % int(resid))
        return self.__residues[residue_id]

    def set_amino(self, resid, aa_type='A'):
        residue_id = int(resid)
        if residue_id not in self.__residues:
            return False
        self.__residues[residue_id] = get_amino(aa_type)
        return True

    def write(self, fh):
        if not hasattr(fh, 'write'):
            raise Exception('Invalid object type [expects file object]')
        res_ids = self.residue_ids()
        for i, r in enumerate(res_ids):
            line = "ATOM  %5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f " % (i,
                                                                                 'CA',
                                                                                 self.__residues[r],
                                                                                 'A',
                                                                                 int(r),
                                                                                 float(self.__structure[r]['x']),
                                                                                 float(self.__structure[r]['y']),
                                                                                 float(self.__structure[r]['z']),
                                                                                 1.0,
                                                                                 float(self.__bfactor[r]))
            fh.write("%s\n" % line)

    def __len__(self):
        return len(self.__structure)

    def __str__(self):
        s = ''
        for aa in self.sequence():
            s = s + aa.name(one_letter_code=True)
        return s


def do_mutation(ca_trace, res_id, to_residue):
    assert isinstance(ca_trace, CaTrace)
    residue_ids = ca_trace.residue_ids()
    assert res_id in residue_ids
    ca_trace.set_amino(res_id, to_residue)
    return ca_trace

