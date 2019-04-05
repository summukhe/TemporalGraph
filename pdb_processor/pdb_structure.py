"""
This program file defines the basic class object
for storing pdb 3D structure. The basic class structure
only allows to read and access the coordinates and related
information in a systematic way.
"""

import os
from amino_acids import get_amino


__author__  = "Sumanta Mukherjee"
__version__ = "1.0"

__all__ = ['CaTrace']


class CaTrace():
    """
    This is a basic pdb structure class. This class reads only
    CA trace atom trace from a given pdb file. This class relies
    on the ATOM tags of the record. In case of modified residues
    with HETATM tag, this class may results in missing residues
    in the trace.
    """
    def __init__(self, name, chainId, entry):
        self.__name = name
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
        resId = int(resid)
        if resId not in self.__structure:
            raise KeyError('Invalid residue id: %d' % int(resid))
        return (self.__structure[resId]['x'],
                self.__structure[resId]['y'],
                self.__structure[resId]['z'])

    def get_amino(self, resid):
        resId = int(resid)
        if resId not in self.__residues:
            raise KeyError('Invalid residue id: %d' % int(resid))
        return self.__residues[resId]

    def __len__(self):
        return len(self.__structure)

    def __str__(self):
        s = ''
        for aa in self.sequence():
            s = s + aa.name(one_letter_code=True)
        return s



