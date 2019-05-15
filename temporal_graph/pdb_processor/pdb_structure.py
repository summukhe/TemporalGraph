"""
This program file defines the basic class object
for storing pdb 3D structure. The basic class structure
only allows to read and access the coordinates and related
information in a systematic way.
"""

import logging
from .amino_acids import get_amino
from temporal_graph.force_field import FFManager


__version__ = "1.0"

__all__ = ['CaTrace', 'do_mutation', 'PDBStructure', 'pdb_to_catrace']


class CaTrace:
    """
    This is a basic pdb structure class. This class reads only
    CA trace atom trace from a given pdb file. This class relies
    on the ATOM tags of the record. In case of modified residues
    with HETATM tag, this class may results in missing residues
    in the trace.
    """
    def __init__(self, name, entry, chainId='A'):
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

    @property
    def name(self):
        return self.__name

    @property
    def chain(self):
        return self.__chain_id

    @property
    def size(self):
        return len(self.__structure)

    @property
    def sequence(self):
        s = list()
        for k in sorted([int(r) for r in self.__residues.keys()]):
            s.append(self.__residues[k])
        return s

    @property
    def residue_ids(self):
        return sorted([int(k) for k in self.__structure.keys()])

    def key(self, residue_id):
        if residue_id in self.__residues:
            return '%s%d' % (self.__residues[residue_id], residue_id)
        return ''

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
        res_ids = self.residue_ids
        for i, r in enumerate(res_ids):
            line = "ATOM  %5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f " % (i,
                                                                                 'CA',
                                                                                 self.__residues[r],
                                                                                 self.__chain_id,
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
        for aa in self.sequence:
            s = s + aa.name(one_letter_code=True)
        return s


def do_mutation(ca_trace, res_id, to_residue):
    assert isinstance(ca_trace, CaTrace)
    residue_ids = ca_trace.residue_ids
    assert res_id in residue_ids
    ca_trace.set_amino(res_id, to_residue)
    return ca_trace


class PDBStructure:
    def __init__(self, name, entry, chainId='A'):
        self.__name = name
        self.__chain_id = chainId
        self.__structure = dict()
        self.__residues = dict()
        self.__bfactor = dict()

        logger = logging.getLogger(name='pdb_processor.PDBStructure')
        assert isinstance(entry, list)
        ff = FFManager()
        keys = ['resid', 'resname', 'atomname', 'atomid', 'x', 'y', 'z']
        for item in entry:
            assert isinstance(item, dict)
            for k in keys:
                assert k in item
            residue_id, residue_name, atom_name, atom_id = item['resid'], item['resname'], item['atomname'], item['atomid']
            if residue_id not in self.__structure:
                self.__structure[residue_id] = dict()
            if not ff.is_valid_atom(residue_name, atom_name):
                logger.error('Error validating residue, atom pair (%s, %s)' % (residue_name, atom_name))
                raise Exception('Error: validating residue, atom pair (%s, %s)' % (residue_name, atom_name))
            self.__structure[residue_id][atom_name] = {'x': item['x'],
                                                       'y': item['y'],
                                                       'z': item['z'],
                                                       'id': atom_id}
            self.__residues[residue_id] = get_amino(residue_name)
            if residue_id not in self.__bfactor:
                self.__bfactor[residue_id] = {}
            if 'bfactor' in item:
                self.__bfactor[residue_id][atom_name] = item['bfactor'] if 'bfactor' in item else 100

    @property
    def name(self):
        return self.__name

    @property
    def chain(self):
        return self.__chain_id

    @property
    def size(self):
        return len(self.__structure)

    @property
    def sequence(self):
        s = list()
        for k in sorted([int(r) for r in self.__residues.keys()]):
            s.append(self.__residues[k])
        return s

    @property
    def residue_ids(self):
        return sorted([int(k) for k in self.__structure.keys()])

    def find_residue(self, aa_type):
        return [rid for rid in self.__residues.keys() if self.__residues[rid] == aa_type]

    def atom_list(self, residue_id):
        if residue_id in self.__structure:
            return [self.__structure[residue_id][atm] for atm in self.atom_names(residue_id)]
        return []

    def atom_names(self, residue_id):
        if residue_id in self.__structure:
            return [ atm for atm in sorted(self.__structure[residue_id].keys(), key=lambda x: self.__structure[residue_id][x]['id'])]
        return []

    def key(self, residue_id):
        if residue_id in self.__residues:
            return '%s%d' % (self.__residues[residue_id], residue_id)
        return ''

    def write(self, fh):
        if not hasattr(fh, 'write'):
            raise Exception('Invalid object type [expects file object]')
        res_ids = self.residue_ids
        for i, r in enumerate(res_ids):
            for j, aa in enumerate(self.atom_list(r)):
                line = "ATOM  %5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f " % (self.__structure[r][aa]['id'],
                                                                                     aa,
                                                                                     self.__residues[r],
                                                                                     'A',
                                                                                     int(r),
                                                                                     float(self.__structure[r][aa]['x']),
                                                                                     float(self.__structure[r][aa]['y']),
                                                                                     float(self.__structure[r][aa]['z']),
                                                                                     1.0,
                                                                                     float(self.__bfactor[r][aa]))
            fh.write("%s\n" % line)

    def __len__(self):
        return len(self.__structure)

    def __str__(self):
        s = ''
        for aa in self.sequence:
            s = s + aa.name(one_letter_code=True)
        return s

    def residue_name(self, residue_id, one_letter_code=False):
        if residue_id in self.__residues:
            return self.__residues[residue_id].name(one_letter_code=one_letter_code)
        return None

    def xyz(self, residue_id, atom_name):
        if (residue_id in self.__structure) and (atom_name in self.__structure[residue_id]):
            atom_details = self.__structure[residue_id][atom_name]
            return atom_details['x'], atom_details['y'], atom_details['z']
        return None, None, None

    def bfactor(self, residue_id, atom_name):
        if (residue_id in self.__bfactor) and (atom_name in self.__bfactor[residue_id]):
            return self.__bfactor[residue_id][atom_name]
        return None


def pdb_to_catrace(pdb_struct):
    assert isinstance(pdb_struct, PDBStructure)
    residue_ids = pdb_struct.residue_ids
    entries = list()
    logger = logging.getLogger(name='pdb_processor.pdb_to_catrace')
    for r in residue_ids:
        resname = pdb_struct.residue_name(r)
        x, y, z = pdb_struct.xyz(r, 'CA')
        bfactor = pdb_struct.bfactor(r, 'CA')
        if (x is None) or (y is None) or (z is None):
            logger.error('Missing CA in residue (%d)' % r)
            raise Exception('Missing CA in the protein for residue (%d)' % r)
        entries.append({'resname': resname,
                        'resid': r,
                        'x': x,
                        'y': y,
                        'z': z,
                        'bfactor':bfactor})
    return CaTrace(pdb_struct.name, entries, pdb_struct.chain)
