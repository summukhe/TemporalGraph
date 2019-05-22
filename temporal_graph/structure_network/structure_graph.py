import numpy as np
from copy import deepcopy
from temporal_graph.network_analysis import GeometricGraph3d
from temporal_graph.spatial_ds import *
from temporal_graph.pdb_processor import *
from temporal_graph.force_field import *


__version__ = "1.0"
__all__ = ['contact_graph',
           'potential_contact_graph',
           'contact_energy_graph']


def contact_graph(pdb_structure,
                  cutoff=12,
                  potential='charmm',
                  weight_normalizer=FFNormalizer()):
    assert isinstance(pdb_structure, PDBStructure) or \
           isinstance(pdb_structure, CaTrace)
    if isinstance(cutoff, DistanceCutoff):
        cutoff = cutoff.cutoff
    assert cutoff > 0
    assert potential in ['energy', 'charmm', 'mj']
    if potential == 'energy':
        if isinstance(pdb_structure, CaTrace):
            return potential_contact_graph(pdb_structure,
                                           cutoff=DistanceCutoff(def_cutoff=cutoff),
                                           potential='charmm')
        else:
            g = contact_energy_graph(pdb_structure,
                                     contact_radius=cutoff,
                                     energy_score=weight_normalizer)
    else:
        if isinstance(pdb_structure, PDBStructure):
            structure = pdb_to_catrace(pdb_structure)
        else:
            structure = deepcopy(pdb_structure)
        g = potential_contact_graph(structure,
                                    cutoff=DistanceCutoff(def_cutoff=cutoff),
                                    potential=potential)
    return g


def potential_contact_graph(ca_trace, cutoff=DistanceCutoff(), potential='mj'):
    assert isinstance(ca_trace, CaTrace)
    assert isinstance(cutoff, DistanceCutoff)
    assert potential in ['mj', 'charmm']
    res_ids = ca_trace.residue_ids
    c_graph = GeometricGraph3d(directed=False)
    for r in res_ids:
        amino_key = ca_trace.key(r)
        amino_crd = Coordinate3d(*ca_trace.xyz(r))
        c_graph.add_vertex(amino_key, attribute=amino_crd)

    for ri in res_ids:
        amino_i = ca_trace.get_amino(ri)
        x_i, y_i, z_i = ca_trace.xyz(ri)
        for rj in res_ids:
            if ri < rj:
                amino_j = ca_trace.get_amino(rj)
                x_j, y_j, z_j = ca_trace.xyz(rj)
                c = cutoff(amino_i, amino_j)
                d = np.sqrt((x_i-x_j)**2 + (y_i-y_j)**2 + (z_i-z_j)**2)
                if d <= c:
                    p = get_pair_potential(amino_i, amino_j, d, pot_type=potential)
                    c_graph.add_edge('%s%d' % (amino_i, ri),
                                     '%s%d' % (amino_j, rj),
                                     weight=p)
    return c_graph


def contact_energy_graph(pdb_struct,
                         contact_radius=12,
                         epsilon=1.,
                         elec_only=False,
                         summed=True,
                         energy_score=FFNormalizer()):
    assert isinstance(pdb_struct, PDBStructure)
    assert isinstance(energy_score, FFNormalizer)
    ca_trace = pdb_to_catrace(pdb_struct)
    residues = ca_trace.residue_ids
    x_lst, y_lst, z_lst = [], [], []
    for r in residues:
        x, y, z = ca_trace.xyz(r)
        x_lst.append(x)
        y_lst.append(y)
        z_lst.append(z)
    grid = Grid3D(max_coord=Coordinate3d(np.max(x_lst), np.max(y_lst), np.max(z_lst)),
                  min_coord=Coordinate3d(np.min(x_lst), np.min(y_lst), np.min(z_lst)),
                  spacing=2)
    for r in residues:
        grid.register_obj(r, Coordinate3d(*ca_trace.xyz(r)))
    neighbors = dict()
    for r1 in residues:
        neighbors[r1] = {r2: 0 for r2 in grid.neighbors(r1, contact_radius)}
    ff = FFManager()
    for r1 in neighbors:
        residue_name1 = pdb_struct.residue_name(r1)
        atom_names1 = pdb_struct.atom_names(r1)
        for r2 in neighbors[r1]:
            residue_name2 = pdb_struct.residue_name(r2)
            atom_names2 = pdb_struct.atom_names(r2)
            for atom1 in atom_names1:
                for atom2 in atom_names2:
                    d = distance(Coordinate3d(*pdb_struct.xyz(r1, atom1)),
                                 Coordinate3d(*pdb_struct.xyz(r2, atom2)))
                    neighbors[r1][r2] += ff.energy(residue_name1,
                                                   atom1,
                                                   residue_name2,
                                                   atom2,
                                                   distance=d,
                                                   epsilon=epsilon,
                                                   elec_only=elec_only,
                                                   summed=summed)
    c_graph = GeometricGraph3d(directed=False)
    for r in residues:
        c_graph.add_vertex(pdb_struct.key(r),
                           attribute=Coordinate3d(*pdb_struct.xyz(r, 'CA')))

    for r1 in neighbors:
        for r2 in neighbors[r1]:
            c_graph.add_edge(pdb_struct.key(r1),
                             pdb_struct.key(r2),
                             weight=energy_score(neighbors[r1][r2]))
    return c_graph

