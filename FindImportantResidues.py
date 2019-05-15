import os
import sys
import logging
import argparse
from CLIUtility import *
from temporal_graph import *

sys.path.append( os.path.dirname(os.path.abspath(__file__)))


def read_site_residue(site_file, structure):
    assert os.path.isfile(site_file)
    assert isinstance(structure, PDBStructure) or is_valid_file(structure, CaTrace)
    with open(site_file, "r") as f:
        site_residues = [int(r) for r in f.readlines()]
    residue_ids = structure.residue_ids
    site_reskey = [structure.key(r) for r in site_residues if r in residue_ids]
    return site_reskey


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mark important residues between sites in a pdb structure")
    parser.add_argument('--pdb', dest='pdb',
                        help="pdb file missing residues must be resolved.",
                        required=True,
                        metavar='FILE', type=lambda x: is_valid_file(parser, x))

    parser.add_argument('--chain', dest='chain',
                        help="pdb chain to consider",
                        required=False,
                        metavar='STRING', type=str)

    parser.add_argument('--site1', dest='site1',
                        help="text file containing list of residue ids @ site1",
                        required=True,
                        metavar='FILE', type=lambda x: is_valid_file(parser, x))

    parser.add_argument('--site2', dest='site2',
                        help="text file containing list of residue ids @ site2",
                        required=True,
                        metavar='FILE', type=lambda x: is_valid_file(parser, x))

    parser.add_argument('--method', dest='method',
                        help="string suggesting type of analysis to be performed ('centrality', 'mincut')",
                        required=False,
                        metavar="STRING", type=str)

    parser.add_argument('--potential', dest='potential',
                        help="string suggesting potential to be used for the analysis ('mj', 'energy')",
                        required=False,
                        metavar="STRING", type=str)

    parser.add_argument('--cutoff', dest='cutoff',
                        help="integer value defining the radius to use for contact cutoff (default: 12)",
                        required=False,
                        metavar='INTEGER', type=int)

    parser.add_argument('--debug', dest='debug',
                        help="detailed logging",
                        action='store_true')

    parser.set_defaults(debug=False)
    parser.set_defaults(chain="A")
    parser.set_defaults(method='centrality')
    parser.set_defaults(potential='mj')
    parser.set_defaults(cutoff=12)
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    logger = logging.getLogger(name="MAIN")
    pdbs = read_pdb(args.pdb)

    if len(pdbs) == 0:
        logger.error('No pdb found!!')
        exit(1)

    if len(pdbs) == 1:
        logger.debug('One instance read!!')

    instance = pdbs[0]

    assert isinstance(instance, dict)

    if args.chain not in instance:
        logger.error("Can not find chain (%s)" % args.chain)
        exit(1)

    structure = instance[args.chain]
    site1 = read_site_residue(args.site1, structure)
    site2 = read_site_residue(args.site2, structure)

    if args.method not in ['mincut', 'centrality']:
        logger.error('Invalid method (%s)' % args.method)
        exit(1)

    if args.potential not in ['mj', 'energy']:
        logger.error('Invalid potential (%s)' % args.potential)
        exit(1)

    if args.method == 'mincut':
        __, residues = between_site_residues_by_mincut(pdb_structure=structure,
                                                       site1=site1,
                                                       site2=site2,
                                                       type=args.potential,
                                                       contact_radius=args.cutoff)
    elif args.method == "centrality":
        residues = between_site_residues_by_stpath(pdb_structure=structure,
                                                   site1=site1,
                                                   site2=site2,
                                                   type=args.potential,
                                                   contact_radius=args.cutoff)
    ordered_residues = list(reversed(sorted([(k, v) for k, v in residues.items()], key=lambda x: x[1])))
    for k, v in ordered_residues:
        print('%s %.3f' % (k,v))
