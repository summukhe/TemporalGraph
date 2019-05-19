import os
import sys
import logging
import argparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from CLIUtility import *
import temporal_graph

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

    parser.add_argument('--potential', dest='potential',
                        help="string suggesting potential to be used for the analysis ('mj', 'energy', 'charmm')",
                        required=False,
                        metavar="STRING", type=str)

    parser.add_argument('--debug', dest='debug',
                        help="detailed logging",
                        action='store_true')

    parser.set_defaults(debug=False)
    parser.set_defaults(chain="A")
    parser.set_defaults(potential='charmm')
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    logger = logging.getLogger(name="MAIN")
    pdbs = temporal_graph.read_pdb(args.pdb)
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

    g = temporal_graph.contact_graph(structure,
                                     cutoff=12,
                                     potential=args.potential)

    print('%d x %d' % (g.order, g.size))
    path_filter = temporal_graph.GeometricPathFilter(forward=True,
                                                     distance_cutoff=6)
    g_it = temporal_graph.AllPathIterator(g, path_filter)
    for s in site1:
        paths = g_it.all_path(s, stop_vertex=site2, min_path_length=4)
        print(paths)
