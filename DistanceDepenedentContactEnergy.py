import os
import sys
import logging
import argparse
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from CLIUtility import *
from temporal_graph import *


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch distance dependent contact energy "
                                                 "in residue residue interaction structure graph")
    parser.add_argument('--pdb', dest='pdb',
                        help="pdb file missing residues must be resolved.",
                        required=True,
                        metavar='FILE', type=lambda x: is_valid_file(parser, x))

    parser.add_argument('--chain', dest='chain',
                        help="pdb chain to consider",
                        required=False,
                        metavar='STRING', type=str)

    parser.add_argument('--cutoff', dest='cutoff',
                        help="integer value defining the radius to use for contact cutoff (default: 12)",
                        required=False,
                        metavar='INTEGER', type=int)

    parser.add_argument('--debug', dest='debug',
                        help="detailed logging",
                        action='store_true')

    parser.set_defaults(debug=False)
    parser.set_defaults(chain="A")
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
    g = contact_energy_graph(structure, contact_radius=args.cutoff)
    for u, v in g.edges:
        u_resid = int(u[3:])
        v_resid = int(v[3:])
        d = distance(Coordinate3d(*structure.xyz(u_resid, "CA")),
                     Coordinate3d(*structure.xyz(v_resid, "CA")))
        print("%s %s %.4f %.3f" % (u, v, g.weight(u, v), d))

