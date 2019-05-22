import os
import sys
import json
import logging
import argparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from CLIUtility import *
import temporal_graph

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rank mutants by it effect in network " +
                                                 "communication between sites in a pdb structure")
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

    parser.add_argument('--residue', dest='residue',
                        help='residue id on which the mutation selection will be performed!',
                        required=True,
                        metavar='INTEGER', type=int)

    parser.add_argument('--amino', dest='amino',
                        help='mutated residue type',
                        required=True,
                        metavar='AMINO', type=str)

    parser.add_argument('--method', dest='method',
                        help="string suggesting type of analysis to be performed ('centrality', 'mincut')",
                        required=False,
                        metavar="STRING", type=str)

    parser.add_argument('--min-energy', dest='min_energy',
                        help="minimum energy in the contact trail (default: 1.5)",
                        required=False,
                        metavar="FLOAT", type=float)

    parser.add_argument('--cutoff', dest='cutoff',
                        help="integer value defining the radius to use for contact cutoff (default: 12)",
                        required=False,
                        metavar='INTEGER', type=int)

    parser.add_argument('--debug', dest='debug',
                        help="detailed logging",
                        action='store_true')

    parser.add_argument('--backward-flow', dest='flow_forward',
                        help="allow geometric retrace of trail",
                        action='store_false')

    parser.add_argument('--min-trail-length', dest='min_path_length',
                        required=False,
                        help="Minimum length of the trail to be considered in the analysis (default: 3)",
                        type=int)

    parser.add_argument('--max-trail-length', dest='max_path_length',
                        required=False,
                        help="Maximum length of the trail to be considered "
                             "in the analysis (default: 10)",
                        type=int)

    parser.set_defaults(debug=False)
    parser.set_defaults(flow_forward=True)
    parser.set_defaults(chain="A")
    parser.set_defaults(method='centrality')
    parser.set_defaults(min_energy=1.5)
    parser.set_defaults(cutoff=12)
    parser.set_defaults(min_path_length=3)
    parser.set_defaults(max_path_length=10)

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
    site1 = read_site_residue(args.site1, structure, key=False)
    site2 = read_site_residue(args.site2, structure, key=False)
    residue_id = args.residue
    assert residue_id in structure.residue_ids

    if args.method not in ['mincut', 'centrality']:
        logger.error('Invalid method (%s)' % args.method)
        exit(1)

    if residue_id in site1 + site2:
        logger.warning('Residue selected should not be part of the site')

    logger.debug("Starting mutation analysis!!")

    ca_backbone = temporal_graph.pdb_to_catrace(structure)
    ca_backbone.set_amino(args.residue, args.amino)
    g = temporal_graph.contact_graph(ca_backbone, cutoff=args.cutoff, potential='charmm')

    site1_key = [ca_backbone.key(r) for r in site1]
    site2_key = [ca_backbone.key(r) for r in site2]

    __, scores = temporal_graph.all_path_centrality(g,
                                                    source=site1_key,
                                                    target=site2_key,
                                                    forward_path=args.flow_forward,
                                                    maximum_distance=args.cutoff,
                                                    minimum_weight=args.min_energy,
                                                    min_path_length=args.min_path_length,
                                                    max_path_length=args.max_path_length)

    print(json.dumps(scores, indent=2))
