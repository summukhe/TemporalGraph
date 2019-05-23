import os
import sys
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

    parser.add_argument('--method', dest='method',
                        help="string suggesting type of analysis to be performed ('centrality', 'mincut')",
                        required=False,
                        metavar="STRING", type=str)

    parser.add_argument('--min-energy', dest='min_energy',
                        help="minimum energy in the contact trail (default: 1.5)",
                        required=False,
                        metavar="STRING", type=float)

    parser.add_argument('--cutoff', dest='cutoff',
                        help="integer value defining the radius to "
                             "use for contact cutoff (default: 12)",
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

    parser.add_argument('--ntrails', dest='ntrail',
                        required=False,
                        help="Maximum number of path to be used for each site for the "
                             "centrality analysis (default: 25)",
                        type=int)

    parser.add_argument('--path-overlap', dest='overlap',
                        required=False,
                        help="Allowed overlap distance between paths originated from "
                             "the same vertex (default: 0.5)",
                        type=float)

    parser.add_argument('--backtrack', dest='backtrack',
                        required=False,
                        help="Program hack for fast stack unrolling to find non overlapping "
                             "trails. (Default: 0.5)",
                        type=float)

    parser.set_defaults(debug=False)
    parser.set_defaults(flow_forward=True)
    parser.set_defaults(chain="A")
    parser.set_defaults(method='centrality')
    parser.set_defaults(min_energy=2.5)
    parser.set_defaults(cutoff=12)
    parser.set_defaults(min_path_length=3)
    parser.set_defaults(max_path_length=10)
    parser.set_defaults(ntrail=25)
    parser.set_defaults(overlap=0.5)
    parser.set_defaults(backtrack=0.5)

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
    curr_residue, result = temporal_graph.all_path_mutation_evaluation(pdb_structure=structure,
                                                                       residue_id=residue_id,
                                                                       site1=site1,
                                                                       site2=site2,
                                                                       method=args.method,
                                                                       flow_forward=args.flow_forward,
                                                                       minimum_energy=args.min_energy,
                                                                       maximum_distance=args.cutoff,
                                                                       ntrails=args.ntrail,
                                                                       min_path_length=args.min_path_length,
                                                                       max_path_length=args.max_path_length,
                                                                       allowed_path_overlap=args.overlap,
                                                                       backtrack_fraction=args.backtrack)

    fscore = result.normalized_score

    for aa, value in fscore.items():
        print("%1s %1s %10.5f" % (result.ref_amino, aa, value))
