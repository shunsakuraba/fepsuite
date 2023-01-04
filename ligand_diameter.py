import argparse
import itertools

import mdtraj
import numpy

import common_gmx_files

def main(args: argparse.Namespace):
    ndx = common_gmx_files.parse_index(args.index)
    ligand = ndx[args.ligand_mol]

    trajfile = args.trajectory
    if trajfile is None:
        trajfile = args.structure
    distpair = list(itertools.product(ligand, ligand))
    diameters = numpy.zeros((0,))
    for chunk in mdtraj.iterload(trajfile, top=args.structure):
        distances = mdtraj.compute_distances(chunk, distpair) # sized [fr, pair]
        diameter_per_frame = numpy.amax(distances, axis=1)
        diameters = numpy.concatenate((diameters, diameter_per_frame))
    
    if len(diameters) == 1:
        safeval = diameters[0] * 1.5
    else:
        maxd = numpy.amax(diameters)
        safeval = max(maxd + numpy.std(diameters) * 6, maxd * 1.4)
    
    print("avg", numpy.mean(diameters))
    print("max", numpy.amax(diameters))
    print("safe", safeval)

def init_args():
    parser = argparse.ArgumentParser(description="Get ligand diameter (necessary for determining enough rlist)",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--structure", type=str, required=True, help="Structure")
    parser.add_argument("--trajectory", type=str, default=None, help="Trajectory")
    parser.add_argument("--index", type=str, required=True, help="Gromacs index file name")
    parser.add_argument("--ligand-mol", type=str, default="Ligand", help="Index name for ligand selection")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()

    main(args)
    
