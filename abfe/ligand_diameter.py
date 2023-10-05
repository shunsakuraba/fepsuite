import argparse
import itertools

import mdtraj
import numpy

import common_gmx_files

def main(args: argparse.Namespace):
    trajfile = args.trajectory
    if trajfile is None:
        trajfile = args.structure
    if args.index:
        ndx = common_gmx_files.parse_index(args.index)
        ligand = ndx[args.ligand_mol]  
    else:
        # all atoms in the system
        ligandtraj = mdtraj.load(args.structure)
        ligand = list(range(ligandtraj.n_atoms))

    distpair = list(itertools.product(ligand, ligand))
    diameters = numpy.zeros((0,))
    for chunk in mdtraj.iterload(trajfile, top=args.structure):
        distances = mdtraj.compute_distances(chunk, distpair) # sized [fr, pair]
        diameter_per_frame = numpy.amax(distances, axis=1)
        diameters = numpy.concatenate((diameters, diameter_per_frame))
    
    if len(diameters) == 1:
        safeval = diameters[0] * 1.5
    else:
        # Currently we are using regression assuming the tail behave similar to the normal function.
        obs_points = numpy.array([0.90, 0.93, 0.95, 0.97, 0.99])
        invsf = [1.2815515655446004, 1.4757910281791706, 1.6448536269514729, 1.8807936081512511, 2.3263478740408408] # scipy.stats.norm.isf(1-obs_points), hardcoded to remove scipy package dependence
        target_invsf = 3.7190164854556804 # isf(1e-4)
        quantiles = numpy.quantile(diameters, obs_points)
        # solve least square to get the distance.
        # A = [log1m 1]
        A = numpy.vstack([invsf, numpy.ones_like(invsf)]).T
        slope, intercept = numpy.linalg.lstsq(A, quantiles, rcond=None)[0]
        safeval = slope * target_invsf + intercept

    print("avg", numpy.mean(diameters))
    print("max", numpy.amax(diameters))
    print("safe", safeval)

def init_args():
    parser = argparse.ArgumentParser(description="Get ligand diameter (necessary for determining enough rlist)",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--structure", type=str, required=True, help="Structure")
    parser.add_argument("--trajectory", type=str, default=None, help="Trajectory")
    parser.add_argument("--index", type=str, default=None, help="Gromacs index file name")
    parser.add_argument("--ligand-mol", type=str, default="Ligand", help="Index name for ligand selection")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()

    main(args)
    
