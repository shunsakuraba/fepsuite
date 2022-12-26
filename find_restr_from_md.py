import argparse
import math
import common_gmx_files
import numpy
import mdtraj
import itertools

def find_restraints(args):
    refstructure = mdtraj.load(args.topology)
    topology = refstructure.topology
    def iterbegin_trajectory():
        return mdtraj.iterload(args.trajectory, top=topology)

    indices = common_gmx_files.parse_index(args.index)
    protix = indices[args.prot_sel]
    ligix = indices[args.lig_sel]
    ligheavyix = list(set(topology.select("not type H")).intersection(set(ligix)))
    # before scanning the trajectory, check whether --anchor-atoms are sane
    anchor_names = args.anchor_atoms.split(',')
    anchors = list(set(topology.select(f'name {" ".join(anchor_names)}')).intersection(protix))
    if len(anchors) == 0:
        raise RuntimeError("No match for anchor atoms in system")
    if len(set(ligheavyix).intersection(anchors)) > 0:
        raise RuntimeError("Ligand and receptor atoms are overlapping")

    # find neighboring atoms in both state
    neighbors_lig = mdtraj.compute_neighbors(refstructure, args.search_dist, anchors, haystack_indices=ligheavyix)[0]
    neighbors_anchor = mdtraj.compute_neighbors(refstructure, args.search_dist, ligheavyix, haystack_indices=anchors)[0]

    if len(neighbors_lig) == 0 or len(neighbors_anchor) == 0:
        raise RuntimeError(f"Number of neighbor atoms in ligand was {len(neighbors_lig)}, in anchor was {len(neighbors_anchor)}")

    atom_pairs = list(itertools.product(neighbors_anchor, neighbors_lig))
    distances = mdtraj.compute_distances(refstructure, atom_pairs)[0, :] # [1, num_pairs] -> [num_pairs]
    filtered_pairs = [p for p, d in zip(atom_pairs, distances) if d <= args.search_dist]

    print(f"There are {len(filtered_pairs)} pairs of atoms considered")
    if len(filtered_pairs) == 0:
        raise RuntimeError(f"Could not find pairs within {args.search_dist} nm")
    # for each combination find a possible combination of 6 atoms

    def find_bonded(a, searchindex):
        threshold = 0.22 # 0.22 nm, should be sufficiently long for any atoms
        found = mdtraj.compute_neighbors(refstructure, threshold, [a], [x for x in searchindex if x != a])[0]
        return found
    
    def find_two_bonds(a, searchindex):
        bonded = find_bonded(a, searchindex)
        ret = []
        for b in bonded:
            angled = find_bonded(b, searchindex)
            angled = [x for x in angled if x != a]
            for c in angled:
                ret.append((a, b, c))
        return ret

    # enumerate evey possible 6-atom combinations
    anclig_list = []
    for (anc, lig) in filtered_pairs:
        anc_angles = find_two_bonds(anc, anchors)
        lig_angles = find_two_bonds(lig, ligheavyix)
        anc_rev_angles = [(c, b, a) for (a, b, c) in anc_angles]
        anclig_list.extend(itertools.product(anc_rev_angles, lig_angles))

    print(f"There are {len(anclig_list)} sets of 6-atom tuples considered")

    # determine the best residue and its ordering
    dihed_indices_a = [[anc[0], anc[1], anc[2], lig[0]] for anc, lig in anclig_list]
    dihed_indices_b = [[anc[1], anc[2], lig[0], lig[1]] for anc, lig in anclig_list]
    dihed_indices_c = [[anc[2], lig[0], lig[1], lig[2]] for anc, lig in anclig_list]
    angle_indices_a = [[anc[1], anc[2], lig[0]] for anc, lig in anclig_list]
    angle_indices_b = [[anc[2], lig[0], lig[1]] for anc, lig in anclig_list]
    dist_indices = [[anc[2], lig[0]] for anc, lig in anclig_list]

    npairs = len(dist_indices)
    # all these values have the same dimension of [frame, pairs]
    diheds_a = numpy.zeros((0, npairs))
    diheds_b = numpy.zeros((0, npairs))
    diheds_c = numpy.zeros((0, npairs))
    angles_a = numpy.zeros((0, npairs))
    angles_b = numpy.zeros((0, npairs))
    dists = numpy.zeros((0, npairs))
    for chunk in iterbegin_trajectory():
        # calculate actual distances/angles/dihedrals
        dist_chunk = mdtraj.compute_distances(chunk, dist_indices)
        dists = numpy.concatenate((dists, dist_chunk), axis=0)

        angle_chunk_a = mdtraj.compute_angles(chunk, angle_indices_a)
        angles_a = numpy.concatenate((angles_a, angle_chunk_a), axis=0)

        angle_chunk_b = mdtraj.compute_angles(chunk, angle_indices_b)
        angles_b = numpy.concatenate((angles_b, angle_chunk_b), axis=0)

        dihed_chunk_a = mdtraj.compute_dihedrals(chunk, dihed_indices_a)
        diheds_a = numpy.concatenate((diheds_a, dihed_chunk_a), axis=0)

        dihed_chunk_b = mdtraj.compute_dihedrals(chunk, dihed_indices_b)
        diheds_b = numpy.concatenate((diheds_b, dihed_chunk_b), axis=0)

        dihed_chunk_c = mdtraj.compute_dihedrals(chunk, dihed_indices_c)
        diheds_c = numpy.concatenate((diheds_c, dihed_chunk_c), axis=0)

    def periodic_diheds(dh):
        def periodic_impl(dh, baseval):
            diff = dh[:, :] - baseval[numpy.newaxis, :] # baseval have [pairs] dimension
            dh[:, :] -= (2 * math.pi) * numpy.round(diff / (2 * math.pi))
        periodic_impl(dh, dh[0, :])
        xmean = numpy.mean(dh, axis=0)
        xmean -= (2 * math.pi) * numpy.round(xmean / (2 * math.pi)) # shift the average to be around [-pi, pi]
        periodic_impl(dh, xmean)

    periodic_diheds(diheds_a)
    periodic_diheds(diheds_b)
    periodic_diheds(diheds_c)

    total_weights = (
            args.distance_weight * numpy.var(dists, axis=0) +
            args.angle_weight * (numpy.var(angles_a, axis=0) + numpy.var(angles_b, axis=0)) + 
            args.dihedral_weight * (numpy.var(diheds_a, axis=0) + numpy.var(diheds_b, axis=0) + numpy.var(diheds_c, axis=0))
            )
    # prevent ang < 45 or ang > 135 (log sin theta part gets unstable)
    avgangle_a_all = numpy.mean(angles_a, axis=0)
    avgangle_b_all = numpy.mean(angles_b, axis=0)
    banned_a = numpy.logical_or(avgangle_a_all < (math.pi * 45.0 / 180.0), avgangle_a_all > (math.pi * 135.0 / 180.0))
    banned_b = numpy.logical_or(avgangle_b_all < (math.pi * 45.0 / 180.0), avgangle_b_all > (math.pi * 135.0 / 180.0))
    banned = numpy.logical_or(banned_a, banned_b)
    print(f"Removed {numpy.sum(banned)} tuples due to bad angles")
    total_weights[banned] = 1e+6

    bestcand = numpy.argmin(total_weights)
    avgdist = numpy.mean(dists[:, bestcand])
    avgangle_a = numpy.mean(angles_a[:, bestcand])
    avgangle_b = numpy.mean(angles_b[:, bestcand])
    avgdihed_a = numpy.mean(diheds_a[:, bestcand])
    avgdihed_b = numpy.mean(diheds_b[:, bestcand])
    avgdihed_c = numpy.mean(diheds_c[:, bestcand])
    bestanclig = anclig_list[bestcand]
    with open(args.output, "w") as ofh:
        print("# ancA ancB ancC ligA ligB ligC", file=ofh)
        print(bestanclig[0][0], bestanclig[0][1], bestanclig[0][2],
              bestanclig[1][0], bestanclig[1][1], bestanclig[1][2], file=ofh)
        print("#", 
                topology.atom(bestanclig[0][0]),
                topology.atom(bestanclig[0][1]),
                topology.atom(bestanclig[0][2]),
                topology.atom(bestanclig[1][0]),
                topology.atom(bestanclig[1][1]),
                topology.atom(bestanclig[1][2]),
                file=ofh)
        print("# avg dist / angle-a,b / dihed a,b,c (nm or radian)", file=ofh)
        print(avgdist, avgangle_a, avgangle_b, avgdihed_a, avgdihed_b, avgdihed_c, file=ofh)
        print('#', numpy.std(dists[:, bestcand]), numpy.std(angles_a[:, bestcand]), numpy.std(angles_b[:, bestcand]),
                numpy.std(diheds_a[:, bestcand]), numpy.std(diheds_b[:, bestcand]), numpy.std(diheds_c[:, bestcand]), file=ofh)

def init_args():
    parser = argparse.ArgumentParser(description="Generate restraint key atoms in protein and key atoms in ligand",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--prot-sel", type=str, default="Receptor", help="Index name for receptor selection")
    parser.add_argument("--lig-sel", type=str, default="Ligand", help="Index name for ligand selection")
    parser.add_argument("--index", type=str, default="index.ndx", help="Index file name")
    parser.add_argument("--search-dist", type=float, default=0.5, help="Distance as an upper limit to search contacting atoms. [nm]")
    parser.add_argument("--anchor-atoms", type=str, default="CB,CA,C,N,O", help="Anchor atom names.")
    parser.add_argument("--topology", type=str, help="Topology file, possibly .pdb or .gro file.", required=True)
    parser.add_argument("--trajectory", type=str, help="Trajectory file.", required=True)
    parser.add_argument("--output", type=str, help="Output configuration file", required=True)
    parser.add_argument("--distance-weight", type=float, default=4184.0, help="Weight to the distance stdev. Default weight values comes from restraint energy constants.")
    parser.add_argument("--angle-weight", type=float, default=41.84, help="Weight to the angle stdev.")
    parser.add_argument("--dihedral-weight", type=float, default=41.84, help="Weight to the dihedral stdev.")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()
    
    find_restraints(args)
