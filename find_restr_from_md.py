import argparse
import math
import warnings

def find_restraints(args):
    # importing at this timing is discouraged in python, but it's too heavy to just show --help
    import numpy
    import mdtraj 

    nm_of_angstrom = 0.1

    refstructure = mdtraj.load(args.topology)
    topology = refstructure.topology
    iterbegin_trajectory = lambda: mdtraj.iterload(args.trajectory, top=topology)

    protix = topology.select(args.prot_sel)
    ligix = topology.select(args.lig_sel)
    ligheavyix = topology.select("(%s) and not type H" % args.lig_sel)
    # before scanning the trajectory, check whether --anchor-atoms are sane
    anchor_names = args.anchor_atoms.split(',')
    [anc1, anc2, anc3] = anchor_names
    _anchors = topology.select("(%s) and (name %s %s %s)" % (args.prot_sel, anc1, anc2, anc3))

    # pass 0: determine ligand representative atoms
    # first atom: center of the ligand. Assumes that refstructure makes ligand whole.
    center = mdtraj.compute_center_of_mass(refstructure.atom_slice(ligheavyix))[0, :]
    d2s = numpy.sum((refstructure.xyz[0, ligheavyix, :] - center[numpy.newaxis, :])**2, axis=1)
    lig_a_ix = ligheavyix[numpy.argmin(d2s)]
    print("ligand first atom:", lig_a_ix, topology.atom(lig_a_ix))

    # second atom: farthest from the first
    lig_a = refstructure.xyz[0, lig_a_ix, :]
    d2s = numpy.sum((refstructure.xyz[0, ligheavyix, :] - center[numpy.newaxis, :])**2, axis=1)
    lig_c_ix = ligheavyix[numpy.argmax(d2s)]
    assert lig_a_ix != lig_c_ix
    print("ligand second atom:", lig_c_ix, topology.atom(lig_c_ix))

    # TODO FIXME(shun): this condition is irrational, needs a fix
    # third atom: a-b-c angle > 90 deg, farthest from the first
    angle_indices = [[lig_a_ix, ix, lig_c_ix] for ix in ligheavyix]
    angle_indices2 = [[ix, lig_a_ix, lig_c_ix] for ix in ligheavyix]
    angle_indices3 = [[ix, lig_c_ix, lig_a_ix] for ix in ligheavyix]
    assert numpy.shape(numpy.array(angle_indices))[1] == 3

    def calc_angle(indices):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            angles = mdtraj.compute_angles(refstructure, indices)[0]
            angles[numpy.isnan(angles)] = 0.0
        return angles
    angles = calc_angle(angle_indices)
    angles2 = calc_angle(angle_indices2)
    angles3 = calc_angle(angle_indices3)
    rangle = (90.0 * math.pi / 180.0) 
    #mask = (angles < rangle) & (angles2 < rangle) & (angles3 < rangle)
    mask = (angles2 < rangle) & (angles3 < rangle)
    d2s[mask] = 0.0
    lig_b_ix = ligheavyix[numpy.argmax(d2s)]
    assert lig_a_ix != lig_b_ix
    assert lig_c_ix != lig_b_ix
    print("ligand third atom:", lig_b_ix, topology.atom(lig_b_ix))


    # pass 1: find neighboring protein atoms
    neighbor_sets = set()
    for chunk in iterbegin_trajectory():
        # update neighboring atoms list
        matcheslist = mdtraj.compute_neighbors(chunk, nm_of_angstrom * args.search_dist,
                ligix, protix, periodic=True)
        for m in matcheslist:
            neighbor_sets |= set(m)

    # pass 1.5: determine residues
    residues = set()
    for a in neighbor_sets:
        residues.add(topology.atom(a).residue.index)
    residues = list(residues)
    residues.sort()
    print("Nearby residue ids (0-origin):", residues)
    residues_filtered = []
    anchors = []
    # I do believe this is an extreme case, but just in case...
    for r in residues:
        anchors_cur = [topology.select("(%s) and resid %d and name %s" % (args.prot_sel, r, n)) for n in anchor_names]
        fail = False
        for al in anchors_cur:
            if len(al) != 1:
                fail = True
                break
        if fail:
            continue
        anchors_cur = [x for [x] in anchors_cur]
        residues_filtered.append(r)
        anchors.append(numpy.array(anchors_cur, dtype=int))
    print("Filtered out %d residues" % (len(residues) - len(residues_filtered)))

    # generate pass
    lig_indices = numpy.array([lig_a_ix, lig_b_ix, lig_c_ix], dtype=int)
    ind_perm = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
    anclig_list = []
    for (i, r) in enumerate(residues_filtered):
        for ancperm in ind_perm:
            for ligperm in ind_perm:
                anclig_list.append((anchors[i][ancperm], lig_indices[ligperm]))

    # pass 2: determine the best residue and its ordering
    dihed_indices_a = [[anc[0], anc[1], anc[2], lig[0]] for anc, lig in anclig_list]
    dihed_indices_b = [[anc[1], anc[2], lig[0], lig[1]] for anc, lig in anclig_list]
    dihed_indices_c = [[anc[2], lig[0], lig[1], lig[2]] for anc, lig in anclig_list]
    angle_indices_a = [[anc[1], anc[2], lig[0]] for anc, lig in anclig_list]
    angle_indices_b = [[anc[2], lig[0], lig[1]] for anc, lig in anclig_list]
    dist_indices = [[anc[2], lig[0]] for anc, lig in anclig_list]

    npairs = len(dist_indices)
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
        for i in range(1, dh.shape[1]):
            diff = dh[:, i] - dh[:, i - 1]
            dh[:, i] -= (2 * math.pi) * numpy.round(diff / (2 * math.pi))

    periodic_diheds(diheds_a)
    periodic_diheds(diheds_b)
    periodic_diheds(diheds_c)

    print(dists.shape)
    print(angles_a.shape)
    print(angles_b.shape)
    print(diheds_a.shape)
    print(diheds_b.shape)
    print(diheds_c.shape)

    total_weights = (
            args.distance_weight * numpy.var(dists, axis=0) +
            args.angle_weight * (numpy.var(angles_a, axis=0) + numpy.var(angles_b, axis=0)) + 
            args.dihedral_weight * (numpy.var(diheds_a, axis=0) + numpy.var(diheds_b, axis=0) + numpy.var(diheds_c, axis=0))
            )
    # prevent ang < 45 or ang > 135 (log sin theta part gets unstable)
    # FIXME(shun): actually such an angle should not be chosen from the beginning
    avgangle_a_all = numpy.mean(angles_a, axis=0)
    avgangle_b_all = numpy.mean(angles_b, axis=0)
    banned_a = numpy.logical_or(avgangle_a_all < (math.pi * 45.0 / 180.0), avgangle_a_all > (math.pi * 135.0 / 180.0))
    banned_b = numpy.logical_or(avgangle_b_all < (math.pi * 45.0 / 180.0), avgangle_b_all > (math.pi * 135.0 / 180.0))
    banned = numpy.logical_or(banned_a, banned_b)
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
    parser.add_argument("--prot-sel", type=str, default="protein", help="Selection text for protein")
    parser.add_argument("--lig-sel", type=str, default="not (protein or water or resname NA CL K)", help="Selection text for ligand")
    parser.add_argument("--search-dist", type=float, default=5.0, help="Distance as an upper limit to search contacting atoms. [angstrom]")
    parser.add_argument("--anchor-atoms", type=str, default="CA,C,N", help="Anchor atom names. Atoms are chosen to be in the same residue.")
    parser.add_argument("--topology", type=str, help="Topology file, possibly .pdb or .gro file.", required=True)
    parser.add_argument("--trajectory", type=str, help="Trajectory file.", required=True)
    parser.add_argument("--output", type=str, help="Output configuration file", required=True)
    parser.add_argument("--distance-weight", type=float, default=1.0, help="Weight to the distance stdev")
    parser.add_argument("--angle-weight", type=float, default=1.0, help="Weight to the angle stdev")
    parser.add_argument("--dihedral-weight", type=float, default=1.0, help="Weight to the dihedral stdev")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()
    
    find_restraints(args)
