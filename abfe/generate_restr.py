import argparse
import math
import sys

def abc(x):
    return chr(ord("A") + x)
def generate(args):
    with open(args.restrinfo) as fh:
        ls = fh.readlines()
        ls = [x for x in ls if not x.startswith("#")]
        anchor_atoms = [int(x) + 1 for x in ls[0].split()] # +1 for itp/ndx being 1-origin
        avgs = [float(x) for x in ls[1].split()]
        for i in range(1,6):
            avgs[i] = avgs[i] * 180.0 / math.pi # all format used degrees, convert it
            avgs[i] -= 360.0 * round(avgs[i] / 360.0) # normalize to -180 to 180
    if args.itp is not None:
        with open(args.itp, "w") as ofh:
            print("[ intermolecular_interactions ]", file=ofh)
            print("[ bonds ]", file=ofh)
            if args.decouple_B:
                print("%d %d 6 %.3f %.3f %.3f 0.000" % (anchor_atoms[2], anchor_atoms[3], avgs[0], args.distance_spring, avgs[0]), file=ofh)
            else:
                print("%d %d 6 %.3f %.3f" % (anchor_atoms[2], anchor_atoms[3], avgs[0], args.distance_spring), file=ofh)
            print(file=ofh)
            print("[ angles ]", file=ofh)
            if args.decouple_B:
                print("%d %d %d 1 %.3f %.3f %.3f 0.000" % (anchor_atoms[1], anchor_atoms[2], anchor_atoms[3], avgs[1], args.angle_spring, avgs[1]), file=ofh)
                print("%d %d %d 1 %.3f %.3f %.3f 0.000" % (anchor_atoms[2], anchor_atoms[3], anchor_atoms[4], avgs[2], args.angle_spring, avgs[2]), file=ofh)
            else:
                print("%d %d %d 1 %.3f %.3f" % (anchor_atoms[1], anchor_atoms[2], anchor_atoms[3], avgs[1], args.angle_spring), file=ofh)
                print("%d %d %d 1 %.3f %.3f" % (anchor_atoms[2], anchor_atoms[3], anchor_atoms[4], avgs[2], args.angle_spring), file=ofh)
            print(file=ofh)
            print("[ dihedrals ]", file=ofh)
            # uses functype=2 for restraining
            if args.decouple_B:
                print("%d %d %d %d 2 %.3f %.3f %.3f 0.000" % (anchor_atoms[0], anchor_atoms[1], anchor_atoms[2], anchor_atoms[3], avgs[3], args.dihedral_spring, avgs[3]), file=ofh)
                print("%d %d %d %d 2 %.3f %.3f %.3f 0.000" % (anchor_atoms[1], anchor_atoms[2], anchor_atoms[3], anchor_atoms[4], avgs[4], args.dihedral_spring, avgs[4]), file=ofh)
                print("%d %d %d %d 2 %.3f %.3f %.3f 0.000" % (anchor_atoms[2], anchor_atoms[3], anchor_atoms[4], anchor_atoms[5], avgs[5], args.dihedral_spring, avgs[5]), file=ofh)
            else:
                print("%d %d %d %d 2 %.3f %.3f" % (anchor_atoms[0], anchor_atoms[1], anchor_atoms[2], anchor_atoms[3], avgs[3], args.dihedral_spring), file=ofh)
                print("%d %d %d %d 2 %.3f %.3f" % (anchor_atoms[1], anchor_atoms[2], anchor_atoms[3], anchor_atoms[4], avgs[4], args.dihedral_spring), file=ofh)
                print("%d %d %d %d 2 %.3f %.3f" % (anchor_atoms[2], anchor_atoms[3], anchor_atoms[4], anchor_atoms[5], avgs[5], args.dihedral_spring), file=ofh)
            print(file=ofh)
    else:
        with open(args.ndx, "w") as ofh:
            # Write down anchor molecules
            print("".join(["[ anchor%s ]\n%d\n" 
                             % (abc(i), a) 
                           for (i, a) in enumerate(anchor_atoms)]),
                  file=ofh)
        with open(args.mdp, "w") as ofh:
            print("""
            pull = yes
            pull-nstxout = 0
            pull-nstfout = 0
            pull-ngroups = 6 ; A-F
            pull-ncoords = 6 ; 1 bond, 2 angles, 3 dihedrals
            pull-pbc-ref-prev-step-com = yes
            """, file=ofh)

            for i in range(6):
                print("pull-group%d-name    = anchor%s" % (i + 1, abc(i)), file=ofh)
                print("pull-group%d-pbcatom = %d" % (i + 1, anchor_atoms[i]), file=ofh)

            dref = {
                "distref_23": avgs[0],
                "distK_23_A": args.distance_spring,
                "distK_23_B": args.distance_spring,
                "angleref_123": avgs[1],
                "angleK_123_A": args.angle_spring,
                "angleK_123_B": args.angle_spring,
                "angleref_234": avgs[2],
                "angleK_234_A": args.angle_spring,
                "angleK_234_B": args.angle_spring,
                "dihedralref_0123": avgs[3],
                "dihedralK_0123_A": args.dihedral_spring,
                "dihedralK_0123_B": args.dihedral_spring,
                "dihedralref_1234": avgs[4],
                "dihedralK_1234_A": args.dihedral_spring,
                "dihedralK_1234_B": args.dihedral_spring,
                "dihedralref_2345": avgs[5],
                "dihedralK_2345_A": args.dihedral_spring,
                "dihedralK_2345_B": args.dihedral_spring
                }
            if args.decouple_B:
                dref["distK_23_B"] = 0.
                dref["angleK_123_B"] = 0.
                dref["angleK_234_B"] = 0.
                dref["dihedralK_0123_B"] = 0.
                dref["dihedralK_1234_B"] = 0.
                dref["dihedralK_2345_B"] = 0.

            print("""
            pull-coord1-type = umbrella
            pull-coord1-geometry = distance
            pull-coord1-groups = 3 4
            pull-coord1-dim = Y Y Y
            pull-coord1-k = {distK_23_A}
            pull-coord1-kB = {distK_23_B}
            pull-coord1-init = {distref_23}

            pull-coord2-type = umbrella
            pull-coord2-geometry = angle
            pull-coord2-groups = 3 2 3 4
            pull-coord2-dim = Y Y Y
            pull-coord2-k = {angleK_123_A}
            pull-coord2-kB = {angleK_123_B}
            pull-coord2-init = {angleref_123}

            pull-coord3-type = umbrella
            pull-coord3-geometry = angle
            pull-coord3-groups = 4 3 4 5
            pull-coord3-dim = Y Y Y
            pull-coord3-k = {angleK_234_A}
            pull-coord3-kB = {angleK_234_B}
            pull-coord3-init = {angleref_234}

            pull-coord4-type = umbrella
            pull-coord4-geometry = dihedral
            pull-coord4-groups = 1 2 2 3 3 4
            pull-coord4-dim = Y Y Y
            pull-coord4-k = {dihedralK_0123_A}
            pull-coord4-kB = {dihedralK_0123_B}
            pull-coord4-init = {dihedralref_0123}
            
            pull-coord5-type = umbrella
            pull-coord5-geometry = dihedral
            pull-coord5-groups = 2 3 3 4 4 5
            pull-coord5-dim = Y Y Y
            pull-coord5-k = {dihedralK_1234_A}
            pull-coord5-kB = {dihedralK_1234_B}
            pull-coord5-init = {dihedralref_1234}
           
            pull-coord6-type = umbrella
            pull-coord6-geometry = dihedral
            pull-coord6-groups = 3 4 4 5 5 6
            pull-coord6-dim = Y Y Y
            pull-coord6-k = {dihedralK_2345_A}
            pull-coord6-kB = {dihedralK_2345_B}
            pull-coord6-init = {dihedralref_2345}
           
            """.format(**dref), file=ofh)

def init_args():
    parser = argparse.ArgumentParser(description="Generate restraint key atoms in protein and key atoms in ligand",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--restrinfo", type=str, required=True, help="Restraint information file")
    parser.add_argument("--itp", type=str, default=None, help="Output restraint file")
    parser.add_argument("--mdp", type=str, default=None, help="Output pull-code file")
    parser.add_argument("--ndx", type=str, default=None, help="Output index file")
    parser.add_argument("--distance-spring", type=float, default=4184, help="Spring constant to apply (kJ/mol/nm^2)")
    parser.add_argument("--angle-spring", type=float, default=41.84, help="Weight to the angle stdev (kJ/mol/rad^2)")
    parser.add_argument("--dihedral-spring", type=float, default=41.84, help="Weight to the dihedral stdev (kJ/mol/rad^2)")
    parser.add_argument("--decouple-B", action="store_true",help="Set decoupling mode")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()

    [bitp, bmdp, bndx] = [x is None for x in [args.itp, args.mdp, args.ndx]]
    if not ((bitp, bmdp, bndx) == (True, False, False) or
            (bitp, bmdp, bndx) == (False, True, True)):
        print("specify either: (1) itp file, or (2) mdp and ndx files\n")
        sys.exit(1)
    generate(args)

