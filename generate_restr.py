import argparse
import math

def generate_itp(args):
    with open(args.restrinfo) as fh, open(args.itp, "w") as ofh:
        ls = fh.readlines()
        ls = [x for x in ls if not x.startswith("#")]
        anchor_atoms = [int(x) + 1 for x in ls[0].split()] # +1 for itp being 1-origin
        avgs = [float(x) for x in ls[1].split()]
        for i in range(1,6):
            avgs[i] = avgs[i] * 180.0 / math.pi # all format used degrees, convert it
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


def init_args():
    parser = argparse.ArgumentParser(description="Generate restraint key atoms in protein and key atoms in ligand",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--restrinfo", type=str, required=True, help="Restraint information file")
    parser.add_argument("--itp", type=str, required=True, help="Output restraint file")
    parser.add_argument("--distance-spring", type=float, default=418.68, help="Spring constant to apply (kJ/mol/nm^2)")
    parser.add_argument("--angle-spring", type=float, default=4.1868, help="Weight to the angle stdev (kJ/mol/rad^2)")
    parser.add_argument("--dihedral-spring", type=float, default=4.1868, help="Weight to the dihedral stdev (kJ/mol/rad^2)")
    parser.add_argument("--decouple-B", action="store_true",help="Set decoupling mode")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()

    generate_itp(args)
