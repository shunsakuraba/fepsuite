import re
import numpy
import os.path
import argparse
import math

modes = ["charging", "restrain", "annihilation-lig", "annihilation-complex"]

def lambda_schedule(mode, n):
    # return lambda values for each mode
    # actually for all cases we just return with equal distribution
    if mode == "charging":
        return numpy.linspace(0., 1., n, endpoint=True)
    if mode == "restrain":
        return numpy.linspace(0., 1., n, endpoint=True)
    if mode == "annihilation-lig":
        if False:
            base = numpy.linspace(0., math.pi / 2, n, endpoint=True)
            # state B: annihilated, dense lambda needed
            ret = numpy.sin(base)
            ret[0] = 0.
            ret[-1] = 1.
            return ret
        return numpy.linspace(0., 1., n, endpoint=True)
    if mode == "annihilation-complex":
        if False:
            base = numpy.linspace(0., math.pi / 2, n, endpoint=True)
            # state B: annihilated, dense lambda needed
            ret = numpy.sin(base)
            ret[0] = 0.
            ret[-1] = 1.
            return ret
        return numpy.linspace(0., 1., n, endpoint=True)
    raise RuntimeError("Unsupported mode")

def generate_decoupling(args):
    mode = args.mode

    with open(args.cominfo) as fh:
        lines = [l for l in fh if not l.startswith("#")]
        complex_com = [float(x) for x in lines[0].split()]
        lig_com = [float(x) for x in lines[1].split()]
        anchors = [int(x) for x in lines[2].split()]

    # read template file, small enough
    with open(os.path.join(args.template_dir, mode + ".mdp")) as fh:
        template = fh.read()
    if args.additional_mdp is not None:
        if mode == "charging":
            removelines = re.compile("\\s+(pull[-_]ncoords|pull[-_]ngroups|pull)\\s*=")
        else:
            removelines = re.compile("dummy-never-match")
        with open(args.additional_mdp) as fh:
            for l in fh:
                if not removelines.match(l):
                    template += l
    # just feed into format function

    lambda_values = lambda_schedule(mode, args.N)
    lambda_values_str = " ".join(["%.4f" % l for l in lambda_values])
    params = {
            "lambdas_formatted": lambda_values_str,
            "group_mol": args.mol_name,
            "group_ligand": "grp-lig",
            "group_complex": "grp-complex",
            "anchor_ligand": anchors[1] + 1,
            "anchor_complex": anchors[0] + 1,
            "complex": complex_com,
            "lig": lig_com
            }
    for i in range(args.N):
        params["lambda_state"] = i
        with open(os.path.join(args.output_mdp, "%s-%d.mdp" % (mode, i)), "w") as ofh:
            ofh.write(template.format(**params))

def init_args():
    parser = argparse.ArgumentParser(description="Generate restraint key atoms in protein and key atoms in ligand",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cominfo", type=str, default="cominfo", help="Center of mass information file")
    parser.add_argument("--mol-name", type=str, default="MOL", help="target molecul name")
    parser.add_argument("--output-mdp", type=str, default='.', help="mdp output directory")
    parser.add_argument("--additional-mdp", type=str, default=None, help="Additional mdp file (typically pull)")
    parser.add_argument("--mode", type=str, required=True, help="Stage (charging, restrain, annihilation-lig, or annihilation-xomplex)")
    parser.add_argument("-N", type=int, required=True, help="Number of states")
    parser.add_argument("--template-dir", type=str, required=True, help="Template directory")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()
    generate_decoupling(args)



