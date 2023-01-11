import argparse

def resurrect(args):
    flexsol = []
    with open(args.flexible) as fh:
        section = None
        in_solvent = False
        for lraw in fh:
            if section is None:
                if lraw.startswith("*"):
                    continue
            l = lraw.rstrip()
            lcomment = l.split(';', 1)
            l = lcomment[0]
            l = l.strip()
            ls = l.split()
            if l.startswith("["):
                sectionstr = l.lstrip("[")
                sectionstr = sectionstr.rstrip("]")
                sectionstr = sectionstr.strip()
                section = sectionstr
                if in_solvent and section == "moleculetype":
                    in_solvent = False
                elif in_solvent:
                    flexsol.append(lraw)
                continue
            if in_solvent:
                flexsol.append(lraw)
            if l == "":
                pass
            elif l.startswith("#"):
                raise RuntimeError("topology is not preprocessed")
            elif section == "moleculetype":
                molname = ls[0]
                if args.solvent == molname:
                    in_solvent = True
    with open(args.topology) as fh, open(args.output, "w") as ofh:
        in_solvent = False
        for lraw in fh:
            if section is None:
                if lraw.startswith("*"):
                    continue
            l = lraw.rstrip()
            lcomment = l.split(';', 1)
            l = lcomment[0]
            l = l.strip()
            ls = l.split()
            if l.startswith("["):
                sectionstr = l.lstrip("[")
                sectionstr = sectionstr.rstrip("]")
                sectionstr = sectionstr.strip()
                section = sectionstr
                if in_solvent and section == "moleculetype":
                    in_solvent = False
                    ofh.write("#else\n")
                    for lflex in flexsol:
                        ofh.write(lflex)
                    ofh.write("#endif\n")
                ofh.write(lraw)
                continue
            ofh.write(lraw)
            if l == "":
                pass
            elif l.startswith("#"):
                raise RuntimeError("topology is not preprocessed")
            elif section == "moleculetype":
                molname = ls[0]
                if args.solvent == molname:
                    in_solvent = True
                    ofh.write("#ifndef FLEXIBLE\n")

def init_args():
    parser = argparse.ArgumentParser(description="Copy different SOL definition to allow flexible solvent during minimization",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--flexible", type=str, help="Input topology file with flexible version (must be preprocessed)", required=True)
    parser.add_argument("--topology", type=str, help="Input topology file (must be preprocessed)", required=True)
    parser.add_argument("--output", type=str, help="Output topology file", required=True)
    parser.add_argument("--solvent", type=str, default="SOL", help="Solvent molecule name")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()
    
    resurrect(args)

