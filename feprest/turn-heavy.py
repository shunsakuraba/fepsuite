import argparse

def main(args):
    section = None
    molecule = None
    ignored_molecules = set(args.ignore_moleculetype.split())
    outbuf = []
    with open(args.topology) as fh:
        for l in fh:
            outbuf.append(l)
            if section is None:
                if l.startswith("*"):
                    continue
            l = l.rstrip()
            lcomment = l.split(';', 1)
            l = lcomment[0]
            l = l.strip()
            if l == "":
                continue
            elif l.startswith("#"):
                continue
            elif l.startswith("["):
                sectionstr = l.lstrip("[")
                sectionstr = sectionstr.rstrip("]")
                sectionstr = sectionstr.strip()
                section = sectionstr
                continue
            ls = l.split()
            if ls == []:
                continue
            if section == "moleculetype":
                molecule = l.split()[0]
            elif section == "atoms" and molecule not in ignored_molecules:
                if len(ls) > 7:
                    # only when mass exists. Otherwise we just ignore.
                    # FIXME: perhaps we should check atomtype?
                    to_be_replaced = False
                    for replace_at in (7, 10): # state A, state B
                        if len(ls) > replace_at:
                            mass = float(ls[replace_at])
                            if mass < args.threshold:
                                ls[replace_at] = str(args.weight)
                                to_be_replaced = True
                    if to_be_replaced:
                        del outbuf[-1] # replace the line
                        outbuf.append(" ".join(ls) + "\n")
    with open(args.output, "w") as ofh:
        for l in outbuf:
            ofh.write(l)
                    
def parse_args():
    parser = argparse.ArgumentParser(description="""Turn hydrogens heavy""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--topology', '-p', action='store', type=str, required=True,
                        help="Input topology file (must be preprocessed)")
    parser.add_argument('--output', '-o', action='store', type=str, required=True,
                        help="Base directory to save state (default: current dir)")
    parser.add_argument('--ignore-moleculetype', action='store', type=str, default="SOL",
                        help="Moleculetype listed in this will be ignored")
    parser.add_argument('--threshold', action='store', type=float, default=3.5,
                        help="Atoms lighter than this value (amu) will get its mass buffed")
    parser.add_argument('--weight', action='store', type=float, default=8.0,
                        help="Atom masses will be changed to this value")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)
