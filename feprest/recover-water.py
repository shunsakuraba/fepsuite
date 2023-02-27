import argparse
import os.path
import sys

# After preprocessing, waters are all rigidified by settle, because #ifdef FLEXIBLE is removed. We want to set it back to flexible because the optimization gets unstable.

def main(args):
    topname = args.topology
    writeto = args.output

    with open(topname) as fh:
        section = None
        molname = None
        outbuf = []
        ignore = False
        for l in fh:
            if not ignore:
                outbuf.append(l)
            if l.startswith("["):
                section = l.split()[1]
                if section == "moleculetype":
                    if ignore:
                        outbuf.append(l)
                    ignore = False # end ignoring 
                continue
            if ';' in l:
                l = l.split(';')[0]
            if l.strip() == "":
                continue
            if section == "moleculetype":
                molname = l.split()[0]
                if molname == args.water_moltype:
                    del outbuf[-1] # prevent dup
                    ignore = True
                    # Copy contents of *.water.itp
                    with open(f"{args.water_dir}/{args.ff}.water.itp") as itpfh:
                        for li in itpfh:
                            if li.startswith('[') and li.split()[1] == "moleculetype":
                                pass # prevent dup
                            else:
                                outbuf.append(li)

    with open(writeto, "w") as ofh:
        for l in outbuf:
            ofh.write(l)
                    
def parse_args():
    parser = argparse.ArgumentParser(description="""Turn preprocessed solvent-containing force field back to un-preprocessed topology""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--topology', '-p', action='store', type=str, required=True,
                        help="Input topology file (must be preprocessed)")
    parser.add_argument('--output', '-o', action='store', type=str, required=True,
                        help="Base directory to save state (default: current dir)")
    parser.add_argument('--water-moltype', action='store', type=str, default="SOL",
                        help="This molecule type will be converted back to that of (--ff).water.itp")
    parser.add_argument('--water-dir', action='store', type=str, default=os.path.join(os.path.split(sys.argv[0])[0], "water_ion_models"),
                        help=".itp files under this directory is used to generate updated topology")
    parser.add_argument('--ff', action='store', type=str, required=True,
                        help='Force field type, (this name).water.itp will be appended and used')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)
