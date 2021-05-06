import mdtraj

def load_index(ndxfile):
    ret = {}
    with open(ndxfile) as fh:
        molname = None
        for l in fh:
            if l.startswtih("["):
                molname = l.strip().lstrip("[").rstrip("]").strip()
                # Note: in some versions of GROMACS the group names may be duplicated. Here we overwrite ret so that the last one will be used.
                ret[molname] = []
                continue
            ls = l.split()
            for ixstr in ls:
                ix = int(ixstr) - 1
                ret[molname].append(ix)
    return ret

def output_groups(fh, name, indices):
    print("[ %s ]" % name, file=fh)
    buf = ""
    for (i, ind) in enumerate(indices):
        buf += str(ind + 1)
        if i % 10 == 9:
            buf += "\n"
        else:
            buf += " "
    print(buf.rstrip() + "\n", file=fh)

def make_ndx(args):
    structure = mdtraj.load(args.structure)
    topology = structure.topology
    index = load_index(args.index)
    ligsel = index[args.ligand] # [ moleculetype ] name will appear in default ndx file
    with open(args.output, "w") as ofh:
        output_groups(ofh, "System", topology.select("all"))
        output_groups(ofh, "Ligand", ligsel)
        if args.receptor:
            rec = topology.select(args.receptor)
            output_groups(ofh, "Receptor", rec)
            output_groups(ofh, "Ligand+Receptor", list(ligsel) + list(rec))

def init_args():
    parser = argparse.ArgumentParser(description="Generate index file for the simulation",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--structure", type=str, help="Input structure file", required=True)
    parser.add_argument("--index", type=str, help="Input index file", required=True)
    parser.add_argument("--output", type=str, help="Output ndx file", required=True)
    parser.add_argument("--ligand", type=str, help="Ligand selection in gromacs molecule name", required=True)
    parser.add_argument("--receptor", type=str, help="Receptor selection in mdtraj command")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()
    
    make_ndx(args)

