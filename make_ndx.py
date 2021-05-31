import mdtraj
import argparse
import common_gmx_files

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
    gmxtop = common_gmx_files.parse_top(args.topology)
    structure = mdtraj.load(args.structure)
    topology = structure.topology
    atomndx = 0
    found = False
    ligndx = None
    for (molname, molcount) in gmxtop["system"]:
        natom = len(gmxtop["moleculetypes"][molname])
        if molname == args.ligand:
            if molcount == 0:
                continue
            if molcount != 1:
                raise RuntimeError("More than one ligands are not supported")
            if found:
                raise RuntimeError("Ligand molecules appeared more than once")
            found = True
            ligndx = range(atomndx, atomndx + natom)
        atomndx += molcount * natom
    if ligndx is None:
        raise RuntimeError("Ligand not found: %s (candidates are: %s)" % (args.ligand, " ".join(gmxtop["moleculetypes"].keys())))
    with open(args.output, "w") as ofh:
        output_groups(ofh, "System", topology.select("all"))
        output_groups(ofh, "Ligand", ligndx)
        if args.receptor is not None:
            rec = topology.select(args.receptor)
            output_groups(ofh, "Receptor", rec)
            output_groups(ofh, "Ligand+Receptor", sorted(list(ligndx) + list(rec)))

def init_args():
    parser = argparse.ArgumentParser(description="Generate index file for the simulation",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--structure", type=str, help="Input structure file", required=True)
    parser.add_argument("--topology", type=str, help="Input topology file (must be preprocessed)", required=True)
    parser.add_argument("--output", type=str, help="Output ndx file", required=True)
    parser.add_argument("--ligand", type=str, help="Ligand selection in gromacs molecule name", required=True)
    parser.add_argument("--receptor", type=str, help="Receptor selection in mdtraj command")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()
    
    make_ndx(args)

