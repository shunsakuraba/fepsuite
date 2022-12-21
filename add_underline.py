
import mdtraj
import sys
import argparse


def main(args):
    confname = args.structure
    topname = args.topology
    writeto = args.output
    target_molecule = args.target_molecule

    modified = []
    with open(topname) as fh:
        section = None
        molname = None
        for l in fh:
            if l.startswith("["):
                section = l.split()[1]
                continue
            if ';' in l:
                l = l.split(';')[0]
            if l.strip() == "":
                continue
            if section == "moleculetype":
                molname = l.split()[0]
                continue
            if section == "atoms" and molname == "merged":
                ls = l.split()
                atomno = int(ls[0]) - 1
                atype = ls[1]
                btype = ls[8]
                acharge = float(ls[6])
                bcharge = float(ls[9])
                if atype != btype or acharge != bcharge: # charge should be literally equal
                    modified.append(atomno)

    #print(modified)
    gro = mdtraj.load(confname)
    protein = gro.topology.select(target_molecule)
    neighbors = mdtraj.compute_neighbors(gro, args.distance, query_indices=modified, haystack_indices=protein)
    neighbors = neighbors[0]
    resids = set()
    for n in neighbors:
        resids.add(gro.topology.atom(n).residue.resSeq)

    with open(writeto, "w") as ofh:
        print("; undelined resids = ", sorted(resids), file=ofh)
        # reprocess topology
        with open(topname) as fh:
            section = None
            molname = None
            for l in fh:
                comment = ""
                if l.startswith("["):
                    section = l.split()[1]
                elif ';' in l:
                    (l, comment) = l.split(';', 1)
                    comment = "; " + comment
                elif l.strip() == "":
                    pass
                elif section == "moleculetype":
                    molname = l.split()[0]
                elif section == "atoms" and molname == "merged":
                    ls = l.split()
                    resno = int(ls[2])
                    if resno in resids:
                        l = " ".join([ls[0], ls[1]+"_"] + ls[2:]) + "\n"
                print(l.rstrip() + comment, file=ofh)


def parse_args():
    parser = argparse.ArgumentParser(description="""Add underline to atoms in topology file for the REST2 calculation""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--structure', '-c', action='store', type=str, required=True,
                        help="Input structure file (used in mdtraj)")
    parser.add_argument('--topology', '-t', action='store', type=str, required=True,
                        help="Input topology file (must be preprocessed)")
    parser.add_argument('--output', '-o', action='store', type=str, required=True,
                        help="Output topology file")
    parser.add_argument('--target-moleucle', action='store', type=str, default="not (resname HOH SOL NA CL Na Cl K SOD CLA)",
                        help="Target molecules (in mdtraj syntax)")
    parser.add_argument('--distance', action='store', type=float,  default=0.4,
                        help='Residues having atoms within (this value) nm from perturbed atoms are "hot" region in REST2')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)

