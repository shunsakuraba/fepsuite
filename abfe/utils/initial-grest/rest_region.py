#!/usr/bin/env python3

import sys
import mdtraj
import numpy
import argparse
import common_gmx_files

def main(args):
    system = mdtraj.load(args.structure)
    top = system.topology

    top_contents = common_gmx_files.parse_top(args.input)
    atomcnt = 0
    atom_start_table = {}
    for molname, molcount in top_contents["system"]:
        if molcount <= 0:
            continue
        if molname in atom_start_table:
            print("Warning: molecule %s appeared more than once" % molname)
        atom_start_table[molname] = atomcnt
        atomcnt += molcount * len(top_contents["moleculetypes"][molname])
    print(atom_start_table)

    prot = top.select(args.receptor)
    if args.target_mdtraj:
        lig = top.select(args.target)
    else:
        molname = args.target_gmx
        assert molname is not None
        lig = [atom_start_table[molname] + i for i in range(len(top_contents["moleculetypes"][molname]))]

    cutoff = args.range # nm
    prot_int = sorted(list(mdtraj.compute_neighbors(system, cutoff, lig, prot)[0])) # query = lig, haystack = prot, searches receptor resids near ligand
    if prot_int == []:
        raise RuntimeError("Nearby atoms not found")

    prot_resids = set([top.atom(i).residue.index for i in prot_int])

    prot_selstr = " or ".join(["resid %d" % r for r in prot_resids]) # use residue id to combat the multiemer case

    #print(prot_selstr, file=sys.stderr)

    target_indices = list(lig) + list(top.select("protein and (%s)" % prot_selstr))

    #print(target_indices, file=sys.stderr)

    target_indices = set(target_indices)

    with open(args.input) as fh, open(args.output, "w") as ofh:
        molname = None
        state = None
        for lraw in fh:
            l = lraw.split(";")[0].strip()
            if l.startswith("["):
                state = l.strip("[] \n")
            elif l == "":
                pass
            else:
                if state == "moleculetype":
                    ls = l.split()
                    molname = ls[0]
                elif state == "atoms":
                    ls = l.split()
                    assert molname in top_contents["moleculetypes"]
                    if molname in atom_start_table:
                        offset = atom_start_table[molname]
                        atomid_of_mol = int(ls[0]) - 1
                        if offset + atomid_of_mol in target_indices:
                            ls[1] += "_"
                            print(atomid_of_mol)
                        ofh.write(" ".join(ls) + "\n")
                        continue
            ofh.write(lraw)

        for lraw in fh:
            ofh.write(lraw)

def init_args():
    parser = argparse.ArgumentParser(description="Add underline to atom types so that REST / gREST program can work",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--structure", type=str, required=True, help="Structure, used for mdtraj selection")
    parser.add_argument("--target-mdtraj", type=str, help="Target molecule, this molecule will have scaled parameters. MDtraj selection syntax.")
    parser.add_argument("--target-gmx", type=str, help="Target molecule, this molecule will have scaled parameters. Gromacs moleculetype.")
    parser.add_argument("--receptor", type=str, required=True, help="Receptor molecules near the target will have sacled parameters. MDtraj selection syntax.")
    parser.add_argument("--range", type=float, required=True, help="Distance of \"near the target\" (unit: nm)")
    parser.add_argument("--input", type=str, required=True, help="input GROMACS topology file")
    parser.add_argument("--output", type=str, required=True, help="output GROMACS topology file")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()

    main(args)

