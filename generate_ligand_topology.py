import argparse
import math
import sys
import mdtraj
import common_gmx_files

def output_topology_with_only_ligand(topology_in, molname, topology_out):
    section = None
    with open(topology_in) as fh, open(topology_out, "w") as ofh:
        for lraw in fh:
            if section is None:
                if lraw.startswith("*"):
                    ofh.write(lraw)
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
                ofh.write(lraw)
                continue
            if l == "":
                pass
            elif l.startswith("#"):
                raise RuntimeError("topology is not preprocessed")
            elif section == "molecules":
                if ls[0] == molname:
                    ofh.write(lraw)
                continue
            ofh.write(lraw)

def generate(args):
    ptop = common_gmx_files.parse_top(args.topology)
    molcomposition = ptop['system']
    moleculetypes = ptop['moleculetypes']
    ligatoms = moleculetypes[args.mol]

    # validation 1: total ligand numbers
    totcount = 0
    for (m, count) in molcomposition:
        if m == args.mol:
            totcount += count
    if totcount == 0:
        raise RuntimeError("Ligand did not appear in topology")
    if totcount > 1:
        raise RuntimeError("Ligand appeared more than onece in topology")
    
    # validation 2: check index mask
    atomptr = 0
    selected_from_top = None
    for (m, c) in molcomposition:
        mollist = moleculetypes[m]
        natom = len(mollist)
        if m == args.mol:
            selected_from_top = [atomptr + i for i in range(natom)]
        atomptr += natom * c
    if args.index:
        indexinfo = common_gmx_files.parse_index(args.index)
        selected = indexinfo[args.ligand_group]
        if selected_from_top != selected:
            raise RuntimeError("Index file mismatched with topology")

    # validation completed
    structure = mdtraj.load(args.structure)
    outstr = structure.atom_slice(selected_from_top)
    outstr.save(args.output_ligand_structure)
    
    # output topology
    output_topology_with_only_ligand(args.topology, args.mol, args.output_ligand_topology)

    # compute total charge
    totalcharge = 0.
    moleculeinfo = moleculetypes[args.mol]
    for a in moleculeinfo:
        totalcharge += a[4]

    with open(args.total_charge, "w") as ofh:
        print("%.5f" % totalcharge, file=ofh)

def init_args():
    parser = argparse.ArgumentParser(description="Generate ligand topology",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mol", type=str, required=True, help="Moleculetype name for the target ligand")
    parser.add_argument("--topology", type=str, required=True, help="GROMACS topology file")
    parser.add_argument("--structure", type=str, required=True, help="Input (complex) structure")
    parser.add_argument("--index", type=str, help="Input index file")
    parser.add_argument("--ligand-group", type=str, default="Ligand", help="Input group")
    parser.add_argument("--output-ligand-structure", type=str, required=True, help="Ligand structure output (PDB recommended)")
    parser.add_argument("--output-ligand-topology", type=str, required=True, help="Ligand topology output (top file)")
    parser.add_argument("--total-charge", type=str, required=True, help="Total charge of the ligand will be output to this file")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()

    generate(args)

