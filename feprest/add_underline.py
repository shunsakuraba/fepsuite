# FIXME(shun): this script needs overhaul
import mdtraj
import argparse
from collections import defaultdict

def main(args):
    confname = args.structure
    topname = args.topology
    writeto = args.output
    target_molecule = args.target_molecule
    ignored_sections = set(args.non_perturbed_moleculetype.split())

    perturbed = []
    with open(topname) as fh:
        section = None
        molecule = None
        molcomposition = []
        natom_moleculetypes = {}
        natom = 0
        for l in fh:
            if l.startswith("["):
                section = l.split()[1]
                if section in ["molecules", "system", "moleculetype"]:
                    # here we do a cleanup, as "molecules" / "system" must be after the end of molecular section
                    if molecule is not None:
                        natom_moleculetypes[molecule] = natom
                        molecule = None
                        natom = 0
                continue
            if ';' in l:
                l = l.split(';')[0]
            if l.strip() == "":
                continue
            ls = l.split()
            if section == "moleculetype":
                molecule = ls[0]
                continue
            if section == "atoms":
                atomno = int(ls[0]) - 1
                atype = ls[1]
                btype = ls[8] if len(ls) > 8 else atype
                acharge = float(ls[6])
                bcharge = float(ls[9]) if len(ls) > 9 else acharge
                if (atype != btype or acharge != bcharge) and molecule not in ignored_sections:
                    perturbed.append((molecule, atomno))
                natom += 1
            elif section == "molecules":
                molecule = ls[0]
                nmol = int(ls[1])
                molcomposition.append((molecule, nmol))
    
    # get starting numbers for each molecule
    atomno_starts = []
    curatomno = 0
    starts_by_name = defaultdict(list)
    for (m, n) in molcomposition:
        atomno_starts.append(curatomno)
        natom_mol = natom_moleculetypes[m]
        for k in range(n):
            starts_by_name[m].append(curatomno + natom_mol * k)
        curatomno += natom_mol * n

    # get indices of perturbed atoms in the whole molecule
    perturbed_atomno = []
    for (m, a) in perturbed:
        if len(starts_by_name[m]) > 1 and not args.ignore_perturbing_multiple_molecules:
            raise RuntimeError(f"Molecule named {m} are perturbed but there are multiple ({len(starts_by_name[m])}) molecules in the system")
        for ba in starts_by_name[m]:
            perturbed_atomno.append(ba + a)
    perturbed_atomno.sort()

    # load conformation
    gro = mdtraj.load(confname)
    protein = gro.topology.select(target_molecule)
    neighbors = mdtraj.compute_neighbors(gro, args.distance, query_indices=perturbed_atomno, haystack_indices=protein)
    neighbors = neighbors[0]
    rest_residues = set()
    rest_atoms = set()

    # list residues that need to be RESTed
    for n in neighbors:
        res = gro.topology.atom(n).residue
        rest_residues.add(res)
    
    # list atoms that need to be perturbed
    for r in rest_residues:
        for a in res.atoms:
            rest_atoms.add(a.index)

    with open(writeto, "w") as ofh:
        # residue ids
        residue_ids = [r.resSeq for r in rest_residues]
        print("; undelined resids = ", sorted(residue_ids), file=ofh)
        # reprocess topology
        with open(topname) as fh:
            section = None
            molecule = None
            molecule_offset = None
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
                    molecule = l.split()[0]
                    if molecule in starts_by_name and len(starts_by_name[molecule]) >= 1:
                        molecule_offset = starts_by_name[molecule][0]
                    else:
                        molecule_offset = None
                elif section == "atoms" and molecule not in ignored_sections and molecule_offset is not None:
                    # 3rd condition == molecule is actually listed in [ molecules ]
                    ls = l.split()
                    atomno = int(ls[0]) - 1
                    if molecule_offset + atomno in rest_atoms:
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
    parser.add_argument('--non-perturbed-moleculetype', action='store', type=str, default="SOL SOL2pos SOL2neg",
                        help="Molecule names which is ignored when determining residues to perturb (in GROMACS moleculetype name)")
    parser.add_argument('--target-molecule', action='store', type=str, default="not (resname HOH SOL NA CL Na Cl K SOD CLA)",
                        help="Target molecules (in mdtraj syntax)")
    parser.add_argument('--distance', action='store', type=float,  default=0.4,
                        help='Residues having atoms within (this value) nm from perturbed atoms are "hot" region in REST2')
    parser.add_argument('--ignore-perturbing-multiple-molecules', action='store_true', default=False,
                        help='Perturbing multiple molecules typically mean something is wrong and forbidden by default. With this flag this script ignores the case.')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)

