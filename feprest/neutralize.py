#!/usr/bin/python3

import sys
import random
import numpy
import mdtraj
import argparse
import os.path

class WaterModel:
    def __init__(self, itpfile):
        self.natom = None
        self.posresname = None
        self.negresname = None
        self.atomnames = None
        self.coord = None
        self.contents = []

        self.parse_itpfile(itpfile)
    
    def parse_itpfile(self, file):
        self.contents = []
        with open(file) as fh:
            section = None
            molecule = None
            for l in fh:
                self.contents.append(l.rstrip())
                if l.startswith("; SOLINFO"):
                    ls = l.split()
                    assert ls[0] == ";"
                    assert ls[1] == "SOLINFO", "Space/tab must be inserted after SOLINFO"
                    assert len(ls) == 5, "wrong SOLINFO fomat"
                    self.natom = int(ls[2])
                    self.posresname = ls[3]
                    self.negresname = ls[4]
                elif l.startswith("; SOLCOORD"):
                    ls = l.split()
                    assert ls[0] == ";"
                    assert ls[1] == "SOLCOORD", "Space/tab must be inserted after SOLCOORD"
                    assert self.natom is not None, "SOLCOORD must be after SOLINFO"
                    assert len(ls) == 2 + 3 * self.natom, "SOLCOORD atom number mismatch with SOLINFO"
                    self.coord = []
                    for ix in range(2, 2 + 3 * self.natom, 3):
                        self.coord.append(numpy.array([float(ls[ix]), float(ls[ix + 1]), float(ls[ix + 2])]))
                else:
                    lsc = l.split(';', 1)
                    lcmd = lsc[0]
                    if lcmd.strip() == "":
                        continue
                    lcmd = lcmd.lstrip()
                    if lcmd.startswith('['):
                        ls = lcmd.split()
                        section = ls[1]
                        continue
                    if section == "moleculetype":
                        ls = lcmd.split()
                        molecule = ls[0]
                        if molecule == "SOL2pos":
                            self.atomnames = []
                        continue
                    if section == "atoms" and molecule == "SOL2pos":
                        ls = lcmd.split()
                        assert int(ls[0]) == len(self.atomnames) + 1, "[ atoms ] section atom number mismatch"
                        self.atomnames.append(ls[4])


def sign(x):
    if x < 0:
        return -1
    elif x > 0:
        return 1
    else:
        return 0

def main(args):
    watermodel = WaterModel(os.path.join(args.feprest_dir, args.ff + ".ion.itp"))

    with open(args.topology) as fh, \
         open(args.output_topology, 'w') as ofh:
        sectionname = None
        moleculetype = None
        dmol = {}
        natoms = {}
        molecules = []
        remain = []
        for lraw in fh:
            l = lraw.split(";")[0].strip()
            ls = l.split()
            if l == "":
                pass
            elif l.startswith("["):
                if sectionname == "molecules":
                    # section name after molecules
                    remain.append(lraw)
                    for lraw in fh:
                        remain.append(lraw)
                    break
                sectionname = l.split()[1]
                if sectionname == "system":
                    for l in watermodel.contents:
                        ofh.write(l + "\n")
            elif l.startswith("#"):
                sys.stderr.write("Warning: input file must be preprocessed, but seems not\n")
            elif sectionname == "moleculetype":
                moleculetype = ls[0]
                _nrexcl = int(ls[1])
                dmol[moleculetype] = 0.
                natoms[moleculetype] = 0
            elif sectionname == "atoms":
                # sum up charges
                natoms[moleculetype] += 1
                if len(ls) >= 10:
                    # no charge difference unless atoms are perturbed
                    chga = float(ls[6])
                    chgb = float(ls[9])
                    dmol[moleculetype] += (chgb - chga)
            elif sectionname == "molecules":
                mol = ls[0]
                nmol = int(ls[1])
                molecules.append((mol, nmol))
                # do not write [ molecules ] yet, because we need to modify them
                continue
            ofh.write(lraw)
        
        dtot = 0
        mol_begin = [0]
        for (mol, nmol) in molecules:
            dtot += dmol[mol] * nmol
            mol_begin.append(mol_begin[-1] + natoms[mol] * nmol)
            print(f"DEBUG mol {mol} nmol {nmol} dmol[m] {dmol[mol]} dtot {dtot}")
        print(f"Total charge change: {dtot:.4f}")
        nchg = abs(round(dtot))
        if abs(dtot) - nchg > 1e-3:
            raise RuntimeError(f"Noninteger charge perturbation")
        sgn = sign(dtot)
        if sgn == -1:
            print("Perturbing water(s) into positive ion(s)")
            # easy case: auto and posonly both use "SOL2pos"
            # residue name & molecule name
            frommol = "SOL" 
            # number of atoms 
            nfrom = watermodel.natom
            fepmol = "SOL2pos"
            tomol = watermodel.posresname
        elif sgn == 1:
            if args.mode == "auto":
                print("Perturbing water(s) into negative ion(s)")
                # residue name & molecule name
                frommol = "SOL" 
                # number of atoms 
                nfrom = watermodel.natom
                fepmol = "SOL2neg"
                tomol = watermodel.negresname
            elif args.mode == "posonly":
                print("Perturbing positive ion(s) into water(s)")
                # residue name & molecule name
                frommol = watermodel.posresname
                # number of atoms 
                nfrom = 1
                fepmol = "pos2SOL"
                tomol = "SOL"
            else:
                raise RuntimeError(f"Unsupported mode: {args.mode}")
        else:
            print("No perturbation required")
            frommol = None
            nfrom = watermodel.natom
        nto = watermodel.natom # always perturbed molecules take this size
        
        output_mol = []
        exchange_mol = None
        # Find the last (frommol) molecule section. This peculiar feature is to support the following case:
        # [ molecules ]
        # Protein 1
        # SOL 35
        # SOL 13367
        # In the above case, SOL 35 is likely to be crystalline water and may not be preferable to be converted nor shuffled in the further processing.
        for (i, (mol, nmol)) in list(enumerate(molecules))[::-1]:
            if exchange_mol is not None:
                output_mol.append((mol, nmol))
            else:
                if mol == frommol and nmol >= nchg:
                    exchange_mol = i
                    if nmol > nchg:
                        output_mol.append((mol, nmol - nchg))
                    output_mol.append((fepmol, nchg))
                else:
                    output_mol.append((mol, nmol))
        if exchange_mol is None and frommol is not None:
            raise RuntimeError(f"Unable to find {nchg} consecutive {frommol} molecules in the system ")

        output_mol = output_mol[::-1]
        for (mol, nmol) in output_mol:
            ofh.write(f"{mol}     {nmol}\n")
        for l in remain:
            ofh.write(l)
    # end of topology read/write fh/ofh

    structure = mdtraj.load(args.gro)
    nonsolvent = structure.topology.select(f"not (resname SOL {watermodel.posresname} {watermodel.negresname})")
    if exchange_mol is None:
        bsolv = 0
        esolv = 0
        far = []
        far_chosen = []
        target_unchosen = []
    else:
        bsolv = mol_begin[exchange_mol]
        esolv = mol_begin[exchange_mol + 1]
        targets = list(range(bsolv, esolv, nfrom))
        far = mdtraj.compute_neighbors(structure, args.exclude_distance, nonsolvent, targets)[0]
        far_chosen = random.sample(list(far), nchg)
        target_unchosen = sorted(set(targets).difference(set(far_chosen)))
    with open(args.gro) as fh, \
        open(args.output_gro, 'w') as ofh:
        title = next(fh)
        ofh.write(title)
        lraw = next(fh)
        natom = int(lraw.strip())
        nshift = nchg * (nto - nfrom)
        new_natom = natom + nshift # may or may not change
        ofh.write(f"{new_natom:5d}\n") # note natom may exceed 5 digits but still this is OK
        outputbuf = []
        coords = []
        for lraw in fh:
            coords.append(lraw)
        box = coords[-1]
        del coords[-1]

        # output as-is until we hit the "solvent" region
        for i in range(0, bsolv):
            outputbuf.append(coords[i])
        
        # output sampled ions (perturbed)
        for c in far_chosen:
            if nfrom < nto:
                assert nfrom == 1 # other cases are not considered in the code below
                l = coords[c]
                # need to complement atoms not existing in the input file
                _resid = l[0:5]
                _resname = l[5:10]
                _atomname = l[10:15]
                atomno = l[15:20].strip()
                x = float(l[20:28].strip())
                y = float(l[28:36].strip())
                z = float(l[36:44].strip())
                basecrd = numpy.array([x, y, z])
                remain = "\n"
                if len(l) > 44:
                    remain = l[44:] # may have velocity info
                for i in range(nto):
                    crd = basecrd + watermodel.coord[i]
                    # atomname is the same but gromacs will not complain
                    newline = l[0:5] + "%-5s" % "SOL" + "%5s" % watermodel.atomnames[i] + l[15:20] + "%8.3f%8.3f%8.3f" % tuple(crd) + remain
                    outputbuf.append(newline)
            else:
                assert nfrom == nto
                for i in range(nfrom):
                    outputbuf.append(coords[c + i])
        
        # output ions not chosen (as-is)
        for c in target_unchosen:
            for i in range(nfrom):
                outputbuf.append(coords[c + i])            
        
        # to the end
        for i in range(esolv, len(coords)):
            outputbuf.append(coords[i])

        for lw in outputbuf:
            ofh.write(lw)
        ofh.write(box)

def parse_args():
    parser = argparse.ArgumentParser(description="Neutralize ions in topology as well as gro file", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--gro', '-c', action='store', type=str, required=True,
                        help="Input structure file (used in both internal parser and mdtraj)")
    parser.add_argument('--topology', '-t', action='store', type=str, required=True,
                        help="Input topology file (must be preprocessed)")
    parser.add_argument('--output-topology', '-o', action='store', type=str, required=True,
                        help="Output topology file")
    parser.add_argument('--output-gro', '-g', action='store', type=str, required=True,
                        help="Output gro file")
    parser.add_argument('--exclude-distance', action='store', type=float, default=0.4,
                        help='Water / ions having atoms within (this value) nm from non-solvent atoms are chosen at random and mutated')
    parser.add_argument('--feprest-dir', action='store', type=str, default=os.path.join(os.path.split(sys.argv[0])[0], "water_ion_models"),
                        help=".itp files under this directory is used to generate updated topology")
    parser.add_argument('--mode', action='store', type=str, default="auto",
                        help='Processing mode. Either "auto" or "posonly"')
    parser.add_argument('--ff', action='store', type=str, required=True,
                        help='Force field type, (this name).itp will be appended and used')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)

