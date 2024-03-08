import sys
import os

# FIXME: I know this script's quality is abysmal and needs updating

[_, topfrom, pdbfrom, targetresid, topto, pdbto] = sys.argv
targetresid = int(targetresid)
targetresids = [targetresid-1, targetresid, targetresid+1] # FIXME: when mutation resides at the termini this is apparently wrong. The current workaround is to make residue number gapped.


curatom = 0
atomidmap = {}
numatoms = {
        "bonds": 2,
        "angles": 3,
        "dihedrals": 4,
        "pairs": 2,
        "cmap": 5,
        "position_restraints": 1
        }
with open(topfrom) as fh, open(topto, "w") as ofh:
    def writels():
        print(" ".join([str(s) for s in ls]) + comment, file=ofh)
    sectionname = None
    moleculename = None
    for l in fh:
        ltmp = l.split(';', 1)
        if len(ltmp) > 1:
            comment = ";" + ltmp[1]
        else:
            comment = ""
        ls = ltmp[0].split()
        if len(ls) == 0:
            pass
        elif ls[0] == '[':
            sectionname = ls[1]
        elif l[0] == "#":
            pass
        elif sectionname == "moleculetype":
            moleculename = ls[0]
        elif sectionname == "atoms" and moleculename == "merged":
            assert len(ls) >= 9, "moleculetype \"merged\" must have all atoms and residues perturbed"
            resid = int(ls[2])
            atomid = int(ls[0])
            if resid in targetresids:
                curatom += 1
                atomidmap[atomid] = curatom
                ls[0] = str(curatom)
                ls[5] = str(curatom)
                writels()
            perturbed = False
            for i, j, converter in [(1, 8, str), (6, 9, float), (7, 10, float)]:
                if converter(ls[i]) != converter(ls[j]):
                    perturbed = True
            if perturbed and resid != targetresid:
                message = f"""Residue {resid}-{ls[3]} atom {atomid} name {ls[4]} is perturbed but this residue is not the target of mutation.
This is presumably error, and this is likely to happen when His resiudes have different states before and after the mutation and subsequent structure optimization.
If you want to proceed, change all His residues in your system into appropriate names (HID/HIE for amber and HSD/HSE for charmm) to forcefully fix the protonation state."""
                raise RuntimeError(message)
            continue
        elif moleculename == "merged" and sectionname in numatoms:
            na = numatoms[sectionname]
            violate = False
            for i in range(na):
                aid = int(ls[i])
                if aid not in atomidmap:
                    violate = True
                    break
                ls[i] = atomidmap[aid]
            if violate:
                continue
            writels()
            continue
        ofh.write(l)

with open(pdbfrom) as fh, open(pdbto, "w") as ofh:
    for l in fh:
        if l.startswith("ATOM  "):
            resid = l[22:26].strip()
            if int(resid) in targetresids:
                ofh.write(l)




