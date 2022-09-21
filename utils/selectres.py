import sys



[_, topfrom, pdbfrom, targetresid, topto, pdbto] = sys.argv
targetresid = int(targetresid)
targetresids = [targetresid-1, targetresid, targetresid+1] # FIXME: termini


curatom = 0
atomidmap = {}
numatoms = {
        "bonds": 2,
        "angles": 3,
        "dihedrals": 4,
        "pairs": 2,
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
            resid = int(ls[2])
            atomid = int(ls[0])
            if resid in targetresids:
                curatom += 1
                atomidmap[atomid] = curatom
                ls[0] = curatom
                ls[5] = curatom
                writels()
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
            ls = l.split()
            if int(ls[5]) in targetresids:
                ofh.write(l)




