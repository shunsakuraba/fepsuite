import sys

[_, topname, ndxname] = sys.argv

moltable = {}
scaled_indices = []
tot_natom = 0
with open(topname) as fh, open(ndxname, "w") as ofh:
    section = None
    molname = None
    scaledatoms = []
    natom = 0
    for l in fh:
        if l.startswith("["):
            if section == "atoms":
                moltable[molname] = (natom, scaledatoms)
            section = l.split()[1]
            if section == "atoms":
                scaledatoms = []
                natom = 0
            continue
        if ';' in l:
            l = l.split(';')[0]
        if l.strip() == "":
            continue
        if section == "moleculetype":
            molname = l.split()[0]
            continue
        if section == "atoms":
            ls = l.split()
            aindex = int(ls[0])
            atype = ls[1]
            if atype.endswith("_") or (len(ls) >= 10 and (ls[8] != atype or float(ls[9]) != float(ls[6]))): # ls[8]: B at, ls[6]/[9]: chg
                scaledatoms.append(aindex - 1)
            natom += 1
        if section == "molecules":
            ls = l.split()
            molname = ls[0]
            molcount = int(ls[1])
            (natom, scaledatoms) = moltable[molname]
            for i in range(molcount):
                for c in scaledatoms:
                    scaled_indices.append(c + i * natom + tot_natom)
            tot_natom += natom * molcount
            continue
    print("[ System ]", file=ofh)
    for i in range(tot_natom):
        print(i + 1, file=ofh)
    print(file=ofh)
    print("[ hot ]", file=ofh)
    for i in scaled_indices:
        print(i + 1, file=ofh)
    print(file=ofh)


