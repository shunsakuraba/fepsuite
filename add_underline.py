
import mdtraj
import sys

[_, confname, topname, writeto] = sys.argv

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
protein = gro.topology.select("not (resname HOH NA CL Na Cl K)")
# neighbors are defined as within 4 aa
neighbors = mdtraj.compute_neighbors(gro, 0.4, query_indices=modified, haystack_indices=protein)
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
                print(l.rstrip(), file=ofh)
                continue
            elif ';' in l:
                (l, comment) = l.split(';', 1)
                comment = "; " + comment
            if l.strip() == "":
                pass
            elif section == "moleculetype":
                molname = l.split()[0]
            elif section == "atoms" and molname == "merged":
                ls = l.split()
                resno = int(ls[2])
                if resno in resids:
                    l = " ".join([ls[0], ls[1]+"_"] + ls[2:]) + "\n"
            print(l.rstrip() + comment, file=ofh)
