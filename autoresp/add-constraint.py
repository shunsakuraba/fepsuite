import sys
import pybel
import openbabel as ob
import os.path
import re

if len(sys.argv) != 4:
    print >> sys.stderr, "Usage: %s [Gaussian(Z-matrix)] [pdb (or mol2)] [Constraint file]" % sys.argv[0]
    sys.exit(1)

gaufile = sys.argv[1]
namefile = sys.argv[2]
constraintfile = sys.argv[3]
nametype = os.path.splitext(namefile)[1]
if nametype == "":
    print >> sys.stderr, "File extension must be set: \"%s\"" % namefile
    sys.exit(1)
if nametype[0] == '.':
    nametype = nametype[1:]

m = pybel.readfile(nametype, namefile).next()
if m == None:
    print >> sys.stderr, "Failed to read from %s" % namefile
    sys.exit(1)
obmol = m.OBMol

name2idx = {}
print >> sys.stderr, "Index, AtomID:"
for a in ob.OBMolAtomIter(obmol):
    atomid = a.GetResidue().GetAtomID(a).strip()
    index = a.GetIdx()
    print >> sys.stderr, index, atomid
    name2idx[atomid] = index

constraints = []
with open(constraintfile, "rt") as fh:
    for l in fh:
        if l == "":
            continue
        ls = l.split()
        if ls[0] == "dihedral":
            constraints.append((ls[0], 
                                name2idx[ls[1]],
                                name2idx[ls[2]],
                                name2idx[ls[3]],
                                name2idx[ls[4]],
                                float(ls[5])))
        else:
            print >> sys.stderr, "Unknown constraint type: \"%s\"" % l
            sys.exit(1)

with open(gaufile, "rt") as fh:
    for l in fh:
        if len(l) == 0 or l[0] != "#":
            sys.stdout.write(l)
            continue
        # "#" line
        l = re.sub(r"\bopt\b", "popt", l, flags = re.I)
        sys.stdout.write(l)
        break

    for l in fh:
        if l.strip() == "":
            sys.stdout.write("Geom=ModRedundant\n")
            sys.stdout.write("\n")
            break
        else:
            l = re.sub(r"\bopt\b", "popt", l, flags = re.I)
            sys.stdout.write(l)
            continue
    # comment line
    l = fh.next()
    sys.stdout.write(l)
    # blank
    l = fh.next()
    assert(l.strip() == "")
    sys.stdout.write(l)

    # charge and duplication
    l = fh.next()
    sys.stdout.write(l)
    
    # coordinates
    for l in fh:
        if l.strip() == "":
            sys.stdout.write(l)
            break
        if len(l.split()) != 4:
            print >> sys.stderr, "This script does not support Z-Matrix format"
            print >> sys.stderr, "caught problem at: \"%s\"" % l.strip()
            sys.exit(1)
        sys.stdout.write(l)

    # output target angle first 
    for c in constraints:
        if c[0] == "dihedral":
            sys.stdout.write("D %d %d %d %d %f\n" % 
                             (c[1], c[2], c[3], c[4], c[5]))

    # then constraints
    for c in constraints:
        if c[0] == "dihedral":
            sys.stdout.write("D %d %d %d %d F\n" % 
                             (c[1], c[2], c[3], c[4]))
    
    for l in fh:
        sys.stdout.write(l)


