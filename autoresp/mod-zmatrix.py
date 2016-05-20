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
            angle = ls[5]
            if angle == "asis":
                pass
            else:
                angle = float(angle)
            constraints.append((ls[0], 
                                name2idx[ls[1]],
                                name2idx[ls[2]],
                                name2idx[ls[3]],
                                name2idx[ls[4]],
                                angle))
        else:
            print >> sys.stderr, "Unknown constraint type: \"%s\"" % l
            sys.exit(1)

with open(gaufile, "rt") as fh:
    for l in fh:
        if len(l) == 0 or l[0] != "#":
            sys.stdout.write(l)
            continue
        # "#" line
        l = re.sub(r"\bopt\b", "popt(Z-matrix)", l, flags = re.I)
        sys.stdout.write(l)
        break

    for l in fh:
        if l.strip() == "":
            sys.stdout.write("\n")
            break
        else:
            l = re.sub(r"\bopt\b", "popt(Z-matrix)", l, flags = re.I)
            sys.stdout.write(l)
            continue
    # comment line
    l = fh.next()
    sys.stdout.write(l)
    # blank
    l = fh.next()
    assert(l.strip() == "")
    sys.stdout.write(l)

    # charge and multiplicity
    l = fh.next()
    sys.stdout.write(l)
    
    # coordinates
    dihmaps = {}
    aid = 1
    for l in fh:
        ls = l.split()
        if len(ls[0]) > 2 or l.strip() == "":
            sys.stdout.write(l)
            break
        if len(ls) not in [1, 3, 5, 7]:
            print >> sys.stderr, "This script does not support Cartesian format"
            print >> sys.stderr, "caught problem at: \"%s\"" % l.strip()
            sys.exit(1)
        if len(ls) == 7:
            dihmaps[(int(ls[5]), int(ls[3]), int(ls[1]), aid)] = ls[6]
        aid += 1
        sys.stdout.write(l)

    fixvars = {}
    ndihconstraints = 0
    for c in constraints:
        if c[0] == "dihedral":
            if (c[1], c[2], c[3], c[4]) in dihmaps:
                fixvars[dihmaps[(c[1], c[2], c[3], c[4])]] = c[5]
            if (c[4], c[3], c[2], c[1]) in dihmaps:
                if c[5] == "asis":
                    fixvars[dihmaps[(c[1], c[2], c[3], c[4])]] = c[5]
                else:
                    fixvars[dihmaps[(c[4], c[3], c[2], c[1])]] = -c[5]
            ndihconstraints += 1
    
    if len(fixvars) != ndihconstraints:
        print >> sys.stderr, "Not all constraints are covered!"
        sys.exit(1)

    for l in fh:
        ls = l.split('=')
        if ls == []:
            sys.stdout.write(l)
            break
        if ls[0].strip() in fixvars:
            if fixvars == "asis":
                sys.stdout.write("%s= %12.5f F\n" % (ls[0], ls[1]))
            else:
                sys.stdout.write("%s= %12.5f F\n" % (ls[0], fixvars[ls[0]]))
            continue
        else:
            sys.stdout.write(l)

    for l in fh:
        sys.stdout.write(l)


