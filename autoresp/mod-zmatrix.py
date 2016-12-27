import sys
import pybel
import openbabel as ob
import os.path
import re
from collections import OrderedDict

# XXX: quick hack
redundant = False
if len(sys.argv) >= 2 and sys.argv[1] == '--redundant':
    redundant = True
    # whoa
    del sys.argv[1]

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

def getangle(i, j, k):
    return obmol.GetAngle(obmol.GetAtom(i), obmol.GetAtom(j), obmol.GetAtom(k))

def getlength(i, j):
    return obmol.GetBond(obmol.GetAtom(i), obmol.GetAtom(j)).GetLength()

constraints = []
with open(constraintfile, "rt") as fh:
    for l in fh:
        if l == "":
            continue
        ls = l.split()
        if ls[0] == "dihedral":
            # normalize to x-y-z-w where ind(w) > ind(x)
            [a1, a2, a3, a4] = [name2idx[ls[i]] for i in range(1, 5)]
            sign = 1
            if a1 > a4:
                sign = -1
                ls[1:5] = ls[4:0:-1]
            if (a2 > a4 or a3 > a4):
                print >> sys.stderr, "Cannot assign dihedral to the current constraint (%s)" % l.rstrip()
                sys.exit(1)
            angle = ls[5]
            if angle == "asis":
                angle = obmol.GetTorsion(a1, a2, a3, a4)
                print >> sys.stderr, "ASIS angle of %s-%s-%s-%s from model: %f " % (ls[1], ls[2], ls[3], ls[4], angle)
                pass
            else:
                angle = float(angle) * sign
            constraints.append((ls[0],
                                a1, a2, a3, a4,
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
        if not redundant:
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
    dihmaps = OrderedDict()
    fixvars = OrderedDict()
    modvars = {}
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
            # reorder it if necessary
            ds = tuple([int(ls[x]) for x in [5, 3, 1]] + [aid])
            reorder_dihedral = False
            for c in constraints:
                if c[0] == "dihedral":
                    if c[4] == aid:
                        newds = (c[1], c[2], c[3], c[4])
                        if newds == ds:
                            continue
                        newl = ("%s  %d  %s  %d  %s  %d  %s\n" % 
                             (ls[0], c[3], ls[2], c[2], ls[4], c[1], ls[6]))
                        fixvars[ls[6]] = obmol.GetTorsion(c[1], c[2], c[3], c[4])
                        modvars[ls[4]] = getangle(c[2], c[3], c[4])
                        modvars[ls[2]] = getlength(c[3], c[4])
                        reorder_dihedral = True
            if reorder_dihedral:
                print >> sys.stderr, "Dihedral reordered: \"%s\" to \"%s\"" % (l.strip(), newl.strip())
                l = newl
            dihmaps[ds] = ls[6]

        aid += 1
        sys.stdout.write(l)

    ndihconstraints = 0
    for c in constraints:
        if c[0] == "dihedral":
            if (c[1], c[2], c[3], c[4]) in dihmaps:
                fixvars[dihmaps[(c[1], c[2], c[3], c[4])]] = c[5]
            ndihconstraints += 1
    
    if len(fixvars) != ndihconstraints:
        print >> sys.stderr, "Not all constraints are covered!"
        sys.exit(1)

    lastislf = False
    for l in fh:
        ls = l.split('=')
        if ls == []:
            sys.stdout.write(l)
            lastislf = True
            break
        variable = ls[0].strip()
        if variable in fixvars:
            sys.stdout.write("%s= %12.5f%s\n" % (ls[0], fixvars[variable], [" F", ""][redundant]))
            continue
        elif variable in modvars:
            sys.stdout.write("%s= %12.5f\n" % (ls[0], modvars[variable]))
        else:
            sys.stdout.write(l)
            lastislf = (l.strip() == "")

    for l in fh:
        sys.stdout.write(l)
        lastislf = (l.strip() == "")

    if redundant:
        if (not lastislf):
            sys.stdout.write("\n") 
        for (k, v) in dihmaps.iteritems():
            if v in fixvars:
                (a1, a2, a3, a4) = k
                sys.stdout.write("D %d %d %d %d F\n" % (a1, a2, a3, a4))
