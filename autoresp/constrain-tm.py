import sys
import pybel
import openbabel as ob
import os.path
import re
import math

# constrain angles for turbomole

if len(sys.argv) != 6:
    print >> sys.stderr, "Usage: %s [structure] [pdb (or mol2) for naming] [Constraint file] [tmole output] [dihedral spec output]" % sys.argv[0]
    sys.exit(1)

basefile = sys.argv[1]
namefile = sys.argv[2]
constraintfile = sys.argv[3]
tmoleout = sys.argv[4]
constraintout = sys.argv[5]
def typeoffn(f):
    ftype = os.path.splitext(f)[1]
    if ftype == "":
        print >> sys.stderr, "File extension must be set: \"%s\"" % f
        sys.exit(1)
    if ftype[0] == '.':
        ftype = ftype[1:]
    if ftype == "log":
        ftype = "g09"
    return ftype

basetype = typeoffn(basefile)
nametype = typeoffn(namefile)

m = pybel.readfile(basetype, basefile).next()
if m == None:
    print >> sys.stderr, "Failed to read from %s" % basefile
    sys.exit(1)
obmol = m.OBMol
m_name = pybel.readfile(nametype, namefile).next()
if m == None:
    print >> sys.stderr, "Failed to read from %s" % namefile
    sys.exit(1)
obmol_name = m_name.OBMol

name2idx = {}
print >> sys.stderr, "Index, AtomID:"
for a in ob.OBMolAtomIter(obmol_name):
    atomid = a.GetResidue().GetAtomID(a).strip()
    index = a.GetIdx()
    print >> sys.stderr, index, atomid
    # 1-origin
    name2idx[atomid] = index

def getangle(i, j, k):
    return obmol.GetAngle(obmol.GetAtom(i), obmol.GetAtom(j), obmol.GetAtom(k))

def getlength(i, j):
    return obmol.GetBond(obmol.GetAtom(i), obmol.GetAtom(j)).GetLength()

# Read torsion constraints and modify structure if necessary

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
                angle = float(angle)
                obmol.SetTorsion(a1, a2, a3, a4, angle * 2 * math.pi / 360)
            constraints.append((ls[0],
                                a1, a2, a3, a4,
                                angle))
        else:
            print >> sys.stderr, "Unknown constraint type: \"%s\"" % l
            sys.exit(1)

# output tmol format
obconversion = ob.OBConversion()
obconversion.SetInAndOutFormats("tmol", "tmol")
tmolstr = obconversion.WriteString(obmol)
tmolstr = tmolstr.replace("ho", "h")
with open(tmoleout, 'w') as fh:
    fh.write(tmolstr)

# output dihedral constraint
with open(constraintout, 'w') as fh:
    for c in constraints:
        if c[0] == "dihedral":
            fh.write("f tors %3d %3d %3d %3d\n" %
                     (c[1], c[2], c[3], c[4]))
