import sys
import pybel
import openbabel as ob
import os.path
import re
import math

# constrain angles for turbomole

if len(sys.argv) not in [3,4]:
    print >> sys.stderr, "Usage: %s (structure) [pdb (or mol2) for naming] (Dihedral connected by '-'s)" % sys.argv[0]
    sys.exit(1)

basefile = sys.argv[1]
dihedral = sys.argv[-1]
if len(sys.argv) == 3:
    namefile = basefile
else:
    namefile = sys.argv[2]

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
#print >> sys.stderr, "Index, AtomID:"
for a in ob.OBMolAtomIter(obmol_name):
    atomid = a.GetResidue().GetAtomID(a).strip()
    index = a.GetIdx()
    #print >> sys.stderr, index, atomid
    # 1-origin
    name2idx[atomid] = index

def getangle(i, j, k):
    return obmol.GetAngle(obmol.GetAtom(i), obmol.GetAtom(j), obmol.GetAtom(k))

def getlength(i, j):
    return obmol.GetBond(obmol.GetAtom(i), obmol.GetAtom(j)).GetLength()

# Read torsion constraints and modify structure if necessary
[a1,a2,a3,a4] = [name2idx[x] for x in dihedral.split('-')]
angle = obmol.GetTorsion(a1, a2, a3, a4)
print "%.8f" % angle
