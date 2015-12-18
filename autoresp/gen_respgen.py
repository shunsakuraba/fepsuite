import pybel
import openbabel as ob
import sys
import os.path
import math

if len(sys.argv) != 3:
    print >> sys.stderr, "Usage: %s [Base mol] [RESP condition]" % sys.argv[0]
    sys.exit(1)

basefile = sys.argv[1]
basetype = os.path.splitext(basefile)[1]
basetype = basetype.lstrip('.')
respcondfile = sys.argv[2]

if basetype == "":
    print >> sys.stderr, "File extension must be set: \"%s\"" % basefile
    sys.exit(1)

b = pybel.readfile(basetype, basefile).next()
if b == None:
    print >> sys.stderr, "Failed to read from %s" % basefile
    sys.exit(1)
bmol = b.OBMol

# load restraint conditions
restraints = []
with open(respcondfile, "rt") as fh:
    for l in fh:
        l = l.strip()
        if len(l) == 0 or l[0] == '#':
            continue
        ls = l.split()
        charge = float(ls[0])
        atoms = ls[1:]
        restraints.append((charge, atoms))

def atomid(atom):
    return atom.GetResidue().GetAtomID(atom).strip()

name2idx = {}
for a in ob.OBMolAtomIter(backmol):
    name2idx[atomid(a)] = a.GetIdx()

for (charge, atoms) in restraints:
    if len(atoms) == 1:
        atom = atoms[0]
        print "CHARGE %12.6f %4d %5s" % (charge, name2idx[atom], atom)
    else:
        print "GROUP %4d %12.6f" % (len(atom), charge)
        for atom in atoms:
            print "ATOM %4d %5s" % (name2idx[atom], atom)

