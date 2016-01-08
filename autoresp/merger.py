import pybel
import openbabel as ob
import sys
import os.path
import math

if len(sys.argv) != 4:
    print >> sys.stderr, "Usage: %s [Base mol] [Backbone mol] [output file]" % sys.argv[0]
    sys.exit(1)

basefile = sys.argv[1]
basetype = os.path.splitext(basefile)[1]
backbonefile = sys.argv[2]
backbonetype = os.path.splitext(backbonefile)[1]
basetype = basetype.lstrip('.')
backbonetype = backbonetype.lstrip('.')

if basetype == "":
    print >> sys.stderr, "File extension must be set: \"%s\"" % basefile
    sys.exit(1)
if backbonetype == "":
    print >> sys.stderr, "File extension must be set: \"%s\"" % backbonefile
    sys.exit(1)

b = pybel.readfile(basetype, basefile).next()
if b == None:
    print >> sys.stderr, "Failed to read from %s" % basefile
    sys.exit(1)
bmol = b.OBMol

back = pybel.readfile(backbonetype, backbonefile).next()
if back == None:
    print >> sys.stderr, "Failed to read from %s" % backbonefile
    sys.exit(1)
backmol = back.OBMol

def atomid(atom):
    return atom.GetResidue().GetAtomID(atom).strip()
def find_D_atoms(mol):
    """Find D atoms or D-starting atoms from mol"""

    for a in ob.OBMolAtomIter(mol):
        aid = atomid(a)
        if a.GetIsotope() == 2 or aid[0].upper() == 'D':
            d = a
            print >> sys.stderr, "Connector atom: ", a.GetIdx(), aid
            p = ob.OBAtomAtomIter(a).next()
            print >> sys.stderr, " Parent atom: ", p.GetIdx(), atomid(p)
            return (d, p)

(base_d, base_p) = find_D_atoms(bmol)
(backbone_d, backbone_p) = find_D_atoms(backmol)

dpos = backbone_d.GetVector()
ppos = backbone_p.GetVector()

# align base structure so that dpos matches base_p and
# ppos matches base_b. Note b<=>p are inverted.
bmol.Align(base_d, base_p, ppos, dpos)

# remove D atoms
bmol.DeleteAtom(base_d)
backmol.DeleteAtom(backbone_d)

# Connect two molecules
n = backmol.NumAtoms()
ap = backbone_p.GetIdx()
bp = base_p.GetIdx() + n
backmol += bmol

builder = ob.OBBuilder()
builder.Connect(backmol, ap, bp)

# Fix "HO" get misunderstood as Holmium (atom number 67)
for a in ob.OBMolAtomIter(backmol):
    if a.GetAtomicNum() == 67:
        a.SetAtomicNum(1)
    if a.GetAtomicNum() == 0:
        aid = atomid(a)
        if len(aid) >= 2 and aid[0:2] == "OP":
            # OP1 OP2
            a.SetAtomicNum(8)

# Get residue name
def get_residue_name(mol):
    a = ob.OBMolAtomIter(bmol).next()
    return a.GetResidue().GetName()

# Reset residue name
rn = get_residue_name(bmol)
for r in ob.OBResidueIter(backmol):
    r.SetName(rn)

# Fix atom names
for a in ob.OBMolAtomIter(backmol):
    a.GetResidue().SetAtomID(a, "%-4s" % a.GetResidue().GetAtomID(a).strip())

# Topology completed. Evaluate the vdW
ff = ob.OBForceField.FindForceField("UFF")
def eval_energy(mol):
    ff.Setup(mol)
    return ff.E_VDW()
def find_torsion_of_bond(mol, a, b):
    for dh in ob.OBMolTorsionIter(mol):
        if ((dh[1] + 1, dh[2] + 1) == (a, b) or
            (dh[1] + 1, dh[2] + 1) == (b, a)):
            return tuple([mol.GetAtom(dh[i] + 1) for i in range(4)])
            break
torsion = find_torsion_of_bond(backmol, ap, bp)
print >> sys.stderr, ("Scanning around %s-%s-%s-%s" %
                      tuple([atomid(a) for a in torsion]))
stepsize = 10
nround = 360 / stepsize
minenergy = 1e+9
minangle = -1
for i in range(nround):
    newangle = float(i * stepsize)
    backmol.SetTorsion(torsion[0],
                       torsion[1],
                       torsion[2],
                       torsion[3],
                       newangle * 2 * math.pi / 360)
    energy = eval_energy(backmol)
    print >> sys.stderr, "angle = %4f VDW = %f" % (newangle, energy)
    if energy < minenergy:
        minenergy = energy
        minangle = newangle
print >> sys.stderr, "--------"
print >> sys.stderr, "Minimum collision angle = %f" % minangle
backmol.SetTorsion(torsion[0],
                   torsion[1],
                   torsion[2],
                   torsion[3],
                   minangle * 2 * math.pi / 360)

outputfile = sys.argv[3]
outputtype = os.path.splitext(backbonefile)[1]
outputtype = outputtype.lstrip('.')

output = pybel.Outputfile(outputtype, outputfile, overwrite=True)
output.write(pybel.Molecule(backmol))





