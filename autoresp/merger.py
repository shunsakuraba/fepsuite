import pybel
import openbabel as ob
import sys
import os.path

if len(sys.argv) != 3:
    print >> sys.stderr, "Usage: %s [Base mol] [Backbone mol]" % sys.argv[0]
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


# This seems not necessary at all?
def repair_sprious_bond(mol):
    """Repair topology by removing sprious bond b/w H-C of H-O-C"""
    
    for a in ob.OBMolAtomIter(mol):
        if a.GetAtomicNum() == 1:
            # Hydrogen
            acrd = a.GetVector()
            lmin = 1e+5
            for b in ob.OBAtomAtomIter(a):
                bcrd = b.GetVector()
                diff = ob.vector3(acrd)
                diff -= bcrd
                lmin = min(lmin, diff.length_2())
            # This logic may look stupid...
            for bb in ob.OBAtomBondIter(a):
                if a.GetIdx() == bb.GetBeginAtomIdx():
                    bid = bb.GetEndAtomIdx()
                else:
                    bid = bb.GetBeginAtomIdx()
                b = mol.GetAtom(bid)
                bcrd = b.GetVector()
                diff = ob.vector3(acrd)
                diff -= bcrd
                if diff.length_2() != lmin:
                    print >> sys.stderr, "Removing bond ", a.GetResidue().GetAtomID(a).strip(), b.GetResidue().GetAtomID(b).strip()
                    mol.DeleteBond(bb)
                    break

repair_sprious_bond(bmol)
repair_sprious_bond(backmol)

def find_D_atoms(mol):
    """Find D atoms or D-starting atoms from mol"""

    for a in ob.OBMolAtomIter(mol):
        atomid = a.GetResidue().GetAtomID(a).strip()
        if a.GetIsotope() == 2 or atomid[0].upper() == 'D':
            d = a
            print >> sys.stderr, "Connector atom: ", a.GetIdx(), atomid
            p = ob.OBAtomAtomIter(a).next()
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


n = backmol.NumAtoms()
ap = backbone_p.GetIdx()
bp = base_p.GetIdx() + n
backmol += bmol

builder = ob.OBBuilder()
builder.Connect(backmol, ap, bp)

output = pybel.Outputfile('pdb', "debug.pdb", overwrite=True)
output.write(pybel.Molecule(backmol))





