#!/usr/bin/env python

import sys
import openbabel as ob
import pybel
import sys
import os.path
import math

if len(sys.argv) not in [5, 6]:
    print >> sys.stderr, "Usage: %s [uncharged mol] [base mol] [known partial charges] [output mol2]" % sys.argv[0]
    sys.exit(1)

integral = True
if sys.argv[1] == "-nonint":
    integral = False
    del(sys.argv[1])

bbfile = sys.argv[1]
bbtype = os.path.splitext(bbfile)[1]
bbtype = bbtype.lstrip('.')

basefile = sys.argv[2]
basetype = os.path.splitext(basefile)[1]
basetype = basetype.lstrip('.')

chargefile = sys.argv[3]

outputfile = sys.argv[4]
outputtype = os.path.splitext(outputfile)[1]
outputtype = outputtype.lstrip('.')

if bbtype == "":
    print >> sys.stderr, "File extension must be set: \"%s\"" % bbfile
    sys.exit(1)

if basetype == "":
    print >> sys.stderr, "File extension must be set: \"%s\"" % basefile
    sys.exit(1)

if bbtype == "pdb":
    print >> sys.stderr, "Warning: Uncharged mol file must not be pdb (must have charges, mol2 recommended); This feature is only useful to generate the charge list; continue running"
    #sys.exit(1)

if basetype == "pdb":
    print >> sys.stderr, "Base file must not be pdb (must have charges, mol2 recommended)"
    sys.exit(1)

bb = pybel.readfile(bbtype, bbfile).next()
if bb == None:
    print >> sys.stderr, "Failed to read from %s" % bbfile
    sys.exit(1)
bbmol = bb.OBMol

base = pybel.readfile(basetype, basefile).next()
if bb == None:
    print >> sys.stderr, "Failed to read from %s" % basefile
    sys.exit(1)
basemol = base.OBMol

chargemap = {}
with open(chargefile, "rt") as fh:
    for l in fh:
        ls = l.split()
        chargemap[ls[0]] = float(ls[1])

def atomid(atom):
    return atom.GetResidue().GetAtomID(atom).strip()

chargemap_base = {}
for a in ob.OBMolAtomIter(basemol):
    aid = atomid(a)
    chargemap_base[aid] = a.GetPartialCharge()

chargemap_used = []
bbmol.SetAutomaticPartialCharge(False)
for a in ob.OBMolAtomIter(bbmol):
    aid = atomid(a)
    if aid in chargemap:
        if aid in chargemap_used:
            print >> sys.stderr, "The same atom name (%s) used more than once" % aid
            sys.exit(1)
        chargemap_used.append(aid)
        c = chargemap[aid]
        a.SetPartialCharge(c)
    else:
        # copy from chargemap_base
        if aid not in chargemap_base:
            print >> sys.stderr, "Atom %s does not found in base mol" % aid
            sys.exit(1)
        c = chargemap_base[aid]
        a.SetPartialCharge(c)
for aid in chargemap:
    if aid not in chargemap_used:
        print >> sys.stderr, "Predefined atom %s does not appear" % aid
        sys.exit(1)

if outputtype == "charge":
    with open(outputfile, "wt") as fh:
        for a in ob.OBMolAtomIter(bbmol):
            print >> fh, a.GetPartialCharge()
else:
    output = pybel.Outputfile(outputtype, outputfile, overwrite=True)
    output.write(pybel.Molecule(bbmol))

total = 0.0
for a in ob.OBMolAtomIter(bbmol):
    total += a.GetPartialCharge()
print >> sys.stderr, "Total charges = %.6f" % total

if integral and abs(total - round(total)) > 0.05:
    print >> sys.stderr, "Total charge is not integral"
    sys.exit(1)




