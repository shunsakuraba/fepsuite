#!/usr/bin/env python

import sys
import prepi

if len(sys.argv) < 4:
    print "Usage: %s  (prep) (dihedralspec) (newatomtype)" % sys.argv[0]
    print "example: %s ionsine.prep O4\'-C1\'-N9-C8 CP" % sys.argv[0]
    sys.exit(1)

prep = prepi.prepi(sys.argv[1])

diheds = sys.argv[2].split('-')
replace_atom = diheds[3]
diheds_inds = [prep.atommap[a] for a in diheds]
diheds_ats = [prep.atoms[i][1] for i in diheds_inds]
original_atomtype = diheds_ats[3]
new_atomtype = sys.argv[3]

with open(sys.argv[1], "r") as fh:
    for l in fh:
        l = l.rstrip("\n")
        if l[6:10] == "%-4s" % diheds[3]:
            print "%s%-2s%s" % (l[0:12], new_atomtype, l[14:])
        else:
            print l
