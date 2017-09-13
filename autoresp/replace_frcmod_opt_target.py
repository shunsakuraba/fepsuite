#!/usr/bin/env python


import sys
import re



[_exe, opt, frcmod] = sys.argv

res = {}
with open(opt) as fh:
    for l in fh:
        [dup, ene, phase] = l.split()
        dup = int(dup)
        ene = float(ene)
        phase = float(phase)
        res[dup] = (ene, phase)

pat = re.compile(r"Opt target (\d)")
with open(frcmod) as fh:
    for l in fh:
        m = pat.search(l)
        if m:
            dup = int(m.group(1))
            (ene, phase) = res[dup]
            print ("%15s%9.3lf%14.3lf%16.3lf         %s" % 
                   (l[0:15], ene, phase, float(l[38:54].strip()), ""))
        else:
            sys.stdout.write(l)



