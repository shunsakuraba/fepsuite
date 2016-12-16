#!/usr/bin/python

# CCSD(T)/CBS calculation to reproduce Zgarbova 

import re
import sys
import numpy as np
import math

fortran_float_re = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?'

def extract_multi(fname, patstr):
    pat = re.compile(patstr)
    with open(fname, 'r') as fh:
        for l in fh:
            m = pat.search(l)
            if m:
                matches = m.groups()
                rets = [float(matched.replace('D', 'E')) for matched in matches]
                return rets

def extract(fname, patstr):
    return extract_multi(fname, patstr)[0]

# get various data
trunk = sys.argv[1]
if trunk.find("%s") == -1:
    print >> sys.stderr, "Example: {} foo.%s.log"
    sys.exit(1)

# um, after processing file I realized I don't need mp2pvdz ...
fccsd = trunk % "ccsdpvdz"

eccsd = extract(fccsd, r"CCSD\(T\)\s*=\s*({})".format(fortran_float_re))
emp2 = extract(fccsd, r"EUMP2\s*=\s*({})".format(fortran_float_re))
deccsd = eccsd - emp2

# Get RHF and E2 for cbs
fpvtz = trunk % "mp2pvtz"
fpvqz = trunk % "mp2pvqz"

rhf_pat = r"E\(RHF\)\s*=\s*(%s)" % fortran_float_re
erhf3 = extract(fpvtz, rhf_pat)
erhf4 = extract(fpvqz, rhf_pat)

e2_pat = r"E2\s*=\s*({0})\s+EUMP2\s*=\s*({0})".format(fortran_float_re)
ecorr3 = extract_multi(fpvtz, e2_pat)[0]
ecorr4 = extract_multi(fpvqz, e2_pat)[0]

# RHF CBS: Halkier et al. CPL 302 437 (1999)
# Uses exponential form:
# E_x = E_inf + B exp (- alpha x)
# alpha = 1.54 is recommended for pvtz-pvqz two-point fitting
alpha = 1.54
m = np.array([[1, math.exp(-3 * alpha)],[1, math.exp(-4*alpha)]])
solution = np.linalg.solve(m, np.array([erhf3, erhf4]))
erhf_inf = solution[0]
B = solution[1]

# Corr CBS: Halkier et al. CPL 286 243 (1998)
# Uses polynomial form:
# E_x = E_inf + A X^-3
m = np.array([[1, math.exp(-3 * alpha)],[1, math.exp(-4*alpha)]])
solution = np.linalg.solve(m, np.array([ecorr3, ecorr4]))
ecorr_inf = solution[0]
A = solution[1]

#print deccsd
#print erhf3, erhf4, erhf_inf
#print ecorr3, ecorr4, ecorr_inf
print erhf_inf + ecorr_inf + deccsd
