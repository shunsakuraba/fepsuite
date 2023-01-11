#!/usr/bin/env python

import pymbar
import sys
import re
import argparse
import numpy as np
import os.path

# in kJ/mol/K (to fit GROMACSy output)
gasconstant = 0.008314472

parser = argparse.ArgumentParser(description = 'Debugged caller for pymbar')
parser.add_argument('xvgs', metavar=".xvg", type = str, nargs = '+',
                    help = 'xvg files')
#parser.add_argument('--temp', dest="temperature", type = float)
parser.add_argument('--save-dir', help="save result to this directory", type = str, default = None)
parser.add_argument('--begin', help="begin time (ps)", type = float, default = 0.0)
parser.add_argument('--end', help="end time (ps)", type = float, default = 0.0)

opts = parser.parse_args()

# L: observation states
# K: evaluation states (in this case, and usually, L == K)
# N: samples (N_K)
betas = []
# energies will be represented as energies[K][L][N]
energies = []
floatpat = r'[+-]?(?:\d+(\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
nsamples = []
begintime = opts.begin
endtime = opts.end
if endtime == 0.0:
    endtime = 1e+100
states = None

for f in opts.xvgs:
    headers = []
    # ordered in data[N][K]
    data = []
    with open(f, "rt") as fh:
        for l in fh:
            if len(l) == 0:
                continue
            if l[0] in ['#', '@']:
                headers.append(l)
            else:
                ls = l.split()
                tt = float(ls[0])
                if tt < begintime:
                    continue
                if tt >= endtime:
                    # quick hack
                    break
                data.append([float(x) for x in ls])

    # get temperature
    pat = re.compile("T = (%s) \(K\)" % floatpat)
    beta = None
    for l in headers:
        m = pat.search(l)
        if m:
            temperature = float(m.group(1))
            beta = 1.0 / (gasconstant * temperature)
            betas.append(beta)
            assert(betas[0] == betas[-1])
    # convert potentials into u_k
    enetemp = []
    nsample = 0
    for vs in data:
        t = vs[0]
        if t < begintime or t >= endtime:
            continue
        pv = vs[-1]
        es = vs[1:-1]
        es = [beta * (e + pv) for e in es]
        enetemp.append(es)
        nsample += 1
        if states == None:
            states = len(es)
        else:
            assert(states == len(es))
    nsamples.append(nsample)
    if energies == []:
        energies = [[] for _ in range(states)]
    for k in range(states):
        energies[k].append([])
        for i in range(nsample):
            energies[k][-1].append(enetemp[i][k])
    print "Finished loading %s" % f
    sys.stdout.flush()

print "Converting from list to np array"
sys.stdout.flush()
u = np.array(energies, np.float64)
nk = np.array(nsamples, np.int)
print "Initializing MBAR"
sys.stdout.flush()
mb = pymbar.MBAR(u, nk, verbose=True)
print "Computing dFE"
sys.stdout.flush()
Deltaf_ij, dDeltaf_ij, Theta_ij = mb.getFreeEnergyDifferences()
print "Completed."
sys.stdout.flush()

kbt = 1.0 / betas[0]
print "Delta F (kJ/mol) = "
for i in range(states):
    print Deltaf_ij[i, 0] * kbt

# FIXME: this might be squared value, so kbt**2 might be appropriate.
# it is undocumented in pymbar and thus we need to read the source
print "dDelta F matrix (kJ/mol) = "
for i in range(states):
    for j in range(states):
        print dDeltaf_ij[i, j] * kbt,
    print
print "Theta matrix = "
print Theta_ij * kbt


if opts.save_dir:
    np.save(os.path.join(opts.save_dir, "Deltaf_ij.npy"), Deltaf_ij)
    np.save(os.path.join(opts.save_dir, "dDeltaf_ij.npy"), dDeltaf_ij)
    np.save(os.path.join(opts.save_dir, "Theta_ij.npy"), Theta_ij)

[Delta_f_ij, dDelta_f_ij, Delta_u_ij, dDelta_u_ij, Delta_s_ij, dDelta_s_ij] = mb.computeEntropyAndEnthalpy(verbose=True)

print "Delta H (kJ/mol) = "
for i in range(states):
    print Delta_u_ij[i, 0] * kbt
print "Delta S (kJ/mol/K) = "
for i in range(states):
    print Delta_s_ij[i, 0] * gasconstant

if opts.save_dir:
    np.save(os.path.join(opts.save_dir, "Deltau_ij.npy"), Delta_u_ij)
    np.save(os.path.join(opts.save_dir, "dDeltau_ij.npy"), dDelta_u_ij)
    np.save(os.path.join(opts.save_dir, "Deltas_ij.npy"), Delta_s_ij)
    np.save(os.path.join(opts.save_dir, "dDeltas_ij.npy"), dDelta_s_ij)

