#!/usr/bin/env python

import pymbar
import sys
import re
import argparse
import numpy as np
import os.path
import pickle

# in kJ/mol/K (to fit GROMACSy output)
gasconstant = 0.008314472
# in kcal/mol/K 
gasconstant_kcal = 0.0019872036


parser = argparse.ArgumentParser(description = 'Convergence tester')
parser.add_argument('--xvgs', metavar=".xvg", type = str, required=True,
                    help = 'xvg file path, "%sim", "%eval", and "%part" will be replaced by appropriate numbers')
#parser.add_argument('--pvs', metavar=".xvg", type = str, required=True,
#                    help = 'pV value files. "%sim" and "%part" will be replaced by appropriate numbers')
parser.add_argument('--nsim', metavar="N", type = int, required=True,
                    help = 'number of simulations')
parser.add_argument('--neval', metavar="N", type = int, required=True,
                    help = 'number of evaluation potentials')
parser.add_argument('--minpart', metavar="N", type = int, required=True,
                    help = 'part number begin')
parser.add_argument('--maxpart', metavar="N", type = int, required=True,
                    help = 'part number end')
parser.add_argument('--temp', help="Temperature (K)", type = float, default=300.0)
parser.add_argument('--save-dir', help="save result to this directory", type = str, default = None)
parser.add_argument('--subsample', help="subsample interval", type = int, default = 1)
parser.add_argument('--split', help="Number of chunks", type = int, default = 10)

opts = parser.parse_args()

# FIXME TODO: mbar.py seems to store the energy by [L][K][N]???
# ==> see mbar.py, correct representation is indeed [L][K][N].
# L: observation states
# K: evaluation states (in this case, and usually, L == K)
# N: samples (N_K)
betas = []
# energies will be represented as energies[K][L][N] <= no, due to mbar's rotten documentation
energies = {}
floatpat = r'[+-]?(?:\d+(\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
nsamples = {}
time_all = {}

def parse_xvgs(fs, subsample):
    data = []
    tprev = -1
    times = []
    for f in fs:
        with open(f) as fh:
            samplecount = 0
            for l in fh:
                if len(l) == 0:
                    continue
                if l[0] in ['#', '@']:
                    pass
                else:
                    ls = l.split()
                    tt = float(ls[0])
                    if tprev == tt:
                        # prevent double counting
                        continue
                    tprev = tt
                    times.append(tt)
                    if(samplecount % subsample == 0):
                        data.append(float(ls[1]))
                    samplecount += 1
    return np.array(times), np.array(data)

tprev = None
for isim in range(opts.nsim):
    #pvfiles = []
    #for part in range(opts.minpart, opts.maxpart + 1):
        #pvfiles.append(opts.pvs.replace("%sim", str(isim)).replace("%part", "%04d"%part))
    #times_, pvs = parse_xvgs(pvfiles, opts.subsample)
   
    beta = 1. / (gasconstant * opts.temp)
    for ieval in range(opts.neval):
        if abs(ieval - isim) > 1:
            continue
        files = []
        for part in range(opts.minpart, opts.maxpart + 1):
            f = opts.xvgs.replace("%sim", str(isim)).replace("%eval", str(ieval)).replace("%part", "%04d"%part)
            files.append(f)
        times, data = parse_xvgs(files, opts.subsample)
        energies[(isim, ieval)] = beta * data
        print("Finished loading %s with %d frames" % (", ".join(files), len(data)))
        if ieval == isim:
            time_all[isim] = times
            nsamples[isim] = len(data)
        sys.stdout.flush()

# Do simple bar with 2x2 matrix
print("Performing sectioned BAR")
def bar(emat, time_all, btime, etime):
    dgtot = 0.0
    for isim in range(opts.nsim - 1):
        mask_isim = np.logical_and(time_all[isim] > btime, time_all[isim] <= etime)
        mask_isim1 = np.logical_and(time_all[isim + 1] > btime, time_all[isim + 1] <= etime)
        #print("DEBUG", emat[(isim, isim)].shape)
        u = np.array([[emat[(isim, isim)][mask_isim], emat[(isim, isim + 1)][mask_isim]],
                      [emat[(isim + 1, isim)][mask_isim1], emat[(isim + 1, isim + 1)][mask_isim1]]], np.float64)
        nk = np.array([np.sum(mask_isim), np.sum(mask_isim1)])
        #print(u.shape, nk.shape)
        mb = pymbar.MBAR(u, nk, verbose=False)
        rettuple = mb.getFreeEnergyDifferences()
        #print("DEBUG", type(rettuple))
        #print("DEBUG len", len(rettuple))
        Deltaf = rettuple[0]
        #print("DEBUG", type(Deltaf))
        #print("DEBUG shape", Deltaf.shape)
        dgtot += Deltaf[0, 1] # F[isim + 1] - F[isim], fixed comment 2021 May 20th
    return dgtot # F[opts.nsim - 1] - F[0]

#print("BAR:", bar(energies, time_all, 2100, 4100) * gasconstant* opts.temp)
#print("BAR (100-4100):", bar(energies, time_all, 100, 4100) * gasconstant* opts.temp)
results = []
for i in range(opts.split):
    tmin = time_all[0][0]
    tmax = time_all[0][-1]
    twidth = (tmax - tmin) / opts.split
    tspan = (i + 1) * twidth
    t0 = tmin + 0.5 * tspan
    t1 = tmin + tspan
    barres = bar(energies, time_all, t0, t1) 
    print("BAR %f-%f:" % (t0, t1), barres * gasconstant* opts.temp)
    results.append((t0, t1, barres))

pickle.dump(results, open("%s/results.pickle" % opts.save_dir, "wb"))


