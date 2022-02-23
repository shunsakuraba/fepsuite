#!/usr/bin/env python

import pymbar
import sys
import re
import argparse
import numpy
import os.path
import pickle
import collections

SimEval = collections.namedtuple("SimEval", ["sim", "eval"])

# in kJ/mol/K (to fit GROMACSy output)
gasconstant = 0.008314472
# in kcal/mol/K 
gasconstant_kcal = 0.0019872036


def parse_args():
    parser = argparse.ArgumentParser(description = 'Convergence tester')
    parser.add_argument('--xvgs', metavar=".xvg", type = str, required=True,
                        help = 'xvg file path, "%sim", "%eval", and "%part" will be replaced by appropriate numbers')
    #parser.add_argument('--pvs', metavar=".xvg", type = str, required=True,
    #                    help = 'pV value files. "%sim" and "%part" will be replaced by appropriate numbers')
    parser.add_argument('--nsim', metavar="N", type = int, required=True,
                        help = 'number of simulations')
    parser.add_argument('--minpart', metavar="N", type = int, required=True,
                        help = 'part number begin')
    parser.add_argument('--maxpart', metavar="N", type = int, required=True,
                        help = 'part number end')
    parser.add_argument('--temp', help="Temperature (K)", type = float, default=300.0)
    parser.add_argument('--save-dir', help="save result to this directory", type = str, default = os.getcwd())
    parser.add_argument('--subsample', help="subsample interval", type = int, default = 1)
    parser.add_argument('--split', help="Number of chunks", type = int, default = 10)
    parser.add_argument('--show-intermediate', help="Show intermediate cumsum", action="store_true")

    opts = parser.parse_args()
    return opts

floatpat = r'[+-]?(?:\d+(\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'

def parse_deltae(fs, subsample, simindex):
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
                    eval_pot_pair = [(int(ls[i]), float(ls[1 + i])) for i in range(1, len(ls), 2)]
                    if(samplecount % subsample == 0):
                        for (evix, evpot) in eval_pot_pair:
                            data.append((tt, evix, evpot))
                    samplecount += 1
    return data

def bar(emat, time_all, nsim, btime, etime, show_intermediate):
    dgtot = 0.0
    for isim in range(nsim - 1):
        basestate = SimEval(sim=isim, eval=isim)
        #print(time_all[basestate])
        mask_isim = numpy.logical_and(time_all[basestate] > btime, time_all[basestate] <= etime)
        nmasked = numpy.sum(mask_isim)
        assert nmasked > 0
        assert len(emat[basestate]) == len(emat[SimEval(sim=isim, eval=isim + 1)])
        assert len(emat[basestate]) == len(emat[SimEval(sim=isim + 1, eval=isim + 1)])
        assert len(emat[basestate]) == len(emat[SimEval(sim=isim + 1, eval=isim + 1)])
        # uses u_kn representation
        # K: evaluation states
        # N: samples (N_K)
        u = numpy.empty((2, nmasked * 2))
        # evaluated by isim
        u[0, 0:nmasked] = emat[basestate][mask_isim]
        u[0, nmasked:nmasked*2] = emat[SimEval(sim=isim+1, eval=isim)][mask_isim]
        u[1, 0:nmasked] = emat[SimEval(sim=isim, eval=isim+1)][mask_isim]
        u[1, nmasked:nmasked*2] = emat[SimEval(sim=isim+1, eval=isim+1)][mask_isim]

        nk = numpy.array([nmasked, nmasked])


        # print("debug shape:", u.shape, nk.shape)
        mb = pymbar.MBAR(u, nk, verbose=False)
        rettuple = mb.getFreeEnergyDifferences(compute_uncertainty=False, warning_cutoff=1)
        #print("DEBUG", type(rettuple))
        #print("DEBUG len", len(rettuple))
        Deltaf = rettuple[0]
        #print("DEBUG", type(Deltaf))
        #print("DEBUG shape", Deltaf.shape)
        dgtot += Deltaf[0, 1] # F[isim + 1] - F[isim], fixed comment 2021 May 20th
        if show_intermediate:
            print("INT", isim, isim+1, Deltaf[0,1], dgtot)
    return dgtot # F[opts.nsim - 1] - F[0]

def main():
    opts = parse_args()

    betas = []
    energies = {}
    nsamples = {}
    time_all = {}
    beta = 1. / (gasconstant * opts.temp)

    for isim in range(opts.nsim):
        files = []
        for part in range(opts.minpart, opts.maxpart + 1):
            f = opts.xvgs.replace("%sim", str(isim)).replace("%part", "%04d"%part)
            files.append(f)
        data = parse_deltae(files, opts.subsample, isim)
        
        for (t, st, energy) in data:
            se = SimEval(sim=isim, eval=st)
            if se not in energies:
                energies[se] = []
                time_all[se] = []
            if len(time_all[se]) > 0 and time_all[se][-1] == t:
                # dup frame
                continue
            energies[se].append(energy)
            time_all[se].append(t)

        nsamples[isim] = len(time_all[SimEval(sim=isim, eval=isim)])
        print("Finished loading %s with %d points %d frames" % (", ".join(files), len(data), len(time_all[se]))) # time_all[se] shall exist
        sys.stdout.flush()

    # is this legit?
    #for k in energies:
    for k in list(energies.keys()):
        energies[k] = numpy.array(energies[k]) * beta
        time_all[k] = numpy.array(time_all[k])

    # Do simple bar with 2x2 matrix
    print("Performing sectioned BAR")

    results = []
    for i in range(opts.split):
        tmin = time_all[SimEval(0,0)][0]
        tmax = time_all[SimEval(0,0)][-1]
        twidth = (tmax - tmin) / opts.split
        tspan = (i + 1) * twidth
        t0 = tmin + 0.5 * tspan
        t1 = tmin + tspan
        barres = bar(energies, time_all, opts.nsim, t0, t1, opts.show_intermediate) 
        print("BAR %f-%f:" % (t0, t1), barres * gasconstant * opts.temp)
        results.append((t0, t1, barres))

    pickle.dump(results, open("%s/results.pickle" % opts.save_dir, "wb"))

main()

