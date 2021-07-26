#!/usr/bin/env python

import pymbar
import sys
import re
import argparse
import numpy
import os.path
import pickle

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
                    basepot = None
                    for (evix, evpot) in eval_pot_pair:
                        if evix == simindex:
                            basepot = evpot
                    if(samplecount % subsample == 0):
                        for (evix, evpot) in eval_pot_pair:
                            if evix != simindex:
                                #data.append((tt, eval_pot_pair))
                                data.append((tt, evix, basepot, evpot))
                    samplecount += 1
    return data

def bar(emat, time_all, btime, etime):
    nsim = len(emat)
    dgtot = 0.0
    for isim in range(nsim - 1):
        mask_isim = numpy.logical_and(time_all[isim] > btime, time_all[isim] <= etime)
        #print("DEBUG", emat[(isim, isim)].shape)
        u = numpy.array([[emat[isim][(0, 0)][mask_isim], emat[isim][(0, 1)][mask_isim]],
                      [emat[isim][(1, 0)][mask_isim], emat[isim][(1, 1)][mask_isim]]], numpy.float64)
        nk = numpy.array([numpy.sum(mask_isim), numpy.sum(mask_isim)])

        #print(u.shape, nk.shape)
        mb = pymbar.MBAR(u, nk, verbose=False)
        rettuple = mb.getFreeEnergyDifferences(compute_uncertainty=False, warning_cutoff=1)
        #print("DEBUG", type(rettuple))
        #print("DEBUG len", len(rettuple))
        Deltaf = rettuple[0]
        #print("DEBUG", type(Deltaf))
        #print("DEBUG shape", Deltaf.shape)
        dgtot += Deltaf[0, 1] # F[isim + 1] - F[isim], fixed comment 2021 May 20th
    return dgtot # F[opts.nsim - 1] - F[0]

def main():
    opts = parse_args()

    # Correct representation is indeed [L][K][N]. (note some previous mbar documentation was wrong)
    # L: observation states
    # K: evaluation states (in this case, and usually, L == K)
    # N: samples (N_K)
    betas = []
    energies2x2 = {} # maps isim_base -> (isim - isimbase, ieval - isimbase) -> array(E)
    nsamples = {}
    time_all = {}

    for isim in range(opts.nsim):
        beta = 1. / (gasconstant * opts.temp)
        files = []
        for part in range(opts.minpart, opts.maxpart + 1):
            f = opts.xvgs.replace("%sim", str(isim)).replace("%part", "%04d"%part)
            files.append(f)
        data = parse_deltae(files, opts.subsample, isim)
        for ievalshift in [-1, 1]:
            ieval = isim + ievalshift
            selected = list(filter(lambda x: x[1] == ieval, data))
            if ievalshift == -1:
                eneix = isim - 1
                if eneix < 0: 
                    continue
            else:
                eneix = isim
                assert eneix not in energies2x2
                energies2x2[eneix] = {}
            # Note the order of ieval and isim, because this value is opponent's coordinate eval'd by isim
            energies2x2[eneix][(ieval - eneix, isim - eneix)] = beta * numpy.array([opponent - mye for (_t, _, mye, opponent) in selected])
            energies2x2[eneix][(isim - eneix, isim - eneix)] = beta * numpy.array([mye for (_t, _, mye, _opponent) in selected])

            if ievalshift == 1:
                time_all[isim] = numpy.array([t for (t, _, _, _) in selected])
        nsamples[isim] = len(time_all[isim])
        print("Finished loading %s with %d frames" % (", ".join(files), len(data)))
        sys.stdout.flush()
    # At this point energies2x2 contains all information to calculate dE

    # sanity checks
    for isim_base in range(opts.nsim - 1):
        assert isim_base in energies2x2
        assert (0, 0) in energies2x2[isim_base]
        assert (0, 1) in energies2x2[isim_base]
        assert (1, 0) in energies2x2[isim_base]
        assert (1, 1) in energies2x2[isim_base]
        #print("DEBUG:", len(time_all[isim_base]))
        #print(len(energies2x2[isim_base][(0,0)]))
        #print(len(energies2x2[isim_base][(0,1)]))
        #print(len(energies2x2[isim_base][(1,0)]))
        #print(len(energies2x2[isim_base][(1,1)]))

    # reconstruct matrix elements
    for isim_base in range(opts.nsim - 1):
        # Delta-E added to base evaluation
        energies2x2[isim_base][0, 1] += energies2x2[isim_base][0, 0]
        energies2x2[isim_base][1, 0] += energies2x2[isim_base][1, 1]
    
    # Do simple bar with 2x2 matrix
    print("Performing sectioned BAR")

    results = []
    for i in range(opts.split):
        tmin = time_all[0][0]
        tmax = time_all[0][-1]
        twidth = (tmax - tmin) / opts.split
        tspan = (i + 1) * twidth
        t0 = tmin + 0.5 * tspan
        t1 = tmin + tspan
        barres = bar(energies2x2, time_all, t0, t1) 
        print("BAR %f-%f:" % (t0, t1), barres * gasconstant* opts.temp)
        results.append((t0, t1, barres))

    pickle.dump(results, open("%s/results.pickle" % opts.save_dir, "wb"))

main()

