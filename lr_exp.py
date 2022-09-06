import pyedr
import argparse
import numpy
import scipy.special

sum_factors = ["LJ (SR)", "LJ-14", "LJ recip.", "LJC-14 q", "LJC Pairs NB", "Disper. corr."]

gasconstant = 0.0083144599 # kJ/mol/K

def expcal(du, kbt):
    # ln<exp(X)> = ln (1/N sum[exp(x)]) = logsumexp(x) - log(N)
    lnexp = scipy.special.logsumexp(du) - numpy.log(len(du))
    return - 1./kbt * lnexp

def lr_exp(args):
    original = pyedr.pyedr.EDRFile(args.short)
    reeval = pyedr.pyedr.EDRFile(args.long)
    original_e_ixs = []
    reeval_e_ixs = []

    for (i, enx) in enumerate(original.nms):
        if enx.name in sum_factors:
            print("Loading '%s' from short ranged edr" % enx.name)
            original_e_ixs.append(i)
    for (i, enx) in enumerate(reeval.nms):
        if enx.name in sum_factors:
            print("Loading '%s' from long ranged edr" % enx.name)
            reeval_e_ixs.append(i)

    kbt = gasconstant * args.temp
    ts = []
    deltas = []
    with open(args.output, "w") as ofh:
        for (origfr, reevalfr) in zip(iter(original), iter(reeval)):
            if origfr.t != reevalfr.t:
                raise RuntimeError("Two time frames did not match")
            if origfr.t <= args.time_begin:
                continue
            sum_orig = 0.
            sum_reeval = 0.
            for i in original_e_ixs:
                sum_orig += origfr.ener[i].e
            for i in reeval_e_ixs:
                sum_reeval += reevalfr.ener[i].e
            deltas.append(- (sum_reeval - sum_orig) / kbt)
            ts.append(origfr.t)
    # Perform B-fold blockwise error estimation
    deltas = numpy.array(deltas)

    with open(args.output, "w") as ofh:
        estimates = []
        for b in range(args.block):
            ixb = (b + 0) * len(deltas) // args.block
            ixe = (b + 1) * len(deltas) // args.block
            estimates.append(expcal(deltas[ixb:ixe], kbt))
            print(ts[ixb], ts[ixe - 1], "%.4f" % estimates[-1], file=ofh)

        estmean = numpy.mean(estimates)
        eststd = numpy.std(estimates, ddof=1)

        print("%.4f\t%.4f" % (estmean, eststd), file=ofh)

def init_args():
    parser = argparse.ArgumentParser(description="Calculate long-range correction",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--long", type=str, help="Input edr file, re-evaled with long range", required=True)
    parser.add_argument("--short", type=str, help="Input edr file, short range evaluation", required=True)
    parser.add_argument("--output", type=str, help="Output ndx file", required=True)
    parser.add_argument("--temp", type=float, help="Temperature (K)", required=True)
    parser.add_argument("--time-begin", type=float, default=0., help="Starting time (ps)")
    parser.add_argument("--block", type=int, default=5, help="Numbe of blocks to split (to esimate the error)")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()
    
    lr_exp(args)

    #python3 $ABFE_ROOT/lr-exp.py --long $ID/$lr_output.edr --short $ID/$non_lr_output.$non_lr_repl/$non_lr_output.edr --output $ID/$lr_output.lrc.txt

