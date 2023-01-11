import argparse
import mdtraj


def check(args):
    with open(args.rms) as fh:
        nsum = 0
        xsum = 0.
        for l in fh:
            if l.startswith("@") or l.startswith("#"):
                continue
            ls = l.split()
            if ls == []:
                continue
            ls = [float(x) for x in ls]
            t = ls[0]
            v = ls[1]
            xsum += v
            nsum += 1
        average = xsum / nsum
        if average > args.threshold:
            raise RuntimeError("RMS was too large (average %f, threshold %f)" % (average, args.threshold))

def init_args():
    parser = argparse.ArgumentParser(description="Check RMS of the trajectory",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--rms", type=str, required=True, help="xvg file representing the rms of the unrestrained run")
    parser.add_argument("--threshold", type=float, default=0.4, help="average rms thresholding (nm)")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()

    check(args)

