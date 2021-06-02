#!/usr/bin/env python3

import sys
import mdtraj
import re
import argparse

def main(args):
    topology = mdtraj.load(args.base_structure).topology

    nst = []
    middles = []
    with open(args.cluster_log) as fh:
        pat = re.compile(r"^cl. \|")
        for l in fh:
            if pat.match(l):
                break
        cluster_id = 0
        for l in fh:
            cols = [s.strip() for s in l.split("|")] 
            if cols[0] != "":
                assert cluster_id + 1 == int(cols[0])
                cluster_id = int(cols[0])
                nst.append(int(cols[1].split()[0]))
                middles.append(float(cols[2].split()[0]))
  
    totalframe = sum(nst)
    selected = [(nst[i], middles[i]) for i in range(len(nst)) if nst[i] >= args.threshold * totalframe]
    
    with open("%s.txt" % args.output_prefix, "w") as ofh:
        print(len(selected), file=ofh)
        for (isel, (ns, s)) in enumerate(selected):
            print(isel, ns / totalframe, s, file=ofh)

    for chunk in mdtraj.iterload(args.traj, top=topology):
        for (ifr, t) in enumerate(chunk.time):
            for (isel, (_ns, s)) in enumerate(selected):
                if abs(s - t) < 1e-3:
                    chunk[ifr].save_pdb("%s%d.pdb" % (args.output_prefix, isel))

def init_args():
    parser = argparse.ArgumentParser(description="Pick major cluster and write corresponding structure files",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cluster-log", type=str, required=True, help="Cluster log file")
    parser.add_argument("--base-structure", type=str, required=True, help="Base structure file")
    parser.add_argument("--traj", type=str, required=True, help="Trajectory file (without fitting recommended)")
    parser.add_argument("--threshold", type=float, required=True, help="Fraction of snapshots to output cluster")
    parser.add_argument("--output-prefix", type=str, required=True, help="output file prefix")

    return parser.parse_args()

if __name__ == "__main__":
    args = init_args()

    main(args)

