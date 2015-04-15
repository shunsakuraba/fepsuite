#!/usr/bin/python

import argparse
import copy

parser = argparse.ArgumentParser(description='Split topology file into 3 states.')
parser.add_argument('--top', metavar='TOP', type=str,
                    help='topology file', required = True)
parser.add_argument('--output-A', dest='outputA', type=str,
                    help='Output in-state-A spring morphing', required = True)
parser.add_argument('--output-B', dest='outputB', type=str,
                    help='Output in-state-B spring morphing', required = True)
parser.add_argument('--output-fep', dest='outputFEP', type=str,
                    help='Output A-to-B morphing', required = True)

args = parser.parse_args()

with (open(args.outputA, "w")) as Afh, (
      open(args.outputB, "w")) as Bfh, (
      open(args.outputFEP, "w")) as FEPfh, (
      open(args.top, "r")) as fh:
    state = ""
    phantomInA = set()
    phantomInB = set()
    for lraw in fh:
        l = lraw.strip()
        if l == "" or l[0] == ';':
            Afh.write(lraw)
            Bfh.write(lraw)
            FEPfh.write(lraw)
            continue
        if l[0] == '[':
            state = l.split()[1]
            Afh.write(lraw)
            Bfh.write(lraw)
            FEPfh.write(lraw)
            continue
        if state == 'atoms':
            st = l.split()
            Afh.write(" ".join(["%5s" % st[0],
                                "%4s" % st[1],
                                "%3s" % st[2],
                                "%4s" % st[3],
                                "%5s" % st[4],
                                "%4s" % st[5],
                                "%12s" % st[6],
                                "%12s" % st[7],
                                "%4s" % st[1],
                                "%12s" % st[6], 
                                "%12s" % st[7],
                                "\n"]))
            Bfh.write(" ".join(["%5s" % st[0],
                                "%4s" % st[8],
                                "%3s" % st[2],
                                "%4s" % st[3],
                                "%5s" % st[4],
                                "%4s" % st[5],
                                "%12s" % st[9],
                                "%12s" % st[10],
                                "%4s" % st[8],
                                "%12s" % st[9], 
                                "%12s" % st[10],
                                "\n"]))
            FEPfh.write(lraw)
            if st[1] == "PHA":
                phantomInA.add(int(st[0]))
            if st[8] == "PHA":
                phantomInB.add(int(st[0]))
        elif state in ["bonds", "angles", "dihedrals"]:
            l = l.split(";")[0]
            st = l.split()
            natom = { "bonds": 2,
                      "angles": 3,
                      "dihedrals": 4 }[state]
            atoms = [int(st[i]) for i in range(natom)]
            func = int(st[natom])
            npar = (len(st) - (natom + 1)) / 2
            Avars = st[natom+1:natom+1+npar]
            Bvars = st[natom+1+npar:natom+1+npar*2]
            A0vars = copy.copy(Avars)
            B1vars = copy.copy(Bvars)

            isphA = False
            for a in atoms:
                if a in phantomInA:
                    isphA = True
            isphB = False
            for a in atoms:
                if a in phantomInB:
                    isphB = True

            if state == "bonds" and func == 6:
                # phantom constraint
                if isphA:
                    assert(float(Bvars[1]) == 0.0)
                    Avars[1] = "%12.5e" % 0.0
                if isphB:
                    assert(float(Avars[1]) == 0.0)
                    Bvars[1] = "%12.5e" % 0.0
            else:
                # standard bonds / angles
                if isphA:
                    # phantom bonds in A
                    assert(float(Avars[1]) == 0.0)
                    assert(state == "dihedrals" or
                           float(Avars[0]) == float(Bvars[0]))
                    Avars = Bvars
                if isphB:
                    # phantom bonds in B
                    assert(float(Bvars[1]) == 0.0)
                    assert(state == "dihedrals" or
                           float(Avars[0]) == float(Bvars[0]))
                    Bvars = Avars

            Afh.write(" ".join(["%4d" % x for x in atoms] +
                               ["%4d" % func] + 
                               ["%12s" % v for v in A0vars] + 
                               ["%12s" % v for v in Avars] +
                               ["\n"]))
            Bfh.write(" ".join(["%4d" % x for x in atoms] +
                               ["%4d" % func] + 
                               ["%12s" % v for v in Bvars] + 
                               ["%12s" % v for v in B1vars] +
                               ["\n"]))
            FEPfh.write(" ".join(["%4d" % x for x in atoms] +
                                 ["%4d" % func] + 
                                 ["%12s" % v for v in Avars] + 
                                 ["%12s" % v for v in Bvars] +
                                 ["\n"]))                
        else:
            Afh.write(lraw)
            Bfh.write(lraw)
            FEPfh.write(lraw)
