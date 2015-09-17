#!/usr/bin/python

import argparse
import copy
import os.path

parser = argparse.ArgumentParser(description='Split topology file into 2 molecules.')
parser.add_argument('--resid', metavar='resid', type=int,
                    help='residue number to include (less than or equal to this number is included in the 1st molecule)', required = True)
parser.add_argument('--top', metavar='TOP', type=str,
                    help='topology file', required = True)
parser.add_argument('--output-1', dest='output1', type=str,
                    help='Output first half molecule', required = True)
parser.add_argument('--output-2', dest='output2', type=str,
                    help='Output second half molecule', required = True)
parser.add_argument('--pdb', metavar='pdb', type=str,
                    help='pdb (or gro) file for coordinate', required = True)
parser.add_argument('--pdbout1', metavar='pdbout1', type=str,
                    help='pdb (or gro) file for coordinate', required = True)
parser.add_argument('--pdbout2', metavar='pdbout2', type=str,
                    help='pdb (or gro) file for coordinate', required = True)

args = parser.parse_args()

class union_find:
    """Tarjan's famous disjoint-set data structure.

    """
    def __init__(self, n):
        self.uftable_ = [-1] * n
    
    def union(self, i, j):
        pi = self.find(i)
        pj = self.find(j)
        if pi != pj:
            if self.uftable_[pj] < self.uftable_[pi]:
                temp = pi
                pi = pj
                pj = temp
            self.uftable_[pi] += self.uftable_[pj];
            self.uftable_[pj] = pi
    
    def find(self, i):
        if self.uftable_[i] < 0:
            return i
        else:
            newroot = self.find(self.uftable_[i])
            self.uftable_[i] = newroot
            return newroot

    def is_connected(self, i, j):
        pi = self.find(i)
        pj = self.find(j)
        return pi == pj

split_of_atoms = {}
charges = [None, 0.0, 0.0]
atomptrs = [None, 0, 0]
with (open(args.top, "r")) as fh, (
      open(args.output1, "w")) as fh1, (
      open(args.output2, "w")) as fh2:
    writeselector = [None, fh1, fh2]
    state = ""
    atomid = 0
    for lraw in fh:
        l = lraw.strip()
        if l == "" or l[0] == ';' or l[0] == '#':
            fh1.write(lraw)
            fh2.write(lraw)
            continue
        if l[0] == '[':
            state = l.split()[1]
            fh1.write(lraw)
            fh2.write(lraw)
            continue
        if state == 'atoms':
            sp = l.split()
            resnum = int(sp[2])
            charge = float(sp[6])
            if resnum <= args.resid:
                st = 1
            else:
                st = 2
            split_of_atoms[atomid] = (st, atomptrs[st])
            sp1 = l.split(None, 1)
            lout = "%5d %s\n" % (atomptrs[st] + 1, sp1[1])
            writeselector[st].write(lout)
            atomptrs[st] += 1
            charges[st] += charge
            atomid += 1
        elif state in ["bonds", "angles", "dihedrals", "pairs", "exclusions"]:
            sp = l.split()
            a = int(sp[0]) - 1
            (st, _) = split_of_atoms[a]
            natomids = {"bonds": 2,
                        "angles": 3,
                        "dihedrals": 4,
                        "pairs": 2,
                        "exclusions": 2 }[state]
            qs = l.split(None, natomids)
            lout = ""
            for i in range(natomids):
                (st2, na) = split_of_atoms[int(qs[i]) - 1]
                if st != st2:
                    if state == "exclusions":
                        continue
                    else:
                        assert False, "state %s, atoms %s spans multiple sequences. Perhaps you specified wrong resid?" % (state, repr(qs))
                lout = "%s%5d" % (lout, na + 1)
            lout = "%s %s\n" % (lout, qs[natomids])
            writeselector[st].write(lout)
        elif state == "molecules":
            sp = l.split()
            if sp[0] == "NA":
                # assumes negatively charged
                fh1.write("%s %d\n" % (sp[0], int(sp[1]) + round(charges[2])))
                fh2.write("%s %d\n" % (sp[0], int(sp[1]) + round(charges[1])))
            else:
                fh1.write(lraw)
                fh2.write(lraw)
        else:
            fh1.write(lraw)
            fh2.write(lraw)

ispdb = (os.path.splitext(args.pdb)[1] == '.pdb')
with (open(args.pdb, "r")) as pdbfh, (
      open(args.pdbout1, "w")) as pdb1, (
      open(args.pdbout2, "w")) as pdb2:
    if not ispdb:
        header = pdbfh.next()
        pdb1.write(header)
        pdb2.write(header)
        natom = int(pdbfh.next().strip())
        pdb1.write("%5d\n" % (natom - atomptrs[2] + round(charges[2])))
        pdb2.write("%5d\n" % (natom - atomptrs[1] + round(charges[1])))
    atomid = 0
    nacount = 0
    for l in pdbfh:
        if ispdb and (l[0:6] not in ["ATOM  ", "HETATM"]):
            continue
        if (not ispdb) and len(l) > 50:
            pdb1.write(l)
            pdb2.write(l)
            break
        if len(split_of_atoms) <= atomid:
            if ispdb:
                resname = l[17:20]
            else:
                resname = l[5:10]
            if resname.strip() == "NA":
                if -round(charges[2]) <= nacount:
                    pdb1.write(l)
                if -round(charges[1]) <= nacount:
                    pdb2.write(l)
                nacount = nacount + 1
                continue
            pdb1.write(l)
            pdb2.write(l)
            continue
        (st, _) = split_of_atoms[atomid]
        atomid += 1
        ([None, pdb1, pdb2][st]).write(l)

