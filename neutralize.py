#!/usr/bin/python

import sys
import random

natompersol = 3 # FIXME for TIP4P

[_, topin, groin, topout, groout, positive_at, negative_at] = sys.argv[0:7]
if len(sys.argv) > 7:
    otype = sys.argv[7]
    htype = sys.argv[8]
else:
    otype = "OW"
    htype = "HW"

dtotal = 0.0

with open(topin) as fh, \
     open(topout, 'w') as ofh:
    for lraw in fh:
        l = lraw.split(";")[0].strip()
        ls = l.split()
        if l == "":
            pass
        elif l.startswith("["):
            sectionname = l.split()[1]
            if sectionname == "system":
                # TODO: what to do with other water / ion models?
                # ^ Just prepare itp file and it will be done =D
                ofh.write("""[ moleculetype ]
;  molname       nrexcl ; TIP3P model

  SOLchg             2

[ atoms ]
;    nr   type  resnr residue  atom   cgnr     charge       mass    atomB
""")
                if dtotal < -0.5:
                    # increase charge
                    ofh.write(f"""
     1     {otype}      1     SOL     OW      1     -0.834    8.00000    {positive_at}      1.000  8.00000
     2     {htype}      1     SOL    HW1      1      0.417    4.00000    PHA     0.000  4.00000
     3     {htype}      1     SOL    HW2      1      0.417    4.00000    PHA     0.000  4.00000
""")
                elif dtotal > 0.5:
                    # decrease charge
                    ofh.write(f"""
     1     {otype}      1     SOL     OW      1     -0.834    8.00000    {negative_at}     -1.000  8.00000
     2     {htype}      1     SOL    HW1      1      0.417    4.00000    PHA     0.000  4.00000
     3     {htype}      1     SOL    HW2      1      0.417    4.00000    PHA     0.000  4.00000
""")
                else:
                    ofh.write(f"""
     1     {otype}      1     SOL     OW      1     -0.834    8.00000    {otype}     -0.834  8.00000
     2     {htype}      1     SOL    HW1      1      0.417    4.00000    {htype}     0.417  4.00000
     3     {htype}      1     SOL    HW2      1      0.417    4.00000    {htype}     0.417  4.00000
""")
                    

                ofh.write("""

[ constraints ]
;  i j   funct   length
1 2 1 0.09572
1 3 1 0.09572
2 3 1 0.15139

[ exclusions ]
1   2   3
2   1   3
3   1   2
""")

        elif l.startswith("#"):
            sys.stderr.write("Warning: input file must be preprocessed, but seems not\n")
        elif sectionname == "moleculetype":
            moleculetype = ls[0]
            nrexcl = int(ls[1])
        elif sectionname == "atoms" and moleculetype == "merged":
            # sum up charges
            chga = float(ls[6]) # FIXME: actually there are 4 types, see nucfep to see how to fix
            chgb = float(ls[9])
            dtotal += (chgb - chga)
        elif sectionname == "molecules" and ls[0] == "SOL":
            # convert one water into ion
            nsol = int(ls[1])
            nchg = abs(int(round(dtotal)))
            ofh.write("SOLchg\t%d\n" % nchg)
            ofh.write("SOL\t%d\n" % (nsol - nchg))
            continue
        ofh.write(lraw)

with open(groin) as fh, \
     open(groout, 'w') as ofh:
    title = True
    cursol = []
    sols = []
    for lraw in fh:
        if not title and len(lraw) > 10 and lraw[5:10].strip() == "SOL":
            cursol.append(lraw)
            if len(cursol) == 3:
                sols.append(cursol)
                cursol = []
                if len(sols) == nsol:
                    # print in a shuffled order
                    random.shuffle(sols)
                    for s in sols:
                        for lraw_ in s:
                            ofh.write(lraw_)
            continue
        ofh.write(lraw)
        title = False

