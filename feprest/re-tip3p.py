import sys


[_prog, inputfile, outputfile] = sys.argv


with open(inputfile) as fh, open(outputfile, "w") as ofh:
    buf = []
    state = None
    molname = None
    for lraw in fh:
        if lraw.startswith('['):
            state = lraw.split()[1]
            if state == "moleculetype" or state == "system":
                # print buffer
                if molname == "SOL":
                    ofh.write('#include "../tip3p-molecule.itp"\n')
                    buf = []
                else:
                    for lb in buf:
                        ofh.write(lb)
                    buf = []
        elif state == "moleculetype":
            ltmp = lraw.split(';', 1)[0].strip()
            if ltmp != '':
                molname = ltmp.split()[0]
        buf.append(lraw)
    for lb in buf:
        ofh.write(lb)

