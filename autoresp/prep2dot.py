#!/usr/bin/python

import sys

# skip 4 lines
for _i in range(4):
    l = sys.stdin.next()

(rn, _, _) = sys.stdin.next().split()

# 
for _i in range(2):
    l = sys.stdin.next()

bonds=[]
atoms=[]
atommap={}
for l in sys.stdin:
    if l.strip() == "":
        break
    an = l[6:10].strip()
    at = l[10:14].strip()
    if at == "DU":
        continue
    bt = int(l[20:24].strip())
    atommap[an] = len(atoms)
    atoms.append((an, at))
    if bt - 4 >= 0:
        bonds.append((an, atoms[bt - 4][0]))

while True:
    l = sys.stdin.next()
    lc = l.strip()
    if lc == "":
        continue
    if lc == "LOOP":
        for l in sys.stdin:
            if l.strip() == "":
                break
            a1 = l[0:5].strip()
            a2 = l[5:10].strip()
            bonds.append((a1, a2))
    if lc == "STOP" or lc == "DONE":
        break
    else:
        for l in sys.stdin:
            if l.strip() == "":
                break

def name(n):
    i = atommap[n]
    (an, at) = atoms[i]
    return '"%s(%s)"' % (an, at)

def color(t):
    return { "H": "\"#606060\"",
             "O": "red",
             "C": "black",
             "N": "blue",
             "P": "orange" }.get(t[0], "black")

print "graph %s {" % rn
for (n, t) in atoms:
    print "    %s [fontcolor=%s];" % (name(n), color(t))
for (a, b) in bonds:
    print "    %s -- %s;" % (name(a), name(b))
print "}"


