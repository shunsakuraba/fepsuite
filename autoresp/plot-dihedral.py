#/usr/bin/python
import matplotlib
matplotlib.use('Agg')

import sys
import pandas
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

nx = 2
ny = 2
(fig, axs) = plt.subplots(ny, nx, figsize=(10,10), dpi=120)

names = "c2c3 c1c2 o4c1 c4o4".split()

ftrunk = sys.argv[1]

ix = 0
for n in names:
    t = pandas.read_table("%s.dih.%s.txt" % (ftrunk, n), header=None, sep=' ', names=["ty", "angle", "dihed"])
    ax = axs[ix // 2, ix % 2]
    ax.plot(t.angle, t.dihed, 'o')
    ax.set_xlim([0, 36])
    ax.set_title(n)
    ix += 1
plt.tight_layout()
plt.savefig("dihedrals.png", format="png")






