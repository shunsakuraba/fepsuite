#/usr/bin/python
import matplotlib
matplotlib.use('Agg')

import sys
import pandas
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

inpfile = sys.argv[1]
outfile = sys.argv[2]

nx = 1
ny = 2
(fig, axs) = plt.subplots(ny, nx, figsize=(10,10), dpi=120)

t = pandas.read_table(inpfile, header=None, sep=' ', names=["ty", "angle", "energy"])

for (kind, ix) in [("RNA", 0), ("DNA", 1)]:
    ax = axs[ix]
    dat = t[t.ty == kind]
    ax.plot(dat.angle, dat.energy, 'o')
    ax.set_xlim([0, 36])
plt.tight_layout()
plt.savefig(outfile, format="png")






