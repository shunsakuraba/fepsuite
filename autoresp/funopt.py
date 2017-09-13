# function optimization
# determines both angle and coefficients
import matplotlib
matplotlib.use('Agg')

import sys
import numpy as np
import cma
import math
import copy
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

Nfit = 4

def basefuncs(x, n, p):
    return 1.0 + math.cos((n * x - p) * math.pi / 180.0)
vecbases = np.vectorize(basefuncs)

def getlinear(angles, delta_es, weights, phases):
    allevals = []
    allevals.append(np.array([1.0] * len(angles)))
    for i in range(len(phases)):
        n = i + 1
        p = phases[i]
        bases = vecbases(angles, n, p)
        allevals.append(bases)
    allevals = np.array(allevals)
    sqweights = np.sqrt(weights)
    weighted_delta = copy.deepcopy(delta_es)
    for j in range(len(sqweights)):
        allevals[:, j] *= sqweights[j]
        weighted_delta[j] *= sqweights[j]
    (x, residuals, rank, sing) = np.linalg.lstsq(allevals.T, weighted_delta)
    residsum = np.linalg.norm(residuals)
    return (residsum, x)

def geteval(angles, facs, phases):
    allevals = np.zeros_like(np.array(angles))
    allevals += facs[0] # constants
    for i in range(len(phases)):
        n = i + 1
        p = phases[i]
        bases = vecbases(angles, n, p)
        allevals += facs[i + 1] * bases
    return allevals

def test_linear():
    angles = np.array(range(36)) * 10.0
    weights = [1.0] * 36
    delta = (100 + 
             0.9656 * vecbases(angles, 1,  68.79) +
             1.0740 * vecbases(angles, 2,  15.64) +
             0.4575 * vecbases(angles, 3, 171.58) +
             0.3092 * vecbases(angles, 4,  19.09))
    print getlinear(angles, delta, weights, [68.79, 15.64, 171.58, 19.09])
    print getlinear(angles, delta, weights, [70, 15, 170, 19])
    print delta

#test_linear()

if len(sys.argv) <= 4:
    print >> sys.stderr, "Usage: %s (qm) (mm) (weight) (output) [--hartree]" % sys.argv[0]
    sys.exit(1)

def read_values(fn):
    energies = {}
    with open(fn) as fh:
        for l in fh:
            (state, angle, val) = l.split()
            angle = int(angle)
            val = float(val)
            energies[(state, angle)] = val

    ret = []
    adjust_dna = energies[("RNA", 3)] - energies[("DNA", 3)]
    for i in range(36):
        if i >= 3 and i <= 15:
            energy = adjust_dna + energies[("DNA", i)]
        else:
            energy = energies[("RNA", i)]
        ret.append(energy)
    return ret

qms = read_values(sys.argv[1])
mms = read_values(sys.argv[2])
weights = file(sys.argv[3]).readlines()

weights = [float(x.strip()) for x in weights]

if len(sys.argv) > 5:
    if sys.argv[5] == "--hartree":
        # You have: (1 hartree * avogadro) / (kcal_th/mol)
        # You want: 
        #        Definition: 627.50947
        conversion = 627.50947
        qms = [conversion * x for x in qms]
        mms = [conversion * x for x in mms]
    else:
        print >> sys.stderr, "Unknown option: \"%s\"" % sys.argv[5]
        sys.exit(1)


qms = np.array(qms)
qmsmin = np.min(qms)
qms -= qmsmin

mms = np.array(mms)
mmsmin = np.min(mms)
mms -= mmsmin

deltaes = qms - mms

# assumes equal angle partition
dangle = 360.0 / len(deltaes)
angles = np.array(range(len(deltaes))) * dangle

# find initial guess of phases
phaseguess = np.array([0.0] * Nfit)
for i in range(Nfit):
    k = i + 1
    bestguess = 0.
    bestscore = 1e+9
    for a in angles:
        if a >= (359.9 / k):
            break
        phaseguess[i] = a
        (sc, _coeffs) = getlinear(angles, deltaes, weights, phaseguess)
        # print i, a, sc
        if sc < bestscore:
            bestscore = sc
            bestguess = a
    phaseguess[i] = bestguess

def ffun(phases):
    (r, _x) = getlinear(angles, deltaes, weights, phases)
    return r
cmaret = cma.fmin(ffun, phaseguess, dangle)
phaseopt = cmaret[0]

print >> sys.stderr, "Initial phase guess: ", phaseguess
print >> sys.stderr, "Final phase result: ", phaseopt

(r, facs) = getlinear(angles, deltaes, weights, phaseopt)

# final tweak: 
# if factor is negative, just shift the phase by pi and flip sign
# V cos(nx) = -V cos(nx + pi)
# and round to [-pi, pi]
phaseopt_tweak = copy.deepcopy(phaseopt)
for i in range(Nfit):
    if facs[i + 1] < 0:
        phaseopt_tweak[i] += 180.0
        phaseopt_tweak[i] -= 360.0 * round(phaseopt_tweak[i] / 360.0)

(r2, facs2) = getlinear(angles, deltaes, weights, phaseopt_tweak)

print >> sys.stderr, r, facs
print >> sys.stderr, r2, facs2

with open(sys.argv[4], "w") as fh:
    for i in range(Nfit):
        print >> fh, i + 1, facs2[i + 1], phaseopt_tweak[i]
        print i + 1, facs2[i + 1], phaseopt_tweak[i]

# Generate plot
(fig, axs) = plt.subplots(2, 1, figsize=(10,10), dpi=120)

fev = geteval(angles, facs2, phaseopt_tweak)

ax = axs[0]
ax.plot(angles, deltaes, 'o')
ax.plot(angles, fev, '-')
ax.set_xlim([0, 360])

ax = axs[1]
ax.plot(angles, qms, 'or')
ax.plot(angles, mms, 'xb')
ax.plot(angles, mms + fev, '-b')
ax.set_xlim([0, 360])

plt.tight_layout()
plt.savefig("funopt.png", format="png")

        
