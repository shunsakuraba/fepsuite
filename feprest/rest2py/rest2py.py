#!/usr/bin/python
# Copyright 2018 Shun Sakuraba
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



# This script runs on python 2.7. Also tested on python3 but not sure whether it runs perfectly.
# Why not python 3 on 2018? Because 3 may not be installed on some supercomputer centers. Heck.

# I made this script because
# (1) sed-based REST2 provided by PLUMED2-hrex was slow (due to O(N^2) loop at dihedrals)
# (2) Wanted a capability to also convert topologies with FEP safely
# (3) Wanted to speed the calculation up by removing charge FEP
# Kudos to partial_tempering.sh

import argparse
import sys
import math
import copy

# I love this in C++.
def next_permutation(ls):
    if len(ls) <= 1:
        return (False, ls)
    i = len(ls) - 1
    while True:
        j = i
        i -= 1
        if ls[i] < ls[j]:
            k = len(ls) - 1
            while not (ls[i] < ls[k]):
                k -= 1
            tmp = ls[i]
            ls[i] = ls[k]
            ls[k] = tmp
            ls[j:] = ls[j:][::-1]
            return (True, ls)
        if i == 0:
            ls = ls[::-1]
            return (False, ls)

# based on gromacs 5.1.3+ behaviour, best match first.
# (see https://redmine.gromacs.org/issues/1901 )
def find_matching_dihedral(dihtype, ai, aj, ak, al, dfun):
    origindex = [ai, aj, ak, al]
    for nmatch in [4,3,2,1,0]:
        wildcards = [i >= nmatch for i in range(4)]
        while True:
            for fwddir in [True, False]:
                if fwddir:
                    ixs = copy.deepcopy(origindex)
                else:
                    ixs = origindex[::-1]
                for i in range(4):
                    if wildcards[i]:
                        ixs[i] = "X"
                key = tuple(ixs + [dfun])
                if key in dihtype:
                    return dihtype[key]
            (ret, wildcards) = next_permutation(wildcards)
            if not ret:
                break
    raise RuntimeError("Could not find dihedral for %s-%s-%s-%d" % ai, aj, ak, al, dfun)

def print_dihed_parameters(args, ofh, ai, aj, ak, al, dihfun, sscale, end_restrain, params_list, is_dih_restraints=False):
    for params in params_list:
        ofh.write("%5d %5d %5d %5d %2d" % (ai, aj, ak, al, dihfun))
        # This table slightly eases the mess for dihfun / fep
        if is_dih_restraints:
            (n_nonfep, to_scale, intind) = {
                1: (3, [2], [])
            }[dihfun]
        else:
            (n_nonfep, to_scale, intind) = {
                1: (3, [1], [2]),
                2: (2, [1], []),
                3: (6, [0, 1, 2, 3, 4, 5], []),
                4: (3, [1], [2]),
                5: (4, [0, 1, 2, 3], []),
                9: (3, [1], [2])
            }[dihfun]
        if len(params) not in [n_nonfep, 2 * n_nonfep]:
            raise RuntimeError("Number of args in dihedrals: expected %d or %d, but was %d (%d-%d-%d-%d:%d)" %
                                               (n_nonfep, 2 * n_nonfep, len(params),
                                                ai, aj, ak, al, dihfun))
        if (end_restrain and
                            len(params) == 2 * n_nonfep):
                            # calculate parameters
            for i in range(n_nonfep):
                if i in to_scale:
                    vA = float(params[i])
                    vB = float(params[i + n_nonfep])
                    if vA == 0.:
                                        # A real(no dih restr) -> B phantom (restr)
                        lm = args.end_restrain_dihedralA
                    elif vB == 0.:
                                        # A phantom -> B real
                        lm = args.end_restrain_dihedralB
                    else:
                        raise RuntimeError("Unexpected input @ end_restrain processing")
                    v = vA * (1. - lm) + vB * lm
                    params[i] = v
                    params[i + n_nonfep] = v
                                    
        for j in range(len(params) // n_nonfep):
            for i in range(n_nonfep):
                v = params[i + j * n_nonfep]
                if i in intind:
                    try:
                        v = int(v)
                    except ValueError:
                        if args.ignore_noninteger_periodicity:
                            vf = float(v)
                            vi = int(vf)
                            if abs(vi - vf) > 1e-2:
                                raise RuntimeError("Periodicity should be integer but was %f" % vf)
                            v = vi
                        else:
                            raise RuntimeError("Periodicity should be integer")
                    fmts = "%1d" % v
                else:
                    v = float(v)
                    if i in to_scale:
                        v *= sscale
                    fmts = "%16.8e" % v
                ofh.write(" %s" % fmts)
        ofh.write("\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Converts topology with REST2 protocol.
    Input topology must be preprocessed topology (generated by gmx grompp -pp). Also the atom types you want to scale with REST2 must be appended with "_". See PLUMED2's partial_tempering module for details.
    """)
    parser.add_argument('topology', action='store', type=str, 
                        help='Topology file (must be preprocessed and underlined)')
    parser.add_argument('output', action='store', type=str, 
                        help='Topology file output')
    parser.add_argument('--temp0', action='store', type=float, default=300.0, 
                        help='Simulation temperature')
    parser.add_argument('--temp', action='store', type=float, required=True,
                        help='Requesting temperature (lambda=T0/T)')
    parser.add_argument('--unify-charge', action='store_true', default=False,
                        help='Turn charge FEP to single-charge')
    parser.add_argument('--charge-lambda', action='store', type=float, default=None, 
                        help='Lambda value used in --unify-charge')
    parser.add_argument('--charge-lambda-A', action='store', type=float, default=None,
                        help='Lambda value used for A(real) -> B(phantom) atoms')
    parser.add_argument('--charge-lambda-B', action='store', type=float, default=None,
                        help='Lambda value used for A(phantom) -> B(real) atoms')
    parser.add_argument('--end-restrain-dihedralA', action='store', type=float, default=None,
                        help='Diherdal restraints of A(real) -> B(phantom) atoms use different lambda to prevent a sampling issue. Recommended to have values quickly growing to 1 near B state')
    parser.add_argument('--end-restrain-dihedralB', action='store', type=float, default=None,
                        help='Diherdal restraints of A(phantom) -> B(real) atoms use different lambda to prevent a sampling issue. Recommended to have values quickly growing to 1 near A state')
    parser.add_argument('--end-restrain-comment', action='store', type=str, default="dummy conn.",
                        help='Dihed parameters with this keyword in the comment is recognized as phantom-restraining dihedrals.')
    parser.add_argument('--exclude-peptide', action='store_true', default=True,
                        help='Exclude dihedral angle scaling for non-prolyl peptide bonds. This is necessary to prevent cis peptide formation (see also Neale et al., JCTC 2016, 12, 1989; doi: 10.1021/acs.jctc.5b01022).')
    parser.add_argument('--no-exclude-peptide', dest='exclude_peptide', action='store_false', default=False,
                        help='Stop excluding dihedral angle scaling for non-prolyl peptide bonds.')
    
    parser.add_argument('--suffix', action='store', type=str, default="_", 
                        help='Suffix character')
    parser.add_argument('--disable-mass-unification', dest="unify_mass", action='store_false', 
                        help='Hrex may not work correctly if the mass is exchanged, so without this option, the masses are changed to the maximum of the two states. Default is unifying them. Use this option only when you know what you are doing.')
    parser.add_argument('--ignore-noninteger-periodicity', dest='ignore_noninteger_periodicity', action='store_true',
                        help='Some dihedral parameters require periodicity sections, which should typically be integer. Instead of aborting on float, this option ignores non-integer input.')

    args = parser.parse_args()

    # args consistency check
    if args.unify_charge and args.charge_lambda is None:
        sys.stderr.write("--unify-charge requires --charge-lambda\n")
        parser.print_help()
        sys.exit(1)

    scale = args.temp0 / args.temp
    if scale > 1.0:
        sys.stderr.write("Warning: scaling factor exceeds 1.0 (temp > temp0)\n")
    sqrtscale = math.sqrt(scale)

    with open(args.topology) as fh, open(args.output, "w") as ofh:
        bondtype_of_atomtype = {}
        atomtype_info = {}
        dummy_atomtypes = set()
        dihtype = {}
        molecule = None
        sectiontype = None
        coulombrule = None
        residues = None
        atomnames = None

        def is_mainchain_nonpro(ai, aj, ak, al):
            ai -= 1
            aj -= 1
            ak -= 1
            al -= 1
            if (atomnames[ai], atomnames[aj], atomnames[ak], atomnames[al])\
               in [('O', 'C', 'N', 'H'), ('H', 'N', 'C', 'O')]: # AMBER nomenclature
                if atomnames[ai] == 'O' and residues[al] != 'PRO':
                    return True
                if atomnames[ai] == 'H' and residues[ai] != 'PRO':
                    return True
            if (atomnames[ai], atomnames[aj], atomnames[ak], atomnames[al])\
               in [('O', 'C', 'N', 'HN'), ('HN', 'N', 'C', 'O')]: # CHARMM nomenclature
                if atomnames[ai] == 'O' and residues[al] != 'PRO':
                    return True
                if atomnames[ai] == 'HN' and residues[ai] != 'PRO':
                    return True
            return False

        for lraw in fh:
            while lraw.endswith("\\\n"):
                lraw = "%s %s" % (lraw[:-2], next(fh)) # this is an undocumented feature in grompp, space is implicit
            ltmp = lraw.split(';', 1)
            if len(ltmp) == 1:
                l = ltmp[0]
                comment = ""
            else:
                l = ltmp[0]
                comment = ";" + ltmp[1]
            l = l.strip()
            ls = l.split()
            if l.startswith('#'):
                sys.stderr.write("The topology file is not preprocessed")
                sys.exit(1)
            if l.startswith('['):
                sectiontype = ls[1]
                ofh.write(lraw)
                continue

            # blank line
            if len(ls) == 0:
                ofh.write(lraw)
                continue

            if sectiontype is None:
                pass
            elif sectiontype == 'defaults':
                combrule = int(ls[1])
            elif sectiontype == 'atomtypes':
                if len(ls) < 6:
                    raise RuntimeError("Atomtype contains < 6 fields")
                # here everything is mess but toppush.cpp is actually super mess 
                if len(ls[5]) == 1 and ls[5][0].isalpha():
                    have_bonded_type = True
                    have_atomic_number = True
                elif len(ls[3]) == 1 and ls[3][0].isalpha():
                    have_bonded_type = False
                    have_atomic_number = False
                else:
                    have_bonded_type = ls[1][0].isalpha()
                    have_atomic_number = not have_bonded_type

                atomtype = ls[0]
                (mass, charge, particle, sigc6, epsc12) = ls[1 + int(have_bonded_type) + int(have_atomic_number):]

                if have_bonded_type:
                    bondtype = ls[1]
                else:
                    bondtype = atomtype
                if have_atomic_number:
                    atomic_ix = 1 + int(have_bonded_type)
                    atomic_number = int(ls[atomic_ix])
                else:
                    atomic_number = 0 # ??
                    
                # store this because we use in [ atoms ] section
                bondtype_of_atomtype[atomtype] = bondtype

                mass = float(mass)
                charge = float(charge)
                sigc6 = float(sigc6)
                epsc12 = float(epsc12)
                is_dummy = epsc12 == 0.

                atomtype_info[atomtype] = (charge, mass)
                if is_dummy:
                    dummy_atomtypes.add(atomtype)

                # for each atomtype print scaled one
                if combrule == 1: # C6 and C12
                    scale1 = scale
                    scale2 = scale
                else:
                    scale1 = 1.0
                    scale2 = scale

                # also scales the charge just in case (which doesn't exist in partial_tempring.sh)
                ofh.write('%4s %4s %4d %16.8e %16.8e %1s %16.8e %16.8e ; scaled\n' %
                          (atomtype + args.suffix,
                           bondtype + args.suffix, # bond type is also suffixed, this is almost never used except cmaptypes
                           atomic_number,
                           mass,
                           charge * sqrtscale, # charge is scaled by sqrt(c)
                           particle,
                           sigc6 * scale1,
                           epsc12 * scale2))
                # and the original line is also appended after this big if-elif-else block
            elif sectiontype == 'pairtypes' or sectiontype == 'nonbond_params':
                if combrule == 1: # C6 and C12
                    scale1 = scale
                    sqrtscale1 = sqrtscale
                    scale2 = scale
                    sqrtscale2 = sqrtscale
                else:
                    scale1 = 1.0
                    sqrtscale1 = 1.0
                    scale2 = scale
                    sqrtscale2 = sqrtscale
                sigc6 = float(ls[3])
                epsc12 = float(ls[4])
                ofh.write("%s%s %s %s %f %f ; scaled\n" %
                          (ls[0], args.suffix, ls[1], ls[2], sigc6 * sqrtscale1, epsc12 * sqrtscale2))
                ofh.write("%s %s%s %s %f %f ; scaled\n" %
                          (ls[0], ls[1], args.suffix, ls[2], sigc6 * sqrtscale1, epsc12 * sqrtscale2))
                ofh.write("%s%s %s%s %s %f %f ; scaled\n" %
                          (ls[0], args.suffix, ls[1], args.suffix,ls[2], sigc6 * scale1, epsc12 * scale2))
            elif sectiontype in ['bondtypes', 'angletypes', 'constrainttypes']:
                if sectiontype in ['bondtypes', 'constrainttypes']:
                    natomtype = 2
                elif sectiontype == 'angletypes':
                    natomtype = 3
                atomtypes = ls[0:natomtype]
                bondfun = int(ls[natomtype])
                # bond func is not scaled at all
                values = ls[natomtype + 1:]
                for bitmask in range(1, 2**natomtype):
                    newatomtypes = []
                    for (i, a) in enumerate(atomtypes):
                        if ((bitmask >> i) & 1) == 1:
                            newatomtypes.append(a + args.suffix)
                        else:
                            newatomtypes.append(a)
                    ofh.write("%s %d %s\n" % (" ".join(newatomtypes),
                                              bondfun,
                                              " ".join(values)))

            elif sectiontype == 'angletypes':
                atomtypes = ls[0:2]
                bondfun = int(ls[2])
                # bond func is not scaled at all
                values = ls[5:]
                for bitmask in range(1, 2**2):
                    newatomtypes = []
                    for (i, a) in enumerate(atomtypes):
                        if ((bitmask >> i) & 1) == 1:
                            newatomtypes.append(a + args.suffix)
                        else:
                            newatomtypes.append(a)
                    ofh.write("%s %s %d %s\n" % tuple(newatomtypes + [bondfun] + [" ".join(values)]))
            elif sectiontype == 'dihedraltypes':
                (ai, aj, ak, al) = ls[0:4]
                dihfun = int(ls[4])
                values = ls[5:]
                key = (ai, aj, ak, al, dihfun)
                if dihfun == 9:
                    # allows multiple dihedraltype for fn = 9
                    if key not in dihtype:
                        dihtype[key] = []
                    dihtype[key].append(values)
                else:
                    if key in dihtype:
                        for (i, e) in enumerate(dihtype[key]):
                            d = abs(float(values[i]) - float(e))
                            if d > 1e-20:
                                raise RuntimeError("Multiple dihedral for dihfun = %d, %s-%s-%s-%s"
                                                   % (dihfun, ai, aj, ak, al))
                    else:
                        dihtype[key] = values
                continue # suppress printing, we won't use dihedraltypes.
            elif sectiontype == 'cmaptypes':
                # For CMAP, we have no way to directly specify [ cmap ] to GROMACS.
                # Thus all conversion must be specified here.
                atomtypes = ls[0:5]
                cmapfun = int(ls[5])
                assert cmapfun == 1
                grid1 = int(ls[6])
                grid2 = int(ls[7])
                ngrid = grid1 * grid2 # current cmap should have this number
                for bitmask in range(2 ** 5):
                    if bitmask == 0:
                        cmapscale = 1.0
                        commentstr = ""
                    else:
                        cmapscale = scale
                        commentstr = " ; scaled"
                    newatomtypes = []
                    for (i, a) in enumerate(atomtypes):
                        if ((bitmask >> i) & 1) == 1:
                            newatomtypes.append(a + args.suffix)
                        else:
                            newatomtypes.append(a)
                    ofh.write("%s 1 %d %d" % (" ".join(newatomtypes), grid1, grid2))
                    assert len(ls[8:]) == ngrid
                    for (i, v) in enumerate(ls[8:]):
                        if i % grid2 == 0:
                            ofh.write("\\\n") # Hack: prevent 4095 chars line limit
                        else:
                            ofh.write(" ")
                        # AFAIK, all cmap dihedrals have nothing to do with omega peptide angle, thus scaling is justified here.
                        ofh.write("%.8f" % (float(v) * cmapscale))
                    ofh.write("%s\n\n" % commentstr) # extra eol, just in case, I felt a bad smell in the cmap code...
                # indeed we can print lraw as-is, but due to gromacs' 4095 char limit, here we complete raw values as well
                continue
            elif sectiontype == 'moleculetype':
                molecule = ls[0]
                # These None are sentinels for 1-origin access
                bondtype_list = [None]
                scaled = [None]
                atomnames = []
                residues = []
            elif sectiontype == 'atoms':
                aindex = int(ls[0])
                atomtype = ls[1]
                # remove suffix to store atomname
                is_scaled = atomtype.endswith(args.suffix)
                if is_scaled:
                    canonical_atomtype = atomtype[:-len(args.suffix)]
                else:
                    canonical_atomtype = atomtype
                assert(aindex == len(scaled))
                scaled.append(is_scaled)
                bondtype_list.append(bondtype_of_atomtype[canonical_atomtype])

                # charge & mass is optional parameters, oof...
                if len(ls) > 6:
                    charge = float(ls[6])
                else:
                    (charge, _) = atomtype_info[canonical_atomtype]
                if len(ls) > 7:
                    mass = float(ls[7])
                else:
                    (_, mass) = atomtype_info[canonical_atomtype]
                if len(ls) > 8:
                    atomtypeB = ls[8]
                    if atomtypeB.endswith(args.suffix):
                        sys.stderr.write("Warning: typeB in atom line should not be suffixed (molecule %s, atom index %d)\n" % (molecule, aindex))
                    (chargeB, massB) = atomtype_info[atomtypeB]
                    fep = True
                    if is_scaled and not atomtypeB.endswith(args.suffix):
                        atomtypeB += args.suffix
                else:
                    fep = False
                    chargeB = 0.
                    massB = 0.
                    atomtypeB = None
                if len(ls) > 9:
                    chargeB = float(ls[9])
                if len(ls) > 10:
                    massB = float(ls[10])

                if fep and args.unify_charge:
                    # lambda = 0: 100% A, lambda = 1: 100% B
                    real_atomtypeA = atomtype
                    if atomtype.endswith(args.suffix):
                        real_atomtypeA = atomtype[:-1]
                    real_atomtypeB = atomtypeB
                    if atomtypeB.endswith(args.suffix):
                        real_atomtypeB = atomtypeB[:-1]
                        
                    phantomA = (charge == 0. and real_atomtypeA in dummy_atomtypes)
                    phantomB = (chargeB == 0. and real_atomtypeB in dummy_atomtypes)
                    # charge_lambda_A is used when A(real) => B(phantom)
                    if (args.charge_lambda_A is not None and
                        not phantomA and phantomB):
                        clambda = args.charge_lambda_A
                    # charge_lambda_B is used when A(phantom) => B(real)
                    elif (args.charge_lambda_B is not None and
                          phantomA and not phantomB):
                        clambda = args.charge_lambda_B
                    else:
                        clambda = args.charge_lambda

                    newcharge = charge * (1. - clambda) + chargeB * clambda
                    charge = newcharge
                    chargeB = newcharge

                if fep and args.unify_mass:
                    massL = max(mass, massB)
                    mass = massL
                    massB = massL

                if is_scaled:
                    charge *= sqrtscale
                    chargeB *= sqrtscale

                _resnr = ls[2]
                resname = ls[3]
                residues.append(resname)
                atomname = ls[4]
                atomnames.append(atomname)
                    
                ofh.write("%5d %4s %4s %4s %4s %5s %16.8e %16.8e" % (aindex, atomtype, ls[2], ls[3], ls[4], ls[5], charge, mass))
                if fep:
                    ofh.write(" %4s %16.8e %16.8e%s\n" % (atomtypeB, chargeB, massB, comment.rstrip()))
                else:
                    ofh.write(" %s\n" % comment.rstrip())
                continue
            elif sectiontype == 'dihedrals':
                dihfun = int(ls[4])
                if dihfun in [1,2,3,4,5,9]:
                    (ai, aj, ak, al) = [int(x) for x in ls[0:4]]
                    sscale = 1.0
                    if scaled[ai]:
                        sscale *= sqrtscale
                    if scaled[al]:
                        sscale *= sqrtscale
                    if len(ls) == 5:
                        # must load dihedral table
                        (ti, tj, tk, tl) = [bondtype_list[x] for x in [ai, aj, ak, al]]
                        ofh.write("; parameters for %s-%s-%s-%s, fn=%d\n" % (ti, tj, tk, tl, dihfun))
                        params_tmp = find_matching_dihedral(dihtype, ti, tj, tk, tl, dihfun)
                        matched = True
                    else:
                        params_tmp = ls[5:]
                        matched = False
                    if args.end_restrain_dihedralA is not None and args.end_restrain_dihedralB is not None:
                        end_restrain = (comment.find(args.end_restrain_comment) != -1)
                    else:
                        end_restrain = False
                    if dihfun == 9 and matched:
                        params_list = params_tmp
                    else:
                        params_list = [params_tmp]

                    if args.exclude_peptide and is_mainchain_nonpro(ai, aj, ak, al) and sscale != 1.0:
                        sys.stderr.write("Prevented peptide bond scaling %s %d-%d-%d-%d (%s:%s-%s:%s-%s:%s-%s:%s)\n" %
                                         (molecule,
                                          ai, aj, ak, al,
                                          residues[ai-1], atomnames[ai-1],
                                          residues[aj-1], atomnames[aj-1],
                                          residues[ak-1], atomnames[ak-1],
                                          residues[al-1], atomnames[al-1]))
                        sscale = 1.0

                    print_dihed_parameters(args, ofh, ai, aj, ak, al, dihfun, sscale, end_restrain, params_list)
                    continue
            elif sectiontype == 'dihedral_restraints':
                dihfun = int(ls[4])
                if dihfun in [1]:
                    (ai, aj, ak, al) = [int(x) for x in ls[0:4]]
                    if args.end_restrain_dihedralA is not None and args.end_restrain_dihedralB is not None:
                        end_restrain = True
                    else:
                        end_restrain = False
                    params_list = [ls[5:]]
                else:
                    raise RuntimeError("Unsupported dihedral restraints function")
                sscale = 1.0
                print_dihed_parameters(args, ofh, ai, aj, ak, al, dihfun, sscale, end_restrain, params_list, is_dih_restraints=True)
                continue

            # With a few exceptions, we just print as-is
            ofh.write(lraw)






