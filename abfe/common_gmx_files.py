import sys

def parse_top(topfile):
    atomtypes = {}
    molcomposition = []
    moleculetypes = {}
    section = None
    molecule = None
    atomlist = None
    defaults = None
    with open(topfile) as fh:
        for l in fh:
            if section is None:
                if l.startswith("*"):
                    continue
            l = l.rstrip()
            lcomment = l.split(';', 1)
            l = lcomment[0]
            l = l.strip()
            if l == "":
                continue
            elif l.startswith("#"):
                raise RuntimeError("topology is not preprocessed")
            elif l.startswith("["):
                sectionstr = l.lstrip("[")
                sectionstr = sectionstr.rstrip("]")
                sectionstr = sectionstr.strip()
                section = sectionstr
                if section in ["molecules", "system", "moleculetype"]:
                    # here we do a cleanup, as "molecules" / "system" must be after the end of molecular section
                    if molecule is not None:
                        moleculetypes[molecule] = atomlist
                        molecule = None
                        atomlist = []
                continue
            ls = l.split()
            if ls == []:
                continue
            if section == "atomtypes":
                # here everything is mess but toppush.cpp is actually super mess 
                if len(ls[5]) == 1 and ls[5].isalpha():
                    # "If field 5 is a single char we have both."
                    have_bonded_type = True
                    have_atomic_number = True
                elif len(ls[3]) == 1 and ls[3].isalpha():
                    # "If field 3 (starting from 0) is a single char, 
                    #  we have neither bonded_type or atomic numbers."
                    have_bonded_type = False
                    have_atomic_number = False
                else:
                    # The logic on this part was changed after GROMACS resolved issue 4120.
                    # "Attempt parsing field 1 to integer. If conversion fails, we do not have an atomic number but a bonded type.
                    # Unfortunately, int() in Python is more permissive (e.g. "3_10" is 310) than C atoi(), so I change this part to be more restrictive
                    have_atomic_number = all((c.isdigit() for c in ls[1]))
                    have_bonded_type = not have_atomic_number

                atomtype = ls[0]
                (mass, charge, particle, sigc6, epsc12) = ls[1 + int(have_bonded_type) + int(have_atomic_number):]
                mass = float(mass)
                charge = float(charge)
                sigc6 = float(sigc6)
                epsc12 = float(epsc12)

                if have_bonded_type:
                    bondtype = ls[1]
                    if all((c.isdigit() for c in ls[1])):
                        raise RuntimeError(f"""[ atomtypes ] contains bondtype "{atomtype}", but the bondtype for this atomtype consists of all digits.
This is considered invalid atomtype in GROMACS.""")
                else:
                    bondtype = atomtype
                if have_atomic_number:
                    atomic_ix = 1 + int(have_bonded_type)
                    atomic_number = int(ls[atomic_ix])
                else:
                    atomic_number = 0 # ??
                
                if defaults[1] == 1:
                    # geometric average with C6/C12
                    if sigc6 == 0.0 or epsc12 == 0.0:
                        eps = 0.0
                        sigma = 0.0
                    else:
                        sigma = (epsc12 / sigc6) ** (1. / 6.)
                        eps = sigc6 / (4.0 * (sigma ** 6))
                elif defaults[1] == 2 or defaults[1] == 3:
                    # arithmetic mean with sig/eps (2)
                    # geometric mean with sig/eps (3)
                    eps = epsc12
                    sigma = sigc6
                else:
                    raise RuntimeError(f"Unknown [ defaults ] cmb-rule: {defaults[1]}")
                atomtypes[atomtype] = (bondtype, particle, atomic_number, mass, charge, sigma, eps)
            elif section == "moleculetype":
                molecule = ls[0]
                atomlist = []
            elif section == "atoms":
                atomtype = ls[1]
                resnr = int(ls[2])
                resname = ls[3]
                atomname = ls[4]
                charge = float(ls[6]) if len(ls) >= 7 else None
                atomlist.append((atomtype, resnr, resname, atomname, charge))
            elif section == "molecules":
                molname = ls[0]
                nmol = int(ls[1])
                molcomposition.append((molname, nmol))
            elif section == "defaults":
                defaults = ls
                defaults[0] = int(defaults[0]) # nbfunc
                defaults[1] = int(defaults[1]) # cmbrule
                defaults[3] = float(defaults[3]) # fudgeLJ
                defaults[4] = float(defaults[4]) # fudgeQQ

    return { 'atomtypes': atomtypes, 'system': molcomposition, 'moleculetypes': moleculetypes }


def parse_index(ndx):
    group = None
    ret = {}
    with open(ndx) as fh:
        for l in fh:
            ls = l.split()
            if l.startswith('['):
                group = ls[1]
                if group in ret:
                    raise RuntimeError("Same group appeared more than once in the index file")
                else:
                    ret[group] = []
                continue
            else:
                for x in ls:
                    ret[group].append(int(x) - 1)
    return ret

