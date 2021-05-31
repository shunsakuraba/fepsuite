
def parse_top(topfile):
    molcomposition = []
    moleculetypes = {}
    section = None
    molecule = None
    atomlist = None
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
            if section == "moleculetype":
                molecule = ls[0]
                atomlist = []
            elif section == "atoms":
                atomtype = ls[1]
                resnr = int(ls[2])
                resname = ls[3]
                atomname = ls[4]
                charge = float(ls[6])
                atomlist.append((atomtype, resnr, resname, atomname, charge))
            elif section == "molecules":
                molname = ls[0]
                nmol = int(ls[1])
                molcomposition.append((molname, nmol))
    return { 'system': molcomposition, 'moleculetypes': moleculetypes }


def parse_index(ndx):
    group = None
    ret = {}
    with open(ndx) as fh:
        for l in fh:
            ls = l.split()
            if l.startswith('['):
                group = ls[1]
                if group in ret:
                    raise RuntimeError("Multiple groups in index file")
                else:
                    ret[group] = []
                continue
            else:
                for x in ls:
                    ret[group].append(int(x) - 1)
    return ret

