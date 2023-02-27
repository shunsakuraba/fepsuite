import copy
import encodings
import shutil
import subprocess
import os.path
import os
import sys
from tokenize import maybe
from typing import List
import Bio.PDB
import Bio.SeqUtils
import argparse
import re
import warnings

def log_with_color(message):
    col = 33 # yellow?
    prefix = ">> "
    suffix = ""
    if os.isatty(sys.stdout.fileno()):
        prefix = "\x1b[" + str(col) + ";1m>> "
        suffix = "\x1b[0m"
    print(prefix + message + suffix, file=sys.stdout)

def check_call_verbose(cmdline: List[str]):
    message = " ".join(cmdline)
    log_with_color(message)
    subprocess.check_call(cmdline)

def seq1(aa):
    a = Bio.SeqUtils.seq1(aa)
    if a != "X":
        return a
    nonstandard_map = {
        'HID': 'H', 'HIE': 'H', 'HIP': 'H',
        'HSD': 'H', 'HSE': 'H', 'HSP': 'H',
        'GLH': 'E', 'ASH': 'D', # FIXME: charmm uses ASPP and GLUP but what should I do here?
        'LYN': 'K', 'LSN': 'K',
        'CYX': 'C', 'CYM': 'C'
    }
    if aa in nonstandard_map:
        return nonstandard_map[aa]
    raise RuntimeError("Unknown residue name " + aa)

def parse_pdb(pdb):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Bio.PDB.PDBExceptions.PDBConstructionWarning)
        p = Bio.PDB.PDBParser()
        structure = p.get_structure('X', pdb)
    for model in structure:
        seq = ''
        reschainmap = {}
        order = []
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                aaa = residue.get_resname()
                a = seq1(aaa)
                seq += str(a)
                resid = residue.get_id()[1] # get_id returns triplet, but only 2nd (resseq) is used
                if resid not in reschainmap:
                    reschainmap[resid] = {}
                if chain_id in reschainmap[resid]:
                    raise RuntimeError(f"Duplicated Chain-resid combination at {chain_id}-{resid}")
                reschainmap[resid][chain_id] = a
                order.append((chain_id, resid))
        return (seq, order, reschainmap) # returns at the first structure
    raise RuntimeError("PDB ended before the first structure appears")

def parse_mut(mutstr: str):
    # accepted syntax example:
    # 22A
    # A:V23P
    # g:S152K
    patseq = re.compile(r"(?:(?P<chain>[A-Za-z]):)?(?P<before>[A-Z])?(?P<resid>\d+)(?P<after>[A-Z])$")
    muts = mutstr.split("_")
    res = []
    for mu in muts:
        matcher = patseq.match(mu)
        d = matcher.groupdict()
        d['resid'] = int(d['resid'])
        res.append(d)
    
    return res

# TODO: this function is defined but unused
def update_mutinfo(mutations, basepdbinfo):
    ret_mutations = []
    _pdbseq, pdborder, reschainmap = basepdbinfo
    for (c, r) in pdborder:
        res = reschainmap[r][c]
        for m in mutations:
            if m["resid"] == r and (m["chain"] is None or m["chain"] == c):
                # matched mutation
                if m["before"] is not None and m["before"] != res:
                    raise RuntimeError(f"Mutation ({m}) is not compatible with chain '{c}' resid {r}")
                newmut = copy.copy(m)
                newmut["before"] = res
                ret_mutations.append(newmut)
    if len(ret_mutations) != len(mutations):
        raise RuntimeError(f"Not all mutations are consumed, requested mutations: {mutations}, but found: {ret_mutations}")
    return ret_mutations

def generate_gmx(args, basepdbinfo, ddir, mutations):
    if not os.path.exists(ddir):
        os.mkdir(ddir)
    else:
        if os.path.exists(f"{args.wtdir}/conf.pdb") and os.path.exists(f"{args.wtdir}/topol_can.top"):
            return # already done
    if mutations == []:
        # wild-type, just copy
        shutil.copy(args.pdb, f"{ddir}/completed.pdb")
    else:
        with open(f"{ddir}/seq.txt", "w") as ofh:
            _pdbseq, pdborder, reschainmap = basepdbinfo
            curchain = None
            totmut = 0
            for (c, r) in pdborder:
                if curchain != c:
                    if curchain is not None:
                        print("", file=ofh)
                    chainname = c if c.strip() != "" else "X"
                    print(f">{chainname}", file=ofh)
                    curchain = c
                res = reschainmap[r][c]
                mutated = False
                for m in mutations:
                    if m["resid"] == r and (m["chain"] is None or m["chain"] == c):
                        if m["before"] is not None and m["before"] != res: # this check should not be necessary here
                            raise RuntimeError(f"Mutation ({m}) is not compatible with chain '{c}' resid {r}")
                        res = m["after"]
                        mutated = True
                        totmut += 1
                        break

                # Hack for FASPR: keep Cys rotamer intact
                if res == "C" and not mutated:
                    res = "c"

                print(res, end="", file=ofh) # keep Cys rotamer so that disulfide bonds are intact
            print("", file=ofh)
            if totmut != len(mutations):
                raise RuntimeError(f"Not all mutations are consumed, requested mutations: {mutations}")

        check_call_verbose([args.faspr, "-i", args.pdb, "-o", f"{ddir}/completed.pdb", "-s", f"{ddir}/seq.txt"])
        # faspr often returns 0 even if there is an error
        if not os.path.exists(f"{ddir}/completed.pdb") or \
            os.path.getmtime(f"{ddir}/completed.pdb") < os.path.getmtime(f"{ddir}/seq.txt"):
            raise RuntimeError("FASPR run failed    ")

    curdir = os.getcwd()
    os.chdir(ddir)
    if os.path.exists(f"../{args.ff}.ff") and not os.path.exists(f"{args.ff}.ff"):
        # use local force field file
        os.symlink(f"../{args.ff}.ff", f"{args.ff}.ff")

    # Generate GMX topology as well as structures
    check_call_verbose([args.gmx] + f"pdb2gmx -f completed.pdb -o conf.pdb -water {args.water_model} -ignh -ff {args.ff} -merge all".split())


    # "canonicalize" the topology
    with open("topol.top") as ifh, open("topol_m.top", "w") as ofh:
        out = True
        for l in ifh:
            # remove #ifdef POSRES ... #include "ions.itp"
            # TODO FIXME: this part should be fixed by updating nucfepgen
            if re.search("POSRES", l):
                out = False
            if out:
                ofh.write(l)
            if re.search("ions.itp", l):
                out = True
    check_call_verbose([args.python3, f"{args.feprest}/utils/pp.py", "-o", "topol_pp.top", "topol_m.top"])
    check_call_verbose([args.python3, f"{args.feprest}/rest2py/canonicalize_top.py", "topol_pp.top", "topol_can.top"])

    os.chdir(curdir)
    return

def difficulty(muts):

    def charge(r):
        if r in "DE":
            return -1
        elif r in "RK":
            return 1
        else:
            return 0

    difficult = False
    totcharge = 0
    for m in muts:
        mutfrom = m["before"]
        mutto = m["after"]
        totcharge -= charge(mutfrom)
        totcharge += charge(mutto)
        if mutto in "PFYW" or mutfrom in "PFYW":
            difficult = True
    return (difficult, totcharge != 0)

def generate_para_conf(args, muts):
    (difficult, charged) = difficulty(muts)
    if charged:
        difficult = True
    if difficult:
        with open("para_conf.zsh", "w") as ofh:
            print(f"NREP={args.difficult_nrep}", file=ofh)
            if charged:
                print(f"CHARGED=yes", file=ofh)

def mut_name(mutstr: str):
    return mutstr.replace(":", "_")

def solvate_conf_topology(radius):
    check_call_verbose([args.gmx] + f"editconf -f fepbase.pdb -d {radius} -bt dodecahedron -o conf_box".split())
    maybe_relative = ""
    if os.path.exists(f"{args.ff}.ff"):
        maybe_relative = "./"

    with open("fepbase.top") as fh, open("topol_solvated.top", "w") as ofh:
        for l in fh:
            ls = l.split()
            if ls == ['[', 'system', ']']:
                print(f'#include "{maybe_relative}{args.ff}.ff/{args.water_model}.itp"', file=ofh)
                print(f'#include "{maybe_relative}{args.ff}.ff/ions.itp"', file=ofh)
            ofh.write(l)
    waterbox = "spc216.gro"
    if args.water_model in ["tip4p", "tip4pew"]:
        waterbox = "tip4p.gro"
    check_call_verbose([args.gmx] + f"solvate -cp conf_box -p topol_solvated -cs {waterbox} -o conf_solvated.pdb".split())
    with open("dummy.mdp", "w") as ofh:
        pass # clear dummy file
    log_with_color("Note: if error about \"no default Bond types\" occurs on next grompp, remove molecules containing [ bonds ] from corresponding \"ions.itp\""
        "or add bond, angle, dihedral parameters. This happens only with recent CHARMM36 parameters")
    check_call_verbose([args.gmx] +  "grompp -f dummy.mdp -p topol_solvated.top -c conf_solvated.pdb -po dummy_out -o topol_solvated -maxwarn 1".split())
    shutil.copy("topol_solvated.top", "topol_ionized.top")
    cmds = [args.gmx] + f"genion -s topol_solvated -o conf_ionized.pdb -p topol_ionized.top -pname {args.ion_positive} -nname {args.ion_negative} -conc {args.ion} -neutral".split()
    log_with_color(" ".join(cmds))
    genion_pipe = subprocess.Popen(cmds, stdin=subprocess.PIPE)
    genion_pipe.communicate(bytes("SOL\n", "utf-8"))

def main(args):
    #pdbseq, pdborder, reschainmap = parse_pdb(basestr)
    basepdbinfo = parse_pdb(args.pdb)
    mutname = mut_name(args.mutation)

    mutation_list = parse_mut(args.mutation)
    mutation_list = update_mutinfo(mutation_list, basepdbinfo)
    generate_gmx(args, basepdbinfo, args.wtdir, [])
    generate_gmx(args, basepdbinfo, mutname, mutation_list)

    # TODO: use better string (especially with chain information)
    fepdir = f"{args.wtdir}_{mutname}"
    refdirs = [f"{fepdir}_ref"]
    if len(mutation_list) > 1:
        refdirs = []
        for i, m in enumerate(mutation_list):
            refdirs.append(f"{fepdir}_ref{i+1}")
    mutdir = f"{mutname}"
    if not os.path.exists(fepdir):
        log_with_color(f"mkdir {fepdir}")
        os.mkdir(fepdir)
    curdir = os.getcwd()
    log_with_color(f"cd {fepdir}")
    os.chdir(fepdir)
    if os.path.exists(f"../{args.ff}.ff") and not os.path.exists(f"{args.ff}.ff"):
        # use local force field file
        log_with_color(f"ln -s ../{args.ff}.ff {args.ff}.ff")
        os.symlink(f"../{args.ff}.ff", f"{args.ff}.ff")
    check_call_verbose([f"{args.nucfepgen}/nucfepgen"] + f"-A ../wt/conf.pdb -B ../{mutdir}/conf.pdb -a ../wt/topol_can.top -b ../{mutdir}/topol_can.top -O fepbase.pdb -o fepbase.top --structureOA fepbase_A.pdb --structureOB fepbase_B.pdb --protein --honor-resnames --generate-restraint CA --disable-cmap-error".split())
    solvate_conf_topology(args.solv)
    generate_para_conf(args, mutation_list)
    log_with_color(f"cd ..")
    os.chdir(curdir)

    for refdir, m in zip(refdirs, mutation_list):
        # FIXME: we do not consider the case that mutations are continuous like 25A_26A
        if not os.path.exists(refdir):
            log_with_color(f"mkdir {refdir}")
            os.mkdir(refdir)
        log_with_color(f"cd {refdir}")
        os.chdir(refdir)
        if os.path.exists(f"../{args.ff}.ff") and not os.path.exists(f"{args.ff}.ff"):
            # use local force field file
            log_with_color(f"ln -s ../{args.ff}.ff {args.ff}.ff")
            os.symlink(f"../{args.ff}.ff", f"{args.ff}.ff")
        check_call_verbose([args.python3, f"{args.feprest}/utils/selectres.py"] + f"../{fepdir}/fepbase.top ../{fepdir}/fepbase.pdb {m['resid']} fepbase.top fepbase.pdb".split())
        solvate_conf_topology(args.solv_ref)
        generate_para_conf(args, [m])
        log_with_color(f"cd ..")
        os.chdir(curdir)


class ArgumentDefaultsHelpFormatterWithRawFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _fill_text(self, text, width, indent):
        return ''.join(indent + line for line in text.splitlines(keepends=True))


def argparse_options():
    parser = argparse.ArgumentParser(description="""Automatic mutation introduction and calculation preparation
    
    The mutation information in --mutation accepts following formats:
    22V      Change 22 to Val
    A22V     Check 22 is Ala and change that to Val
    H:15S    Chain H resid 15 is changed to Ser
    56Y_31F  Two mutations
    """, formatter_class=ArgumentDefaultsHelpFormatterWithRawFormatter)
    gmxbin = None
    if "GMXBIN" in os.environ:
        gmxbin = os.environ["GMXBIN"] + "/gmx"

    fepgenpath = os.path.realpath(os.path.split(os.path.realpath(__file__))[0] + "/../fepgen/fepgen")

    parser.add_argument("--gmx", default=gmxbin, required=(gmxbin is None or not os.path.exists(gmxbin)), help="path to GROMACS gmx")
    parser.add_argument("--faspr", required=True, help="Location of FASPR binary")
    parser.add_argument("--feprest", required=True, help="Path to feprest pipeline directory")
    parser.add_argument("--fepgen", default=(fepgenpath if os.path.exists(fepgenpath) else None), help="Path to fepgen directory")
    parser.add_argument("--python3", default=sys.executable, help="Path to python3")
    
    parser.add_argument("--wtdir", default="wt", help="Where we store wild-type residue informations")
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--mutation", required=True, help="Mutation string")
    parser.add_argument("--solv", default=0.5, type=float, help="Water thickness around the whole system (unit: nm)")
    parser.add_argument("--solv-ref", default=1.2, type=float, help="Water thickness around the reference system (unit: nm)")
    parser.add_argument("--ion", default=0.15, type=float, help="Ion concentration (mol/L)")
    parser.add_argument("--ion-positive", default="NA", type=str, help="Positive Ion name")
    parser.add_argument("--ion-negative", default="CL", type=str, help="Negative Ion name")
    parser.add_argument("--ff", required=True, type=str, help="Force field to use")
    parser.add_argument("--water-model", default="tip3p", type=str, help="Water model")

    parser.add_argument("--difficult-nrep", default=64, help="Num of replica on difficult (involving YFWP or charge-changing) cases")

    return parser.parse_args()

if __name__ == "__main__":
    args = argparse_options()
    main(args)
