import argparse
import copy
from dataclasses import dataclass
import os
import os.path
import re
import shutil
import subprocess
import sys
import warnings
from typing import Dict, List, Optional, Sequence, Tuple, Union

from mutation import HashableMutations, Mutation, parse_mutations

@dataclass
class PDBInfoSummary:
    seq: str
    order: List[Tuple[str, int]]
    res_chain_dic: Dict[int, Dict[str, str]] # accessible through [resid][chain], allowing wildcard chain search

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

def seq1(aa: str):
    aa = aa.upper()
    aa_map = {
        'ALA': 'A',
        'CYS': 'C',
        'ASP': 'D',
        'GLU': 'E',
        'PHE': 'F',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LYS': 'K',
        'LEU': 'L',
        'MET': 'M',
        'ASN': 'N',
        'PRO': 'P',
        'GLN': 'Q',
        'ARG': 'R',
        'SER': 'S',
        'THR': 'T',
        'VAL': 'V',
        'TRP': 'W',
        'TYR': 'Y',
        'HID': 'H', 'HIE': 'H', 'HIP': 'H',
        'HSD': 'H', 'HSE': 'H', 'HSP': 'H',
        'GLH': 'E', 'ASH': 'D', # FIXME: charmm uses ASPP and GLUP but what should I do here?
        'LYN': 'K', 'LSN': 'K',
        'CYX': 'C', 'CYM': 'C'
    }
    if aa in aa_map:
        return aa_map[aa]
    raise RuntimeError("Unknown residue name " + aa)

def parse_pdb(pdb) -> PDBInfoSummary:
    with open(pdb) as fh:
        seq = ''
        reschainmap = {}
        order = []
        for l in fh:
            # note this also ignores HETATM
            if not l.startswith('ATOM  '):
                continue
            # From PDB specification: (note columns starts from 1 in spec)
            # COLUMNS        DATA  TYPE    FIELD        DEFINITION
            # -------------------------------------------------------------------------------------
            #  1 -  6        Record name   "ATOM  "
            #  7 - 11        Integer       serial       Atom  serial number.
            # 13 - 16        Atom          name         Atom name.
            # 17             Character     altLoc       Alternate location indicator.
            # 18 - 20        Residue name  resName      Residue name.
            # 22             Character     chainID      Chain identifier.
            # 23 - 26        Integer       resSeq       Residue sequence number.
            # 27             AChar         iCode        Code for insertion of residues.
            name = l[12:16].strip()
            if name != "CA":
                continue
            chain_id = l[21]
            aaa = l[17:20].strip()
            if l[20] != " ":
                warnings.warn("Residue name may be four-lettered (nonstandard extension of PDB) but the program currently does not support")
            a = seq1(aaa)
            seq += str(a)
            resid = int(l[22:26].strip())
            if resid not in reschainmap:
                reschainmap[resid] = {}
            if chain_id in reschainmap[resid]:
                raise RuntimeError(f"Duplicated Chain-resid combination at {chain_id}-{resid}")
            reschainmap[resid][chain_id] = a
            order.append((chain_id, resid))
    return PDBInfoSummary(seq=seq, order=order, res_chain_dic=reschainmap) # returns at the first model

def update_mutinfo(mutations: Sequence[Mutation], basepdbinfo: PDBInfoSummary) -> List[Mutation]:
    ret_mutations = []
    pdborder = basepdbinfo.order
    reschainmap = basepdbinfo.res_chain_dic
    for (c, r) in pdborder:
        res = reschainmap[r][c]
        for m in mutations:
            if m.resid == r and (m.chain is None or m.chain == c):
                # matched mutation
                if m.before_res is not None and m.before_res != res:
                    raise RuntimeError(f"Mutation ({m}) is not compatible with chain '{c}' resid {r}")
                newmut = copy.copy(m)
                newmut.before_res = res
                ret_mutations.append(newmut)
    if len(ret_mutations) != len(mutations):
        raise RuntimeError(f"Not all mutations are consumed, requested mutations: {mutations}, but found: {ret_mutations}")
    return ret_mutations

class DefaultProteinMutationGenerator:
    pdb: str

    # software paths
    gmx: str
    python3: str
    faspr: str
    fepsuite: str

    # ff, water and ions
    ff: str
    water_model: str
    solv: float
    solv_ref: float
    ion: float
    ion_positive: str
    ion_negative: str

    # tunable parameters
    difficult_nrep: int = 64

    def __init__(self, pdb: str, gmx: str, python3: str, faspr: str, fepsuite: str,
                 ff: str, water_model: str,
                 solv: float, solv_ref: float,
                 ion: float, ion_positive: str, ion_negative: str,
                 difficult_nrep: int):
        self.pdb = pdb
        self.gmx = gmx
        self.python3 = python3
        self.faspr = faspr
        self.fepsuite = fepsuite

        self.ff = ff
        self.water_model = water_model
        self.solv = solv
        self.solv_ref = solv_ref
        self.ion = ion
        self.ion_positive = ion_positive
        self.ion_negative = ion_negative

        self.difficult_nrep = difficult_nrep

    def generate_gmx(self, basepdbinfo, ddir, mutations: Sequence[Mutation]):
        if not os.path.exists(ddir):
            os.mkdir(ddir)
        else:
            if os.path.exists(f"{ddir}/conf.pdb") and os.path.exists(f"{ddir}/topol_can.top") and os.path.exists(f"{ddir}/gmx_success"):
                if os.path.getmtime(f"{ddir}/gmx_success") > os.path.getmtime(self.pdb):
                    print(f"Directory {ddir} has been generated successfully in previous run, skipping the structure/topology generation")
                    return # already done
                else:
                    print(f"Directory {ddir} has been generated successfully in the previous run, but base PDB's access time was newer, thus we restart the structure/topoology generation")
        if mutations == []:
            # wild-type, just copy
            self.generate_wt(self.pdb, f"{ddir}/completed.pdb")
        else:
            self.generate_mutant(basepdbinfo, ddir, mutations)

        curdir = os.getcwd()
        os.chdir(ddir)
        if os.path.exists(f"../{self.ff}.ff") and not os.path.exists(f"{self.ff}.ff"):
            # use local force field file if necessary
            os.symlink(f"../{self.ff}.ff", f"{self.ff}.ff")

        # Generate GMX topology as well as structures
        self.generate_gmx_top()

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
        check_call_verbose([self.python3, f"{self.fepsuite}/feprest/tools/pp.py", "-o", "topol_pp.top", "topol_m.top"])
        check_call_verbose([self.python3, f"{self.fepsuite}/feprest/rest2py/canonicalize_top.py", "topol_pp.top", "topol_can.top"])

        # write a blank file
        with open("gmx_success", "w") as _ofh:
            pass

        os.chdir(curdir)
        return

    def generate_gmx_top(self):
        """From conf.pdb, generate completed.pdb and topol.top in the current directory"""
        check_call_verbose([self.gmx] + f"pdb2gmx -f completed.pdb -o conf.pdb -water {self.water_model} -ignh -ff {self.ff} -merge all".split())

    def generate_mutant(self, basepdbinfo: PDBInfoSummary, ddir: str, mutations: Sequence[Mutation]):
        with open(f"{ddir}/seq.txt", "w") as ofh:
            pdborder = basepdbinfo.order
            reschainmap = basepdbinfo.res_chain_dic
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
                    print(repr((c,r)), repr(m))
                    if m.resid == r and (m.chain is None or m.chain == c):
                        if m.before_res is not None and m.before_res != res: # this check should not be necessary here
                            raise RuntimeError(f"Mutation ({m}) is not compatible with chain '{c}' resid {r}")
                        res = m.after_res
                        mutated = True
                        totmut += 1
                        break

                # Hack for FASPR: keep Cys rotamer intact
                if res == "C" and not mutated:
                    res = "c" # keep Cys rotamer so that disulfide bonds are intact

                print(res, end="", file=ofh)
            print("", file=ofh)
            if totmut != len(mutations):
                raise RuntimeError(f"Not all mutations are consumed, requested mutations: {mutations}")

        check_call_verbose([self.faspr, "-i", self.pdb, "-o", f"{ddir}/completed.pdb", "-s", f"{ddir}/seq.txt"])
        # faspr often returns 0 even if there is an error
        if not os.path.exists(f"{ddir}/completed.pdb") or \
                os.path.getmtime(f"{ddir}/completed.pdb") < os.path.getmtime(f"{ddir}/seq.txt"):
            raise RuntimeError("FASPR run failed")

    def generate_wt(self, input, output):
        shutil.copy(input, output)

    @staticmethod
    def difficulty(muts: Sequence[Mutation]):

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
            mutfrom = m.before_res
            mutto = m.after_res
            totcharge -= charge(mutfrom)
            totcharge += charge(mutto)
            if mutto in "PFYW" or mutfrom in list("PFYW"):
                difficult = True
                if mutto in "FY" and mutfrom in list("FY"):
                    difficult = False
        return (difficult, totcharge != 0)

    def generate_para_conf(self, muts: Sequence[Mutation]):
        (difficult, charged) = self.difficulty(muts)
        if charged:
            difficult = True
        if difficult:
            with open("para_conf.zsh", "w") as ofh:
                print(f"NREP={self.difficult_nrep}", file=ofh)

    @staticmethod
    def mut_name(mutstr: str):
        return mutstr.replace(":", "_")

    def solvate_curdir_conf_topology(self, radius: float):
        """Convert conf.pdb + topol.top into conf_ionized.pdb and topol_ionized.top"""
        check_call_verbose([self.gmx] + f"editconf -f fepbase.pdb -d {radius} -bt dodecahedron -o conf_box".split())
        maybe_relative = ""
        if os.path.exists(f"{self.ff}.ff"):
            maybe_relative = "./"

        with open("fepbase.top") as fh, open("topol_solvated.top", "w") as ofh:
            for l in fh:
                ls = l.split()
                if ls == ['[', 'system', ']']:
                    print(f'#include "{maybe_relative}{self.ff}.ff/{self.water_model}.itp"', file=ofh)
                    print(f'#include "{maybe_relative}{self.ff}.ff/ions.itp"', file=ofh)
                ofh.write(l)
        waterbox = "spc216.gro"
        if self.water_model in ["tip4p", "tip4pew"]:
            waterbox = "tip4p.gro"
        check_call_verbose([self.gmx] + f"solvate -cp conf_box -p topol_solvated -cs {waterbox} -o conf_solvated.pdb".split())
        with open("dummy.mdp", "w") as ofh:
            pass # clear dummy file
        log_with_color("Note: if error about \"no default Bond types\" occurs on next grompp, remove molecules containing [ bonds ] from corresponding \"ions.itp\""
            "or add bond, angle, dihedral parameters. This happens only with recent CHARMM36 parameters")
        check_call_verbose([self.gmx] +  "grompp -f dummy.mdp -p topol_solvated.top -c conf_solvated.pdb -po dummy_out -o topol_solvated -maxwarn 1".split())
        shutil.copy("topol_solvated.top", "topol_ionized.top")
        cmds = [self.gmx] + f"genion -s topol_solvated -o conf_ionized.pdb -p topol_ionized.top -pname {self.ion_positive} -nname {self.ion_negative} -conc {self.ion} -neutral".split()
        log_with_color(" ".join(cmds))
        genion_pipe = subprocess.Popen(cmds, stdin=subprocess.PIPE)
        genion_pipe.communicate(bytes("SOL\n", "utf-8"))

    def generate(self, muts:Union[str, HashableMutations], muts_name: Optional[str], base_muts: Union[str, HashableMutations] = "", base_name: Optional[str] = None, base_mut_separator: str = "_") -> Tuple[str, List[str]]:
        #pdbseq, pdborder, reschainmap = parse_pdb(basestr)
        basepdbinfo = parse_pdb(self.pdb)

        fepdir = []
        mutation_lists = []
        for (_mode, dirname, mutstr_or_data) in \
            [ ("base", base_name, base_muts),
              ("mutant", muts_name, muts) ]:
            if type(mutstr_or_data) == str:
                if mutstr_or_data == "":
                    mutation_list = []
                else:
                    mutation_list = parse_mutations(mutstr_or_data)
            else:
                mutation_list = mutstr_or_data.to_mutations()
            mutation_list = update_mutinfo(mutation_list, basepdbinfo)
            if dirname is None:
                dirname = "_".join([str(m) for m in mutation_list])
                if mutation_list == []:
                    dirname = "wt"
            if dirname == "":
                dirname = "wt"
            self.generate_gmx(basepdbinfo, dirname, mutation_list)
            fepdir.append(dirname)
            mutation_lists.append(mutation_list)
        
        mutation_list_diff = list(set(mutation_lists[1]) - set(mutation_lists[0]))
        inv_mutation_list = list(set(mutation_lists[0]) - set(mutation_lists[1]))
        if inv_mutation_list != []:
            raise RuntimeError("Implement the case when reverse mutations must be considered")
            # TODO FIXME: in addition to the reverse mutation, also consider the case like base: 1A mut: 1L

        # At this moment two directories with conf/top were generated. Start generating FEP'ed structure
        fepdirname = base_mut_separator.join(fepdir)
        refdirs = [f"{fepdirname}_ref"]
        if len(mutation_list_diff) > 1:
            refdirs = []
            for i, m in enumerate(mutation_list_diff):
                refdirs.append(f"{fepdirname}_ref{i+1}")
        log_with_color(f"mkdir {fepdirname}")
        os.makedirs(fepdirname, exist_ok=True)
        curdir = os.getcwd()
        log_with_color(f"cd {fepdirname}")
        os.chdir(fepdirname)
        if os.path.exists(f"../{self.ff}.ff") and not os.path.exists(f"{self.ff}.ff"):
            # use local force field file
            log_with_color(f"ln -s ../{self.ff}.ff {self.ff}.ff")
            os.symlink(f"../{self.ff}.ff", f"{self.ff}.ff")
        check_call_verbose([f"{self.fepsuite}/feprest/fepgen/fepgen"] + f"-A ../{fepdir[0]}/conf.pdb -B ../{fepdir[1]}/conf.pdb -a ../{fepdir[0]}/topol_can.top -b ../{fepdir[1]}/topol_can.top -O fepbase.pdb -o fepbase.top --structureOA fepbase_A.pdb --structureOB fepbase_B.pdb --protein --honor-resnames --generate-restraint CA --disable-cmap-error".split())
        self.solvate_curdir_conf_topology(self.solv)
        self.generate_para_conf(mutation_list_diff)
        log_with_color(f"cd ..")
        os.chdir(curdir)

        # Do the same for the reference systems
        for refdir, m in zip(refdirs, mutation_list_diff):
            # FIXME: we currently do not consider the case that mutations are continuous like 25A_26A
            log_with_color(f"mkdir {refdir}")
            os.makedirs(refdir, exist_ok=True)
            log_with_color(f"cd {refdir}")
            os.chdir(refdir)
            if os.path.exists(f"../{self.ff}.ff") and not os.path.exists(f"{self.ff}.ff"):
                # use local force field file
                log_with_color(f"ln -s ../{self.ff}.ff {self.ff}.ff")
                os.symlink(f"../{self.ff}.ff", f"{self.ff}.ff")
            check_call_verbose([self.python3, f"{self.fepsuite}/feprest/tools/selectres.py"] + f"../{fepdirname}/fepbase.top ../{fepdirname}/fepbase.pdb {m.resid} fepbase.top fepbase.pdb".split())
            self.solvate_curdir_conf_topology(self.solv_ref)
            self.generate_para_conf([m])
            log_with_color(f"cd ..")
            os.chdir(curdir)
        
        return fepdirname, refdirs

class ArgumentDefaultsHelpFormatterWithRawFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _fill_text(self, text, _width, indent):
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

    fepsuitepath = os.path.realpath(os.path.split(os.path.realpath(__file__))[0] + "/../..")

    parser.add_argument("--gmx", default=gmxbin, required=(gmxbin is None or not os.path.exists(gmxbin)), help="path to GROMACS gmx")
    parser.add_argument("--faspr", required=True, help="Location of FASPR binary")
    parser.add_argument("--fepsuite", default=os.path.realpath(fepsuitepath), help="Path to fepsuite directory")
    parser.add_argument("--python3", default=sys.executable, help="Path to python3")
    
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--base-mutation", default="", help="Mutation string for base mutation, default=\"\" (wild-type)")
    parser.add_argument("--mutation", required=True, help="Mutation string")
    parser.add_argument("--solv", default=0.5, type=float, help="Water thickness around the whole system (unit: nm)")
    parser.add_argument("--solv-ref", default=1.2, type=float, help="Water thickness around the reference system (unit: nm)")
    parser.add_argument("--ion", default=0.15, type=float, help="Ion concentration (mol/L)")
    parser.add_argument("--ion-positive", default="NA", type=str, help="Positive Ion name")
    parser.add_argument("--ion-negative", default="CL", type=str, help="Negative Ion name")
    parser.add_argument("--ff", required=True, type=str, help="Force field to use")
    parser.add_argument("--water-model", default="tip3p", type=str, help="Water model")
    parser.add_argument("--base-mut-separator", default="_", type=str, help="Separate base mutation and new mutation with this symbol")

    parser.add_argument("--difficult-nrep", default=64, help="Num of replica on difficult (involving YFWP or charge-changing) cases")

    return parser.parse_args()

if __name__ == "__main__":
    args = argparse_options()
    generator = DefaultProteinMutationGenerator(args.pdb, args.gmx, args.python3, args.faspr, args.fepsuite,
                                                args.ff, args.water_model,
                                                args.solv, args.solv_ref, args.ion, args.ion_positive, args.ion_negative,
                                                args.difficult_nrep)
    dirname, refdirs = generator.generate(args.mutation, args.mutation, base_muts=args.base_mutation, base_name=args.base_mutation, base_mut_separator=args.base_mut_separator)
    print(f"Successfully generated FEP-ed structure on \"{dirname}\" and {refdirs}")