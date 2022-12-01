import argparse
import os
import re
from typing import List, Tuple, Set

# filename and line is unnecessary but for debugging...
class IfInfo:
    defname: str = None
    filename: str = None
    line: int = None
    value: bool = None
    state: str = None

    def __str__(self) -> str:
        return f"file {self.filename} line {self.line} if(n)def {self.defname} value {self.value} state {self.state}"

def eval_if(ifinfo: IfInfo):
    if ifinfo.state == "ifdef":
        return ifinfo.value
    elif ifinfo.state == "else":
        return not ifinfo.value
    else:
        raise RuntimeError("Unknown if state")

class State:
    ofh = None
    includes: List[str] = []
    defines: Set[str] = set()
    ifstack: List[IfInfo] = []

def find_include(basepath, incpath, includespec):
    includespec = includespec.strip()
    good = False
    for [l, r] in [['"', '"'], ['<', '>']]:
        if includespec.startswith(l):
            if includespec.endswith(r):
                good = True
                break
    if not good:
        raise RuntimeError(f"Failed to parse '{includespec}'")
    
    path = includespec[1:-1]
    if path.startswith("/"): # absolute
        incpathfind = [""]
    else:
        incpathfind = [os.path.split(basepath)[0] + "/"] 
        for inc in incpath:
            incpathfind.append(inc.rstrip("/") + "/")

    for inc in incpathfind:
        fpath = inc + path
        assert os.path.isabs(fpath)
        if os.path.isfile(fpath):
            return fpath
    raise RuntimeError(f"Failed to find include file: \'{path}\' from " + str(incpath))

def process(inputfile, state: State):             
    with open(inputfile) as fh:
        for (iline, l) in enumerate(fh):
            line = iline + 1 # 1-origin
            ifenable = all([eval_if(x) for x in state.ifstack])
            if l.startswith("#"):
                orders = l[1:]
                orders = orders.lstrip()
                # ifdef/else/endif and #include are processed regardless of enabling state
                if orders.startswith("ifdef") or orders.startswith("ifndef"):
                    if orders.startswith("ifdef"):
                        orders = orders[5:].split()
                        basestate = True
                    elif orders.startswith("ifndef"):
                        orders = orders[6:].split()
                        basestate = False                        
                    if len(orders) > 1:
                        raise RuntimeError("In GROMACS, #ifdef/#ifndef only accepts 1 argument")   
                    ifinfo = IfInfo()
                    ifinfo.defname = orders[0]
                    ifinfo.filename = inputfile
                    ifinfo.line = line
                    ifinfo.value = (orders[0] in state.defines) == basestate 
                    ifinfo.state = "ifdef"
                    state.ifstack.append(ifinfo)
                    continue
                elif orders.startswith("else"):
                    stacktop = state.ifstack.pop()
                    stacktop.state = "else"
                    state.ifstack.append(stacktop)
                    continue
                elif orders.startswith("endif"):
                    state.ifstack.pop()
                    continue
                elif orders.startswith("include"):
                    target = find_include(inputfile, state.includes, orders[7:])
                    process(target, state)
                    continue
                if not ifenable: # define, undef are processed only when enabled
                    continue
                if orders.startswith("define"):
                    orders = orders[6:].split()
                    if len(orders) > 1:
                        raise RuntimeError("In GROMACS, #define only accepts 1 argument, but was: " + str(orders))
                    state.defines.add(orders[0])
                elif orders.startswith("undef"):
                    orders = orders[5:].split()
                    if len(orders) > 1:
                        raise RuntimeError("#undef only accepts 1 argument")
                    state.defines.remove(orders[0])
                else:
                    raise RuntimeError(f"Unknown # directive: {orders.rstrip()}")
            else:
                if ifenable:
                    state.ofh.write(l)

def main(args):
    state = State()
    state.defines = set(args.defines)
    state.includes = args.include
    with open(args.output, "w") as ofh:
        state.ofh = ofh
        process(os.path.abspath(args.input), state)


class ArgumentDefaultsHelpFormatterWithRawFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _fill_text(self, text, width, indent):
        return ''.join(indent + line for line in text.splitlines(keepends=True))

def argparse_options():
    parser = argparse.ArgumentParser(description="""Preprocess top file

    This program preprocesses .top file so that #include and #ifdef are resolved (in mostly GROMACS order)
    """, formatter_class=ArgumentDefaultsHelpFormatterWithRawFormatter)

    parser.add_argument("--include", "-I", default=[], action="append", help="Path to GROMACS topology directory (you should probably source GMXRC instead of specifying this directory)")

    parser.add_argument("-D", dest="defines", default=[], action="append", help="Defined macros")
    parser.add_argument("-o", dest="output", required=True, type=str, help="Output file")
    parser.add_argument("input", metavar="FILE", type=str, help="Input file")

    return parser.parse_args()

if __name__ == "__main__":
    args = argparse_options()

    if args.include == [] and "GMXDATA" in os.environ:
        gmxdatad = os.path.join(os.environ["GMXDATA"], "top")
        if os.path.isdir(gmxdatad):
            args.include = [gmxdatad]
    main(args)
