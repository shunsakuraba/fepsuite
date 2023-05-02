from dataclasses import dataclass
import re
from typing import Dict, List, Optional, Sequence, Tuple

@dataclass(unsafe_hash=True, order=True, frozen=True)
class MutationKey:
    chain: Optional[str]
    resid: int

    def __str__(self):
        if self.chain is None:
            return str(self.resid)
        else:
            return f"{self.chain}:{self.resid}"

@dataclass(unsafe_hash=True, order=True)
class Mutation:
    chain: Optional[str]
    resid: int
    before_res: Optional[str]
    after_res: str

    def to_key_value(self) -> Tuple[MutationKey, str]:
        return (MutationKey(self.chain, self.resid), self.after_res)

    def __str__(self):
        if self.chain is None:
            chain = ""
        else:
            chain = self.chain + ":"
        before_res = self.before_res if self.before_res is not None else ""
        return f"{chain}{before_res}{self.resid}{self.after_res}"
    
MutationsDict = Dict[MutationKey, str]

class HashableMutations:
    d: MutationsDict
    s: str

    def __init__(self, din: MutationsDict):
        self.d = din
        self.s = self.canonical_str(self.d)

    def to_mutations(self) -> List[Mutation]:
        ret: List[Mutation] = []
        for k, v in self.d.items():
            ret.append(Mutation(k.chain, k.resid, None, v))
        ret.sort()
        return ret

    @staticmethod
    def canonical_str(m: MutationsDict) -> str:
        ret = list(m.items())
        ret.sort(key=lambda x: x[0])
        ret = [f"{str(key)}{mut}" for (key, mut) in ret]

        if ret == []:
            return "wt"
        else:
            return "_".join(ret)
    
    def __hash__(self):
        return hash(self.s)
    
    def __eq__(self, other):
        return self.s == other.s
    
    def __str__(self) -> str:
        return self.s
    
    def __repr__(self) -> str:
        return self.s
    
    def __lt__(self, other):
        return self.s < other.s
    
    def __gt__(self, other):
        return self.s > other.s

_mutation_pattern_seq = re.compile(r"(?:(?P<chain>[A-Za-z]):)?(?P<before>[A-Z])?(?P<resid>\d+)(?P<after>[A-Z])$")
_mutation_pattern_seq_nonsingle = re.compile(r"((?P<chain>[A-Za-z]):)?(?P<before>[A-Z])?(?P<resid>\d+)(?P<after>[A-Z]*)$")

def parse_single_mutation(mu: str, allow_nonsingle_after: bool = False) -> Mutation:
    patseq = _mutation_pattern_seq_nonsingle if allow_nonsingle_after else _mutation_pattern_seq
    matcher = patseq.match(mu)
    if matcher is None:
        raise ValueError(f"Mutation string parse failed around \"{mu}\"")
    d = matcher.groupdict()

    chain = None
    if "chain" in d:
        # remove ":"
        chain = d['chain']
    resid = int(d['resid'])
    before_res = None
    if "before" in d:
        before_res = d['before']
    return Mutation(chain, resid, before_res, after_res=d['after'])

def parse_mutations(mutstr: str, allow_nonsingle_after: bool = False) -> List[Mutation]:
    """Parse mutation string into the list of single-point mutations.

    Accepted syntax example include:
    22A
    A:V23P
    g:S152K
    23K_B:V21F
    """

    muts = mutstr.split("_")
    return [parse_single_mutation(m, allow_nonsingle_after) for m in muts]

def to_dict(ms: Sequence[Mutation]) -> MutationsDict:
    return {MutationKey(m.chain, m.resid): m.after_res for m in ms}