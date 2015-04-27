#!/bin/zsh

A=$1
B=$2
DIST=$3

if [[ -z $B ]]; then
    echo "Usage: $0 TRUNK"
    exit 1
fi

if [[ -z $DIST ]]; then
    DIST=1.0
fi

#trap 'echo \"Error at line $LINENO \!\"; exit 1' ERR

../nucfepgen -A ../$A.pdb -a ../${A}_GMX.top -B ../$B.pdb -b ../${B}_GMX.top --maxdist $DIST -O init.pdb -o na.top

solvate(){
    TRUNK=$1
    PDBTRUNK=$2
    gmx_d editconf -f ${PDBTRUNK}.pdb -d 1.5 -bt dodecahedron -o ${PDBTRUNK}_box
    zsh ../solvate.zsh ${TRUNK}
    gmx_d solvate -cp ${PDBTRUNK}_box -p ${TRUNK}_solvated -cs -o ${PDBTRUNK}_solvated
    ln -s ../*.itp .
    touch dummy.mdp
    gmx_d grompp -f dummy.mdp -p ${TRUNK}_solvated.top -c ${PDBTRUNK}_solvated.gro -po dummy_out -o ${TRUNK}_solvated
    cp ${TRUNK}_solvated.top ${TRUNK}_ionized.top
    echo SOL | gmx_d genion -s ${TRUNK}_solvated -o ${TRUNK}_ionized.gro -p ${TRUNK}_ionized.top -pname NA -nname CL -conc 1.0 -neutral || exit 1
    return 0
}

# set up system for paired fep
solvate na init
#gmx_d editconf -f na_ionized.gro -o na_ionized.pdb
python ../split_molecule.py --resid 6 --top na_ionized.top --output-1 na_chain1.top --output-2 na_chain2.top --pdb na_ionized.gro --pdbout1 fep_chain1.gro --pdbout2 fep_chain2.gro
python ../split_states.py --top na_chain1.top --output-A a_chain1.top --output-B b_chain1.top --output-fep fep_chain1.top
python ../split_states.py --top na_chain2.top --output-A a_chain2.top --output-B b_chain2.top --output-fep fep_chain2.top
python ../split_states.py --top na_ionized.top --output-A a_paired.top --output-B b_paired.top --output-fep fep_paired.top

#zsh ../solvate.zsh paired-a
#cat paired-a_solvated.top <(tail -3 fep_ionized.top) > paired-a_ionized.top
#zsh ../solvate.zsh paired-b
#cat paired-b_solvated.top <(tail -3 fep_ionized.top) > paired-b_ionized.top

# split at resid 6

# set up system for chains
#python ../split_states.py --top fep_chain1.top --output-A chain1-a.top --output-B chain1-b.top --output-fep chain1_fep.top
#python ../split_states.py --top fep_chain2.top --output-A chain2-a.top --output-B chain2-b.top --output-fep chain2_fep.top
