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

../nucfepgen -A ../$A.pdb -a ../${A}_GMX.top -B ../$B.pdb -b ../${B}_GMX.top --maxdist $DIST -O init.pdb -o ab.top

# split at resid 6
python ../split_molecule.py --resid 6 --top ab.top --output-1 fep_chain1.top --output-2 fep_chain2.top --pdb init.pdb --pdbout1 init_chain1.pdb --pdbout2 init_chain2.pdb

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
python ../split_states.py --top ab.top --output-A paired-a.top --output-B paired-b.top --output-fep fep.top
solvate fep init
zsh ../solvate.zsh paired-a
cat paired-a_solvated.top <(tail -3 fep_ionized.top) > paired-a_ionized.top
zsh ../solvate.zsh paired-b
cat paired-b_solvated.top <(tail -3 fep_ionized.top) > paired-b_ionized.top

# set up system for chains
python ../split_states.py --top fep_chain1.top --output-A chain1-a.top --output-B chain1-b.top --output-fep chain1_fep.top
solvate fep_chain1 init_chain1
zsh ../solvate.zsh chain1-a
cat chain1-a_solvated.top <(tail -3 fep_ionized.top) > chain1-a_ionized.top
zsh ../solvate.zsh chain1-b
cat chain1-b_solvated.top <(tail -3 fep_ionized.top) > chain1-b_ionized.top

python ../split_states.py --top fep_chain2.top --output-A chain2-a.top --output-B chain2-b.top --output-fep chain2_fep.top
solvate fep_chain2 init_chain2
zsh ../solvate.zsh chain2-a
cat chain2-a_solvated.top <(tail -3 fep_ionized.top) > chain2-a_ionized.top
zsh ../solvate.zsh chain2-b
cat chain2-b_solvated.top <(tail -3 fep_ionized.top) > chain2-b_ionized.top

