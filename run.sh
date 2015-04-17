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

../nucfepgen -A ../$A.pdb -a ../${A}_GMX.top -B ../$B.pdb -b ../${B}_GMX.top --maxdist $DIST -O init.pdb -o ab.top
python ../split_states.py --top ab.top --output-A a.top --output-B b.top --output-fep fep.top || exit 1
gmx_d editconf -f init.pdb -d 1.5 -bt dodecahedron -o init_box
zsh ../solvate.zsh fep
gmx_d solvate -cp init_box -p fep_solvated -cs -o init_solvated
ln -s ../*.itp .
touch dummy.mdp
gmx_d grompp -f dummy.mdp -p fep_solvated.top -c init_solvated.gro -po dummy_out -o fep_solvated
cp fep_solvated.top fep_ionized.top
echo SOL | gmx_d genion -s fep_solvated -o fep_ionized.gro -p fep_ionized.top -pname NA -nname CL -conc 1.0 -neutral
