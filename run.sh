#!/bin/zsh

TRUNK=$1

if [[ -z $TRUNK ]]; then
    echo "Usage: $0 TRUNK"
    exit 1
fi

../nucfepgen -A ../$TRUNK.pdb -a ../${TRUNK}_GMX.top -B ../$TRUNK.pdb -b ../${TRUNK}_GMX.top --maxdist 1.0 -O init.pdb -o ab.top
python ../split_states.py --top ab.top --output-A a.top --output-B b.top --output-fep fep.top
gmx_d editconf -f init.pdb -d 1.5 -bt dodecahedron -o init_box
zsh ../solvate.zsh fep
gmx_d solvate -cp init_box -p fep_solvated -cs -o init_solvated
ln -s ../*.itp .
touch dummy.mdp
gmx_d grompp -f dummy.mdp -p fep_solvated.top -c init_solvated.gro -po dummy_out -o fep_solvated
cp fep_solvated.top fep_ionized.top
echo SOL | gmx_d genion -s fep_solvated -o fep_ionized.gro -p fep_ionized.top -pname NA -nname CL -conc 1.0 -neutral
