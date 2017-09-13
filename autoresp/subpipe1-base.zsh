#!/usr/bin/zsh

if [[ -z $1 ]]; then
    echo "Usage: $0 (molecule trunk)" 1>&2
    echo "Example: $0 inosine" 1>&2
    exit 1
fi

basedir=${0:h}
trunk=$1
basestructure=$trunk.base.pdb
backbone=$basedir/RNA.pdb
backbonesugared=$basedir/RNA_sugar.pdb
backbone_D=$basedir/DNA.pdb
backbonesugared_D=$basedir/DNA_sugar.pdb
basestructurename=$trunk
basestructurename_D=$trunk.DNA

python $basedir/merger.py $basestructure $backbone $basestructurename.monomer.pdb
python $basedir/merger.py $basestructure $backbonesugared $basestructurename.pdb

python $basedir/merger.py $basestructure $backbone_D $basestructurename_D.monomer.pdb
python $basedir/merger.py $basestructure $backbonesugared_D $basestructurename_D.pdb

[[ -e vars.zsh ]] && mv -f vars.zsh vars.zsh.bak
echo "basestructurename=$trunk" > vars.zsh
echo "CHARGE=${CHARGE:-0}" >> vars.zsh
