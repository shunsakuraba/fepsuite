#!/usr/bin/zsh

# Constraint optimization.
# This optimization has two objectives:
# 1. Determine the structure of the molecule for the charge assignment.
# 2. Determine the angle for constraining in later calculations

if [[ -z $2 ]]; then
    echo "Usage: $0 (structure) (constraint file)" 1>&2
    exit 1
fi

basedir=${0:h}
basestructure=$1
baseconstraint=$2

if [[ -z $CHARGE ]]; then
    CHARGE=0
fi

source $basedir/defaults.zsh

trap 'echo "Error returned at previous execution"; exit 1' ZERR
set -x

basestructuretype=${basestructure:e}
if [[ $basestructuretype = "gau" ]]; then
    echo "Base structure file name must not be gau" 1>&2
    exit 1
fi
basestructurename=${basestructure:r}
basename=${basestructure:t:r}

initialgau=${basestructurename}.init.gau

opt1gau=${basestructurename}.opt1.gau
opt1check=${basename}.opt1.chk
opt1log=${basestructurename}.opt1.log

# Use babel, not antechamber (antechamber may swap atoms and emits no information for recovering atom names!)
$OPENBABEL $basestructure -ogzmat -xk "%chk=$opt1check
# HF/6-31G POPT(Z-matrix,maxcycle=5000)" $initialgau

# Fix 'Holmium' problem and broken multiplicity
sed -i "s/Ho/H/;/^0 /c $CHARGE 1" $initialgau

python $basedir/mod-zmatrix.py $initialgau $basestructure $baseconstraint > $opt1gau

# Optimize the structure with constraint
# HF/6-31G level first to remove large crashes
zsh $basedir/g09run.zsh $basename $opt1gau

zsh $basedir/g09fetch.zsh $basename $opt1log


opt2tmol=${basestructurename}.opt2.tmol
opt2dihf=${basestructurename}.opt2.dih
optfintmol=${basestructurename}.opt2fin.tmol
rundir=${basename}/opt2

# Switch to Turbomole
# Define constraints
python $basedir/constrain-tm.py $opt1log $basestructure $constraintfile $opt2tmol $opt2dihf

# PBE/TZVPP (as in Zgarbova et al. 2011)
zsh $basedir/turborun.zsh $rundir TITLE="Optimization" CONSTRAINTS=$opt2dihf BASIS="TZVPP" CHARGE=$CHARGE FUNCTIONAL=pbe GRID=m4 COSMO=78.4 DISP3=bj ENERGY_CONV=7 COORD_CONV=4 CYCLE=1000 $opt2tmol 

zsh $basedir/turbofetch.zsh $rundir $optfintmol

